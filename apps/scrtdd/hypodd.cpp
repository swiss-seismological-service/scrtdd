/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 *   You can redistribute and/or modify this program under the             *
 *   terms of the SeisComP Public License.                                 *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   SeisComP Public License for more details.                             *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#include "hypodd.h"
#include "csvreader.h"

#include <seiscomp3/io/recordinput.h>
#include <seiscomp3/core/typedarray.h>
#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/filter/butterworth.h>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/datamodel/station.h>
#include <seiscomp3/datamodel/event.h>
#include <seiscomp3/datamodel/pick.h>
#include <seiscomp3/datamodel/origin.h>
#include <seiscomp3/datamodel/magnitude.h>
#include <seiscomp3/datamodel/databasequery.h>
#include <seiscomp3/utils/files.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <boost/filesystem.hpp>

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp3/logging/log.h>
 
using namespace std;

namespace {

pid_t startExternalProcess(const vector<string> &cmdparams)
{
	pid_t pid;
	string cmdline;

	pid = fork();

	if ( pid < 0 ) // parent
    {
		return pid;
    }
	else if ( pid == 0 ) // child
    {
		vector<char *> params(cmdparams.size());
		for ( size_t i = 0; i < cmdparams.size(); ++i )
        {
			params[i] = (char*)cmdparams[i].c_str();
			if ( i > 0 ) cmdline += " ";
			cmdline += cmdparams[i];
		}
		params.push_back(NULL);

		SEISCOMP_DEBUG("$ %s", cmdline.c_str());
		execv(params[0], &params[0]);
		exit(1);
	}

	return pid;
}

}


namespace Seiscomp {

BgCatalog::BgCatalog(const string& idFile, DataModel::DatabaseQuery* query)
{
	if ( !Util::fileExists(idFile) )
	{
		string msg = "File " + idFile + " does not exist";
		throw runtime_error(msg);
	}

	vector<string> ids;
	vector< map<string,string> > rows = CSV::readWithHeader(idFile, {"id"});

	for(const auto& row : rows)
	{
		const string& id = row.at("id");
		ids.push_back(id);
	}

	initFromIds(ids, query);
}


BgCatalog::BgCatalog(const vector<string>& ids, DataModel::DatabaseQuery* query)
{
	initFromIds(ids, query);
}


void BgCatalog::initFromIds(const vector<string>& ids, DataModel::DatabaseQuery* query)
{
	int eventId = 1;

	for(const string& id : ids)
	{
		DataModel::OriginPtr org = DataModel::Origin::Cast(query->getObject(DataModel::Origin::TypeInfo(), id));
		if ( !org )
		{
			DataModel::EventPtr ev = DataModel::Event::Cast(query->getObject(DataModel::Event::TypeInfo(), id));
			if ( ev )
			{
				org = DataModel::Origin::Cast(query->getObject(DataModel::Origin::TypeInfo(), (ev->preferredOriginID())));
			}
		}

		if ( !org )
		{
			string msg = "Cannot find origin/event with id " + id;
			throw runtime_error(msg);
		}

		query->loadArrivals(org.get());

		int year, month, day, hour, min, sec, usec;
		if ( ! org->time().value().get(&year, &month, &day, &hour, &min, &sec, &usec) )
		{
			string msg = "Cannot fetch origin time for origin " + org->publicID();
			throw runtime_error(msg);
		}

		// Add event
		Event ev;
		ev.id          = eventId++;
		ev.year        = year;
		ev.month       = month;
		ev.day         = day;
		ev.hour        = hour;
		ev.minute      = min;
		ev.second      = sec + double(usec) / 1e6;
		ev.latitude    = org->latitude();
		ev.longitude   = org->longitude();
		ev.depth       = org->depth(); // it `has to be km
//		ev.magnitude   = org->;
		ev.horiz_err   = 0;
		ev.depth_err   = 0;
		ev.tt_residual = 0;
		ev.originId    = org->publicID();
		ev.eventId     = query->getEvent(ev.originId)->publicID();
		_events[ev.id] = ev;

		// Add Phases
		for ( size_t i = 0; i < org->arrivalCount(); ++i )
		{
			DataModel::Arrival *orgArr = org->arrival(i);

			const DataModel::Phase& orgPh = orgArr->phase();
			if (orgPh.code() != "P" && orgPh.code() != "S")
			{
				SEISCOMP_ERROR("Phase is neither P nor S, skip it (origin %s)", org->publicID().c_str());
				continue;
			}

			DataModel::PickPtr pick = DataModel::Pick::Cast(query->getObject(DataModel::Pick::TypeInfo(), orgArr->pickID()));
			if ( !pick )
			{
				SEISCOMP_ERROR("Cannot load pick (%s) (origin %s)", orgArr->pickID().c_str(), org->publicID().c_str());
				continue;
			}

			double travel_time = pick->time().value() - org->time().value();
			if (travel_time < 0)
			{
				SEISCOMP_ERROR("Ignoring pick (%s) with negative travel time (origin %s)",
				               orgArr->pickID().c_str(), org->publicID().c_str());
				continue; 
			}

			// find the station
			string netCode = pick->waveformID().networkCode();
			string staCode = pick->waveformID().stationCode();
			string stationId = netCode + "." + staCode;

			// add station if not already there
			if (_stations.find(stationId) == _stations.end())
			{
				DataModel::Station* orgArrStation = findStation(netCode, staCode, pick->time(), query);
				if ( !orgArrStation )
				{
					string msg = "Cannot find station for pick " + orgArr->pickID()
					             + ", origin " + org->publicID();
					throw runtime_error(msg);
				}
				Station sta;
				sta.id = stationId;
				sta.latitude = orgArrStation->latitude();
				sta.longitude = orgArrStation->longitude();
				sta.elevation = orgArrStation->elevation();
				sta.networkCode = netCode;
				sta.stationCode = staCode;
				_stations[sta.id] = sta;
			}

			// the station has to be there at this point
			const Station& sta = _stations.find(stationId)->second;

			Phase ph;
			ph.event_id    = ev.id;
			ph.station_id  = sta.id;
			ph.travel_time = travel_time;
			ph.weight      = orgArr->weight();
			ph.type        = orgPh.code();
			ph.networkCode =  pick->waveformID().networkCode();
			ph.stationCode =  pick->waveformID().stationCode();
			ph.locationCode =  pick->waveformID().locationCode();
			ph.channelCode =  pick->waveformID().channelCode();
			_phases.insert( pair<string,Phase>(ph.event_id, ph));
		}
	}
}


DataModel::Station* BgCatalog::findStation(const string& netCode, const string& stationCode,
                                            Core::Time atTime,
                                            DataModel::DatabaseQuery* query)
{
	DataModel::Inventory *inv = Client::Inventory::Instance()->inventory();

	for (size_t n = 0; n < inv->networkCount(); ++n )
	{
		DataModel::Network *net = inv->network(n);
		if ( net->start() > atTime ) continue;
		try { if ( net->end() < atTime ) continue; }
		catch ( ... ) {}

		if ( !Core::wildcmp(netCode, net->code()) ) continue;

		for ( size_t s = 0; s < net->stationCount(); ++s )
		{
			DataModel::Station *sta = net->station(s);
			if ( sta->start() > atTime ) continue;
			try { if ( sta->end() < atTime ) continue; }
			catch ( ... ) {}

			if ( !Core::wildcmp(stationCode, sta->code()) ) continue;

			return sta;
		}
	}

	return nullptr;
}

BgCatalog::BgCatalog(const map<string,Station>& stations,
                     const map<string,Event>& events,
                     const multimap<string,Phase>& phases)
{
	_stations = stations;
	_events = events;
	_phases = phases;
}

BgCatalog::BgCatalog(const string& stationFile, const string& catalogFile, const string& phaFile)
{
	if ( !Util::fileExists(stationFile) )
	{
		string msg = "File " + stationFile + " does not exist";
		throw runtime_error(msg);
	}

	if ( !Util::fileExists(catalogFile) )
	{
		string msg = "File " + catalogFile + " does not exist";
		throw runtime_error(msg);
	}

	if ( !Util::fileExists(phaFile) )
	{
		string msg = "File " + phaFile + " does not exist";
		throw runtime_error(msg);
	}

	// required headers "id" "latitude" "longitude" "elevation"
	vector<map<string,string> > stations = CSV::readWithHeader(stationFile);

	for (const auto& row : stations )
	{
		Station sta;
		sta.id = row.at("id");
		sta.latitude = std::stod(row.at("latitude"));
		sta.longitude = std::stod(row.at("longitude"));
		sta.elevation = std::stod(row.at("elevation"));
		sta.networkCode = row.at("networkCode");
		sta.stationCode = row.at("stationCode");
		_stations[sta.id] = sta;
}

	// required headers year, month, day, hour, minute, second, latitude, 
	// longitude, depth, magnitude, horizontal_error, depth_error, travel_time_residual, id
	vector<map<string,string> >events = CSV::readWithHeader(catalogFile);

	for (const auto& row : events )
	{
		Event ev;
		ev.id          = row.at("id");
		ev.year        = std::stoi(row.at("year"));
		ev.month       = std::stoi(row.at("month"));
		ev.day         = std::stoi(row.at("day"));
		ev.hour        = std::stoi(row.at("hour"));
		ev.minute      = std::stoi(row.at("minute"));
		ev.second      = std::stod(row.at("second"));
		ev.latitude    = std::stod(row.at("latitude"));
		ev.longitude   = std::stod(row.at("longitude"));
		ev.depth       = std::stod(row.at("depth"));
		ev.magnitude   = std::stod(row.at("magnitude"));
		ev.horiz_err   = std::stod(row.at("horiz_err"));
		ev.depth_err   = std::stod(row.at("depth_err"));
		ev.tt_residual = std::stod(row.at("tt_residual"));
		ev.originId    = row.at("originId");
		ev.eventId     = row.at("eventId");
		_events[ev.id] = ev;
	}

	// required headers event_id station_id travel_time weight type
	vector<map<string,string> >phases = CSV::readWithHeader(phaFile);

	for (const auto& row : phases )
	{
		Phase ph;
		ph.event_id    = row.at("event_id");
		ph.station_id  = row.at("station_id");
		ph.travel_time = std::stod(row.at("travel_time"));
		ph.weight      = std::stod(row.at("weight"));
		ph.type        = row.at("type");
		ph.networkCode   = row.at("networkCode");
		ph.stationCode   = row.at("stationCode");
		ph.locationCode  = row.at("locationCode");
		ph.channelCode   = row.at("channelCode");
		_phases.emplace(ph.event_id, ph);
	}
}


BgCatalogPtr BgCatalog::merge(BgCatalogPtr other)
{
	map<string,Station> stations = getStations();
	map<string,Event> events = getEvents();
	multimap<string,Phase> phases = getPhases();

	for (const auto& kv :  other->getStations() )
	{
		const BgCatalog::Station& station = kv.second;
		if ( stations.count(station.id) == 0 )
			stations[station.id] = station;
	}

	for (const auto& kv :  other->getEvents() )
	{
		const BgCatalog::Event& event = kv.second;
		if ( events.count(event.id) == 0 )
			events[event.id] = event;

		auto eqlrng = other->getPhases().equal_range(event.id);
		auto eqlrng2 = phases.equal_range(event.id);
		for (auto it = eqlrng.first; it != eqlrng.second; ++it)
		{
			const BgCatalog::Phase& otherPhase = it->second;
			bool otherPhaseFound = false;
			for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
			{
				const BgCatalog::Phase& phase = it2->second;
				if (phase.id == otherPhase.id)
				{
					otherPhaseFound = true;
					break;
				}
			}
			if ( ! otherPhaseFound )
				phases.emplace(event.id, otherPhase);
		} 
	}

	return new BgCatalog(stations, events, phases);
}


HypoDD::HypoDD(const BgCatalogPtr& input, const HypoDDConfig& cfg, string workingDir)
{
	_ddbgc = input;
	_cfg = cfg;
	_workingDir = workingDir;

	if ( !Util::fileExists(_controlFile) )
	{
		string msg = "File " + _controlFile + " does not exist";
		throw runtime_error(msg);
	}

	if ( !Util::pathExists(_workingDir) )
	{
		if ( !Util::createPath(_workingDir) )
		{
			string msg = "Unable to create working directory: " + _workingDir;
			throw runtime_error(msg);
		}
	}
}


HypoDD::~HypoDD()
{
	if ( _workingDirCleanup )
	{
		boost::filesystem::remove_all(_workingDir);
	}
}


BgCatalogPtr HypoDD::relocateCatalog()
{
	SEISCOMP_DEBUG("Starting HypoDD relocator in multiple events mode");

	// Create working directory 
	string catalogWorkingDir = (boost::filesystem::path(_workingDir)/"catalog").string(); 
	if ( !Util::pathExists(catalogWorkingDir) )
	{
		if ( !Util::createPath(catalogWorkingDir) )
		{
			string msg = "Unable to create working directory: " + catalogWorkingDir;
			throw runtime_error(msg);
		}
	}

	// Create station.dat for ph2dt and hypodd (if not already generated)
	string stationFile = (boost::filesystem::path(catalogWorkingDir)/"station.dat").string();
	if ( ! Util::fileExists(stationFile) )
	{
		createStationDatFile(stationFile);
	}

	// Create phase.dat for ph2dt (if not already generated)
	string phaseFile = (boost::filesystem::path(catalogWorkingDir)/"phase.dat").string();
	if ( ! Util::fileExists(phaseFile) )
	{
		createPhaseDatFile(phaseFile);
	}

	// run ph2dt
	// input files: ph2dt.inp station.dat phase.dat
	// output files: station.sel event.sel event.dat dt.ct
	string dtctFile = (boost::filesystem::path(catalogWorkingDir)/"dt.ct").string();
	string stationSelFile = (boost::filesystem::path(catalogWorkingDir)/"event.sel").string();
	string eventSelfile = (boost::filesystem::path(catalogWorkingDir)/"station.sel").string();
	if ( !Util::fileExists(dtctFile) &&
	     !Util::fileExists(stationSelFile) &&
	     !Util::fileExists(eventSelfile) )
	{
		runPh2dt(catalogWorkingDir, stationFile, phaseFile);
	}

	// Reads the event pairs matched in dt.ct which are selected by ph2dt and
	// calculate cross correlated differential travel_times for every pair.
	// input dt.ct
	// output dt.cc
	string dtccFile = (boost::filesystem::path(catalogWorkingDir)/"dt.cc").string();
	if ( ! Util::fileExists(dtccFile) )
	{ 
		xcorrCatalog(dtctFile, dtccFile);
	}

	// run hypodd
	// input : dt.cc dt.ct event.sel station.sel hypoDD.inp
	// output : hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
	string ddrelocFile = (boost::filesystem::path(catalogWorkingDir)/"hypoDD.reloc").string();
	if ( ! Util::fileExists(ddrelocFile) )
	{
		runHypodd(catalogWorkingDir, dtccFile, dtctFile, eventSelfile, stationSelFile);
	}

	// load a catalog from hypodd output file
	// input: hypoDD.reloc
	return loadRelocatedCatalog(ddrelocFile);
}


BgCatalogPtr HypoDD::relocateSingleEvent(DataModel::Origin *org)
{
	SEISCOMP_DEBUG("Starting HypoDD relocator in single event mode");

	// Create working directory
	string subFolder = generateWorkingSubDir(org);
	if ( Util::pathExists(subFolder) )
	{
		boost::filesystem::remove_all(subFolder);
	}

	string eventWorkingDir = (boost::filesystem::path(_workingDir)/subFolder/"step1").string();
	if ( !Util::createPath(eventWorkingDir) )
	{
		string msg = "Unable to create working directory: " + eventWorkingDir;
		throw runtime_error(msg);
	}

	BgCatalogPtr orgCatalog = new BgCatalog({org->publicID()}, DataModel::DatabaseQuery* query); // TODO pass query to the function?

	// Select neighbouring events
	BgCatalogPtr subCatalog = selectNeighbouringEvents(_ddbgc, orgCatalog, _cfg.ddse.maxIEdis_CT,
	                                                   _cfg.ddse.minNumNeigh_CT, _cfg.ddse.maxNumNeigh_CT,
	                                                   _cfg.ddse.minDTperEvt_CT);
	BgCatalogPtr subCatalogFull = subCatalog->merge(orgCatalog);

	// Create station.dat for hypodd
	string stationFile = (boost::filesystem::path(eventWorkingDir)/"station.dat").string();
	createStationDatFile(stationFile, subCatalogFull);

	// Create event.dat for hypodd
	string eventFile = (boost::filesystem::path(eventWorkingDir)/"event.dat").string();
	createEventDatFile(eventFile, subCatalogFull);

	// Create differential travel times file (dt.ct) for hypodd
	string dtctFile = (boost::filesystem::path(eventWorkingDir)/"dt.ct").string();
	createDtCtFile(subCatalog, orgCatalog, dtctFile);

	// run hypodd
	// input : dt.cc dt.ct event.sel station.sel hypoDD.inp
	// output : hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
	runHypodd(eventWorkingDir, "", dtctFile, eventFile, stationFile);

	// Load the relocated origin from Hypodd
	string ddrelocFile = (boost::filesystem::path(eventWorkingDir)/"hypoDD.reloc").string();
	BgCatalogPtr relocatedCatalog = loadRelocatedCatalog(ddrelocFile);
	BgCatalogPtr relocatedOrigin = extractEvent(relocatedCatalog,  );  // TODO, event id

	if ( _workingDirCleanup )
	{
		boost::filesystem::remove_all(eventWorkingDir);
	}

	eventWorkingDir = (boost::filesystem::path(_workingDir)/subFolder/"step2").string();
	if ( !Util::createPath(eventWorkingDir) )
	{
		string msg = "Unable to create working directory: " + eventWorkingDir;
		throw runtime_error(msg);
	}

	// Select neighbouring events from the relocated origin
	subCatalog = selectNeighbouringEvents(_ddbgc, relocatedOrigin, _cfg.ddse.maxIEdis_CC,
	                                      _cfg.ddse.minNumNeigh_CC, _cfg.ddse.maxNumNeigh_CC,
	                                      _cfg.ddse.minDTperEvt_CC);
	subCatalogFull = subCatalog->merge(relocatedOrigin);

	// Create station.dat for hypodd
	stationFile = (boost::filesystem::path(eventWorkingDir)/"station.dat").string();
	createStationDatFile(stationFile, subCatalogFull);

	// Create event.dat for hypodd
	eventFile = (boost::filesystem::path(eventWorkingDir)/"event.dat").string();
	createEventDatFile(eventFile, subCatalogFull);

	// Create differential travel times file (dt.ct) for hypodd
	dtctFile = (boost::filesystem::path(eventWorkingDir)/"dt.ct").string();
	createDtCtFile(subCatalog, relocatedOrigin, dtctFile);

	// Create cross correlated differential travel times file (dt.cc) for hypodd
	string dtccFile = (boost::filesystem::path(eventWorkingDir)/"dt.cc").string();
	xcorrSingleEvent(subCatalog, relocatedOrigin, dtccFile);

	// run hypodd
	// input : dt.cc dt.ct event.sel station.sel hypoDD.inp
	// output : hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
	runHypodd(eventWorkingDir, dtccFile, dtctFile, eventFile, stationFile);

	// Load the relocated origin from Hypodd
	ddrelocFile = (boost::filesystem::path(eventWorkingDir)/"hypoDD.reloc").string();
	relocatedCatalog = loadRelocatedCatalog(ddrelocFile);
	relocatedOrigin = extractEvent(relocatedCatalog, ); // TODO, event id

	if ( _workingDirCleanup )
	{
		boost::filesystem::remove_all(eventWorkingDir);
	}

	return relocatedOrigin;
}


/*
 *  Write the station.dat input file for ph2dt and hypodd
 *  One station per line:
 *  STA, LAT, LON, ELV, MODID
 *
 *  E.g.
 NCAAS 38.4301 -121.11   12
 NCABA 38.8793 -121.067  25
 NCABJ 39.1658 -121.193  35
 NCABR 39.1381 -121.48   14
 *
 */
void HypoDD::createStationDatFile(string staFileName, BgCatalogPtr catalog)
{
	SEISCOMP_DEBUG("Creating station file %s", staFileName.c_str());

	if ( !catalog )
	{
		catalog = _ddbgc;
	}

	ofstream mystream;
	mystream.open(staFileName);
	if ( !mystream.is_open() ) {
		string msg = "Cannot create file " + staFileName;
		throw runtime_error(msg);
	}
	
	for (const auto& kv :  catalog->getStations() )
	{
		const BgCatalog::Station& station = kv.second;
		mystream << station.id
		         << " " << setprecision(6) << station.latitude
		         << " " << setprecision(6) << station.longitude
		         << " " << station.elevation
		         << endl;
	}
}


/* Write the phase.dat input file for ph2dt
 * ph2dt accepts hypocenter, followed by its travel time data in the following format:
 * #, YR, MO, DY, HR, MN, SC, LAT, LON, DEP, MAG, EH, EZ, RMS, ID
 * followed by nobs lines of observations:
 * STA, TT, WGHT, PHA
 * e.g.
 *
 #  1985  1 24  2 19 58.71  37.8832 -122.2415    9.80 1.40 0.2 0.5 0.0    38542
 NCCSP       2.850  -1.000   P
 NCCSP       2.910   0.016   P
 NCCBW       3.430  -1.000   P
 NCCBW       3.480   0.031   P
 #  1996 11  9  7  8 36.70  37.8810 -122.2457    9.14 1.10 0.6 0.6 0.0   484120
 NCCVP       1.860  -1.000   P
 NCCSP       2.770   0.250   P
 NCCMC       2.810   0.125   P
 NCCBW       3.360   0.250   P
 *
 */
void HypoDD::createPhaseDatFile(string phaseFileName, BgCatalogPtr catalog)
{
	SEISCOMP_DEBUG("Creating phase file %s", phaseFileName.c_str());

	if ( !catalog )
	{
		catalog = _ddbgc;
	}

	ofstream mystream;
	mystream.open(phaseFileName);
	if ( !mystream.is_open() ) {
		string msg = "Cannot create file " + phaseFileName;
		throw runtime_error(msg);
	}

	for (const auto& kv :  catalog->getEvents() )
	{
		const BgCatalog::Event& event = kv.second;
		mystream << "# " << event.year << " " << event.month
		         << " " << event.day   << " " << event.hour
		         << " " << event.minute
		         << " " << setprecision(6) << event.second
		         << " " << setprecision(6) << event.latitude
		         << " " << setprecision(6) << event.longitude
		         << " " << setprecision(4) << event.depth  // must be km.
		         << " " << setprecision(6) << event.magnitude
		         << " " << setprecision(6) << event.horiz_err
		         << " " << setprecision(6) << event.depth_err
		         << " " << setprecision(6) << event.tt_residual
		         << " " << event.id
		         << endl;

		auto eqlrng = catalog->getPhases().equal_range(event.id);
		for (auto it = eqlrng.first; it != eqlrng.second; ++it)
		{
			const BgCatalog::Phase& phase = it->second;
			mystream << phase.station_id
			          << " " << setprecision(6) << phase.travel_time
			          << " " << setprecision(2) << phase.weight
			          << " " << phase.type
			          << endl;
		}
	}
}

/* Write the event.dat input file for hypodd
 * One event per line:
 * DATE, TIME, LAT, LON, DEP, MAG, EH, EV, RMS, ID
 * e.g.
 *
19850124   2195871   37.8832  -122.2415      9.800   1.4    0.15    0.51   0.02      38542
19911126  14274555   37.8738  -122.2432      9.950   1.4    0.22    0.53   0.09     238298
19861019  20503808   37.8802  -122.2405      9.370   1.6    0.16    0.50   0.05      86036
19850814  18015544   37.8828  -122.2497      7.940   2.0    0.13    0.45   0.06      52942
19850527    430907   37.8778  -122.2412      9.050   1.1    0.26    0.75   0.03      48565
19850402   5571645   37.8825  -122.2420      9.440   1.9    0.12    0.30   0.04      45165
 *
 */
void HypoDD::createEventDatFile(string eventFileName, BgCatalogPtr catalog)
{
	SEISCOMP_DEBUG("Creating station file %s", eventFileName.c_str());

	if ( !catalog )
	{
		catalog = _ddbgc;
	}

	ofstream mystream;
	mystream.open(eventFileName);
	if ( !mystream.is_open() ) {
		string msg = "Cannot create file " + eventFileName;
		throw runtime_error(msg);
	}

	for (const auto& kv :  catalog->getEvents() )
	{
		const BgCatalog::Event& event = kv.second;
		mystream << event.year
		         << setfill('0') << setw(2) << event.month
		         << setfill('0') << setw(2) << event.day
		         << " "
		         << setfill('0') << setw(2) << event.hour
		         << setfill('0') << setw(2) << event.minute
		         << setfill('0') << setw(4) << event.second
		         << " " << setprecision(6) << event.latitude
		         << " " << setprecision(6) << event.longitude
		         << " " << setprecision(4) << event.depth  // must be km.
		         << " " << setprecision(6) << event.magnitude
		         << " " << setprecision(6) << event.horiz_err
		         << " " << setprecision(6) << event.depth_err
		         << " " << setprecision(6) << event.tt_residual
		         << " " << event.id
		         << endl;
	}
}

/*
 * run ph2dt
 * input files: ph2dt.inp station.dat phase.dat
 * output files: station.sel event.sel event.dat dt.ct
 */
void HypoDD::runPh2dt(string workingDir, string stationFile, string phaseFile)
{
	SEISCOMP_DEBUG("Running ph2dt...");

	if ( !Util::fileExists(stationFile) )
		throw runtime_error("Unable to run ph2dt, file doesn't exist: " + stationFile);

	if ( !Util::fileExists(phaseFile) )
		throw runtime_error("Unable to run ph2dt, file doesn't exist: " + phaseFile);

	// write ph2dt.inp - input control file for program ph2dt
	ofstream ph2dtinp;
	ph2dtinp.open((boost::filesystem::path(workingDir)/"ph2dt.inp").string());
	if ( !ph2dtinp.is_open() )
		throw runtime_error("Cannot create file ph2dt.inp");
	
	ph2dtinp << stationFile << endl;
	ph2dtinp << phaseFile << endl;
	ph2dtinp << _cfg.ph2dt.minwght << _cfg.ph2dt.maxdist << _cfg.ph2dt.maxsep
	         << _cfg.ph2dt.maxngh  << _cfg.ph2dt.minlnk  << _cfg.ph2dt.minobs
	         << _cfg.ph2dt.maxobs  << endl;
	ph2dtinp.close(); // this flush the data to disk

	// Run ph2dt

	// TODO
}

/*
 * run hypodd executable
 * input files: dt.cc dt.ct event.sel station.sel hypoDD.inp
 * output files: hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
 */
void HypoDD::runHypodd(string workingDir, string dtccFile, string dtctFile,
                       string eventFile, string stationFile)
{
	SEISCOMP_DEBUG("Running hypodd...");

	string fname = (boost::filesystem::path(workingDir)/dtccFile).string();
	if ( !Util::fileExists(fname) )
		throw runtime_error("Unable to run hypodd, file doesn't exist: " + fname);

	fname = (boost::filesystem::path(workingDir)/dtctFile).string();
	if ( !Util::fileExists(fname) )
		throw runtime_error("Unable to run hypodd, file doesn't exist: " + fname);

	fname = (boost::filesystem::path(workingDir)/eventFile).string();
	if ( !Util::fileExists(fname) )
		throw runtime_error("Unable to run hypodd, file doesn't exist: " + fname);

	fname = (boost::filesystem::path(workingDir)/stationFile).string();
	if ( !Util::fileExists(fname) )
		throw runtime_error("Unable to run hypodd, file doesn't exist: " + fname);

	fname = (boost::filesystem::path(workingDir)/_cfg.hypodd.ctrlFile).string();
	if ( !Util::fileExists(fname) )
		throw runtime_error("Unable to run hypodd, file doesn't exist: " + fname);

	// TODO
}


/*
 * Compute distance in km between two points
 */
double HypoDD::computeDistance(double lat1, double lon1, double depth1,
                               double lat2, double lon2, double depth2)
{
	double distance, az, baz;
	Math::Geo::delazi(lat1, lon1, lat2, lon2, &distance, &az, &baz);
	double Hdist = Math::Geo::deg2km(distance);
	double Vdist = abs(depth1 - depth2);
	return std::sqrt( std::pow(Hdist,2), std::pow(Vdist,2) );
}


BgCatalogPtr HypoDD::selectNeighbouringEvents(BgCatalogPtr catalog,
                                              BgCatalogPtr org,
                                              double maxIEdis,
                                              int minNumNeigh,
                                              int maxNumNeigh,
                                              int minDTperEvt)
{
	const BgCatalog::Event& refEv = org->getEvents().begin()->second;

	map<double,string> eventDistances; // distance, eventid

	// loop through every event in the catalog and select the ones within maxIEdis distance
	for (const auto& kv : catalog->getEvents() )
	{
		const BgCatalog::Event& event = kv.second;

		// compute distance between current event and reference origin
		double distance = computeDistance(event.latitude, event.longitude, event.depth,
	                                      refEv.latitude, refEv.longitude, refEv.depth);
		// too far away ?
		if ( distance > maxIEdis )
			continue;

		// not enought phases ?
		auto eqlrng = catalog->getPhases().equal_range(event.id);
		if (minDTperEvt > 0 && std::distance(eqlrng.first, eqlrng.second) < minDTperEvt)
		{
			SEISCOMP_DEBUG("Skipping event %s, not enough phases (minDTperEvt)", event.id);
			continue;
		}

		// keep a list of added events sorted by distance
		eventDistances[distance] = event.id;
	}

	// Limit num neighbors to maxNumNeigh, choose closest ones
	int maxEvents = maxNumNeigh > 0 ? maxNumNeigh : eventDistances.size();
	map<double,string>  selectedEvents(eventDistances.begin(),
	                                   std::next(eventDistances.begin(), maxEvents));

	// Check if enough neighbors were found
	if (minNumNeigh > 0 && selectedEvents.size() < minNumNeigh)
	{
		string msg = "Insufficient number of neighbors in catalog, skip relocation of origin";
		throw runtime_error(msg);
	}

	map<string,BgCatalog::Station> stations;
	map<string,BgCatalog::Event> events;
	multimap<string,BgCatalog::Phase> phases;

	// Add all the selected events within distance
	for (const auto& kv : selectedEvents)
	{
		const BgCatalog::Event& event = catalog->getEvents().at(kv.second);

		// add this event
		events[event.id] = event;

		// add corresponding phases and stations too
		auto eqlrng = catalog->getPhases().equal_range(event.id);
		for (auto it = eqlrng.first; it != eqlrng.second; ++it)
		{
			const BgCatalog::Phase& phase = it->second;
			phases.emplace(event.id, phase);

			auto search = catalog->getStations().find(phase.station_id);
			if (search == catalog->getStations().end())
			{
				string msg = "Malformed catalog: cannot find station " + phase.station_id +
				              " referenced by phase " + phase.id + " for event " + event.id;
				throw runtime_error(msg);
			}

			const BgCatalog::Station& station = search->second;
			stations[station.id] = station;
		}
	}

	return new BgCatalog(stations, events, phases);
}


// load a catalog from hypodd output file
// input: hypoDD.reloc
BgCatalogPtr HypoDD::loadRelocatedCatalog(string ddrelocFile)
{
	SEISCOMP_DEBUG("Loading catalog relocated by hypodd...");

	if ( !Util::fileExists(ddrelocFile) )
		throw runtime_error("Cannot load hypodd relocated catalog file: " + ddrelocFile);

	// TODO
}


BgCatalogPtr HypoDD::extractEvent(BgCatalogPtr catalog, string eventId)
{
	map<string,BgCatalog::Station> stations;
	map<string,BgCatalog::Event> events;
	multimap<string,BgCatalog::Phase> phases;

	for (const auto& kv : catalog->getEvents() )
	{
		const BgCatalog::Event& event = kv.second;

		if ( event.id == eventId )
		{
			events[event.id] = event;
			auto eqlrng = catalog->getPhases().equal_range(event.id);
			for (auto it = eqlrng.first; it != eqlrng.second; ++it)
			{
				const BgCatalog::Phase& phase = it->second;
				phases.emplace(event.id, phase);

				auto search = catalog->getStations().find(phase.station_id);
				if (search == catalog->getStations().end())
					throw runtime_error("Malformed catalog: cannot find station " + phase.station_id
					                    + " referenced by phase " + phase.id + " for event " + event.id);

				const BgCatalog::Station& station = search->second;
				stations[station.id] = station;
			}
			break;
		}
	}

	return new BgCatalog(stations, events, phases);
}


/* 
 * Create differential travel times file (dt.ct) for hypodd
 * This is for single event mode
 *
 * Each event pair is listed by a header line (in free format)
 * #, ID1, ID2
 * followed by nobs lines of observations (in free format):
 * STA, TT1, TT2, WGHT, PHA
 * 
 */
void HypoDD::createDtCtFile(BgCatalogPtr catalog, BgCatalogPtr org, string dtctFile)
{
	SEISCOMP_DEBUG("Creating Catalog travel time file %s", dtctFile.c_str());

	const BgCatalog::Event& refEv = org->getEvents().begin()->second;

	ofstream mystream;
	mystream.open(dtctFile);
	if ( !mystream.is_open() )
		throw runtime_error("Cannot create file " + dtctFile);

	// loop through catalog events
	for (const auto& kv : catalog->getEvents() )
	{
		const BgCatalog::Event& event = kv.second;

		// loop through event phases
		auto eqlrng = catalog->getPhases().equal_range(event.id);
		for (auto it = eqlrng.first; it != eqlrng.second; ++it)
		{
			const BgCatalog::Phase& phase = it->second;

			// loop through reference event (org to relocate) phases
			auto eqlrng2 = org->getPhases().equal_range(refEv.id);
			for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
			{
				const BgCatalog::Phase& phase2 = it2->second;

				if (phase.station_id == phase2.station_id && 
				    phase.type == phase2.type)
				{
					mystream <<  STA, TT1, TT2, WGHT, PHA  
				}
			}
		}
	}
	// TODO
}


/*
 * Reads the event pairs matched in dt.ct which are selected by ph2dt and
 * calculate cross correlated differential travel_times for every pair.
 * input dt.ct
 * output dt.cc
 */
void HypoDD::xcorrCatalog(string dtctFile, string dtccFile)
{
	SEISCOMP_DEBUG("Calculating cross correlated differential travel times...");

	if ( !Util::fileExists(dtctFile) )
		throw runtime_error("Unable to perform cross correlation, cannot find file: " + dtctFile);

	ofstream mystream;
	mystream.open(dtccFile);
	if ( !mystream.is_open() )
		throw runtime_error("Cannot create file " + dtccFile);

	// TODO
}


/*
 * Calculate cross correlated differential travel_times for every event in catalog and org
 * output dt.cc
 */
void HypoDD::xcorrSingleEvent(BgCatalogPtr catalog, BgCatalogPtr org, string dtccFile)
{
	SEISCOMP_DEBUG("Creating Cross correlation differential time file %s", dtccFile.c_str());

	ofstream mystream;
	mystream.open(dtccFile);
	if ( !mystream.is_open() )
		throw runtime_error("Cannot create file " + dtccFile);

	// TODO
}


GenericRecordPtr
HypoDD::loadWaveform(const Core::Time& starttime,
                     double duration,
                     const string& networkCode,
                     const string& stationCode,
                     const string& locationCode,
                     const string& channelCode)
{
	IO::RecordStreamPtr rs = IO::RecordStream::Open( _cfg.xcorr.recordStreamURL.c_str() );
	if ( rs == NULL )
	{
		string msg = "Cannot open RecordStream: " + _cfg.xcorr.recordStreamURL;
		throw runtime_error(msg);
	}

	Core::TimeWindow tw(starttime, duration);
	rs->setTimeWindow(tw);
	rs->addStream(networkCode, stationCode, locationCode, channelCode);

	// Store each record in a RecordSequence
	IO::RecordInput inp(rs.get(), Array::DOUBLE, Record::DATA_ONLY);
	std::shared_ptr<RecordSequence> seq( new TimeWindowBuffer(tw) );
	RecordPtr rec;
	while ( rec = inp.next() )
	{
		seq->feed(rec.get());
	}
	
	if ( seq->empty() )
	{
		ostringstream msg;
		msg << "No data for stream " << networkCode << "." << stationCode << "." 
		    << locationCode << "." << channelCode << " at " << starttime.toString("%FT%T.%fZ");
		throw runtime_error(msg.str());
	}

    GenericRecordPtr trace = new GenericRecord();
    
	if ( !merge(*trace, *seq) )
    {
		ostringstream msg;
		msg << "Data records could not be merged into a single trace:" << networkCode
		      << "." << stationCode << "." << locationCode << "." << channelCode
		      << " at " << starttime.toString("%FT%T.%fZ");
		throw runtime_error(msg.str());
	}

	if ( !trim(*trace, tw) )
    {
		ostringstream msg;
		msg << "Incomplete trace, not enough data for time window: " << networkCode
		      << "." << stationCode << "." << locationCode << "." << channelCode
		      << " at " << starttime.toString("%FT%T.%fZ");
		throw runtime_error(msg.str());
	}

	filter(*trace, true, _cfg.xcorr.filterOrder, _cfg.xcorr.filterFmin,
	       _cfg.xcorr.filterFmax, _cfg.xcorr.filterFsamp);

	return trace;
}

bool HypoDD::merge(GenericRecord &trace, const RecordSequence& seq)
{
	if ( seq.empty() )
	{
		return false;
	}

	RecordCPtr first = seq.front();
	RecordCPtr last;
	double samplingFrequency = first->samplingFrequency();
	Core::TimeSpan maxAllowedGap, maxAllowedOverlap;

	maxAllowedGap = Core::TimeSpan((double)(0.5 / samplingFrequency));
	maxAllowedOverlap = Core::TimeSpan((double)(-0.5 / samplingFrequency));

	trace.setNetworkCode(first->networkCode());
	trace.setStationCode(first->stationCode());
	trace.setLocationCode(first->locationCode());
	trace.setChannelCode(first->channelCode());

	trace.setStartTime(first->startTime());
	trace.setSamplingFrequency(samplingFrequency);

	Array::DataType datatype = first->data()->dataType();
	ArrayPtr arr = ArrayFactory::Create(datatype, datatype, 0, NULL);

	for (const RecordCPtr& rec : seq )
	{
		if ( rec->samplingFrequency() != samplingFrequency ) {
			SEISCOMP_WARNING("%s.%s.%s.%s: record sampling frequencies are not consistent: %f != %f",
			                 trace.networkCode().c_str(),
			                 trace.stationCode().c_str(),
			                 trace.locationCode().c_str(),
			                 trace.channelCode().c_str(),
			                 samplingFrequency, rec->samplingFrequency());
			return false;
		}

		// Check for gaps and overlaps
		if ( last ) {
			Core::TimeSpan diff = rec->startTime()-last->endTime();
			if ( diff > maxAllowedGap ) {
				SEISCOMP_WARNING("%s.%s.%s.%s: gap detected of %d.%06ds",
				                 trace.networkCode().c_str(),
				                 trace.stationCode().c_str(),
				                 trace.locationCode().c_str(),
				                 trace.channelCode().c_str(),
				                 (int)diff.seconds(), (int)diff.microseconds());
				return false;
			}

			if ( diff < maxAllowedOverlap ) {
				SEISCOMP_WARNING("%s.%s.%s.%s: overlap detected of %fs",
				                 trace.networkCode().c_str(),
				                 trace.stationCode().c_str(),
				                 trace.locationCode().c_str(),
				                 trace.channelCode().c_str(),
				                 (double)diff);
				return false;
			}
		}

		arr->append( (Array *)(rec->data()));

		last = rec;
	}

	trace.setData(arr.get());

	return true;
}


bool HypoDD::trim(GenericRecord &trace, const Core::TimeWindow& tw)
{
	int ofs = (int)(double(tw.startTime() - trace.startTime())*trace.samplingFrequency());
	int samples = (int)(tw.length()*trace.samplingFrequency());

	// Not enough data at start of time window
	if ( ofs < 0 )
	{
		SEISCOMP_DEBUG("%s: need %d more samples in past",
		               trace.streamID().c_str(), -ofs);
		return false;
	}

	// Not enough data at end of time window
	if ( ofs+samples > trace.data()->size() )
	{
		SEISCOMP_DEBUG("%s: need %d more samples past the end",
		               trace.streamID().c_str(), trace.data()->size()-samples-ofs);
		return false;
	}

	ArrayPtr sliced = trace.data()->slice(ofs, ofs+samples);

	trace.setStartTime(tw.startTime());
	trace.setData(sliced.get());

	return true;
}


void HypoDD::filter(GenericRecord &trace, bool demeaning,
                    int order, double fmin, double fmax, double fsamp)
{
	DoubleArray *data = DoubleArray::Cast(trace.data());

	if (demeaning)
	{
		double mean = data->mean();
		int cnt = data->size();
		for ( int i = 0; i < cnt; ++i ) (*data)[i] -= mean;
	}

	if ( fmin > 0 && fmax > 0 )
	{
		Math::Filtering::IIR::ButterworthHighLowpass<double> bp(order, fmin, fmax, fsamp);
		bp.setSamplingFrequency(trace.samplingFrequency());
		bp.apply(data->size(), data->typedData());
	}
	else if ( fmin > 0 )
	{
		Math::Filtering::IIR::ButterworthHighpass<double> hp(order, fmin, fsamp);
		hp.setSamplingFrequency(trace.samplingFrequency());
		hp.apply(data->size(), data->typedData());
	}
	else if ( fmax > 0 )
	{
		Math::Filtering::IIR::ButterworthLowpass<double> lp(order, fmax, fsamp);
		lp.setSamplingFrequency(trace.samplingFrequency());
		lp.apply(data->size(), data->typedData());
	}
}


// Creates dir name from event. This id has the following format:
// OriginTime_Lat_Lon_CreationDate
// eg 20111210115715_46343_007519_20111210115740
string HypoDD::generateWorkingSubDir(const DataModel::Origin* org)
{
	char buf[20];

	if ( org == nullptr ) return "";

	string id;
	id = org->time().value().toString("%Y%m%d%H%M%S");

	id += "_";

	// Latitude
	sprintf(buf, "%05d", int(org->latitude().value()*1000));
	id += buf;

	id += "_";

	// Longitude
	sprintf(buf, "%06d", int(org->longitude().value()*1000));
	id += buf;

	id += "_";

	Core::Time t = Core::Time::GMT();

	id += t.toString("%Y%m%d%H%M%S");

	return id;
}


}
