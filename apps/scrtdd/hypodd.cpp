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

#include <seiscomp3/core/datetime.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/core/typedarray.h>
#include <seiscomp3/io/recordinput.h>
#include <seiscomp3/io/records/mseedrecord.h>
#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/filter/butterworth.h>
#include <seiscomp3/client/inventory.h>
#include <seiscomp3/datamodel/station.h>
#include <seiscomp3/datamodel/event.h>
#include <seiscomp3/datamodel/pick.h>
#include <seiscomp3/datamodel/origin.h>
#include <seiscomp3/datamodel/magnitude.h>
#include <seiscomp3/datamodel/amplitude.h>
#include <seiscomp3/utils/files.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <regex>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <sys/wait.h>
#include <unistd.h>

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp3/logging/log.h>
 
using namespace std;
using namespace Seiscomp;
using Seiscomp::Core::stringify;

namespace {

pid_t startExternalProcess(const vector<string> &cmdparams,
                           bool waitChild,
                           const string& workingDir="")
{
	pid_t pid;
	string cmdline;
	vector<char *> params(cmdparams.size());
	for ( size_t i = 0; i < cmdparams.size(); ++i )
	{
		params[i] = (char*)cmdparams[i].c_str();
		if ( i > 0 ) cmdline += " ";
		cmdline += cmdparams[i];
	}
	params.push_back(nullptr);

	SEISCOMP_DEBUG("Executing command: %s", cmdline.c_str());

	pid = fork();

	if ( pid < 0 ) // fork error
	{
		SEISCOMP_ERROR("Error (%d) in fork()", pid);
		return pid;
	}

	if ( pid == 0 ) // child
	{
		if ( ! workingDir.empty() )
		{
			if( chdir(workingDir.c_str()) != 0 )
			{
				exit(1);
			}
		}

		execv(params[0], &params[0]);
		exit(1);
	}
	else // parent
	{
		if (waitChild) // wait for the child to complete
		{
			int   status;
			do {
				pid = waitpid(pid, &status, 0);
			} while (pid == -1 && errno == EINTR);

			if (status != 0)
				SEISCOMP_ERROR("Command exited with non zero value (%d)", status);
		}
	}

	return pid;
}



template <class Map, class Val> typename Map::const_iterator
searchByValue(const Map & SearchMap, const Val & SearchVal)
{
	typename Map::const_iterator iRet = SearchMap.end();
	for (typename Map::const_iterator iTer = SearchMap.begin(); iTer != SearchMap.end(); iTer ++)
	{
		if (iTer->second == SearchVal)
		{
			iRet = iTer;
			break;
		}
	}
	return iRet;
}



template <class T>
T nextPowerOf2(T a, T min=1, T max=1<<31)
{
	int b = min;
	while (b < a)
	{
		b <<= 1;
		if (b > max)
			return -1;
	}
	return b;
}



void copyFileAndReplaceLines(const string& srcFilename,
                             const string& destFilename,
                             map<int,string> linesToReplace,
                             const string& comment="*")
{
	ifstream srcFile(srcFilename);
	ofstream destFile(destFilename);
	if ( ! srcFile.is_open() || ! destFile.is_open() )
	{
		string msg = stringify("Cannot copy %s to %s", srcFilename.c_str(), destFilename.c_str());
		throw runtime_error(msg);
	}

	string line;
	int lineNum = 0;
	while( std::getline(srcFile, line) )
	{
		// increase line number when not a comment
		if ( line.rfind(comment, 0) != 0 )
			lineNum++;

		// replace line
		if ( linesToReplace.find(lineNum) != linesToReplace.end())
		{
			line = linesToReplace[lineNum];
			linesToReplace.erase(lineNum);
		}

		// copy line to output
		destFile << line << std::endl;
	}
}



DataModel::Station* findStation(const string& netCode,
                                const string& stationCode,
                                Core::Time atTime)
{
	DataModel::Inventory *inv = Client::Inventory::Instance()->inventory();
	if ( ! inv )
	{
		SEISCOMP_ERROR("Inventory not available");
		return nullptr;
	}

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

}


namespace Seiscomp {
namespace HDD {



DataModel::PublicObject* DataSource::getObject(const Core::RTTI& classType,
                                                const std::string& publicID)
{
	DataModel::PublicObject* ret = nullptr;

	if ( _eventParameters && ! ret)
	{
		if (classType == DataModel::Pick::TypeInfo() )
			ret = _eventParameters->findPick(publicID);
		else if (classType == DataModel::Amplitude::TypeInfo() )
			ret = _eventParameters->findAmplitude(publicID);
		else if (classType == DataModel::Origin::TypeInfo() )
			ret = _eventParameters->findOrigin(publicID);
		else if (classType == DataModel::Event::TypeInfo() )
			ret = _eventParameters->findEvent(publicID);
	}

 	if ( _cache && ! ret)
	{
		ret = _cache->find(classType, publicID);
	}

   return ret;
}



void DataSource::loadArrivals(DataModel::Origin* org)
{
	bool found = false;

	if ( _eventParameters && ! found)
	{
		DataModel::Origin* epOrg = _eventParameters->findOrigin(org->publicID());
		if ( epOrg )
		{
			found = true;
			for (size_t i = 0; i < epOrg->arrivalCount(); i++)
				org->add(DataModel::Arrival::Cast(epOrg->arrival(i)->clone()));
		}
	}

	if ( _query && ! found )
	{
		found = true;
		_query->loadArrivals(org);
	}
}



DataModel::Event* DataSource::getParentEvent(const std::string& originID)
{
	DataModel::Event* ret = nullptr;

	if ( _eventParameters && ! ret )
	{
		for (size_t i = 0; i <  _eventParameters->eventCount() && !ret; i++)
		{
			DataModel::Event* ev = _eventParameters->event(i);
			for (size_t j = 0; j < ev->originReferenceCount() && !ret; j++)
			{
				DataModel::OriginReference* orgRef = ev->originReference(j);
				if (orgRef->originID() == originID)
					ret = ev;
			}
		}
	}

	if ( _query && ! ret)
	{
		ret = _query->getEvent(originID);
	}

	return ret;
}



void Catalog::initFromIds(Catalog* catalog,
                          const std::vector<std::string>& ids,
                          DataSource& dataSrc)
{
	vector<DataModel::Origin*> origins;

	for(const string& id : ids)
	{
		DataModel::OriginPtr org = dataSrc.get<DataModel::Origin>(id);
		if ( !org )
		{
			DataModel::EventPtr ev = dataSrc.get<DataModel::Event>(id);
			if ( ev )
			{
				org = dataSrc.get<DataModel::Origin>(ev->preferredOriginID());
			}
		}
		if ( !org )
		{
			string msg = "Cannot find origin/event with id " + id;
			throw runtime_error(msg);
		}
		origins.push_back(org.get());
	}

	initFromOrigins(catalog, origins, dataSrc);
}



void Catalog::initFromOrigins(Catalog* catalog,
                              const std::vector<DataModel::Origin*>& orgs,
                              DataSource& dataSrc)
{
	for(DataModel::Origin* org : orgs)
	{
		if ( org->arrivalCount() == 0)
			dataSrc.loadArrivals(org);

		// Add event
		Event ev;
		ev.id          = 0;
		ev.time        = org->time().value().toGMT();
		ev.latitude    = org->latitude();
		ev.longitude   = org->longitude();
		ev.depth       = org->depth(); // km
		ev.horiz_err   = 0;
		ev.depth_err   = 0;
		ev.tt_residual = 0;
		DataModel::MagnitudePtr mag;
		// try to fetch preferred magnitude stored in the event
		DataModel::EventPtr parentEvent = dataSrc.getParentEvent(org->publicID());
		if ( parentEvent )
		{
			mag = dataSrc.get<DataModel::Magnitude>(parentEvent->preferredMagnitudeID());
		}
		if ( mag )
		{
			ev.magnitude = mag->magnitude();
		}
		else
		{
			SEISCOMP_WARNING("Origin %s: cannot load preferred magnitude from parent event, set it to 0",
			                 org->publicID().c_str());
			ev.magnitude = 0.;
		}

		SEISCOMP_DEBUG("Adding origin '%s' to the catalog", org->publicID().c_str());

		catalog->addEvent(ev, false);
		ev = catalog->searchEvent(ev)->second;

		// Add Phases
		for ( size_t i = 0; i < org->arrivalCount(); ++i )
		{
			DataModel::Arrival *orgArr = org->arrival(i);
			const DataModel::Phase& orgPh = orgArr->phase();

			DataModel::PickPtr pick = dataSrc.get<DataModel::Pick>(orgArr->pickID());
			if ( !pick )
			{
				SEISCOMP_ERROR("Cannot load pick '%s' (origin %s)",
				               orgArr->pickID().c_str(), org->publicID().c_str());
				continue;
			}

			// find the station
			Station sta;
			sta.networkCode = pick->waveformID().networkCode();
			sta.stationCode = pick->waveformID().stationCode();

			// add station if not already there
			if (searchStation(sta) == _stations.end())
			{
				DataModel::Station* orgArrStation = findStation(sta.networkCode, sta.stationCode,
				                                                pick->time());
				if ( !orgArrStation )
				{
					string msg = stringify("Cannot find station for arrival '%s' (origin '%s')",
				                           orgArr->pickID().c_str(), org->publicID().c_str());
					throw runtime_error(msg.c_str());
				}
				sta.latitude = orgArrStation->latitude();
				sta.longitude = orgArrStation->longitude();
				sta.elevation = orgArrStation->elevation(); // meter
				catalog->addStation(sta, false);
			}
			// the station has to be there at this point
			sta = catalog->searchStation(sta)->second;

			Phase ph;
			ph.eventId    = ev.id;
			ph.stationId  = sta.id;
			ph.time       = pick->time().value().toGMT();
			try {
				ph.weight = orgArr->weight();
			} catch ( Core::ValueException& ) {
				ph.weight = 1.;
				SEISCOMP_INFO("Pick '%s' (origin %s) has no weight set, use default weight of 1.",
				              orgArr->pickID().c_str(), org->publicID().c_str());
			}

			ph.type        = orgPh.code();
			ph.networkCode =  pick->waveformID().networkCode();
			ph.stationCode =  pick->waveformID().stationCode();
			ph.locationCode =  pick->waveformID().locationCode();
			ph.channelCode =  pick->waveformID().channelCode();
			catalog->addPhase(ph, false);
		}
	}
}



Catalog::Catalog(const std::vector<DataModel::Origin*>& origins,
                 DataSource& dataSrc)
{
	initFromOrigins(this, origins, dataSrc);
}



Catalog::Catalog(const std::vector<std::string>& ids,
                 DataSource& dataSrc)
{

	initFromIds(this, ids, dataSrc);
}



Catalog::Catalog(const std::string& idFile,
                 DataSource& dataSrc)
{
	if ( !Util::fileExists(idFile) )
	{
		string msg = "File " + idFile + " does not exist";
		throw runtime_error(msg);
	}

	vector<string> ids;
	vector< map<string,string> > rows = CSV::readWithHeader(idFile);

	for(const auto& row : rows)
	{
		const string& id = row.at("seiscompId");
		ids.push_back(id);
	}

	initFromIds(this, ids, dataSrc);
}



Catalog::Catalog() : Catalog(map<string,Station>(),
                             map<unsigned,Event>(),
                             multimap<unsigned,Phase>())
{
}



Catalog::Catalog(const map<string,Station>& stations,
                 const map<unsigned,Event>& events,
                 const multimap<unsigned,Phase>& phases)
{
	_stations = stations;
	_events = events;
	_phases = phases;
}



Catalog::Catalog(const string& stationFile, const string& eventFile, const string& phaFile)
{
	if ( !Util::fileExists(stationFile) )
	{
		string msg = "File " + stationFile + " does not exist";
		throw runtime_error(msg);
	}

	if ( !Util::fileExists(eventFile) )
	{
		string msg = "File " + eventFile + " does not exist";
		throw runtime_error(msg);
	}

	if ( !Util::fileExists(phaFile) )
	{
		string msg = "File " + phaFile + " does not exist";
		throw runtime_error(msg);
	}

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

	vector<map<string,string> >events = CSV::readWithHeader(eventFile);

	for (const auto& row : events )
	{
		Event ev;
		ev.id          = std::stoul(row.at("id"));
		ev.time        = Core::Time::FromString(row.at("isotime").c_str(), "%FT%T.%fZ"); //iso format
		ev.latitude    = std::stod(row.at("latitude"));
		ev.longitude   = std::stod(row.at("longitude"));
		ev.depth       = std::stod(row.at("depth"));
		ev.magnitude   = std::stod(row.at("magnitude"));
		ev.horiz_err   = std::stod(row.at("horiz_err"));
		ev.depth_err   = std::stod(row.at("depth_err"));
		ev.tt_residual = std::stod(row.at("tt_residual"));
		_events[ev.id] = ev;
	}

	vector<map<string,string> >phases = CSV::readWithHeader(phaFile);

	for (const auto& row : phases )
	{
		Phase ph;
		ph.eventId    = std::stoul(row.at("eventId"));
		ph.stationId  = row.at("stationId");
		ph.time        = Core::Time::FromString(row.at("isotime").c_str(), "%FT%T.%fZ"); //iso format
		ph.weight      = std::stod(row.at("weight"));
		ph.type        = row.at("type");
		ph.networkCode   = row.at("networkCode");
		ph.stationCode   = row.at("stationCode");
		ph.locationCode  = row.at("locationCode");
		ph.channelCode   = row.at("channelCode");
		_phases.emplace(ph.eventId, ph);
	}
}




CatalogPtr Catalog::merge(const CatalogPtr& other) const
{
	CatalogPtr mergedCatalog = new Catalog(getStations(), getEvents(), getPhases());

	for (const auto& kv :  other->getEvents() )
	{
		const Catalog::Event& event = kv.second;
		mergedCatalog->addEvent(event, true);
		unsigned newEventId = mergedCatalog->searchEvent(event)->first;

		auto eqlrng = other->getPhases().equal_range(event.id);
		for (auto it = eqlrng.first; it != eqlrng.second; ++it)
		{
			Catalog::Phase phase = it->second;

			auto search = other->getStations().find(phase.stationId);
			if (search == other->getStations().end())
			{
				string msg = stringify("Malformed catalog: cannot find station '%s' "
			                       " referenced by phase '%s'",
			                       phase.stationId.c_str(), string(phase).c_str());
				throw runtime_error(msg);
			}
			const Catalog::Station& station = search->second;
			mergedCatalog->addStation(station, true);
			string newStationId = mergedCatalog->searchStation(station)->first;

			phase.eventId = newEventId;
			phase.stationId = newStationId;

			mergedCatalog->addPhase(phase, true);
		} 
	}

	return mergedCatalog;
}


map<string,Catalog::Station>::const_iterator Catalog::searchStation(const Station& station) const
{
	return searchByValue(_stations, station);
}


map<unsigned,Catalog::Event>::const_iterator Catalog::searchEvent(const Event& event) const
{
	return searchByValue(_events, event);
}


map<unsigned,Catalog::Phase>::const_iterator Catalog::searchPhase(const Phase& phase) const
{
	return searchByValue(_phases, phase);
}


bool Catalog::addStation(const Station& station, bool checkDuplicate)
{
	if (checkDuplicate && searchStation(station) != _stations.end())
	{
		return false;
	}
	Station newStation = station;
	newStation.id = newStation.networkCode + newStation.stationCode;
	_stations[newStation.id] = newStation;
	return true;
}


bool Catalog::addEvent(const Event& event, bool checkDuplicate)
{
	decltype(_events)::key_type maxKey = _events.empty() ? 0 : _events.rbegin()->first;
	if (checkDuplicate && searchEvent(event) != _events.end())
	{
		return false;
	}
	Event newEvent = event;
	newEvent.id = maxKey + 1;
	_events[newEvent.id] = newEvent;
	return true;
}


bool Catalog::addPhase(const Phase& phase, bool checkDuplicate)
{
	if (checkDuplicate && searchPhase(phase) != _phases.end())
	{
		return false;
	}
	_phases.emplace(phase.eventId, phase);
	return true;
}

void Catalog::writeToFile(string eventFile, string phaseFile, string stationFile) const
{
	ofstream evStream(eventFile);
	evStream << "id,isotime,latitude,longitude,depth,magnitude,horiz_err,depth_err,tt_residual" << endl;
	for (const auto& kv : _events )
	{
		const Catalog::Event& ev = kv.second;
		evStream << ev.id << "," << ev.time.iso() << "," << ev.latitude << ","
		         << ev.longitude << "," << ev.depth << "," << ev.magnitude << ","
		         << ev.horiz_err << "," << ev.depth_err << "," << ev.tt_residual << endl;
	}

	ofstream phStream(phaseFile);
	phStream << "eventId,stationId,isotime,weight,type,networkCode,stationCode,locationCode,channelCode" << endl;
	for (const auto& kv : _phases )
	{
		const Catalog::Phase& ph = kv.second;
		phStream << ph.eventId << "," << ph.stationId << "," << ph.time.iso() << ","
		         << ph.weight << "," << ph.type << "," << ph.networkCode << ","
		         << ph.stationCode << "," << ph.locationCode << "," << ph.channelCode << endl;
	}

	ofstream staStream(stationFile);
	staStream << "id,latitude,longitude,elevation,networkCode,stationCode" << endl;
	for (const auto& kv : _stations )
	{
		const Catalog::Station& sta = kv.second;
		staStream << sta.id << "," << sta.latitude << "," << sta.longitude << ","
		          << sta.elevation << "," << sta.networkCode << "," << sta.stationCode << endl;
	}
	
}


HypoDD::HypoDD(const CatalogPtr& input, const Config& cfg, const string& workingDir)
{
	_cfg = cfg;
	_workingDir = workingDir;
	_ddbgc = filterOutPhases(input, cfg.validPphases, cfg.validSphases);

	if ( !Util::pathExists(_workingDir) )
	{
		if ( !Util::createPath(_workingDir) )
		{
			string msg = "Unable to create working directory: " + _workingDir;
			throw runtime_error(msg);
		}
	}

	_cacheDir = (boost::filesystem::path(_workingDir)/"wfcache").string();
	if ( !Util::pathExists(_cacheDir) )
	{
		if ( !Util::createPath(_cacheDir) )
		{
			string msg = "Unable to create cache directory: " + _cacheDir;
			throw runtime_error(msg);
		}
	}

	// write filtered catalog for debugging purpose
	_ddbgc->writeToFile((boost::filesystem::path(_workingDir)/"event.csv").string(),
	                    (boost::filesystem::path(_workingDir)/"phase.csv").string(),
	                    (boost::filesystem::path(_workingDir)/"station.csv").string() );
}


HypoDD::~HypoDD()
{
	if ( _workingDirCleanup )
	{
		// delete all expect the cache directory
		for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(_workingDir), {}))
		{       
			if ( ! boost::filesystem::equivalent(entry, _cacheDir) )
				boost::filesystem::remove_all(entry);
		}
	}
}


// Creates dir name from event. This id has the following format:
// OriginTime_Lat_Lon_CreationDate
// eg 20111210115715_46343_007519_20111210115740
string HypoDD::generateWorkingSubDir(const Catalog::Event& ev) const
{
	char buf[20];

	string id;
	id = ev.time.toString("%Y%m%d%H%M%S");

	id += "_";

	// Latitude
	sprintf(buf, "%05d", int(ev.latitude*1000));
	id += buf;

	id += "_";

	// Longitude
	sprintf(buf, "%06d", int(ev.longitude*1000));
	id += buf;

	id += "_";

	Core::Time t = Core::Time::GMT();

	id += t.toString("%Y%m%d%H%M%S");

	return id;
}


/*
 * Build a catalog with requested phases only and for the same event/station pair
 * make sure to have only one phase. If multiple phases are found, keep the first
 * one arrived
 */
CatalogPtr HypoDD::filterOutPhases(const CatalogPtr& catalog,
                                   const std::vector<std::string>& PphaseToKeep,
                                   const std::vector<std::string>& SphaseToKeep) const
{
	multimap<unsigned,Catalog::Phase> filteredPhases;

	// loop through each event
	for (const auto& kv :  catalog->getEvents() )
	{
		const Catalog::Event& event = kv.second;

		// fetch all phases for current event
		auto eqlrng = catalog->getPhases().equal_range(event.id);
		for (auto it = eqlrng.first; it != eqlrng.second; ++it)
		{
			Catalog::Phase phase = it->second; // copying
			if ( find(PphaseToKeep.begin(), PphaseToKeep.end(), phase.type) != PphaseToKeep.end() )
				phase.type = "P";
			else if ( find(SphaseToKeep.begin(), SphaseToKeep.end(), phase.type) != SphaseToKeep.end() )
				phase.type = "S";
			else {
				SEISCOMP_DEBUG("Discard phase (%s), the type is not among the selected ones",
							   string(phase).c_str());
				continue;
			}

			// fetch already selected phases for current event
			bool inserted = false;
			auto eqlrng2 = filteredPhases.equal_range(event.id);
			for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
			{
				Catalog::Phase& existingPhase = it2->second;
				if ( existingPhase.type == phase.type &&
				     existingPhase.stationId == phase.stationId &&
				     existingPhase.time > phase.time)
				{
					existingPhase = phase;
					inserted = true;
					break;
				}
			}
			if ( ! inserted )
				filteredPhases.emplace(phase.eventId, phase);
		}
	}
	return new Catalog(catalog->getStations(), catalog->getEvents(), filteredPhases);
}


CatalogPtr HypoDD::relocateCatalog(bool force)
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
	if ( force || ! Util::fileExists(stationFile) )
	{
		createStationDatFile(stationFile, _ddbgc);
	}

	// Create phase.dat for ph2dt (if not already generated)
	string phaseFile = (boost::filesystem::path(catalogWorkingDir)/"phase.dat").string();
	if ( force || ! Util::fileExists(phaseFile) )
	{
		createPhaseDatFile(phaseFile, _ddbgc);
	}

	// run ph2dt
	// input files: ph2dt.inp station.dat phase.dat
	// output files: station.sel event.sel event.dat dt.ct
	string dtctFile = (boost::filesystem::path(catalogWorkingDir)/"dt.ct").string();
	string stationSelFile = (boost::filesystem::path(catalogWorkingDir)/"station.sel").string();
	string eventSelfile = (boost::filesystem::path(catalogWorkingDir)/"event.sel").string();
	if ( force ||
	     !Util::fileExists(dtctFile) ||
	     !Util::fileExists(stationSelFile) ||
	     !Util::fileExists(eventSelfile) )
	{
		runPh2dt(catalogWorkingDir, stationFile, phaseFile);
	}

	// Reads the event pairs matched in dt.ct which are selected by ph2dt and
	// calculate cross correlated differential travel_times for every pair.
	// input dt.ct
	// output dt.cc
	string dtccFile = (boost::filesystem::path(catalogWorkingDir)/"dt.cc").string();
	if ( force || ! Util::fileExists(dtccFile) )
	{ 
		xcorrCatalog(dtctFile, dtccFile);
	}

	// run hypodd
	// input : dt.cc dt.ct event.sel station.sel hypoDD.inp
	// output : hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
	string ddrelocFile = (boost::filesystem::path(catalogWorkingDir)/"hypoDD.reloc").string();
	if ( force || ! Util::fileExists(ddrelocFile) )
	{
		runHypodd(catalogWorkingDir, dtccFile, dtctFile, eventSelfile, stationSelFile);
	}

	// load a catalog from hypodd output file
	// input: hypoDD.reloc
	return loadRelocatedCatalog(ddrelocFile, _ddbgc);
}


CatalogPtr HypoDD::relocateSingleEvent(const CatalogPtr& singleEvent)
{
	SEISCOMP_DEBUG("Starting HypoDD relocator in single event mode");

	const CatalogPtr evToRelocateCat = filterOutPhases(singleEvent, _cfg.validPphases, _cfg.validSphases);

	// there must be only one event in the catalog, the origin to relocate
	const Catalog::Event& evToRelocate = evToRelocateCat->getEvents().begin()->second;

	// Create working directory
	string subFolder = generateWorkingSubDir(evToRelocate);
	if ( Util::pathExists(subFolder) )
	{
		boost::filesystem::remove_all(subFolder);
	}

	//
	// Step 1: refine location without cross correlation
	//

	string eventWorkingDir = (boost::filesystem::path(_workingDir)/subFolder/"step1").string();
	if ( !Util::createPath(eventWorkingDir) )
	{
		string msg = "Unable to create working directory: " + eventWorkingDir;
		throw runtime_error(msg);
	}

	// Select neighbouring events
	CatalogPtr neighbourCat = selectNeighbouringEvents(_ddbgc, evToRelocate,
	                                                 _cfg.dtt.maxESdist, _cfg.dtt.maxIEdist,
	                                                 _cfg.dtt.minNumNeigh, _cfg.dtt.maxNumNeigh,
	                                                 _cfg.dtt.minDTperEvt);
	// add event to relocate to the neighbour catalog
	neighbourCat = neighbourCat->merge(evToRelocateCat);
	// extract the new id of the event
	unsigned evToRelocateNewId = neighbourCat->searchEvent(evToRelocate)->first;

	// Create station.dat for hypodd
	string stationFile = (boost::filesystem::path(eventWorkingDir)/"station.dat").string();
	createStationDatFile(stationFile, neighbourCat);

	// Create event.dat for hypodd
	string eventFile = (boost::filesystem::path(eventWorkingDir)/"event.dat").string();
	createEventDatFile(eventFile, neighbourCat);

	// Create differential travel times file (dt.ct) for hypodd
	string dtctFile = (boost::filesystem::path(eventWorkingDir)/"dt.ct").string();
	createDtCtFile(neighbourCat, evToRelocateNewId, dtctFile);
 
	// Create an empty cross correlated differential travel times file (dt.cc) for hypodd
	string dtccFile = (boost::filesystem::path(eventWorkingDir)/"dt.cc").string();
	ofstream(dtccFile).close();

	// run hypodd
	// input : dt.cc dt.ct event.sel station.sel hypoDD.inp
	// output : hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
	runHypodd(eventWorkingDir, dtccFile, dtctFile, eventFile, stationFile);

	// Load the relocated origin from Hypodd
	string ddrelocFile = (boost::filesystem::path(eventWorkingDir)/"hypoDD.reloc").string();
	CatalogPtr relocatedCatalog = loadRelocatedCatalog(ddrelocFile, neighbourCat);
	CatalogPtr relocatedEv = extractEvent(relocatedCatalog, evToRelocateNewId);

	if ( _workingDirCleanup )
	{
		boost::filesystem::remove_all(eventWorkingDir);
	}

	//
	// Step 2: relocate the refined location this time with cross correlation
	//

	CatalogPtr relocatedEvWithXcorr;

	try {

		eventWorkingDir = (boost::filesystem::path(_workingDir)/subFolder/"step2").string();
		if ( !Util::createPath(eventWorkingDir) )
		{
			string msg = "Unable to create working directory: " + eventWorkingDir;
			throw runtime_error(msg);
		}

		// Select neighbouring events from the relocated origin
		const Catalog::Event& refinedLoc = relocatedEv->getEvents().begin()->second;
		neighbourCat = selectNeighbouringEvents(_ddbgc, refinedLoc, _cfg.xcorr.maxESdist,
		                                        _cfg.xcorr.maxIEdist,
		                                        _cfg.xcorr.minNumNeigh, _cfg.xcorr.maxNumNeigh,
		                                        _cfg.xcorr.minDTperEvt);
		// add event to relocate to the neighbour catalog
		neighbourCat = neighbourCat->merge(relocatedEv);
		// extract the new id of the event
		unsigned refinedLocNewId = neighbourCat->searchEvent(refinedLoc)->first;

		// Create station.dat for hypodd
		stationFile = (boost::filesystem::path(eventWorkingDir)/"station.dat").string();
		createStationDatFile(stationFile, neighbourCat);

		// Create event.dat for hypodd
		eventFile = (boost::filesystem::path(eventWorkingDir)/"event.dat").string();
		createEventDatFile(eventFile, neighbourCat);

		// Create differential travel times file (dt.ct) for hypodd
		dtctFile = (boost::filesystem::path(eventWorkingDir)/"dt.ct").string();
		createDtCtFile(neighbourCat, refinedLocNewId, dtctFile);

		// Create cross correlated differential travel times file (dt.cc) for hypodd
		dtccFile = (boost::filesystem::path(eventWorkingDir)/"dt.cc").string();
		xcorrSingleEvent(neighbourCat, refinedLocNewId, dtccFile);

		// run hypodd
		// input : dt.cc dt.ct event.sel station.sel hypoDD.inp
		// output : hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
		runHypodd(eventWorkingDir, dtccFile, dtctFile, eventFile, stationFile);

		// Load the relocated origin from Hypodd
		ddrelocFile = (boost::filesystem::path(eventWorkingDir)/"hypoDD.reloc").string();
		relocatedCatalog = loadRelocatedCatalog(ddrelocFile, neighbourCat);
		relocatedEvWithXcorr = extractEvent(relocatedCatalog, refinedLocNewId);

		if ( _workingDirCleanup )
		{
			boost::filesystem::remove_all(eventWorkingDir);
		}

	} catch ( exception &e ) {
		SEISCOMP_ERROR("%s", e.what());
		SEISCOMP_ERROR("It was not possible to use cross correlation when relocating origin");
	}

	return relocatedEvWithXcorr ? relocatedEvWithXcorr : relocatedEv;
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
void HypoDD::createStationDatFile(const string& staFileName, const CatalogPtr& catalog) const
{
	SEISCOMP_DEBUG("Creating station file %s", staFileName.c_str());

	ofstream outStream(staFileName);
	if ( !outStream.is_open() ) {
		string msg = "Cannot create file " + staFileName;
		throw runtime_error(msg);
	}
	
	for (const auto& kv :  catalog->getStations() )
	{
		const Catalog::Station& station = kv.second;
		outStream << stringify("%-12s %12.6f %12.6f %12.f",
		                      station.id.c_str(), station.latitude,
		                      station.longitude, station.elevation);
		outStream << endl;
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
void HypoDD::createPhaseDatFile(const string& phaseFileName, const CatalogPtr& catalog) const
{
	SEISCOMP_DEBUG("Creating phase file %s", phaseFileName.c_str());

	ofstream outStream(phaseFileName);
	if ( !outStream.is_open() ) {
		string msg = "Cannot create file " + phaseFileName;
		throw runtime_error(msg);
	}

	for (const auto& kv :  catalog->getEvents() )
	{
		const Catalog::Event& event = kv.second;

		int year, month, day, hour, min, sec, usec;
		if ( ! event.time.get(&year, &month, &day, &hour, &min, &sec, &usec) )
		{
			SEISCOMP_WARNING("Cannot convert origin time for event '%s'", string(event).c_str());
			continue;
		}

		outStream << stringify("# %d %d %d %d %d %.2f %.6f %.6f %.2f %.4f %.2f %.2f %.2f %u",
		                      year, month, day, hour, min, sec + double(usec)/1.e6,
		                      event.latitude,event.longitude,event.depth,
		                      event.magnitude, event.horiz_err, event.depth_err,
		                      event.tt_residual, event.id);
		outStream << endl;

		auto eqlrng = catalog->getPhases().equal_range(event.id);
		for (auto it = eqlrng.first; it != eqlrng.second; ++it)
		{
			const Catalog::Phase& phase = it->second;

			double travel_time = phase.time - event.time;
			if (travel_time < 0)
			{
				SEISCOMP_ERROR("Ignoring phase '%s' with negative travel time (event '%s')",
				               string(phase).c_str(), string(event).c_str());
				continue; 
			}

			outStream << stringify("%-12s %12.6f %5.2f %4s",
			                      phase.stationId.c_str(), travel_time,
			                      phase.weight, phase.type.c_str());
			outStream << endl;
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
void HypoDD::createEventDatFile(const string& eventFileName, const CatalogPtr& catalog) const
{
	SEISCOMP_DEBUG("Creating station file %s", eventFileName.c_str());

	ofstream outStream(eventFileName);
	if ( !outStream.is_open() )
	{
		string msg = "Cannot create file " + eventFileName;
		throw runtime_error(msg);
	}

	for (const auto& kv :  catalog->getEvents() )
	{
		const Catalog::Event& event = kv.second;

		int year, month, day, hour, min, sec, usec;
		if ( ! event.time.get(&year, &month, &day, &hour, &min, &sec, &usec) )
		{
			SEISCOMP_WARNING("Cannot convert origin time for event '%s'", string(event).c_str());
			continue;
		}

		outStream << stringify("%d%02d%02d  %02d%02d%04d %.4f %.4f %.3f %.2f %.2f %.2f %.2f %u",
		                      year, month, day, hour, min, int(sec * 1e2 + usec / 1e4),
		                      event.latitude, event.longitude, event.depth,
		                      event.magnitude, event.horiz_err, event.depth_err,
		                      event.tt_residual, event.id);
		outStream << endl;
	}
}

/*
 * run ph2dt
 * input files: ph2dt.inp station.dat phase.dat
 * output files: station.sel event.sel event.dat dt.ct
 */
void HypoDD::runPh2dt(const string& workingDir, const string& stationFile, const string& phaseFile) const
{
	SEISCOMP_DEBUG("Running ph2dt...");

	if ( !Util::fileExists(stationFile) )
		throw runtime_error("Unable to run ph2dt, file doesn't exist: " + stationFile);

	if ( !Util::fileExists(phaseFile) )
		throw runtime_error("Unable to run ph2dt, file doesn't exist: " + phaseFile);

	if ( !Util::fileExists(_cfg.ph2dt.ctrlFile) )
		throw runtime_error("Unable to run ph2dt, control file doesn't exist: " + _cfg.ph2dt.ctrlFile);

	// copy control file while replacing input/output file names
	map<int,string> linesToReplace = { {1,stationFile}, {2,phaseFile} };
	copyFileAndReplaceLines(_cfg.ph2dt.ctrlFile,
                            (boost::filesystem::path(workingDir)/"ph2dt.inp").string(),
                            linesToReplace);

	// run ph2dt (use /bin/sh to get stdout/strerr redirection)
	string cmd = stringify("%s %s >ph2dt.out 2>&1",
	                       _cfg.ph2dt.exec.c_str(), "ph2dt.inp");
	::startExternalProcess({"/bin/sh", "-c", cmd}, true, workingDir);
}

/*
 * run hypodd executable
 * input files: dt.cc dt.ct event.sel station.sel hypoDD.inp
 * output files: hypoDD.loc hypoDD.reloc hypoDD.sta hypoDD.res hypoDD.src
 */
void HypoDD::runHypodd(const string& workingDir, const string& dtccFile, const string& dtctFile,
                       const string& eventFile, const string& stationFile) const
{
	SEISCOMP_DEBUG("Running hypodd...");

	if ( !Util::fileExists(dtccFile) )
		throw runtime_error("Unable to run hypodd, file doesn't exist: " + dtccFile);

	if ( !Util::fileExists(dtctFile) )
		throw runtime_error("Unable to run hypodd, file doesn't exist: " + dtctFile);

	if ( !Util::fileExists(eventFile) )
		throw runtime_error("Unable to run hypodd, file doesn't exist: " + eventFile);

	if ( !Util::fileExists(stationFile) )
		throw runtime_error("Unable to run hypodd, file doesn't exist: " + stationFile);

	if ( !Util::fileExists(_cfg.hypodd.ctrlFile) )
		throw runtime_error("Unable to run hypodd, control file doesn't exist: " + _cfg.hypodd.ctrlFile);

	// copy control file while replacing input/output file names
	map<int,string> linesToReplace = {
		{2,dtccFile},
		{3,dtctFile},
		{4,eventFile},
		{5,stationFile},
		{6,"hypoDD.loc"},
		{7,"hypoDD.reloc"},
		{8,"hypoDD.sta"},
		{9,"hypoDD.res"},
		{10,"hypoDD.src"}
	};
	copyFileAndReplaceLines(_cfg.hypodd.ctrlFile,
                            (boost::filesystem::path(workingDir)/"hypodd.inp").string(),
                            linesToReplace);

	// run Hypodd (use /bin/sh to get stdout/strerr redirection)
	string cmd = stringify("%s %s >hypodd.out 2>&1",
	                       _cfg.hypodd.exec.c_str(), "hypodd.inp");
	::startExternalProcess({"/bin/sh", "-c", cmd}, true, workingDir);
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
	// this is an approximation that works when the distance is small
	// and the Earth curvature can be assumed flat
	return std::sqrt( std::pow(Hdist,2) + std::pow(Vdist,2) );
}


CatalogPtr HypoDD::selectNeighbouringEvents(const CatalogPtr& catalog,
                                            const Catalog::Event& refEv,
                                            double maxESdis,
                                            double maxIEdis,
                                            int minNumNeigh,
                                            int maxNumNeigh,
                                            int minDTperEvt)
{
	map<double,unsigned> eventDistances; // distance, eventid

	// loop through every event in the catalog and select the ones within maxIEdis distance
	for (const auto& kv : catalog->getEvents() )
	{
		const Catalog::Event& event = kv.second;

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
			SEISCOMP_DEBUG("Skipping event '%s', not enough phases (minDTperEvt)",
			               string(event).c_str());
			continue;
		}

		// keep a list of added events sorted by distance
		eventDistances[distance] = event.id;
	}

	// Limit num neighbors to maxNumNeigh, choose closest ones
	int maxEvents = maxNumNeigh > 0 ? maxNumNeigh : eventDistances.size();
	map<double,unsigned>  selectedEvents(eventDistances.begin(),
	                                   std::next(eventDistances.begin(), maxEvents));

	// Check if enough neighbors were found
	if (minNumNeigh > 0 && selectedEvents.size() < minNumNeigh)
	{
		string msg = "Insufficient number of neighbors in catalog, skip relocation of origin";
		throw runtime_error(msg);
	}

	map<string,Catalog::Station> stations;
	map<unsigned,Catalog::Event> events;
	multimap<unsigned,Catalog::Phase> phases;

	// Add all the selected events within distance
	for (const auto& kv : selectedEvents)
	{
		const Catalog::Event& event = catalog->getEvents().at(kv.second);

		// add this event
		events[event.id] = event;

		// add corresponding phases and stations too
		auto eqlrng = catalog->getPhases().equal_range(event.id);
		for (auto it = eqlrng.first; it != eqlrng.second; ++it)
		{
			const Catalog::Phase& phase = it->second;
			auto search = catalog->getStations().find(phase.stationId);
			if (search == catalog->getStations().end())
			{
				string msg = stringify("Malformed catalog: cannot find station '%s' "
			                       " referenced by phase '%s' for event '%u'",
			                       phase.stationId.c_str(), string(phase).c_str(),  event.id);
				throw runtime_error(msg);
			}

			const Catalog::Station& station = search->second;
            
			// compute distance between current event and station
			double distance = computeDistance(event.latitude, event.longitude, event.depth,
			                                  station.latitude, station.longitude, -(station.elevation/1000.));
			// too far away ?
			if ( distance > maxESdis )
				continue;

			phases.emplace(event.id, phase);
			stations[station.id] = station;
		}
	}

	return new Catalog(stations, events, phases);
}


/*
 * load a catalog from hypodd output file
 * input: hypoDD.reloc
 *
 * One event per line (written in fixed, but may be read in free format):
 * ID, LAT, LON, DEPTH, X, Y, Z, EX, EY, EZ, YR, MO, DY, HR, MI, SC, MAG, NCCP, NCCS, NCTP,
NCTS, RCC, RCT, CID
 *
 */
CatalogPtr HypoDD::loadRelocatedCatalog(const string& ddrelocFile, const CatalogPtr& originalCatalog)
{
	SEISCOMP_DEBUG("Loading catalog relocated by hypodd...");

	if ( !Util::fileExists(ddrelocFile) )
		throw runtime_error("Cannot load hypodd relocated catalog file: " + ddrelocFile);

	map<string,Catalog::Station> stations = originalCatalog->getStations();
	map<unsigned,Catalog::Event> events = originalCatalog->getEvents();
	multimap<unsigned,Catalog::Phase> phases = originalCatalog->getPhases();

	// read file one line a time
	ifstream in(ddrelocFile);
	while (!in.eof())
	{
		string row;
		std::getline(in, row);
		if (in.bad() || in.fail())
			break;

		// split line on space
		static const std::regex regex(R"([\s]+)", std::regex::optimize);
		std::sregex_token_iterator it{row.begin(), row.end(), regex, -1};
		std::vector<std::string> fields{it, {}};

		// remove the first empty element if the line start with spaces
		if ( !fields.empty() && fields[0] == "")
			fields.erase(fields.begin());

		if (fields.size() != 24)
		{
			SEISCOMP_WARNING("Skipping unrecognized line from '%s' (line='%s')",
			                 ddrelocFile.c_str(), row.c_str());
			continue;
		}

		// load corresponding event and update information
		unsigned eventId = std::stoul(fields[0]);
		auto search = events.find(eventId);
		if (search == events.end())
		{
			string msg = stringify("Malformed catalog: cannot find relocated event %u"
			                       " in the original catalog", eventId);
			throw runtime_error(msg);
		}
		Catalog::Event& event = search->second;
		event.latitude  = std::stod(fields[1]);
		event.longitude = std::stod(fields[2]);
		event.depth     = std::stod(fields[3]);

		int year  = std::stoi(fields[10]);
		int month = std::stoi(fields[11]);
		int day   = std::stoi(fields[12]);
		int hour  = std::stoi(fields[13]);
		int min   = std::stoi(fields[14]);
		double seconds = std::stod(fields[15]);
		int sec  = int(seconds);
		int usec = (seconds - sec) * 1.e6;
	
		event.time = Core::Time(year, month, day, hour, min, sec, usec);
	}

	return new Catalog(stations, events, phases);
}


CatalogPtr HypoDD::extractEvent(const CatalogPtr& catalog, unsigned eventId) const
{
	CatalogPtr eventToExtract = new Catalog();

	auto search = catalog->getEvents().find(eventId);
	if (search == catalog->getEvents().end())
	{
		string msg = stringify("Cannot find event id %u in the catalog.", eventId);
		throw runtime_error(msg);
	}

	const Catalog::Event& event = search->second;
	eventToExtract->addEvent(event, false);

	unsigned newEventId = eventToExtract->searchEvent(event)->first;

	auto eqlrng = catalog->getPhases().equal_range(event.id);
	for (auto it = eqlrng.first; it != eqlrng.second; ++it)
	{
		Catalog::Phase phase = it->second;

		auto search = catalog->getStations().find(phase.stationId);
		if (search == catalog->getStations().end())
		{
			string msg = stringify("Malformed catalog: cannot find station '%s' "
			                       " referenced by phase '%s' for event '%s'",
			                       phase.stationId.c_str(), string(phase).c_str(), string(event).c_str());
			throw runtime_error(msg);
		}
		const Catalog::Station& station = search->second;
		eventToExtract->addStation(station, true);
		string newStationId = eventToExtract->searchStation(station)->first;

		phase.eventId = newEventId;
		phase.stationId = newStationId;

		eventToExtract->addPhase(phase, true);
	} 

	return eventToExtract;
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
void HypoDD::createDtCtFile(const CatalogPtr& catalog,
                            unsigned evToRelocateId,
                            const string& dtctFile) const
{
	SEISCOMP_DEBUG("Creating Catalog travel time file %s", dtctFile.c_str());

	ofstream outStream(dtctFile);
	if ( !outStream.is_open() )
		throw runtime_error("Cannot create file " + dtctFile);

	auto search = catalog->getEvents().find(evToRelocateId);
	if (search == catalog->getEvents().end())
	{
		string msg = stringify("Cannot find event id %u in the catalog.", evToRelocateId);
		throw runtime_error(msg);
	}

	const Catalog::Event& refEv = search->second;

	// loop through catalog events
	for (const auto& kv : catalog->getEvents() )
	{
		const Catalog::Event& event = kv.second;

		if (event.id == evToRelocateId)
			continue;

		int dtCount = 0;
		stringstream evStream;
		evStream << stringify("# %10u %10u", event.id, refEv.id) << endl;

		// loop through event phases
		auto eqlrng = catalog->getPhases().equal_range(event.id);
		for (auto it = eqlrng.first; it != eqlrng.second; ++it)
		{
			const Catalog::Phase& phase = it->second;
			if (phase.weight < _cfg.dtt.minWeight)
				continue;

			// loop through reference event phases
			auto eqlrng2 = catalog->getPhases().equal_range(refEv.id);
			for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
			{
				const Catalog::Phase& refPhase = it2->second;
				if (refPhase.weight < _cfg.dtt.minWeight)
					continue;

				if (phase.stationId == refPhase.stationId && 
				    phase.type == refPhase.type)
				{
					double ref_travel_time = refPhase.time - refEv.time;
					if (ref_travel_time < 0)
					{
						SEISCOMP_WARNING("Ignoring phase '%s' with negative travel time (event '%s')",
						               string(refPhase).c_str(), string(refEv).c_str());
						continue; 
					}
					double travel_time = phase.time - event.time;
					if (travel_time < 0)
					{
						SEISCOMP_WARNING("Ignoring phase '%s' with negative travel time (event '%s')",
						               string(phase).c_str(), string(event).c_str());
						continue; 
					}

					// get common observation weight for pair (FIXME: take the lower one? average?)
					double weight = (refPhase.weight + phase.weight) / 2.0;

					evStream << stringify("%-12s %.6f %.6f %.2f %s",
					                      refPhase.stationId.c_str(), travel_time,
					                      ref_travel_time, weight, refPhase.type.c_str());
					evStream << endl;
					dtCount++;
				}
			}
		}
		if (dtCount >= _cfg.dtt.minDTperEvt)
			outStream << evStream.str();
	}
}


/*
 * Reads the event pairs matched in dt.ct which are selected by ph2dt and
 * calculate cross correlated differential travel_times for every pair.
 * input dt.ct 
 * output dt.cc
 */
void HypoDD::xcorrCatalog(const string& dtctFile, const string& dtccFile)
{
	SEISCOMP_DEBUG("Calculating cross correlated differential travel times...");

	if ( !Util::fileExists(dtctFile) )
		throw runtime_error("Unable to perform cross correlation, cannot find file: " + dtctFile);

	ofstream outStream(dtccFile);
	if ( !outStream.is_open() )
		throw runtime_error("Cannot create file " + dtccFile);

	const std::map<unsigned,Catalog::Event>& events = _ddbgc->getEvents();
	const Catalog::Event *ev1 = nullptr, *ev2 = nullptr;
	int dtCount = 0;
	stringstream evStream;

	// read file one line a time
	ifstream in(dtctFile);
	while (!in.eof())
	{
		string row;
		std::getline(in, row);
		if (in.bad() || in.fail())
			break;

		// split line on space
		static const std::regex regex(R"([\s]+)", std::regex::optimize);
		std::sregex_token_iterator it{row.begin(), row.end(), regex, -1};
		std::vector<std::string> fields{it, {}};

		// remove the first empty element if the line start with spaces
		if ( !fields.empty() && fields[0] == "")
			fields.erase(fields.begin());

		// check beginning of a new event pair line (# ID1 ID2)
		if (fields[0] == "#" && fields.size() == 3)
		{
			unsigned evId1 = std::stoul(fields[1]);
			unsigned evId2 = std::stoul(fields[2]);
			auto search1 = events.find(evId1);
			auto search2 = events.find(evId2);
			if (search1 == events.end() || search2 == events.end())
			{
				string msg = stringify("Relocated catalog contains events ids (%s or %s) "
				                       "that are not present in the original catalog.",
				                       string(*ev1).c_str(), string(*ev2).c_str());
				throw runtime_error(msg.c_str());
			}
			ev1 = &search1->second;
			ev2 = &search2->second;

			// write the pairs has been built up to now
			if (dtCount >= _cfg.xcorr.minDTperEvt)
			{
				outStream << evStream.str();
			}
			evStream = stringstream();
			dtCount = 0;

			evStream << stringify("# %10u %10u       0.0", ev1->id, ev2->id) << endl;
		}
		// observation line (STA, TT1, TT2, WGHT, PHA)
		else if(ev1 != nullptr && ev2 != nullptr && fields.size() == 5)
		{
			string stationId = fields[0];
			string phaseType = fields[4];

			GenericRecordPtr tr1, tr2;

			// loop through event 1 phases
			auto eqlrng = _ddbgc->getPhases().equal_range(ev1->id);
			for (auto it = eqlrng.first; it != eqlrng.second; ++it)
			{
				const Catalog::Phase& phase = it->second;
				if (phase.weight < _cfg.xcorr.minWeight)
					continue;
				if (phase.stationId != stationId ||
				    phase.type != phaseType)
					continue;
				tr1 = getWaveform(*ev1, phase, _wfCache, _wfDiskCache);
			}

			// loop through event 2 phases
			eqlrng = _ddbgc->getPhases().equal_range(ev2->id);
			for (auto it = eqlrng.first; it != eqlrng.second; ++it)
			{
				const Catalog::Phase& phase = it->second;
				if (phase.weight < _cfg.xcorr.minWeight)
					continue;
				if (phase.stationId != stationId ||
				    phase.type != phaseType)
					continue;
				tr2 = getWaveform(*ev2, phase, _wfCache, _wfDiskCache);
			}

			if ( !tr1 || !tr2)
			{
				SEISCOMP_WARNING("Cannot load waveforms. Skipping line '%s' from file '%s'",
				               row.c_str(), dtctFile.c_str());
				continue;
			}

			SEISCOMP_DEBUG("Calculating cross correlation for event pair %u-%u, line '%s' from file '%s'",
			               ev1->id, ev2->id, row.c_str(), dtctFile.c_str());

			double xcorr_coeff, xcorr_dt;
			if ( ! xcorr(tr1, tr2, _cfg.xcorr.maxDelay, xcorr_dt, xcorr_coeff) )
			{
				SEISCOMP_WARNING("Cannot cross correlate traces for events '%s' and '%s', "
				                 "station %s, phase %s. Skipping them.",
				                 string(*ev1).c_str(), string(*ev2).c_str(),
				                 stationId.c_str(), phaseType.c_str() );
				continue;
			}

			if ( xcorr_coeff < _cfg.xcorr.minCoef)
			{
				SEISCOMP_DEBUG("Cross correlation coefficient too low (%f < %f): skip pair",
				               xcorr_coeff, _cfg.xcorr.minCoef);
				continue;
			}

			evStream << stringify("%-12s %.6f %.4f %s",
			                       stationId.c_str(), xcorr_dt, xcorr_coeff, phaseType.c_str());
			evStream << endl;
			dtCount++;
		}
		else
		{
			ev1 = ev2 = nullptr;
			SEISCOMP_WARNING("Skipping unrecognized line from '%s' (line='%s')",
			               dtctFile.c_str(), row.c_str());
		}
	}

	if (dtCount >= _cfg.xcorr.minDTperEvt)
	{
		outStream << evStream.str();
	}
}


/*
 * Compute and store to file differential travel times from cross correlation for pairs
 * of earthquakes.
 *
 * Each event pair is listed by a header line (in free format)
 * #, ID1, ID2, OTC
 * followed by lines with observations (in free format):
 * STA, DT, WGHT, PHA
 *
 */
void HypoDD::xcorrSingleEvent(const CatalogPtr& catalog,
                              unsigned evToRelocateId,
                              const string& dtccFile)
{
	SEISCOMP_DEBUG("Creating Cross correlation differential time file %s", dtccFile.c_str());

	ofstream outStream(dtccFile);
	if ( !outStream.is_open() )
		throw runtime_error("Cannot create file " + dtccFile);

	auto search = catalog->getEvents().find(evToRelocateId);
	if (search == catalog->getEvents().end())
	{
		string msg = stringify("Cannot find event id %u in the catalog.", evToRelocateId);
		throw runtime_error(msg);
	}
	const Catalog::Event& refEv = search->second;

	map<string, GenericRecordPtr> tmpWfCache;

	// loop through catalog events
	for (const auto& kv : catalog->getEvents() )
	{
		const Catalog::Event& event = kv.second;

		if (event.id == evToRelocateId)
			continue;

		int dtCount = 0;
		stringstream evStream;
		evStream << stringify("# %10u %10u       0.0", event.id, refEv.id) << endl;

		// loop through event phases
		auto eqlrng = catalog->getPhases().equal_range(event.id);
		for (auto it = eqlrng.first; it != eqlrng.second; ++it)
		{
			const Catalog::Phase& phase = it->second;
			if (phase.weight < _cfg.xcorr.minWeight)
				continue;

			GenericRecordPtr trace = getWaveform(event, phase, _wfCache, _wfDiskCache);
			if ( !trace )
			{
				SEISCOMP_WARNING("Cannot load phase waveform, skipping xcorr for event '%s' phase '%s')",
								 string(phase).c_str(), string(event).c_str());
				continue;
			}

			// loop through reference event phases
			auto eqlrng2 = catalog->getPhases().equal_range(refEv.id);
			for (auto it2 = eqlrng2.first; it2 != eqlrng2.second; ++it2)
			{
				const Catalog::Phase& refPhase = it2->second;
				if (refPhase.weight < _cfg.xcorr.minWeight)
					continue;

				if (phase.stationId == refPhase.stationId && 
				    phase.type == refPhase.type)
				{
					GenericRecordPtr refTrace = getWaveform(refEv, refPhase, tmpWfCache, false);
					if ( !refTrace )
					{
						SEISCOMP_WARNING("Cannot load phase waveform, skipping xcorr for event '%s' phase '%s')",
						                 string(refPhase).c_str(), string(refEv).c_str());
						continue;
					}

					SEISCOMP_DEBUG("Calculating cross correlation for event pair %u-%u, station %s phase %s",
					               refEv.id, event.id, refPhase.stationId.c_str(), refPhase.type.c_str());

					double xcorr_coeff, xcorr_dt;
					if ( ! xcorr(trace, refTrace, _cfg.xcorr.maxDelay, xcorr_dt, xcorr_coeff) )
					{
						SEISCOMP_WARNING("Cannot cross correlate traces: ev1 '%s' ph1 '%s' "
						               "with ev2 '%s' ph2 '%s', skipping them.",
						               string(refEv).c_str(), string(refPhase).c_str(),
						               string(event).c_str(), string(phase).c_str() );
						continue;
					}

					if ( xcorr_coeff < _cfg.xcorr.minCoef)
					{
						SEISCOMP_DEBUG("Cross correlation coefficient too low (%f < %f): skip pair",
									   xcorr_coeff, _cfg.xcorr.minCoef);
						continue;
					}

					evStream << stringify("%-12s %.6f %.4f %s",
					                      refPhase.stationId.c_str(),
					                      xcorr_dt, xcorr_coeff*xcorr_coeff,
					                      refPhase.type.c_str());
					evStream << endl;
					dtCount++;
				}
			}
		}
		if (dtCount >= _cfg.xcorr.minDTperEvt)
			outStream << evStream.str();
	}
}


bool
HypoDD::xcorr(const GenericRecordPtr& tr1, const GenericRecordPtr& tr2, double maxDelay,
              double& delayOut, double& coeffOut) const
{
	delayOut = 0.;
	coeffOut = -1.;

	if (tr1->samplingFrequency() != tr2->samplingFrequency())
	{
		SEISCOMP_WARNING("Cannot cross correlate traces with different sampling freq (%f!=%f)",
		                 tr1->samplingFrequency(), tr2->samplingFrequency());
		return false;
	}

	double freq = tr1->samplingFrequency();
	int maxDelaySmps = maxDelay * freq; // secs to samples

	// check longest/shortest trace
	bool swap = tr1->data()->size() > tr2->data()->size();
	GenericRecordPtr trShorter = swap ? tr2 : tr1;
	GenericRecordPtr trLonger = swap ? tr1 : tr2;

	const double *smpsS = DoubleArray::ConstCast(trShorter->data())->typedData();
	const double *smpsL = DoubleArray::ConstCast(trLonger->data())->typedData();
	int smpsSsize = trShorter->data()->size();
	int smpsLsize = trLonger->data()->size();

	// Calculate part of the denominator (speed up)
	double denomS;
	for (int idxS = 0; idxS < smpsSsize; idxS++)
	{
		denomS += smpsS[idxS] * smpsS[idxS];
	}

	// Calculate the correlation series (tr1 and tr2 are already demeaned)
	for (int delay = -maxDelaySmps; delay < maxDelaySmps; delay++)
	{
		double numer = 0, denomL = 0;
		for (int idxS = 0; idxS < smpsSsize; idxS++)
		{
			int idxL = idxS + smpsLsize/2 - smpsSsize/2 + delay;
			if (idxL < 0 || idxL >= smpsLsize)
				continue;
			numer  += smpsS[idxS] * smpsL[idxL];
			denomL += smpsL[idxL] * smpsL[idxL];
		}
		double denom = std::sqrt(denomS*denomL);
		double coeff = numer / denom;
		if (coeff > coeffOut)
		{
			coeffOut = coeff;
			delayOut = delay / freq; // samples to secs
		}
	}

	if(swap)
		delayOut = -delayOut;

	return true;
}

/*
 * Return the waveform from the memory cache if present, otherwise load it
 */
GenericRecordPtr
HypoDD::getWaveform(const Catalog::Event& ev,
                    const Catalog::Phase& ph,
                    map<string,GenericRecordPtr>& cache,
                    bool useDiskCache)
{
	Core::Time starttime = ph.time - Core::TimeSpan(_cfg.xcorr.timeBeforePick);
	double duration = _cfg.xcorr.timeBeforePick + _cfg.xcorr.timeAfterPick;

	string wfId = stringify("%s.%s.%s.%s.%s.%.6f",
	                        ph.networkCode.c_str(), ph.stationCode.c_str(),
	                        ph.locationCode.c_str(), ph.channelCode.c_str(),
	                        ph.time.iso().c_str(), duration);

	const auto it = cache.find(wfId);
	if ( it == cache.end() ) // waveform not cached yet
	{
		GenericRecordPtr wf;
		// load waveform
		try {
			wf = loadWaveform(starttime, duration,
			                  ph.networkCode, ph.stationCode,
			                  ph.locationCode, ph.channelCode,
			                  useDiskCache);
		} catch ( ... ) {}
		// cache waveform
		if (wf) 
			cache[wfId] = wf;
		return wf;
	}
	else // the waveform is already in cache
		return it->second;
}


/*
 * Read a waveform from a chached copy on disk if present, otherwise
 * from the configured RecordStream
 */
GenericRecordPtr
HypoDD::loadWaveform(const Core::Time& time,
                     double duration,
                     const string& networkCode,
                     const string& stationCode,
                     const string& locationCode,
                     const string& channelCode,
                     bool useDiskCache) const
{
	Core::TimeWindow tw(time.toLocalTime(), duration);

	string cacheFile = stringify("%s.%s.%s.%s.%s.%.6f.mseed",
	                             networkCode.c_str(), stationCode.c_str(),
	                             locationCode.c_str(), channelCode.c_str(),
	                             tw.startTime().iso().c_str(), tw.length());
	cacheFile = (boost::filesystem::path(_cacheDir)/cacheFile).string();

	GenericRecordPtr trace;
	// First try to read trace from disk cache
	if ( useDiskCache && Util::fileExists(cacheFile) )
	{
		try {
			std::ifstream ifs(cacheFile);
			IO::MSeedRecord cachedRecord(Array::DOUBLE, Record::Hint::DATA_ONLY);
			cachedRecord.read(ifs);
			trace = new GenericRecord(cachedRecord);
			trace->setData(cachedRecord.data()->clone()); // copy data too
		} catch ( ... ) {
			trace = nullptr;
			SEISCOMP_WARNING("Couldn't load cached waveform %s, read it from record stream",
			                 cacheFile.c_str());
		}
	}

	// if the trace is not cached then read it from the configured recordStream
	if ( !trace )
	{
		trace = readWaveformFromRecordStream(tw,networkCode, stationCode, locationCode, channelCode);
		// then save the trace to disk for later usage
		if ( useDiskCache )
		{
			std::ofstream ofs(cacheFile);
			IO::MSeedRecord cachedRecord(*trace);
			int reclen = cachedRecord.data()->size()*cachedRecord.data()->bytes() + 64;
			reclen = nextPowerOf2<int>(reclen, 128, 1048576); // MINRECLEN 128, MAXRECLEN 1048576
			if (reclen > 0)
			{
				try {
					cachedRecord.setOutputRecordLength(reclen);
					cachedRecord.write(ofs);
				} catch ( ... ) {
					SEISCOMP_WARNING("Couldn't write waveform cache to disk %s", cacheFile.c_str());
				}
			}
		}
	}

	if ( !trim(*trace, tw) )
	{
		string msg = stringify("Incomplete trace, not enough data for requested"
		                       " time window (%s.%s.%s.%s from %s length %.2f)",
		                       networkCode.c_str(), stationCode.c_str(),
		                       locationCode.c_str(), channelCode.c_str(),
		                       tw.startTime().iso().c_str(), tw.length());
		throw runtime_error(msg);
	}

	filter(*trace, true, _cfg.xcorr.filterOrder, _cfg.xcorr.filterFmin,
	       _cfg.xcorr.filterFmax, _cfg.xcorr.filterFsamp);

	return trace;
}


GenericRecordPtr
HypoDD::readWaveformFromRecordStream(const Core::TimeWindow& tw,
                                     const string& networkCode,
                                     const string& stationCode,
                                     const string& locationCode,
                                     const string& channelCode) const
{
	IO::RecordStreamPtr rs = IO::RecordStream::Open( _cfg.xcorr.recordStreamURL.c_str() );
	if ( rs == nullptr )
	{
		string msg = "Cannot open RecordStream: " + _cfg.xcorr.recordStreamURL;
		throw runtime_error(msg);
	}

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
		string msg = stringify("Data could not be loaded for stream %s.%s.%s.%s from %s length %.2f)",
		                       networkCode.c_str(), stationCode.c_str(),
		                       locationCode.c_str(), channelCode.c_str(),
		                       tw.startTime().iso().c_str(), tw.length());
		throw runtime_error(msg);
	}

    GenericRecordPtr trace = new GenericRecord();
    
	if ( !merge(*trace, *seq) )
    {
		string msg = stringify("Data records could not be merged into a single trace (%s.%s.%s.%s from %s length %.2f)",
		                       networkCode.c_str(), stationCode.c_str(),
		                       locationCode.c_str(), channelCode.c_str(),
		                       tw.startTime().iso().c_str(), tw.length());
		throw runtime_error(msg);
	}

	return trace;
}


bool HypoDD::merge(GenericRecord &trace, const RecordSequence& seq) const
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
	ArrayPtr arr = ArrayFactory::Create(datatype, datatype, 0, nullptr);

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


bool HypoDD::trim(GenericRecord &trace, const Core::TimeWindow& tw) const
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
                    int order, double fmin, double fmax, double fsamp) const
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


GenericRecordPtr resample(const GenericRecordPtr &trace, int sf, bool average)
{
	if ( sf <= 0 )
		return nullptr;

	GenericRecordPtr newTrace = new GenericRecord(*trace);

	if ( newTrace->samplingFrequency() == sf )
		return newTrace;

	DoubleArray *data = DoubleArray::Cast(newTrace->data());
	double step = newTrace->samplingFrequency() / sf;

	int w = average?step*0.5 + 0.5:0;
	int i = 0;
	double fi = 0.0;
	int cnt = data->size();

	if ( w <= 0 )
	{
		while ( fi < cnt ) {
			(*data)[i++] = (*data)[(int)fi];
			fi += step;
		}
	}
	else
	{
		while ( fi < cnt )
		{
			int ci = (int)fi;
			double scale = 1.0;
			double v = (*data)[ci];

			for ( int g = 1; g < w; ++g )
			{
				if ( ci >= g )
				{
					v += (*data)[ci-g];
					scale += 1.0;
				}

				if ( ci+g < cnt )
				{
					v += (*data)[ci+g];
					scale += 1.0;
				}
			}

			v /= scale;

			(*data)[i++] = v;
			fi += step;
		}
	}

	data->resize(i);
	newTrace->setSamplingFrequency((double)sf);
	newTrace->dataUpdated();

	return newTrace;
}


} // HDD
} // Seiscomp
