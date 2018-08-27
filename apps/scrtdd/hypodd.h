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

#ifndef __SEISCOMP_APPLICATIONS_HYPODD_H__
#define __SEISCOMP_APPLICATIONS_HYPODD_H__

#include <seiscomp3/core/baseobject.h>
#include <seiscomp3/core/genericrecord.h>
#include <seiscomp3/core/recordsequence.h>
#include <seiscomp3/datamodel/databasequery.h>

#include <map>
#include <vector>

#include "hypodd.h"

namespace Seiscomp {

DEFINE_SMARTPOINTER(BgCatalog);

// DD background catalog
class BgCatalog : public Core::BaseObject {
	public:
		struct Station {
			std::string id;
			double latitude;
			double longitude;
			double elevation; // km
			// seiscomp info
			std::string networkCode;
			std::string stationCode;
		};

		struct Event {
			std::string id;
			int year;
			int month;
			int day;
			int hour;
			int minute;
			double second;
			double latitude;
			double longitude;
			double depth;   // km
			double magnitude;
			double horiz_err;
			double depth_err;
			double tt_residual;
			// seiscomp info
			std::string originId;
			std::string eventId;
		};

		struct Phase {
			std::string id;
			std::string event_id;
			std::string station_id;
			double travel_time;  // second
			double weight;       // 0-1 interval
			std::string type;
			// seiscomp info
			std::string networkCode;
			std::string stationCode;
			std::string locationCode;
			std::string channelCode;
		};

		BgCatalog(const std::map<std::string,Station>& stations,
                  const std::map<std::string,Event>& events,
                  const std::multimap<std::string,Phase>& phases);
		BgCatalog(const std::string& stationFile,
		          const std::string& catalogFile,
		          const std::string& phaFile);
		BgCatalog(const std::vector<std::string>& ids, DataModel::DatabaseQuery* query);
		BgCatalog(const std::string& idFile, DataModel::DatabaseQuery* query);
		const std::map<std::string,Station>& getStations() { return _stations;}
		const std::map<std::string,Event>& getEvents() { return _events;}
		const std::multimap<std::string,Phase>& getPhases() { return _phases;}
		BgCatalogPtr merge(BgCatalogPtr other);

	private:
		void initFromIds(const std::vector<std::string>& ids, DataModel::DatabaseQuery* query);
		DataModel::Station* findStation(const std::string& netCode, const std::string& staCode,
		                                Seiscomp::Core::Time, DataModel::DatabaseQuery* query);

		std::map<std::string,Station> _stations; // indexed by station id
		std::map<std::string,Event> _events; //indexed by event id
		std::multimap<std::string,Phase> _phases; //indexed by event id
};

struct HypoDDConfig {
	// catalog relocation specific: ph2dt config
	struct {
		std::string exec = "ph2dt";
		int minwght; // MINWGHT: min. pick weight allowed [-1]
		int maxdist; // MAXDIST: max. distance in km between event pair and stations [200]
		int maxsep;  // MAXSEP: max. hypocentral separation in km [10]
		int maxngh;  // MAXNGH: max. number of neighbors per event [10]
		int minlnk;  // MINLNK: min. number of links required to define a neighbor [8]
		int minobs;  // MINOBS: min. number of links per pair saved [8]
		int maxobs;  // MAXOBS: max. number of links per pair saved [20]
	} ph2dt;
	// single event relocation specific
	struct {
		double maxIEdis_CT = 20.0;   // Max interevent-distance for ct (km)
		double maxIEdis_CC = 10.0;   // Max interevent-distance for cc (km)
		int minNumNeigh_CT = 6;      // Min neighbors in DDBGC for ct (fail if not enough)
		int maxNumNeigh_CT = 20;     // Max neighbors in DDBGC for ct (furthest events are discarded)
		int minNumNeigh_CC = 3;      // Min neighbors in DDBGC for cc (fail if not enough)
		int maxNumNeigh_CC = 20;     // Max neighbors in DDBGC for cc (furthest events are discarded)
		int minDTperEvt_CT = 6;      // Min dt to use an event for ct (Including P+S)
		int minDTperEvt_CC = 4;      // Min pairs to use an event for cc

		double maxESdis_CT = 80.0;   // Max epi-sta epidistance for ct 
		double maxESdis_CC = 50.0;   // Max epi-sta epidistance for cc 
		double minWeight_CT = 0.05;  // Min weight of phases to be considered CT (0-1)
		double minWeight_CC = 0.0;   // Min weight of phases to be considered CC - Not used at the moment, maybe later, at the moment use threshold in pyXCorr.py
		std::vector<std::string> allowedPhases = {"P", "S"};
	} ddse;
	// hypodd executable specific
	struct {
		std::string exec = "hypodd";
		std::string ctrlFile;
	} hypodd;
	// cross correlation specific
	struct {
		std::string recordStreamURL;
		int filterOrder = 3;
		double filterFmin = -1;
		double filterFmax = -1;
		double filterFsamp = 0;
	} xcorr;
};

DEFINE_SMARTPOINTER(HypoDD);

class HypoDD : public Core::BaseObject {
	public:
		HypoDD(const BgCatalogPtr& input, const HypoDDConfig& cfg, std::string workingDir);
		~HypoDD();
		BgCatalogPtr relocateCatalog();
		BgCatalogPtr relocateSingleEvent(DataModel::Origin *org);
		void setWorkingDirCleanup(bool cleanup) { _workingDirCleanup = cleanup; }
	private:
		void createStationDatFile(std::string staFileName, BgCatalogPtr catalog=nullptr);
		void createPhaseDatFile(std::string catFileName, BgCatalogPtr catalog=nullptr);
		void createEventDatFile(std::string eventFileName, BgCatalogPtr catalog=nullptr);
		void createDtCtFile(BgCatalogPtr catalog, BgCatalogPtr org, std::string dtctFile);
		void xcorrCatalog(std::string dtctFile, std::string dtccFile);
		void xcorrSingleEvent(BgCatalogPtr catalog, BgCatalogPtr org, std::string dtccFile);
		void runPh2dt(std::string workingDir, std::string stationFile, std::string phaseFile);
		void runHypodd(std::string workingDir, std::string dtccFile, std::string dtctFile,
		               std::string eventFile, std::string stationFile);
		BgCatalogPtr loadRelocatedCatalog(std::string ddrelocFile);
		double computeDistance(double lat1, double lon1, double depth1,
		                       double lat2, double lon2, double depth2);
		BgCatalogPtr selectNeighbouringEvents(BgCatalogPtr catalog, BgCatalogPtr org,
		                                      double maxIEdis, int minNumNeigh=0,
		                                      int maxNumNeigh=0, int minDTperEvt=0);
		BgCatalogPtr extractEvent(BgCatalogPtr catalog, std::string eventId);
		GenericRecordPtr loadWaveform(const Core::Time& starttime,
		                              double duration,
		                              const std::string& networkCode,
		                              const std::string& stationCode,
		                              const std::string& locationCode,
		                              const std::string& channelCode);
		bool merge(GenericRecord &trace, const RecordSequence& seq);
		bool trim(GenericRecord &trace, const Core::TimeWindow& tw);
		void filter(GenericRecord &trace, bool demeaning=true,
                    int order=3, double fmin=-1, double fmax=-1, double fsamp=0);
		std::string generateWorkingSubDir(const DataModel::Origin *org);
	private:
		std::string _workingDir;
		std::string _controlFile;
		BgCatalogPtr _ddbgc;
		HypoDDConfig _cfg;
		bool _workingDirCleanup = true;
};

}

#endif
