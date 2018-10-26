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
namespace HDD {

DEFINE_SMARTPOINTER(Catalog);

// DD background catalog
class Catalog : public Core::BaseObject {
	public:
		struct Station {
			std::string id;
			double latitude;
			double longitude;
			double elevation; // meter
			std::string networkCode;
			std::string stationCode;
			// this equality works between multiple catalogs
			bool operator==(const Station& other) const
			{
			 return (networkCode == other.networkCode) &&
			        (stationCode == other.stationCode);
			}
			bool operator!=(const Station& other) const
			{
				return !operator==(other);
			}
			operator std::string() const
			{
				return id;
			}
		};

		struct Event {
			unsigned id; // makes it unique in the catalog 
			Core::Time time;
			double latitude;
			double longitude;
			double depth;   // km
			double magnitude;
			double horiz_err;
			double depth_err;
			double tt_residual;
			// equality works between multiple catalogs
			bool operator==(const Event& other) const
			{
			 return (time == other.time) &&
			        (latitude == other.latitude) &&
			        (longitude == other.longitude) &&
			        (depth == other.depth) &&
			        (magnitude == other.magnitude);
			}
			bool operator!=(const Event& other) const
			{
				return !operator==(other);
			}
			operator std::string() const
			{
				return std::to_string(id);
			}
		};

		struct Phase {
			unsigned eventId;
			std::string stationId;
			Core::Time time;
			std::string type;
			double weight;       // 0-1 interval
			std::string networkCode;
			std::string stationCode;
			std::string locationCode;
			std::string channelCode;
			// equality works between multiple catalogs
			bool operator==(const Phase& other) const
			{
				return (time == other.time) &&
				       (type == other.type) &&
				       (networkCode == other.networkCode) &&
				       (stationCode == other.stationCode) &&
				       (locationCode == other.locationCode) &&
				       (channelCode == other.channelCode);
			}
			bool operator!=(const Phase& other) const
			{
				return !operator==(other);
			}
			operator std::string() const
			{
				return type + " " + time.iso() + " " + networkCode + "." +
				       stationCode + "." + locationCode + "." + channelCode + 
				       " evId " + std::to_string(eventId)  + " staId " + stationId;
			}
		};

		Catalog();
		virtual ~Catalog() { }
		Catalog(const std::map<std::string,Station>& stations,
                  const std::map<unsigned,Event>& events,
                  const std::multimap<unsigned,Phase>& phases);
		Catalog(const std::string& stationFile,
		          const std::string& catalogFile,
		          const std::string& phaFile);
		Catalog(const std::vector<std::string>& ids, DataModel::DatabaseQuery* query);
		Catalog(const std::string& idFile, DataModel::DatabaseQuery* query);

		CatalogPtr merge(const CatalogPtr& other) const;

		const std::map<std::string,Station>& getStations() const { return _stations;}
		const std::map<unsigned,Event>& getEvents() const { return _events;}
		const std::multimap<unsigned,Phase>& getPhases() const { return _phases;}

		std::map<std::string,Station>::const_iterator searchStation(const Station&) const;
		std::map<unsigned,Event>::const_iterator searchEvent(const Event&) const;
		std::map<unsigned,Phase>::const_iterator searchPhase(const Phase&) const;

		bool addStation(const Station&, bool checkDuplicate);
		bool addEvent(const Event&, bool checkDuplicate);
		bool addPhase(const Phase&, bool checkDuplicate);

		void writeToFile(std::string eventFile,
		                 std::string phaseFile,
		                 std::string stationFile) const;

	private:
		void initFromIds(const std::vector<std::string>& ids, DataModel::DatabaseQuery* query);
		DataModel::Station* findStation(const std::string& netCode, const std::string& staCode,
		                                Seiscomp::Core::Time, DataModel::DatabaseQuery* query) const;

		std::map<std::string,Station> _stations; // indexed by station id
		std::map<unsigned,Event> _events; //indexed by event id
		std::multimap<unsigned,Phase> _phases; //indexed by event id
};

struct Config {
	
	std::vector<std::string> validPphases = {"P"};
	std::vector<std::string> validSphases = {"S"};

	// ph2dt config specifig (catalog relocation only)
	struct {
		std::string exec = "ph2dt";
		double minwght; // MINWGHT: min. pick weight allowed [-1]
		double maxdist; // MAXDIST: max. distance in km between event pair and stations [200]
		double maxsep;  // MAXSEP: max. hypocentral separation in km [10]
		int maxngh;  // MAXNGH: max. number of neighbors per event [10]
		int minlnk;  // MINLNK: min. number of links required to define a neighbor [8]
		int minobs;  // MINOBS: min. number of links per pair saved [8]
		int maxobs;  // MAXOBS: max. number of links per pair saved [20]
	} ph2dt;

	// hypodd executable specific
	struct {
		std::string exec = "hypodd";
		std::string ctrlFile;
	} hypodd;

	// differential travel time specific
	struct {
		double minWeight;  // Min weight of phases allowed (0-1)
		double maxESdist;   // Max epi-sta epidistance allowed
		double maxIEdist;   // Max interevent-distance allowed (km)
		int minNumNeigh;      // Min neighbors required
		int maxNumNeigh;     // Max neighbors allowed (furthest events are discarded)
		int minDTperEvt;      // Min differential times per event pair required (Including P+S)
	} dtt;

	// cross correlation specific
	struct {
		std::string recordStreamURL;
		double minWeight;  // Min weight of phases allowed (0-1)
		double maxESdist;   // Max epi-sta epidistance allowed
		double maxIEdist;   // Max interevent-distance allowed (km)
		int minNumNeigh;      // Min neighbors required
		int maxNumNeigh;     // Max neighbors allowed (furthest events are discarded)
		int minDTperEvt;      // Min differential times per event pair required (Including P+S)
		double minCoef;    // Min xcorr coefficient required (0-1)
		int filterOrder;
		double filterFmin;
		double filterFmax;
		double filterFsamp;
		double timeBeforePick; // secs
		double timeAfterPick;  // secs
		double maxDelay; //secs
	} xcorr;
};

DEFINE_SMARTPOINTER(HypoDD);

class HypoDD : public Core::BaseObject {
	public:
		HypoDD(const CatalogPtr& input, const Config& cfg, const std::string& workingDir);
		virtual ~HypoDD();
		CatalogPtr relocateCatalog();
		CatalogPtr relocateSingleEvent(const CatalogPtr& orgToRelocate);
		void setWorkingDirCleanup(bool cleanup) { _workingDirCleanup = cleanup; }
		bool workingDirCleanup() { return _workingDirCleanup; }
		void setUseWaveformDiskCache(bool cache) { _wfDiskCache = cache; }
		bool useWaveformDiskCache() { return _wfDiskCache; }
		
	private:
		CatalogPtr filterOutPhases(const CatalogPtr& catalog,
		                           const std::vector<std::string>& PphaseToKeep,
		                           const std::vector<std::string>& SphaseToKeep) const;
		void createStationDatFile(const std::string& staFileName, const CatalogPtr& catalog) const;
		void createPhaseDatFile(const std::string& catFileName, const CatalogPtr& catalog) const;
		void createEventDatFile(const std::string& eventFileName, const CatalogPtr& catalog) const;
		void createDtCtFile(const CatalogPtr& catalog,
		                    unsigned evToRelocateId,
		                    const std::string& dtctFile) const;
		void xcorrCatalog(const std::string& dtctFile, const std::string& dtccFile);
		void xcorrSingleEvent(const CatalogPtr& catalog,
		                      unsigned evToRelocateId,
		                      const std::string& dtccFile);
		bool xcorr(const GenericRecordPtr& tr1, const GenericRecordPtr& tr2, double maxDelay,
              double& delayOut, double& coeffOut) const;
		void runPh2dt(const std::string& workingDir, const std::string& stationFile, const std::string& phaseFile) const;
		void runHypodd(const std::string& workingDir, const std::string& dtccFile, const std::string& dtctFile,
		               const std::string& eventFile, const std::string& stationFile) const;
		CatalogPtr loadRelocatedCatalog(const std::string& ddrelocFile, const CatalogPtr& originalCatalog);
		double computeDistance(double lat1, double lon1, double depth1,
		                       double lat2, double lon2, double depth2);
		CatalogPtr selectNeighbouringEvents(const CatalogPtr& catalog, const Catalog::Event& refEv,
		                                      double maxESdis, double maxIEdis, int minNumNeigh=0,
		                                      int maxNumNeigh=0, int minDTperEvt=0);
		CatalogPtr extractEvent(const CatalogPtr& catalog, unsigned eventId) const;
		GenericRecordPtr getWaveform(const Catalog::Event& ev, const Catalog::Phase& ph,
		                             std::map<std::string,GenericRecordPtr>& cache, bool useDiskCache);
		GenericRecordPtr loadWaveform(const Core::Time& starttime,
		                              double duration,
		                              const std::string& networkCode,
		                              const std::string& stationCode,
		                              const std::string& locationCode,
		                              const std::string& channelCode,
		                              bool useDiskCache) const;
		GenericRecordPtr readWaveformFromRecordStream(const Core::TimeWindow& tw,
		                                              const std::string& networkCode,
		                                              const std::string& stationCode,
		                                              const std::string& locationCode,
		                                              const std::string& channelCode) const;
		bool merge(GenericRecord &trace, const RecordSequence& seq) const;
		bool trim(GenericRecord &trace, const Core::TimeWindow& tw) const;
		void filter(GenericRecord &trace, bool demeaning=true,
                    int order=3, double fmin=-1, double fmax=-1, double fsamp=0) const;
		std::string generateWorkingSubDir(const Catalog::Event& ev) const;
	private:
		std::string _workingDir;
		std::string _cacheDir;
		CatalogPtr _ddbgc;
		Config _cfg;
		bool _workingDirCleanup = true;
		bool _wfDiskCache = false;
		std::map<std::string, GenericRecordPtr> _wfCache;
};

}
}

#endif
