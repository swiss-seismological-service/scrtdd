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


#ifndef __SEISCOMP_APPLICATIONS_RTDD_H__
#define __SEISCOMP_APPLICATIONS_RTDD_H__

#include <seiscomp3/client/application.h>
#include <seiscomp3/processing/amplitudeprocessor.h>
#include <seiscomp3/datamodel/publicobjectcache.h>
#include <seiscomp3/datamodel/eventparameters.h>
#include <seiscomp3/datamodel/origin.h>
#include <seiscomp3/datamodel/amplitude.h>
#include <seiscomp3/datamodel/journaling.h>
#include <seiscomp3/seismology/ttt.h>
#include <seiscomp3/utils/timer.h>

#include "app.h"
#include "hypodd.h"

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp3/logging/log.h>

#include <map>
#include <set>
#include <vector>


namespace Seiscomp {

namespace DataModel {

class Pick;
class Origin;
class Event;

}


class RTDD : public Application {
	public:
		RTDD(int argc, char **argv);
		~RTDD();

		struct Region : public Core::BaseObject {
			virtual bool init(const Application* app, const std::string &prefix) = 0;
			virtual bool isInside(double lat, double lon) const = 0;
		};
		DEFINE_SMARTPOINTER(Region);

	protected:
		void createCommandLineDescription();
		bool validateParameters();

		bool init();
		bool run();
		void done();

		void handleMessage(Core::Message *msg);
		void addObject(const std::string&, DataModel::Object* object);
		void updateObject(const std::string&, DataModel::Object* object);
		void handleRecord(Record *rec) { /* we don't really need this */ }

		void handleTimeout();
		void checkProfileStatus();
		void runNewJobs();

	private:
		DEFINE_SMARTPOINTER(Process);
		DEFINE_SMARTPOINTER(Profile);
		DEFINE_SMARTPOINTER(Cronjob);

		bool addProcess(DataModel::PublicObject* obj);
		bool startProcess(Process *proc);
		void removeProcess(Process *proc);

		bool process(DataModel::Origin *origin);

		void removedFromCache(DataModel::PublicObject *);

		bool send(DataModel::Origin *org);

		DataModel::OriginPtr relocateOrigin(DataModel::Origin *org, ProfilePtr);

		struct Config {
			Config();

			std::string publicIDPattern;
			std::vector<std::string> activeProfiles;
			std::string workingDirectory;
			bool        keepWorkingFiles;
			bool        onlyPreferredOrigin;
			bool        processManualOrigin;
			int         profileTimeAlive; //seconds
			bool        cacheWaveforms;

            // Mode
			bool        forceProcessing;
			bool        testMode;
			double      fExpiry;
			std::string originIDs;
			std::string eventXML;
			std::string forceProfile;
			std::string relocateCatalog;
			std::string dumpCatalog;
			std::string loadCatalog;

            // cron
			int         wakeupInterval;
			bool        logCrontab;
			std::vector<int> delayTimes;

		};

		class Profile : public Core::BaseObject {
			public:
			Profile();
			void load(DataModel::DatabaseQuery* query,
			          DataModel::PublicObjectTimeSpanBuffer* cache,
			          DataModel::EventParameters* eventParameters,
			          const std::string& workingDir,
			          bool cleanupWorkingDir,
			          bool cacheWaveforms,
			          bool preloadData);
			void unload();
			bool isLoaded() { return loaded; }
			Core::TimeSpan inactiveTime() { return Core::Time::GMT() - lastUsage; }
			HDD::CatalogPtr relocateSingleEvent(DataModel::Origin *org);
			HDD::CatalogPtr relocateCatalog(bool force = true);
			bool addIncrementalCatalogEntry(DataModel::Origin *org);

			std::string name;
			std::string earthModelID;
			std::string methodID;
			std::string eventIDFile;
			std::string stationFile;
			std::string eventFile;
			std::string phaFile;
			std::string incrementalCatalogFile;
			RegionPtr   region;
			HDD::Config ddcfg;

			private:
			bool loaded;
			Core::Time lastUsage;
			HDD::HypoDDPtr hypodd;
			DataModel::DatabaseQuery* query;
			DataModel::PublicObjectTimeSpanBuffer* cache;
			DataModel::EventParameters* eventParameters;
		};

		struct Cronjob : public Core::BaseObject {
			std::list<Core::Time> runTimes;
		};

		struct Process : Core::BaseObject {
			Core::Time          created;
			Core::Time          lastRun;
			DataModel::PublicObjectPtr obj;
			CronjobPtr          cronjob;
		};

		typedef std::list<ProcessPtr>              ProcessQueue;
		typedef std::map<std::string, ProcessPtr>  Processes;
		typedef std::set<DataModel::PublicObjectPtr> Todos;

		ProcessQueue               _processQueue;
		Processes                  _processes;
		int                        _cronCounter;
		Todos                      _todos;

		DataModel::PublicObjectTimeSpanBuffer _cache;

		Config                     _config;
		std::list<ProfilePtr>      _profiles;

		DataModel::EventParametersPtr _eventParameters;

		Logging::Channel *_processingInfoChannel;
		Logging::Output  *_processingInfoOutput;

		ObjectLog        *_inputEvts;
		ObjectLog        *_inputOrgs;
		ObjectLog        *_outputOrgs;
};

}

#endif
