/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 *   You can redistribute and/or modify this program under the             *
 *   terms of the "SED Public License for Seiscomp Contributions"          *
 *                                                                         *
 *   You should have received a copy of the "SED Public License for        *
 *   Seiscomp Contributions" with this. If not, you can find it at         *
 *   http://www.seismo.ethz.ch/static/seiscomp_contrib/license.txt         *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   "SED Public License for Seiscomp Contributions" for more details      *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/


#ifndef __RTDD_APPLICATIONS_RTDD_H__
#define __RTDD_APPLICATIONS_RTDD_H__

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

        bool processOrigin(DataModel::Origin *origin,
                           DataModel::OriginPtr& relocatedOrg,
                           std::vector<DataModel::PickPtr>& relocatedOrgPicks,
                           const ProfilePtr& profile, bool forceProcessing=false,
                           bool allowManualOrigin=false, bool doSend=true);

        void relocateOrigin(DataModel::Origin *org, ProfilePtr profile,
                            DataModel::OriginPtr& newOrg,
                            std::vector<DataModel::PickPtr>& newOrgPicks);

        void convertOrigin(const HDD::CatalogCPtr& relocatedOrg,
                           ProfilePtr profile,     // can be nullptr
                           const DataModel::Origin *org, // can be nullptr
                           DataModel::OriginPtr& newOrg,
                           std::vector<DataModel::PickPtr>& newOrgPicks);

        void removedFromCache(DataModel::PublicObject *);

        ProfilePtr getProfile(const DataModel::Origin *origin, const std::string& forceProfile="");
        ProfilePtr getProfile(double latitude, double longitude, const std::string& forceProfile="");

        std::vector<DataModel::OriginPtr> fetchOrigins(const std::string& idFile, std::string options);

        struct Config {
            Config();

            std::string publicIDPattern;
            std::vector<std::string> activeProfiles;
            std::string workingDirectory;
            bool        keepWorkingFiles;
            bool        onlyPreferredOrigin;
            bool        allowManualOrigin;
            int         profileTimeAlive; //seconds
            bool        cacheWaveforms;

            // Mode
            bool        forceProcessing;
            bool        testMode;
            bool        dumpWaveforms;
            double      fExpiry;
            std::string originIDs;
            std::string eventXML;
            std::string forceProfile;
            std::string relocateProfile;
            std::string dumpCatalog;
            std::string mergeCatalogs;
            std::string dumpCatalogXML;
            std::string loadProfile;
            std::string evalXCorr;

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
            void evalXCorr();

            std::string name;
            std::string earthModelID;
            std::string methodID;
            std::string eventIDFile;
            std::string stationFile;
            std::string eventFile;
            std::string phaFile;
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
            unsigned            runCount;
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
