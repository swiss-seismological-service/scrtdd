/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU Affero General Public License as published*
 * by the Free Software Foundation, either version 3 of the License, or    *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This program is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#ifndef __RTDD_APPLICATIONS_RTDD_H__
#define __RTDD_APPLICATIONS_RTDD_H__

#include <map>
#include <set>
#include <vector>

#include "hdd/dd.h"

#include <seiscomp/client/streamapplication.h>
#include <seiscomp/datamodel/eventparameters.h>
#include <seiscomp/datamodel/publicobjectcache.h>
#include <seiscomp/utils/timer.h>

namespace Seiscomp {

namespace DataModel {

class Pick;
class Origin;

} // namespace DataModel

class RTDD : public Client::StreamApplication
{
public:
  RTDD(int argc, char **argv);
  ~RTDD();

  struct Region : public Core::BaseObject
  {
    virtual bool init(const Application *app, const std::string &prefix) = 0;
    virtual bool isInside(double lat, double lon) const                  = 0;
  };
  DEFINE_SMARTPOINTER(Region);

  const char *version() override { return RTDD_VERSION; }

protected:
  void createCommandLineDescription() override;
  bool validateParameters() override;

  bool init() override;
  bool run() override;
  void done() override;

  void handleMessage(Core::Message *msg) override;
  void addObject(const std::string &, DataModel::Object *object) override;
  void updateObject(const std::string &, DataModel::Object *object) override;
  void handleRecord(Record *rec) override
  {                     /* we don't really need this and this is never called */
    RecordPtr tmp(rec); // avoid memory leak
  }

  void handleTimeout() override;
  void checkProfileStatus();
  void runNewJobs();

private:
  DEFINE_SMARTPOINTER(Process);
  DEFINE_SMARTPOINTER(Profile);
  DEFINE_SMARTPOINTER(Cronjob);

  bool addProcess(DataModel::PublicObject *obj);
  bool startProcess(Process *proc);
  void removeProcess(Process *proc);

  bool processOrigin(DataModel::Origin *origin,
                     DataModel::OriginPtr &relocatedOrg,
                     std::vector<DataModel::PickPtr> &relocatedOrgPicks,
                     const ProfilePtr &profile,
                     bool forceProcessing,
                     bool allowAutomaticOrigin,
                     bool allowManualOrigin,
                     bool doSend);

  void relocateOrigin(DataModel::Origin *org,
                      ProfilePtr profile,
                      DataModel::OriginPtr &newOrg,
                      std::vector<DataModel::PickPtr> &newOrgPicks);

  void removedFromCache(DataModel::PublicObject *);

  std::unique_ptr<HDD::Catalog> getCatalog(
      const std::string &catalogPath,
      std::unordered_map<unsigned, DataModel::OriginPtr> *idmap = nullptr);
  ProfilePtr getProfile(const std::string &profile);
  ProfilePtr getProfile(const DataModel::Origin *origin,
                        const std::string &forceProfile = "");
  ProfilePtr getProfile(double latitude,
                        double longitude,
                        const std::string &forceProfile = "");

  void loadProfile(ProfilePtr profile,
                   const HDD::Catalog *alternativeCatalog = nullptr);

  std::vector<DataModel::OriginPtr> fetchOrigins(const std::string &idFile,
                                                 std::string options);
  struct Config
  {
    std::vector<std::string> activeProfiles;
    std::string workingDirectory = "/tmp/rtdd";
    bool saveProcessingFiles     = false;
    bool onlyPreferredOrigin     = true;
    bool allowAutomaticOrigin    = true;
    bool allowManualOrigin       = true;
    int profileTimeAlive         = -1; // seconds
    bool cacheWaveforms          = false;
    bool cacheAllWaveforms       = false;
    bool debugWaveforms          = false;

    // Mode
    bool forceProcessing = false;
    bool testMode        = false;
    double fExpiry       = 1.0;
    std::string originIDs;
    std::string eventXML;
    std::string forceProfile;
    std::string relocateCatalog;
    std::string dumpClusters;
    std::string dumpCatalog;
    std::string dumpCatalogOptions;
    std::string dumpWaveforms;
    std::string mergeCatalogs;
    std::string evalXCorr;
    std::string xcorrCache;
    std::string reloadProfileMsg;
    bool loadProfileWf = false;

    // cron
    int wakeupInterval = 1; // sec
    bool logCrontab    = true;
    std::vector<int> delayTimes;
  };

  class Profile : public Core::BaseObject
  {
  public:
    Profile();
    ~Profile();
    void load(DataModel::DatabaseQuery *query,
              DataModel::PublicObjectTimeSpanBuffer *cache,
              DataModel::EventParameters *eventParameters,
              const std::string &workingDir,
              bool saveProcessingFiles,
              bool cacheWaveforms,
              bool cacheAllWaveforms,
              const HDD::Catalog *alternativeCatalog = nullptr);
    void unload();
    bool isLoaded() { return loaded; }
    void preloadWaveforms();
    void freeResources();
    Core::TimeSpan inactiveTime() { return Core::Time::GMT() - lastUsage; }
    std::unique_ptr<HDD::Catalog> relocateSingleEvent(DataModel::Origin *org);
    std::unique_ptr<HDD::Catalog> relocateCatalog(const std::string &xcorrFile);
    void evalXCorr(const std::string &xcorrFile);
    void dumpWaveforms();
    void dumpClusters();

    std::string name;
    std::string earthModelID;
    std::string methodID;
    std::string eventIDFile;
    std::string stationFile;
    std::string eventFile;
    std::string phaFile;
    std::string tttType;
    std::string tttModel;
    std::string recordStreamURL;
    RegionPtr region;
    HDD::Config ddCfg;
    HDD::ClusteringOptions singleEventClustering;
    HDD::ClusteringOptions multiEventClustering;
    HDD::SolverOptions solverCfg;
    bool detectMissingPhasesAuto;
    bool detectMissingPhasesManual;

  private:
    bool loaded;
    Core::Time lastUsage;
    std::unique_ptr<HDD::DD> dd;
    DataModel::DatabaseQuery *query;
    DataModel::PublicObjectTimeSpanBuffer *cache;
    DataModel::EventParameters *eventParameters;
  };

  struct Cronjob : public Core::BaseObject
  {
    std::list<Core::Time> runTimes;
  };

  struct Process : Core::BaseObject
  {
    Core::Time created;
    unsigned runCount;
    DataModel::PublicObjectPtr obj;
    CronjobPtr cronjob;
  };

  typedef std::list<ProcessPtr> ProcessQueue;
  typedef std::map<std::string, ProcessPtr> Processes;
  typedef std::set<DataModel::PublicObjectPtr> Todos;

  ProcessQueue _processQueue;
  Processes _processes;
  int _cronCounter;
  Todos _todos;

  DataModel::PublicObjectTimeSpanBuffer _cache;

  Config _config;
  std::list<ProfilePtr> _profiles;

  DataModel::EventParametersPtr _eventParameters;

  ObjectLog *_inputEvts;
  ObjectLog *_inputOrgs;
  ObjectLog *_outputOrgs;
};

} // namespace Seiscomp

#endif
