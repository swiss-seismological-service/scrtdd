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

#include "rtddtool.h"
#include "hdd/csvreader.h"
#include "hdd/cvttt.h"
#include "hdd/nllttt.h"
#include "hddsc/sclog.h"
#include "hddsc/scttt.h"
#include "hddsc/scutils.h"
#include "hddsc/scwaveform.h"
#include "msg.h"

#define SEISCOMP_COMPONENT RTDD
#include <seiscomp/logging/log.h>

#include <seiscomp/client/inventory.h>
#include <seiscomp/core/genericrecord.h>
#include <seiscomp/core/strings.h>
#include <seiscomp/core/system.h>
#include <seiscomp/datamodel/event.h>
#include <seiscomp/datamodel/magnitude.h>
#include <seiscomp/datamodel/origin.h>
#include <seiscomp/datamodel/parameter.h>
#include <seiscomp/datamodel/parameterset.h>
#include <seiscomp/datamodel/pick.h>
#include <seiscomp/datamodel/utils.h>
#include <seiscomp/io/archive/xmlarchive.h>
#include <seiscomp/io/records/mseedrecord.h>
#include <seiscomp/logging/channel.h>
#include <seiscomp/logging/filerotator.h>
#include <seiscomp/math/geo.h>
#include <seiscomp/utils/files.h>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace Seiscomp::DataModel;
using HDD::SCAdapter::addToCatalog;
using HDD::SCAdapter::convertOrigin;
using HDD::SCAdapter::printEvalXcorrStats;
using Seiscomp::Core::fromString;
using Seiscomp::Core::stringify;
using PhaseType  = HDD::Catalog::Phase::Type;
using AQ_ACTION  = HDD::SolverOptions::AQ_ACTION;
using DataSource = HDD::SCAdapter::DataSource;

namespace {

using namespace Seiscomp;

template <class T>
bool configGetTypedList(const Client::Application *app,
                        const string &parameter,
                        vector<T> &values,
                        unsigned requiredItems,
                        bool allowEmpty)
{
  vector<string> stringValues;

  try
  {
    stringValues = app->configGetStrings(parameter);
  }
  catch (...)
  {}

  if (stringValues.size() != requiredItems &&
      !(allowEmpty && stringValues.empty()))
  {
    SEISCOMP_ERROR("%s: expected %u values but got %zu", parameter.c_str(),
                   requiredItems, stringValues.size());
    return false;
  }

  for (unsigned i = 0; i < stringValues.size(); i++)
  {
    T val;
    if (!fromString(val, stringValues[i]))
    {
      SEISCOMP_ERROR("%s: invalid value(s)", parameter.c_str());
      return false;
    }
    values.push_back(val);
  }
  return true;
}

// Rectangular region class defining a rectangular region
// by latmin, lonmin, latmax, lonmax.
struct RectangularRegion : public Seiscomp::RTDD::Region
{
  RectangularRegion() {}

  bool init(const Client::Application *app, const string &prefix)
  {

    vector<double> region;
    if (!configGetTypedList(app, prefix + "region", region, 4, true))
      return false;

    if (region.empty())
      isEmpty = true;
    else
    {
      isEmpty = false;
      latMin  = region[0];
      lonMin  = region[1];
      latMax  = region[2];
      lonMax  = region[3];
    }

    return true;
  }

  bool isInside(double lat, double lon) const
  {
    if (isEmpty) return true;

    if (lat < latMin || lat > latMax) return false;

    double lonRange = lonMax - lonMin;
    if (lonRange < 0) lonRange += 360.0;

    double lonDelta = lon - lonMin;
    if (lonDelta < 0) lonDelta += 360.0;

    return lonDelta <= lonRange;
  }

  bool isEmpty;
  double latMin, lonMin;
  double latMax, lonMax;
};

// Rectangular region class defining a circular region
// by lat, lon, radius.
struct CircularRegion : public Seiscomp::RTDD::Region
{
  CircularRegion() {}

  bool init(const Client::Application *app, const string &prefix)
  {

    vector<double> region;
    if (!configGetTypedList(app, prefix + "region", region, 3, true))
      return false;

    if (region.empty())
      isEmpty = true;
    else
    {
      isEmpty = false;
      lat     = region[0];
      lon     = region[1];
      radius  = region[2];
    }

    return true;
  }

  bool isInside(double lat, double lon) const
  {
    if (isEmpty) return true;

    double distance, az, baz;
    Math::Geo::delazi(this->lat, this->lon, lat, lon, &distance, &az, &baz);
    double distKm = Math::Geo::deg2km(distance);
    return distKm <= radius;
  }

  bool isEmpty;
  double lat, lon, radius;
};

string makeUpper(const string &src)
{
  string dest = src;
  for (size_t i = 0; i < src.size(); ++i) dest[i] = toupper(src[i]);
  return dest;
}

bool startsWith(const string &haystack,
                const string &needle,
                bool caseSensitive = true)
{
  string _haystack = haystack;
  string _needle   = needle;
  if (!caseSensitive)
  {
    _haystack = makeUpper(haystack);
    _needle   = makeUpper(needle);
  }
  return _haystack.compare(0, _needle.length(), _needle) == 0;
}

std::vector<std::string> splitString(const std::string &str,
                                     const std::string &split)
{
  std::vector<std::string> tokens;
  boost::split(tokens, str, boost::is_any_of(split), boost::token_compress_on);
  return tokens;
}

Core::Time now; // this is tricky, I don't like globals

} // unnamed namespace

namespace Seiscomp {

RTDD::RTDD(int argc, char **argv) : StreamApplication(argc, argv)
{
  setAutoApplyNotifierEnabled(true);
  setInterpretNotifierEnabled(true);

  setLoadStationsEnabled(true);

  setPrimaryMessagingGroup("LOCATION");

  addMessagingSubscription("EVENT");
  addMessagingSubscription("LOCATION");
  addMessagingSubscription("PICK"); // this is only for caching picks
  addMessagingSubscription("SERVICE_REQUEST");

  setAutoAcquisitionStart(false);
  setAutoCloseOnAcquisitionFinished(false);

  _cache.setPopCallback(boost::bind(&RTDD::removedFromCache, this, _1));
}

RTDD::~RTDD() {}

void RTDD::createCommandLineDescription()
{
  StreamApplication::createCommandLineDescription();

  commandline().addGroup("Mode");
  commandline().addOption(
      "Mode", "reloc-catalog",
      "Relocate the catalog passed as argument in multi-event mode. The "
      "input can be a single file (containing seiscomp origin ids) or a file "
      "triplet (station.csv,event.csv,phase.csv). For events stored "
      "in a XML files add the --ep option. Use in combination with --profile",
      &_config.relocateCatalog, true);
  commandline().addOption(
      "Mode", "origin-id,O",
      "Relocate  the origin (or multiple comma-separated origins) in "
      "signle-event mode and send a message. Each origin will be processed "
      "accordingly to the matching profile region unless the --profile option "
      " is used.",
      &_config.originIDs, true);
  commandline().addOption(
      "Mode", "ep",
      "Event parameters XML file for offline processing of contained origins "
      "(implies --test option). Each contained origin will be processed in "
      "signle-event mode unless --reloc-catalog is provided, which enable "
      "multi-event mode.",
      &_config.eventXML, true);
  commandline().addOption(
      "Mode", "eval-xcorr",
      "Compute cross-correlation statistics for the catalog passed as "
      "argument. The input can be a single file (containing seiscomp origin "
      "ids) or a file triplet (station.csv,event.csv,phase.csv). Use in "
      "combination with --profile",
      &_config.evalXCorr, true);
  commandline().addOption(
      "Mode", "dump-clusters",
      "Find clusters in the catalog passed as argument and save them in "
      "the working directory."
      "The catalog can be a single file (containing seiscomp origin "
      "ids) or a file triplet (station.csv,event.csv,phase.csv). Use "
      "in combination with --profile. The clusters will be saved into "
      "the working directory",
      &_config.dumpClusters, true);
  commandline().addOption(
      "Mode", "dump-wf",
      "Dump processed waveforms of the catalog passed as argument "
      "in the current working directory."
      "The catalog can be a single file (containing seiscomp origin "
      "ids) or a file triplet (station.csv,event.csv,phase.csv). Use "
      "in combination with --profile.",
      &_config.dumpWaveforms, true);
  commandline().addOption(
      "Mode", "load-profile-wf",
      "Load catalog waveforms from the configured recordstream and "
      "save them into the profile working directory. Use in "
      "combination with --profile");
  commandline().addOption(
      "Mode", "send-reload-profile-msg",
      "Send a message to any running scrtdd module requesting to "
      "reload a specific profile passed as argument",
      &_config.reloadProfileMsg, true);
  commandline().addGroup("Catalog");
  commandline().addOption(
      "Catalog", "dump-catalog",
      "Dump the seiscomp event/origin id file passed as argument into "
      "a catalog file triplet (station.csv,event.csv,phase.csv).",
      &_config.dumpCatalog, true);
  commandline().addOption(
      "Catalog", "dump-catalog-options",
      "Allows the --dump-catalog option to accept event ids besides origin "
      "ids. For each event id an origin will be selected following the "
      "provided options whose format is: "
      "'type,evalmode,includeCreator,excludeCreator,region', "
      "where type=preferred|last|first  evalmode=any|onlyManual|onlyAutomatic  "
      "includeCreator=any|author|methodID  excludeCreator=none|author|methodID "
      "region=any|profileName e.g. to select preferred origins of the input"
      "event ids that lie within the region defined for 'myProfile' use "
      "'preferred,any,any,none,myProfile'",
      &_config.dumpCatalogOptions, false);
  commandline().addOption(
      "Catalog", "merge-catalogs",
      "Merge in a single catalog all the catalog file triplets "
      "(station1.csv,event1.csv,phase1.csv,station2.csv,event2.csv,"
      "phase2.csv,...) passed as arguments.",
      &_config.mergeCatalogs, false);
  commandline().addOption(
      "Catalog", "merge-catalogs-keepid",
      "Similar to the --merge-catalogs option but events keep their ids. If "
      "multiple events share the same id, subsequent events will be "
      "discarded.",
      &_config.mergeCatalogs, false);
  commandline().addGroup("ModeOptions");
  commandline().addOption(
      "ModeOptions", "profile",
      "To be used in combination with other options: select the "
      "profile configuration to use",
      &_config.forceProfile, true);
  commandline().addOption(
      "ModeOptions", "xcorr-cache",
      "Specify a file containing precomputed cross-correlation values",
      &_config.xcorrCache, true);
  commandline().addOption(
      "ModeOptions", "expiry,x",
      "Defines the time span in hours after which objects expire.",
      &_config.fExpiry, true);
  commandline().addOption(
      "ModeOptions", "cache-wf-all",
      "All waveforms will be saved to disk cache, even temporarily "
      "ones. Normally only catalog phase waveforms are cached to disk. "
      "This is useful to speed up debugging/testing when the same "
      "origins are repeatedly processed with --origin or --ep options.");
  commandline().addOption(
      "ModeOptions", "test",
      "Test mode, no messages are sent when relocating a single event");
  commandline().addOption("ModeOptions", "xmlout",
                          "Enable XML output when combined with "
                          "--reloc-catalog or --oring-id options");
  commandline().addOption("ModeOptions", "inherit-mag",
                          "Origins inherit the source origin magnitudes "
                          "when combined with --xmlout option");
}

bool RTDD::validateParameters()
{
  if (!StreamApplication::validateParameters()) return false;

  Environment *env = Environment::Instance();

  _config.workingDirectory =
      env->absolutePath(configGetPath("workingDirectory"));

  _config.saveProcessingFiles  = configGetBool("saveProcessingFiles");
  _config.onlyPreferredOrigin  = configGetBool("onlyPreferredOrigins");
  _config.allowAutomaticOrigin = configGetBool("automaticOrigins");
  _config.allowManualOrigin    = configGetBool("manualOrigins");
  _config.activeProfiles       = configGetStrings("activeProfiles");
  _config.logCrontab           = configGetBool("cron.logging");
  _config.delayTimes           = configGetInts("cron.delayTimes");
  _config.profileTimeAlive     = configGetInt("performance.profileTimeAlive");
  _config.cacheWaveforms       = configGetBool("performance.cacheWaveforms");
  _config.testMode             = commandline().hasOption("test");
  _config.loadProfileWf        = commandline().hasOption("load-profile-wf");
  _config.cacheAllWaveforms    = commandline().hasOption("cache-wf-all");

  // disable messaging (offline mode) with certain command line options
  if (!_config.eventXML.empty() || !_config.dumpCatalog.empty() ||
      !_config.mergeCatalogs.empty() || !_config.evalXCorr.empty() ||
      !_config.relocateCatalog.empty() || _config.loadProfileWf ||
      !_config.dumpWaveforms.empty() || !_config.dumpClusters.empty() ||
      (!_config.originIDs.empty() && _config.testMode))
  {
    SEISCOMP_INFO("Disable messaging");
    setMessagingEnabled(false);
    _config.testMode = true; // we won't send any message
  }

  bool profileRequireDB = false;
  if (!_config.evalXCorr.empty() &&
      ::splitString(_config.evalXCorr, ",").size() ==
          1) // single file containing origin ids
  {
    profileRequireDB = true;
  }

  if (!_config.dumpWaveforms.empty() &&
      ::splitString(_config.dumpWaveforms, ",").size() ==
          1) // single file containing origin ids
  {
    profileRequireDB = true;
  }

  if (!_config.dumpClusters.empty() &&
      ::splitString(_config.dumpClusters, ",").size() ==
          1) // single file containing origin ids
  {
    profileRequireDB = true;
  }

  if (!_config.relocateCatalog.empty() &&
      ::splitString(_config.relocateCatalog, ",").size() ==
          1) // single file containing origin ids
  {
    profileRequireDB = true;
  }

  // make sure to load the profile passed via command line too
  std::vector<string> profilesToLoad(_config.activeProfiles);
  if (!_config.forceProfile.empty() &&
      std::find(profilesToLoad.begin(), profilesToLoad.end(),
                _config.forceProfile) == profilesToLoad.end())
  {
    profilesToLoad.push_back(_config.forceProfile);
  }

  bool profilesOK = true;
  for (const string &profileName : profilesToLoad)
  {
    if (profileName.empty()) continue;

    ProfilePtr prof = new Profile;
    prof->name      = profileName;
    string prefix   = string("profile.") + prof->name + ".";

    try
    {
      prof->earthModelID = configGetString(prefix + "earthModelID");
    }
    catch (...)
    {}
    try
    {
      prof->methodID = configGetString(prefix + "methodID");
    }
    catch (...)
    {}
    if (!startsWith(prof->methodID, "RTDD", false))
    {
      prof->methodID = "RTDD" + prof->methodID;
    }
    string regionType;
    try
    {
      regionType = makeUpper(configGetString(prefix + "regionType"));
    }
    catch (...)
    {}
    if (regionType == "RECTANGULAR")
      prof->region = new RectangularRegion;
    else
      prof->region = new CircularRegion;

    if (prof->region == nullptr)
    {
      SEISCOMP_ERROR("profile.%s: invalid region type: %s", prof->name.c_str(),
                     regionType.c_str());
      profilesOK = false;
      continue;
    }

    if (!prof->region->init(this, prefix))
    {
      SEISCOMP_ERROR("profile.%s: invalid region parameters",
                     prof->name.c_str());
      profilesOK = false;
      continue;
    }

    prefix = string("profile.") + prof->name + ".catalog.";

    // For the catalog we can have either a single file (origin ids) or three
    // files ( event, station, phase ). They can also be left empty
    try
    {
      prof->stationFile =
          env->absolutePath(configGetPath(prefix + "stationFile"));
      prof->phaFile   = env->absolutePath(configGetPath(prefix + "phaFile"));
      prof->eventFile = env->absolutePath(configGetPath(prefix + "eventFile"));
    }
    catch (...)
    {
      try
      {
        prof->eventIDFile =
            env->absolutePath(configGetPath(prefix + "eventFile"));
        profileRequireDB = true;
      }
      catch (...)
      {}
    }

    try
    {
      prof->ddCfg.validPphases = configGetStrings(prefix + "P-Phases");
    }
    catch (...)
    {
      prof->ddCfg.validPphases = {"Pg", "P"};
    }
    try
    {
      prof->ddCfg.validSphases = configGetStrings(prefix + "S-Phases");
    }
    catch (...)
    {
      prof->ddCfg.validSphases = {"Sg", "S"};
    }

    prefix = string("profile.") + prof->name +
             ".doubleDifferenceSystem.eventFiltering.";
    try
    {
      prof->singleEventClustering.minDTperEvt =
          configGetInt(prefix + "minNumPhases");
      prof->multiEventClustering.minDTperEvt =
          configGetInt(prefix + "minNumPhases");
    }
    catch (...)
    {
      prof->singleEventClustering.minDTperEvt = 4;
      prof->multiEventClustering.minDTperEvt  = 4;
    }
    try
    {
      prof->singleEventClustering.minNumNeigh =
          configGetInt(prefix + "minNumNeighbours");
      prof->multiEventClustering.minNumNeigh =
          configGetInt(prefix + "minNumNeighbours");
    }
    catch (...)
    {
      prof->singleEventClustering.minNumNeigh = 4;
      prof->multiEventClustering.minNumNeigh  = 4;
    }

    prefix = string("profile.") + prof->name +
             ".doubleDifferenceSystem.phaseFiltering.";

    try
    {
      prof->singleEventClustering.maxDTperEvt =
          configGetInt(prefix + "maxNumPhases");
      prof->multiEventClustering.maxDTperEvt =
          configGetInt(prefix + "maxNumPhases");
    }
    catch (...)
    {
      prof->singleEventClustering.maxDTperEvt = 0;
      prof->multiEventClustering.maxDTperEvt  = 0;
    }

    try
    {
      prof->singleEventClustering.minESdist =
          configGetDouble(prefix + "minStationDistance");
      prof->multiEventClustering.minESdist =
          configGetDouble(prefix + "minStationDistance");
    }
    catch (...)
    {
      prof->singleEventClustering.minESdist = 0;
      prof->multiEventClustering.minESdist  = 0;
    }

    try
    {
      prof->singleEventClustering.maxESdist =
          configGetDouble(prefix + "maxStationDistance");
      prof->multiEventClustering.maxESdist =
          configGetDouble(prefix + "maxStationDistance");
    }
    catch (...)
    {
      prof->singleEventClustering.maxESdist = 0;
      prof->multiEventClustering.maxESdist  = 0;
    }

    try
    {
      prof->singleEventClustering.minEStoIEratio =
          configGetDouble(prefix + "minStationToEventPairDistRatio");
      prof->multiEventClustering.minEStoIEratio =
          configGetDouble(prefix + "minStationToEventPairDistRatio");
    }
    catch (...)
    {
      prof->singleEventClustering.minEStoIEratio = 5;
      prof->multiEventClustering.minEStoIEratio  = 5;
    }

    prefix = string("profile.") + prof->name +
             ".doubleDifferenceSystem.eventPairSelection.singleEvent.";

    try
    {
      prof->singleEventClustering.maxNumNeigh =
          configGetInt(prefix + "maxNumNeighbours");
    }
    catch (...)
    {
      prof->singleEventClustering.maxNumNeigh = 80;
    }

    try
    {
      prof->singleEventClustering.numEllipsoids =
          configGetInt(prefix + "numEllipsoids");
    }
    catch (...)
    {
      prof->singleEventClustering.numEllipsoids = 5;
    }
    try
    {
      prof->singleEventClustering.maxEllipsoidSize =
          configGetDouble(prefix + "maxEllipsoidSize");
    }
    catch (...)
    {
      prof->singleEventClustering.maxEllipsoidSize = 5;
    }

    prefix = string("profile.") + prof->name +
             ".doubleDifferenceSystem.eventPairSelection.multiEvent.";
    try
    {
      prof->multiEventClustering.maxNumNeigh =
          configGetInt(prefix + "maxNumNeighbours");
    }
    catch (...)
    {
      prof->multiEventClustering.maxNumNeigh = 30;
    }

    try
    {
      prof->multiEventClustering.numEllipsoids =
          configGetInt(prefix + "numEllipsoids");
    }
    catch (...)
    {
      prof->multiEventClustering.numEllipsoids = 0;
    }
    try
    {
      prof->multiEventClustering.maxEllipsoidSize =
          configGetDouble(prefix + "maxEllipsoidSize");
    }
    catch (...)
    {
      prof->multiEventClustering.maxEllipsoidSize = 5;
    }
    prefix = string("profile.") + prof->name + ".crossCorrelation.";
    try
    {
      prof->singleEventClustering.xcorrMaxEvStaDist =
          configGetDouble(prefix + "maxStationDistance");
      prof->multiEventClustering.xcorrMaxEvStaDist =
          configGetDouble(prefix + "maxStationDistance");
    }
    catch (...)
    {
      prof->singleEventClustering.xcorrMaxEvStaDist = 0;
      prof->multiEventClustering.xcorrMaxEvStaDist  = 0;
    }
    try
    {
      prof->singleEventClustering.xcorrMaxInterEvDist =
          configGetDouble(prefix + "maxInterEventDistance");
      prof->multiEventClustering.xcorrMaxInterEvDist =
          configGetDouble(prefix + "maxInterEventDistance");
    }
    catch (...)
    {
      prof->singleEventClustering.xcorrMaxInterEvDist = -1;
      prof->multiEventClustering.xcorrMaxInterEvDist  = -1;
    }

    try
    {
      string compatibleChannels =
          configGetString(prefix + "compatibleChannels");
      vector<string> compatibleSets = ::splitString(compatibleChannels, ";");
      for (const auto &cs : compatibleSets)
      {
        vector<string> codes = ::splitString(cs, ",");
        for (size_t i = 0; i < codes.size() - 1; i++)
          for (size_t j = i + 1; j < codes.size(); j++)
          {
            prof->ddCfg.compatibleChannels.push_back({codes[i], codes[j]});
          }
      }
    }
    catch (...)
    {}

    try
    {
      prof->detectMissingPhasesAuto =
          configGetBool(prefix + "detectMissingPhasesAutoOrigin");
    }
    catch (...)
    {
      prof->detectMissingPhasesAuto = true;
    }
    try
    {
      prof->detectMissingPhasesManual =
          configGetBool(prefix + "detectMissingPhasesManualOrigin");
    }
    catch (...)
    {
      prof->detectMissingPhasesManual = false;
    }

    prefix = string("profile.") + prof->name + ".crossCorrelation.p-phase.";
    try
    {
      prof->ddCfg.xcorr[PhaseType::P].startOffset =
          configGetDouble(prefix + "start");
    }
    catch (...)
    {
      prof->ddCfg.xcorr[PhaseType::P].startOffset = -0.50;
    }
    try
    {
      prof->ddCfg.xcorr[PhaseType::P].endOffset =
          configGetDouble(prefix + "end");
    }
    catch (...)
    {
      prof->ddCfg.xcorr[PhaseType::P].endOffset = 0.50;
    }
    try
    {
      prof->ddCfg.xcorr[PhaseType::P].maxDelay =
          configGetDouble(prefix + "maxDelay");
    }
    catch (...)
    {
      prof->ddCfg.xcorr[PhaseType::P].maxDelay = 0.50;
    }
    try
    {
      prof->ddCfg.xcorr[PhaseType::P].minCoef =
          configGetDouble(prefix + "minCCCoef");
    }
    catch (...)
    {
      prof->ddCfg.xcorr[PhaseType::P].minCoef = 0.50;
    }
    try
    {
      prof->ddCfg.xcorr[PhaseType::P].components =
          configGetStrings(prefix + "components");
    }
    catch (...)
    {
      prof->ddCfg.xcorr[PhaseType::P].components = {"Z"};
    }

    prefix = string("profile.") + prof->name + ".crossCorrelation.s-phase.";
    try
    {
      prof->ddCfg.xcorr[PhaseType::S].startOffset =
          configGetDouble(prefix + "start");
    }
    catch (...)
    {
      prof->ddCfg.xcorr[PhaseType::S].startOffset = -0.50;
    }
    try
    {
      prof->ddCfg.xcorr[PhaseType::S].endOffset =
          configGetDouble(prefix + "end");
    }
    catch (...)
    {
      prof->ddCfg.xcorr[PhaseType::S].endOffset = 0.75;
    }
    try
    {
      prof->ddCfg.xcorr[PhaseType::S].maxDelay =
          configGetDouble(prefix + "maxDelay");
    }
    catch (...)
    {
      prof->ddCfg.xcorr[PhaseType::S].maxDelay = 0.50;
    }
    try
    {
      prof->ddCfg.xcorr[PhaseType::S].minCoef =
          configGetDouble(prefix + "minCCCoef");
    }
    catch (...)
    {
      prof->ddCfg.xcorr[PhaseType::S].minCoef = 0.50;
    }
    try
    {
      prof->ddCfg.xcorr[PhaseType::S].components =
          configGetStrings(prefix + "components");
    }
    catch (...)
    {
      prof->ddCfg.xcorr[PhaseType::S].components = {"H"};
    }

    prefix = string("profile.") + prof->name +
             ".crossCorrelation.waveformFiltering.";
    try
    {
      prof->ddCfg.wfFilter.filterStr = configGetString(prefix + "filterString");
    }
    catch (...)
    {
      prof->ddCfg.wfFilter.filterStr = "ITAPER(1)>>BW_HLP(2,1,20)";
    }
    try
    {
      prof->ddCfg.wfFilter.extraTraceLen = configGetDouble(prefix + "margin");
    }
    catch (...)
    {
      prof->ddCfg.wfFilter.extraTraceLen = 1.0;
    }
    try
    {
      prof->ddCfg.wfFilter.resampleFreq =
          configGetDouble(prefix + "resampling");
    }
    catch (...)
    {
      prof->ddCfg.wfFilter.resampleFreq = 0;
    }

    prefix = string("profile.") + prof->name + ".crossCorrelation.snr.";
    try
    {
      prof->ddCfg.snr.minSnr = configGetDouble(prefix + "minSnr");
    }
    catch (...)
    {
      prof->ddCfg.snr.minSnr = 2.;
    }
    try
    {
      prof->ddCfg.snr.noiseStart = configGetDouble(prefix + "noiseStart");
    }
    catch (...)
    {
      prof->ddCfg.snr.noiseStart = -3.0;
    }
    try
    {
      prof->ddCfg.snr.noiseEnd = configGetDouble(prefix + "noiseEnd");
    }
    catch (...)
    {
      prof->ddCfg.snr.noiseEnd = -0.350;
    }
    try
    {
      prof->ddCfg.snr.signalStart = configGetDouble(prefix + "signalStart");
    }
    catch (...)
    {
      prof->ddCfg.snr.signalStart = -0.350;
    }
    try
    {
      prof->ddCfg.snr.signalEnd = configGetDouble(prefix + "signalEnd");
    }
    catch (...)
    {
      prof->ddCfg.snr.signalEnd = 1;
    }

    prefix = string("profile.") + prof->name + ".solver.";
    try
    {
      prof->tttType = configGetString(prefix + "travelTimeTable.tableType");
    }
    catch (...)
    {
      prof->tttType = "libtau";
    }
    try
    {
      prof->tttModel = configGetString(prefix + "travelTimeTable.tableModel");
    }
    catch (...)
    {
      prof->tttModel = "iasp91";
    }

    try
    {
      prof->solverCfg.airQuakes.elevationThreshold =
          configGetDouble(prefix + "airQuakes.elevationThreshold");
    }
    catch (...)
    {
      prof->solverCfg.airQuakes.elevationThreshold = 0;
    }

    try
    {
      string action = makeUpper(configGetString(prefix + "airQuakes.action"));
      if (action == "RESET")
        prof->solverCfg.airQuakes.action = AQ_ACTION::RESET;
      else if (action == "RESET_DEPTH")
        prof->solverCfg.airQuakes.action = AQ_ACTION::RESET_DEPTH;
      else
        prof->solverCfg.airQuakes.action = AQ_ACTION::NONE;
    }
    catch (...)
    {
      prof->solverCfg.airQuakes.action = AQ_ACTION::NONE;
    }

    try
    {
      prof->solverCfg.type = configGetString(prefix + "solverType");
    }
    catch (...)
    {
      prof->solverCfg.type = "LSMR";
    }
    try
    {
      prof->solverCfg.algoIterations = configGetInt(prefix + "algoIterations");
    }
    catch (...)
    {
      prof->solverCfg.algoIterations = 20;
    }
    try
    {
      prof->solverCfg.absLocConstraintStart =
          configGetDouble(prefix + "absoluteLocationConstraint.startingValue");
    }
    catch (...)
    {
      prof->solverCfg.absLocConstraintStart = 0;
    }
    try
    {
      prof->solverCfg.absLocConstraintEnd =
          configGetDouble(prefix + "absoluteLocationConstraint.finalValue");
    }
    catch (...)
    {
      prof->solverCfg.absLocConstraintEnd = 0;
    }
    try
    {
      prof->solverCfg.dampingFactorStart =
          configGetDouble(prefix + "dampingFactor.startingValue");
    }
    catch (...)
    {
      prof->solverCfg.dampingFactorStart = 0.3;
    }
    try
    {
      prof->solverCfg.dampingFactorEnd =
          configGetDouble(prefix + "dampingFactor.finalValue");
    }
    catch (...)
    {
      prof->solverCfg.dampingFactorEnd = 0.3;
    }

    try
    {
      prof->solverCfg.downWeightingByResidualStart =
          configGetDouble(prefix + "downWeightingByResidual.startingValue");
    }
    catch (...)
    {
      prof->solverCfg.downWeightingByResidualStart = 10.;
    }
    try
    {
      prof->solverCfg.downWeightingByResidualEnd =
          configGetDouble(prefix + "downWeightingByResidual.finalValue");
    }
    catch (...)
    {
      prof->solverCfg.downWeightingByResidualEnd = 3.;
    }
    try
    {
      prof->solverCfg.usePickUncertainty =
          configGetBool(prefix + "aPrioriWeights.usePickUncertainties");
    }
    catch (...)
    {
      prof->solverCfg.usePickUncertainty = false;
    }
    try
    {
      prof->solverCfg.absTTDiffObsWeight =
          configGetDouble(prefix + "aPrioriWeights.absoluteTTObsWeight");
    }
    catch (...)
    {
      prof->solverCfg.absTTDiffObsWeight = 1.0;
    }
    try
    {
      prof->solverCfg.xcorrObsWeight =
          configGetDouble(prefix + "aPrioriWeights.xcorrObsWeight");
    }
    catch (...)
    {
      prof->solverCfg.xcorrObsWeight = 1.0;
    }

    prof->ddCfg.diskTraceMinLen =
        configGetDouble("performance.cachedWaveformLength");

    // no reason to make those configurable
    prof->recordStreamURL                 = recordStreamURL();
    prof->singleEventClustering.minWeight = 0;
    prof->multiEventClustering.minWeight  = 0;
    prof->solverCfg.L2normalization       = true;
    prof->solverCfg.solverIterations      = 0;

    _profiles.push_back(prof);
  }

  // disable the database if the inventory is provided by XML file
  // and certain command line options are present too
  if (!isInventoryDatabaseEnabled() && !profileRequireDB &&
      (!_config.eventXML.empty() || !_config.mergeCatalogs.empty() ||
       !_config.evalXCorr.empty() || !_config.relocateCatalog.empty() ||
       !_config.dumpWaveforms.empty() || _config.loadProfileWf ||
       !_config.dumpClusters.empty()))
  {
    SEISCOMP_INFO("Disable database connection");
    setDatabaseEnabled(false, false);
  }

  if (!profilesOK) return false;

  return true;
}

bool RTDD::init()
{

  if (!Application::init()) return false;

  _config.workingDirectory =
      boost::filesystem::path(_config.workingDirectory).string();
  if (!Util::pathExists(_config.workingDirectory))
  {
    if (!Util::createPath(_config.workingDirectory))
    {
      SEISCOMP_ERROR("workingDirectory: failed to create path %s",
                     _config.workingDirectory.c_str());
      return false;
    }
  }

  _inputEvts  = addInputObjectLog("event");
  _inputOrgs  = addInputObjectLog("origin");
  _outputOrgs = addOutputObjectLog("origin", primaryMessagingGroup());

  _cache.setTimeSpan(Core::TimeSpan(_config.fExpiry * 3600.));
  _cache.setDatabaseArchive(query());

  // Enable periodic timer: handleTimeout()
  enableTimer(1);

  // Check each N seconds if a new job needs to be started
  _cronCounter = _config.wakeupInterval;

  HDD::SCAdapter::initLogger();

  return true;
}

bool RTDD::run()
{
  // Send a profile reload request
  if (!_config.reloadProfileMsg.empty())
  {
    RTDDReloadProfileRequestMessage msg;
    msg.setProfile(_config.reloadProfileMsg);
    if (!connection()->send("SERVICE_REQUEST", &msg))
    {
      SEISCOMP_ERROR("Error while sending relocation request message: %s",
                     connection()->lastError().toString());
      return false;
    }
    SEISCOMP_INFO("Relocation request message successfully sent");
    return true;
  }

  // if xml file is provided the load it into _eventParameters
  // this is used by several options
  if (!_config.eventXML.empty())
  {
    IO::XMLArchive ar;
    if (!ar.open(_config.eventXML.c_str()))
    {
      SEISCOMP_ERROR("Unable to open %s", _config.eventXML.c_str());
      return false;
    }
    ar >> _eventParameters;
    ar.close();

    if (!_eventParameters)
    {
      SEISCOMP_ERROR("Event parameters is empty (%s)",
                     _config.eventXML.c_str());
      return false;
    }
  }

  // evaluate cross-correlation and exit
  if (!_config.evalXCorr.empty())
  {
    unique_ptr<HDD::Catalog> catalog = getCatalog(_config.evalXCorr);
    ProfilePtr profile               = getProfile(_config.forceProfile);
    if (!catalog || !profile) return false;
    loadProfile(profile, catalog.get());
    profile->evalXCorr(_config.xcorrCache);
    return true;
  }

  // dump clusters and exit
  if (!_config.dumpClusters.empty())
  {
    unique_ptr<HDD::Catalog> catalog = getCatalog(_config.dumpClusters);
    ProfilePtr profile               = getProfile(_config.forceProfile);
    if (!catalog || !profile) return false;
    loadProfile(profile, catalog.get());
    profile->dumpClusters();
    return true;
  }

  // load catalog waveforms and exit
  if (_config.loadProfileWf)
  {
    ProfilePtr profile = getProfile(_config.forceProfile);
    if (!profile) return false;
    loadProfile(profile);
    profile->preloadWaveforms();
    return true;
  }

  // dump waveforms and exit
  if (!_config.dumpWaveforms.empty())
  {
    unique_ptr<HDD::Catalog> catalog = getCatalog(_config.dumpWaveforms);
    ProfilePtr profile               = getProfile(_config.forceProfile);
    if (!catalog || !profile) return false;
    loadProfile(profile, catalog.get());
    profile->dumpWaveforms();
    return true;
  }

  // dump catalog and exit
  if (!_config.dumpCatalog.empty())
  {
    HDD::Catalog cat;
    DataSource dataSrc(query(), &_cache, _eventParameters.get());
    if (!_config.dumpCatalogOptions.empty())
    {
      vector<DataModel::OriginPtr> origins =
          fetchOrigins(_config.dumpCatalog, _config.dumpCatalogOptions);
      addToCatalog(cat, origins, dataSrc);
    }
    else
    {
      addToCatalog(cat, _config.dumpCatalog, dataSrc);
    }
    cat.writeToFile("event.csv", "phase.csv", "station.csv");
    SEISCOMP_INFO("Wrote files event.csv, phase.csv, station.csv");
    return true;
  }

  // merge catalogs and exit
  if (!_config.mergeCatalogs.empty())
  {
    bool keepEvId = commandline().hasOption("merge-catalogs-keepid");

    std::vector<std::string> tokens = ::splitString(_config.mergeCatalogs, ",");
    if ((tokens.size() % 3) != 0)
    {
      SEISCOMP_ERROR("--merge-catalogs accepts catalog event triplets only");
      return false;
    }

    HDD::Catalog outCat;
    for (size_t i = 0; i < tokens.size(); i += 3)
    {
      SEISCOMP_INFO("Reading and merging %s, %s, %s", tokens[i + 0].c_str(),
                    tokens[i + 1].c_str(), tokens[i + 2].c_str());
      HDD::Catalog cat(tokens[i + 0], tokens[i + 1], tokens[i + 2], true);
      outCat.add(cat, keepEvId);
    }
    outCat.writeToFile("merged-event.csv", "merged-phase.csv",
                       "merged-station.csv");
    SEISCOMP_INFO(
        "Wrote files merged-event.csv, merged-phase.csv, merged-station.csv");
    return true;
  }

  // relocate full catalog and exit
  if (!_config.relocateCatalog.empty())
  {
    std::unordered_map<unsigned, DataModel::OriginPtr> idmap;
    unique_ptr<HDD::Catalog> catalog =
        getCatalog(_config.relocateCatalog, &idmap);
    ProfilePtr profile = getProfile(_config.forceProfile);
    if (!catalog || !profile) return false;

    // if the input catalog is a list of origin ids, then dump its data
    if (!idmap.empty())
    {
      catalog->writeToFile("event.csv", "phase.csv", "station.csv");
      SEISCOMP_INFO(
          "Wrote input catalog files event.csv, phase.csv, station.csv");
    }

    loadProfile(profile, catalog.get());

    unique_ptr<HDD::Catalog> relocatedCat;
    try
    {
      relocatedCat = profile->relocateCatalog(_config.xcorrCache);
    }
    catch (exception &e)
    {
      SEISCOMP_ERROR("Cannot relocate profile catalog: %s", e.what());
      return false;
    }

    // output relocation to xml
    if (commandline().hasOption("xmlout"))
    {
      SEISCOMP_INFO("Converting relocated catalog to XML...");
      DataSource dataSrc(query(), &_cache, _eventParameters.get());
      DataModel::EventParametersPtr evParam = new DataModel::EventParameters();
      evParam->SetRegistrationEnabled(false); // allow existing publicIDs
      for (const auto &kv : relocatedCat->getEvents())
      {
        unique_ptr<HDD::Catalog> ev =
            relocatedCat->extractEvent(kv.second.id, true);
        DataModel::OriginPtr srcOrg;
        try
        {
          srcOrg = idmap.at(kv.second.id); // might not be present
        }
        catch (...)
        {};

        DataModel::OriginPtr newOrg;
        std::vector<DataModel::PickPtr> newOrgPicks;
        bool includeMagnitude = commandline().hasOption("inherit-mag");
        convertOrigin(dataSrc, *ev, srcOrg.get(), author(), agencyID(),
                      profile->methodID, profile->earthModelID,
                      includeMagnitude, true, newOrg, newOrgPicks);

        evParam->add(newOrg.get());
        for (DataModel::PickPtr p : newOrgPicks) evParam->add(p.get());
      }

      IO::XMLArchive ar;
      ar.create("-");
      ar.setFormattedOutput(true);
      ar << evParam;
      ar.close();
    }
    return true;
  }

  // relocate passed origin and exit
  if (!_config.originIDs.empty())
  {
    if (commandline().hasOption("xmlout") && !_eventParameters)
    {
      // no XML input, only XML output
      _eventParameters = new DataModel::EventParameters();
      _config.testMode = true; // we won't send any message
    }

    _config.forceProcessing  = true; // force process of any origin
    _config.profileTimeAlive = 3600; // do not preload profile

    // split multiple origins
    std::vector<std::string> ids = ::splitString(_config.originIDs, ",");
    for (const string &originID : ids)
    {
      OriginPtr org = _cache.get<Origin>(originID);
      if (!org)
      {
        SEISCOMP_ERROR("Origin %s  not found.", originID.c_str());
        continue;
      }

      // Start processing immediately
      _config.delayTimes = {0};
      _cronCounter       = 0;
      addProcess(org.get());
    }

    // output relocation to xml (--ep option provided)
    if (_eventParameters)
    {
      IO::XMLArchive ar;
      ar.create("-");
      ar.setFormattedOutput(true);
      ar << _eventParameters;
      ar.close();
    }

    return true;
  }

  // relocate all origins in xml file and exit
  if (!_config.eventXML.empty())
  {
    _config.forceProcessing  = true; // force process of any origin
    _config.profileTimeAlive = 3600; // do not preload profile

    vector<OriginPtr> origins;
    for (unsigned i = 0; i < _eventParameters->originCount(); i++)
      origins.push_back(_eventParameters->origin(i));

    for (const OriginPtr &org : origins)
    {
      // Start processing immediately
      _config.delayTimes = {0};
      _cronCounter       = 0;

      if (!addProcess(org.get())) return false;
    }

    IO::XMLArchive ar;
    ar.create("-");
    ar.setFormattedOutput(true);
    ar << _eventParameters;
    ar.close();
    return true;
  }

  //
  // real time processing (no other command line options)
  //
  return Application::run();
}

void RTDD::done()
{
  Application::done();

  // Remove crontab log file if exists
  unlink((Environment::Instance()->logDir() + "/" + name() + ".sched").c_str());
}

void RTDD::handleMessage(Core::Message *msg)
{
  Application::handleMessage(msg);

  // Add all events collected by addObject/updateObject
  Todos::iterator it;
  for (it = _todos.begin(); it != _todos.end(); ++it) addProcess(it->get());
  _todos.clear();

  // Reload profile request
  RTDDReloadProfileRequestMessage *reload_req =
      RTDDReloadProfileRequestMessage::Cast(msg);
  if (reload_req)
  {
    SEISCOMP_INFO("Received profile reload request (profile %s)",
                  reload_req->getProfile().c_str());
    RTDDReloadProfileResponseMessage resp;
    ProfilePtr profile = getProfile(reload_req->getProfile());
    if (profile)
    {
      profile->unload();
    }
    else
    {
      SEISCOMP_ERROR("Unknown profile '%s'", reload_req->getProfile().c_str());
      resp.setError(
          stringify("Unknown profile '%s'", reload_req->getProfile().c_str()));
    }
    if (!connection()->send("SERVICE_PROVIDE", &resp))
      SEISCOMP_ERROR("Failed sending profile reload response");
  }
}

void RTDD::addObject(const string &parentID, DataModel::Object *object)
{
  updateObject(parentID, object);
}

void RTDD::updateObject(const string &parentID, Object *object)
{
  Pick *pick = Pick::Cast(object);
  if (pick)
  {
    _cache.feed(pick);
    return;
  }

  Origin *origin = Origin::Cast(object);
  if (origin)
  {
    _cache.feed(origin);
    if (!_config.onlyPreferredOrigin)
    {
      _todos.insert(origin);
      logObject(_inputOrgs, Core::Time::GMT());
    }
    return;
  }

  Event *event = Event::Cast(object);
  if (event)
  {
    _cache.feed(event);
    if (_config.onlyPreferredOrigin)
    {
      _todos.insert(event);
      logObject(_inputEvts, now);
    }
    return;
  }
}

void RTDD::handleTimeout()
{
  checkProfileStatus();
  runNewJobs();
}

/*
 * Periodically clean up profiles unused for some time as they
 * might use lots of memory (waveform data)
 * OR, if the profiles are configured to near expire, make sure
 * they are loaded
 */
void RTDD::checkProfileStatus()
{
  for (ProfilePtr currProfile : _profiles)
  {
    if (!currProfile->isLoaded())
    {
      loadProfile(currProfile);
      if (_config.profileTimeAlive < 0) currProfile->preloadWaveforms();
    }

    // periodic clean up of profiles
    if (_config.profileTimeAlive >= 0)
    {
      Core::TimeSpan expired = Core::TimeSpan(_config.profileTimeAlive, 0);
      if (currProfile->isLoaded() && currProfile->inactiveTime() > expired)
      {
        SEISCOMP_INFO(
            "Profile %s inactive for more than %f seconds: free resources",
            currProfile->name.c_str(), expired.length());
        currProfile->freeResources();
      }
    }
  }
}

void RTDD::runNewJobs()
{
  if (--_cronCounter <= 0)
  {
    // Reset counter
    _cronCounter = _config.wakeupInterval;

    Processes::iterator it;
    now = Core::Time::GMT();

    std::list<ProcessPtr> procToBeRemoved;
    // Update crontab
    for (it = _processes.begin(); it != _processes.end(); it++)
    {
      ProcessPtr proc = it->second;
      CronjobPtr job  = proc->cronjob;

      // Skip processes where nextRun is not set
      if (job->runTimes.empty())
      {
        SEISCOMP_DEBUG("Process %s expired, removing it",
                       proc->obj->publicID().c_str());
        procToBeRemoved.push_back(proc);
        continue;
      }

      Core::Time nextRun = job->runTimes.front();

      // Time of next run in future?
      if (nextRun > now) continue;

      // Remove all times in the past
      while (!job->runTimes.empty() && (job->runTimes.front() <= now))
        job->runTimes.pop_front();

      // Add eventID to processQueue if not already inserted
      if (find(_processQueue.begin(), _processQueue.end(), proc) ==
          _processQueue.end())
      {
        SEISCOMP_DEBUG("Pushing %s to process queue",
                       proc->obj->publicID().c_str());
        _processQueue.push_back(proc);
      }
    }

    for (ProcessPtr &proc : procToBeRemoved) removeProcess(proc.get());

    //
    // Process event queue, but one event only! The next ones will be handled
    // in the next call. This is to avoid being stuck in a long loop where
    // we might miss updateObject/addObject
    //
    if (!_processQueue.empty())
    {
      ProcessPtr proc = _processQueue.front();
      _processQueue.pop_front();
      if (!startProcess(proc.get()))
      {
        SEISCOMP_DEBUG("It is not possible to run job %s: remove it",
                       proc->obj->publicID().c_str());
        // nothing more to do, remove process
        removeProcess(proc.get());
      }
      proc->runCount++;
    }

    // Dump crontab if activated
    if (_config.logCrontab)
    {
      ofstream of((Environment::Instance()->logDir() + "/" + name() + ".sched")
                      .c_str());
      of << "Now: " << now.toString("%F %T") << endl;
      of << "------------------------" << endl;
      of << "[Schedule]" << endl;
      for (it = _processes.begin(); it != _processes.end(); it++)
      {
        ProcessPtr proc    = it->second;
        CronjobPtr cronjob = proc->cronjob;
        if (!cronjob->runTimes.empty())
          of << cronjob->runTimes.front().toString("%F %T") << "\t" << it->first
             << "\t" << (cronjob->runTimes.front() - now).seconds() << endl;
        else
          of << "STOPPED            \t" << it->first << endl;
      }

      // Dump process queue if not empty
      if (!_processQueue.empty())
      {
        of << endl << "[Queue]" << endl;

        ProcessQueue::iterator it;
        for (it = _processQueue.begin(); it != _processQueue.end(); ++it)
          of << "WAITING            \t" << (*it)->obj->publicID() << endl;
      }
    }
  }
}

bool RTDD::addProcess(DataModel::PublicObject *obj)
{
  if (obj == nullptr) return false;

  now = Core::Time::GMT();

  // New process?
  ProcessPtr proc;
  Processes::iterator pit = _processes.find(obj->publicID());
  if (pit == _processes.end())
  {
    SEISCOMP_DEBUG("Adding process [%s]", obj->publicID().c_str());
    // create process
    proc           = new Process;
    proc->created  = now;
    proc->runCount = 0;
    proc->obj      = obj;
    proc->cronjob  = new Cronjob;
    // add process
    _processes[obj->publicID()] = proc;
  }
  else
  {
    SEISCOMP_DEBUG("Update process [%s]: resetting runTimes",
                   obj->publicID().c_str());
    proc = pit->second;
  }

  // populate cronjob
  proc->cronjob->runTimes.clear();
  for (size_t i = 0; i < _config.delayTimes.size(); ++i)
    proc->cronjob->runTimes.push_back(now +
                                      Core::TimeSpan(_config.delayTimes[i], 0));

  SEISCOMP_DEBUG("Update runTimes for [%s]", proc->obj->publicID().c_str());

  handleTimeout();
  return true;
}

// return false when the process cannot run and should not be retried in the
// future
bool RTDD::startProcess(Process *proc)
{
  SEISCOMP_DEBUG("Starting process [%s]", proc->obj->publicID().c_str());

  bool isPreferred = false;
  OriginPtr org;

  // assume process contain an origin (events are relevant only with
  // _config.onlyPreferredOrigin)
  org = Origin::Cast(proc->obj);

  if (!org) // then this must be an event....
  {
    // ...fetch the preferred origin of the event
    EventPtr evt = Event::Cast(proc->obj);
    if (evt)
    {
      org         = _cache.get<Origin>(evt->preferredOriginID());
      isPreferred = true;
    }
  }
  else
  {
    // is 'org'  a preferred origin ?
    DataModel::Event *parentEv = query()->getEvent(org->publicID());
    isPreferred =
        parentEv && (parentEv->preferredOriginID() == org->publicID());
  }

  // If 'onlyPreferredOrigin' is set then make sure we are processing a
  // preferred origin only
  if (_config.onlyPreferredOrigin && !_config.forceProcessing)
  {
    if (!isPreferred)
    {
      SEISCOMP_INFO("Skipping non-preferred origin [%s]",
                    org->publicID().c_str());
      return false;
    }
  }

  if (!org)
  {
    SEISCOMP_DEBUG("Nothing to do for process [%s]",
                   proc->obj->publicID().c_str());
    return false;
  }

  // Find best earth model based on region information and the initial origin
  ProfilePtr currProfile = getProfile(org.get(), _config.forceProfile);

  if (!currProfile)
  {
    SEISCOMP_DEBUG("No profile available, ignoring origin %s",
                   org->publicID().c_str());
    return false;
  }

  // Relocate origin
  OriginPtr relocatedOrg;
  std::vector<DataModel::PickPtr> relocatedOrgPicks;
  return processOrigin(org.get(), relocatedOrg, relocatedOrgPicks, currProfile,
                       _config.forceProcessing, _config.allowAutomaticOrigin,
                       _config.allowManualOrigin, !_config.testMode);
}

void RTDD::removeProcess(Process *proc)
{
  // Remove process from process map
  Processes::iterator pit = _processes.find(proc->obj->publicID());
  if (pit != _processes.end()) _processes.erase(pit);

  // Remove process from queue
  ProcessQueue::iterator qit =
      find(_processQueue.begin(), _processQueue.end(), proc);
  if (qit != _processQueue.end()) _processQueue.erase(qit);
}

// return false when the process cannot run and should not be retried in the
// future
bool RTDD::processOrigin(Origin *origin,
                         OriginPtr &relocatedOrg,
                         std::vector<DataModel::PickPtr> &relocatedOrgPicks,
                         const ProfilePtr &profile,
                         bool forceProcessing,
                         bool allowAutomaticOrigin,
                         bool allowManualOrigin,
                         bool doSend)
{
  relocatedOrg = nullptr;

  if (!origin) return false;

  SEISCOMP_DEBUG("Process origin %s", origin->publicID().c_str());

  // Skip automatic or manul origins if configured so
  if (!forceProcessing)
  {
    bool isManualOrigin = false;
    try
    {
      isManualOrigin =
          origin->evaluationMode() != Seiscomp::DataModel::AUTOMATIC;
    }
    catch (...)
    {
      // origins without an evaluation mode are treated as automatic
    }
    if (isManualOrigin && !allowManualOrigin)
    {
      SEISCOMP_DEBUG("Ignoring manual origin %s", origin->publicID().c_str());
      return false;
    }
    if (!isManualOrigin && !allowAutomaticOrigin)
    {
      SEISCOMP_DEBUG("Ignoring automatic origin %s",
                     origin->publicID().c_str());
      return false;
    }
  }

  if (startsWith(origin->methodID(), profile->methodID, false) &&
      !forceProcessing)
  {
    SEISCOMP_DEBUG("Origin %s was generated by RTDD, skip it",
                   origin->publicID().c_str());
    return false;
  }

  if (isAgencyIDBlocked(objectAgencyID(origin)) && !forceProcessing)
  {
    SEISCOMP_DEBUG("%s: origin's agencyID '%s' is blocked",
                   origin->publicID().c_str(), objectAgencyID(origin).c_str());
    return false;
  }

  SEISCOMP_INFO("Relocating origin %s using profile %s",
                origin->publicID().c_str(), profile->name.c_str());

  try
  {
    relocateOrigin(origin, profile, relocatedOrg, relocatedOrgPicks);
  }
  catch (exception &e)
  {
    SEISCOMP_ERROR("Cannot relocate origin %s (%s)", origin->publicID().c_str(),
                   e.what());
    return true;
  }

  if (!relocatedOrg)
  {
    SEISCOMP_ERROR("processing of origin '%s' failed",
                   origin->publicID().c_str());
    return true;
  }

  SEISCOMP_INFO("Origin %s has been relocated", origin->publicID().c_str());

  //
  // finished processing, send new origin and update journal
  //

  if (_eventParameters)
  {
    // Insert origin to event parameters
    _eventParameters->add(relocatedOrg.get());
    for (DataModel::PickPtr p : relocatedOrgPicks)
      _eventParameters->add(p.get());
  }

  if (connection())
  {
    bool wasEnabled = Notifier::IsEnabled();
    //
    // send origin
    //
    if (doSend)
    {
      SEISCOMP_INFO("Sending origin %s", relocatedOrg->publicID().c_str());

      logObject(_outputOrgs, Core::Time::GMT());

      EventParametersPtr ep = new EventParameters;
      Notifier::Enable();
      ep->add(relocatedOrg.get());
      for (DataModel::PickPtr p : relocatedOrgPicks) ep->add(p.get());
      Notifier::SetEnabled(wasEnabled);

      NotifierMessagePtr msg = Notifier::GetMessage();
      bool result            = false;
      if (msg && connection()) result = connection()->send(msg.get());
      if (!result)
        SEISCOMP_ERROR("%s: sending of relocated origin failed",
                       relocatedOrg->publicID().c_str());
    }
  }

  return true;
}

void RTDD::removedFromCache(Seiscomp::DataModel::PublicObject *po)
{
  // do nothing
}

void RTDD::relocateOrigin(DataModel::Origin *org,
                          ProfilePtr profile,
                          DataModel::OriginPtr &newOrg,
                          std::vector<DataModel::PickPtr> &newOrgPicks)
{
  if (!profile->isLoaded())
  {
    loadProfile(profile);
    if (_config.profileTimeAlive < 0) profile->preloadWaveforms();
  }
  unique_ptr<HDD::Catalog> relocatedOrg = profile->relocateSingleEvent(org);
  bool includeMagnitude = org->evaluationMode() == DataModel::MANUAL;
  DataSource dataSrc(query(), &_cache, _eventParameters.get());
  convertOrigin(dataSrc, *relocatedOrg, org, author(), agencyID(),
                profile->methodID, profile->earthModelID, includeMagnitude,
                false, newOrg, newOrgPicks);
}

std::unique_ptr<HDD::Catalog>
RTDD::getCatalog(const std::string &catalogPath,
                 std::unordered_map<unsigned, DataModel::OriginPtr> *idmap)
{
  std::vector<std::string> tokens = ::splitString(catalogPath, ",");
  try
  {
    if (tokens.size() == 1) // single file containing origin ids
    {
      SEISCOMP_INFO("Loading catalog from origin id file %s",
                    tokens[0].c_str());
      DataSource dataSrc(query(), &_cache, _eventParameters.get());
      unique_ptr<HDD::Catalog> cat(new HDD::Catalog());
      auto _map = addToCatalog(*cat, tokens[0], dataSrc);
      if (idmap) *idmap = _map;
      return cat;
    }
    else if (tokens.size() == 3) // triplet: station.csv,event.csv,phase.csv
    {
      SEISCOMP_INFO(
          "Loading catalog from station file %s, event file %s, phase file %s",
          tokens[0].c_str(), tokens[1].c_str(), tokens[2].c_str());
      return unique_ptr<HDD::Catalog>(
          new HDD::Catalog(tokens[0], tokens[1], tokens[2], true));
    }
    else
    {
      throw runtime_error(
          "The catalog should be a single file containing origin ids or a file "
          "triplet 'station.csv,event.csv,phase.csv'");
    }
  }
  catch (exception &e)
  {
    SEISCOMP_ERROR("Cannot load catalog %s (%s)", catalogPath.c_str(),
                   e.what());
    return nullptr;
  }
}

RTDD::ProfilePtr RTDD::getProfile(const std::string &profile)
{
  if (profile.empty())
  {
    SEISCOMP_ERROR("No profile has been selected");
    return nullptr;
  }
  for (ProfilePtr p : _profiles)
  {
    if (p->name == profile) return p;
  }
  SEISCOMP_ERROR("Profile %s not found", profile.c_str());
  return nullptr;
}

RTDD::ProfilePtr RTDD::getProfile(const DataModel::Origin *origin,
                                  const std::string &forceProfile)
{
  double latitude;
  double longitude;

  try
  {
    latitude  = origin->latitude().value();
    longitude = origin->longitude().value();
  }
  catch (...)
  {
    SEISCOMP_WARNING("Origin %s doesn't have values for lat/lon",
                     origin->publicID().c_str());
    return nullptr;
  }
  return getProfile(latitude, longitude, forceProfile);
}

RTDD::ProfilePtr RTDD::getProfile(double latitude,
                                  double longitude,
                                  const std::string &forceProfile)
{
  ProfilePtr currProfile;

  for (ProfilePtr p : _profiles)
  {
    if (!forceProfile.empty())
    {
      // if user forced a profile, use that
      if (p->name == forceProfile) currProfile = p;
    }
    else
    {
      // if epicenter is inside the configured region, use it
      if (p->region->isInside(latitude, longitude)) currProfile = p;
    }
    if (currProfile) break;
  }

  return currProfile;
}

void RTDD::loadProfile(ProfilePtr profile,
                       const HDD::Catalog *alternativeCatalog)
{
  profile->load(query(), &_cache, _eventParameters.get(),
                _config.workingDirectory, _config.saveProcessingFiles,
                _config.cacheWaveforms, _config.cacheAllWaveforms,
                alternativeCatalog);
}

std::vector<DataModel::OriginPtr> RTDD::fetchOrigins(const std::string &idFile,
                                                     std::string options)
{
  if (!Util::fileExists(idFile))
  {
    throw runtime_error("File " + idFile + " does not exist");
  }

  std::vector<std::string> tokens = ::splitString(options, ",");
  if ((tokens.size() % 5) != 0)
  {
    throw runtime_error("--dump-catalog-options format is: "
                        "type,evalmode,includeCreator,excludeCreator,profile");
  }
  string type           = tokens[0]; // preferred, last, first
  string evalmode       = tokens[1]; // any, onlyManual, onlyAutomatic
  string includeCreator = tokens[2]; // any or a author/methodID
  string excludeCreator = tokens[3]; // none  or a author/methodID
  string profileName    = tokens[4]; // any or profile name

  bool automaticOnly = (evalmode == "onlyAutomatic");
  bool manualOnly    = (evalmode == "onlyManual");

  ProfilePtr profile;
  if (profileName != "any")
  {
    for (ProfilePtr p : _profiles)
    {
      if (p->name == profileName)
      {
        profile = p;
        break;
      }
    }
  }

  SEISCOMP_INFO("Selecting origins with the following characteristics:");
  if (type == "preferred") SEISCOMP_INFO("* PREFERRED only");
  if (type == "first") SEISCOMP_INFO("* arrived FIRST");
  if (type == "last") SEISCOMP_INFO("* arrived LAST");
  if (automaticOnly) SEISCOMP_INFO("* AUTOMATIC only");
  if (manualOnly) SEISCOMP_INFO("* MANUAL only");
  if (includeCreator != "any")
    SEISCOMP_INFO("* whose author or methodID starts with %s",
                  includeCreator.c_str());
  if (excludeCreator != "none")
    SEISCOMP_INFO("* EXCLUDING origins whose author or methodID starts with %s",
                  excludeCreator.c_str());
  if (profile)
    SEISCOMP_INFO("* only origins within %s profile region",
                  profile->name.c_str());

  // fetch origins
  vector<DataModel::OriginPtr> origins;

  for (const auto &row : HDD::CSV::readWithHeader(idFile))
  {
    const string &id = row.at("origin");

    DataModel::OriginPtr org = _cache.get<DataModel::Origin>(id);

    if (!org)
    {
      DataModel::EventPtr ev = _cache.get<DataModel::Event>(id);
      if (ev)
      {
        vector<DataModel::OriginPtr> eventOrigins;

        if (type == "preferred")
        {
          eventOrigins.push_back(
              _cache.get<DataModel::Origin>(ev->preferredOriginID()));
        }
        else
        {
          query()->loadOriginReferences(ev.get());
          for (size_t i = 0; i < ev->originReferenceCount(); i++)
          {
            DataModel::OriginReference *orgRef = ev->originReference(i);
            eventOrigins.push_back(
                _cache.get<DataModel::Origin>(orgRef->originID()));
          }
        }

        Core::Time creationTime;
        for (DataModel::OriginPtr tmpOrg : eventOrigins)
        {
          if (!tmpOrg) continue;

          if (automaticOnly && tmpOrg->evaluationMode() != DataModel::AUTOMATIC)
            continue;

          if (manualOnly && tmpOrg->evaluationMode() != DataModel::MANUAL)
            continue;

          if (includeCreator != "any" &&
              !startsWith(tmpOrg->methodID(), includeCreator, true) &&
              !startsWith(tmpOrg->creationInfo().author(), includeCreator,
                          true))
            continue;

          if (excludeCreator != "none" &&
              (startsWith(tmpOrg->methodID(), excludeCreator, true) ||
               startsWith(tmpOrg->creationInfo().author(), excludeCreator,
                          true)))
            continue;

          if (type == "last" && org &&
              (tmpOrg->creationInfo().creationTime() < creationTime))
            continue;

          if (type == "first" && org &&
              (tmpOrg->creationInfo().creationTime() > creationTime))
            continue;

          if (profile)
          {
            double latitude, longitude;
            try
            {
              latitude  = tmpOrg->latitude().value();
              longitude = tmpOrg->longitude().value();
            }
            catch (...)
            {
              continue;
            }

            if (!profile->region->isInside(latitude, longitude)) continue;
          }

          org          = tmpOrg;
          creationTime = org->creationInfo().creationTime();
        }
      }
    }

    if (!org)
    {
      SEISCOMP_INFO("Cannot find an origin for event %s", id.c_str());
      continue;
    }

    SEISCOMP_INFO("Found origin %s for event %s", org->publicID().c_str(),
                  id.c_str());
    origins.push_back(org);
  }

  return origins;
}

// Profile class

RTDD::Profile::Profile() { loaded = false; }
RTDD::Profile::~Profile() { unload(); }

void RTDD::Profile::load(DatabaseQuery *query,
                         PublicObjectTimeSpanBuffer *cache,
                         EventParameters *eventParameters,
                         const string &workingDir,
                         bool saveProcessingFiles,
                         bool cacheWaveforms,
                         bool cacheAllWaveforms,
                         const HDD::Catalog *alternativeCatalog)
{
  if (loaded) return;

  string pWorkingDir = (boost::filesystem::path(workingDir) / name).string();

  this->query           = query;
  this->cache           = cache;
  this->eventParameters = eventParameters;

  SEISCOMP_INFO("Loading profile %s", name.c_str());

  try
  {
    // load the catalog
    HDD::Catalog ddbgc;
    if (alternativeCatalog) // force this catalog
    {
      ddbgc = *alternativeCatalog;
    }
    else if (!eventIDFile.empty()) // catalog is a list of origin ids
    {
      DataSource dataSrc(query, cache, eventParameters);
      addToCatalog(ddbgc, eventIDFile, dataSrc);
    }
    else // catalog is extended format station.csv,event.csv,phase.csv
    {
      ddbgc = HDD::Catalog(stationFile, eventFile, phaFile);
    }

    // Load the travel time table
    unique_ptr<HDD::TravelTimeTable> ttt;

    if (tttType == "NonLinLoc")
    {
      std::vector<std::string> tokens(::splitString(tttModel, ";"));
      if (tokens.size() != 3 && tokens.size() != 4)
      {
        string msg = stringify(
            "Error while initialzing NLL grids: invalid table model (%s)",
            tttModel.c_str());
        throw runtime_error(msg.c_str());
      }
      string velGridPath   = tokens.at(0);
      string timeGridPath  = tokens.at(1);
      string angleGridPath = tokens.at(2);
      bool swapBytes       = false;
      if (tokens.size() > 3 && tokens.at(3) == "swapBytes")
      {
        swapBytes = true;
      }
      ttt.reset(new HDD::NLL::TravelTimeTable(velGridPath, timeGridPath,
                                              angleGridPath, swapBytes));
    }
    else if (tttType == "ConstVel")
    {
      std::vector<std::string> tokens(::splitString(tttModel, ";"));
      if (tokens.size() != 2)
      {
        string msg = stringify("Error while initialzing ConstVel travel time "
                               "table. Invalid model (%s)",
                               tttModel.c_str());
        throw runtime_error(msg.c_str());
      }
      double pVel = std::stod(tokens.at(0));
      double sVel = std::stod(tokens.at(1));
      ttt.reset(new HDD::ConstantVelocity(pVel, sVel));
    }
    else
    {
      ttt.reset(new HDD::SCAdapter::TravelTimeTable(tttType, tttModel));
    }

    std::unique_ptr<HDD::Waveform::Proxy> wf(
        new HDD::SCAdapter::WaveformProxy(recordStreamURL));

    dd.reset(new HDD::DD(ddbgc, ddCfg, std::move(ttt), std::move(wf)));

    if (saveProcessingFiles)
      dd->enableSaveProcessing(pWorkingDir);
    else
      dd->disableSaveProcessing();

    if (cacheWaveforms)
      dd->enableCatalogWaveformDiskCache(HDD::joinPath(pWorkingDir, "wfcache"));
    else
      dd->disableCatalogWaveformDiskCache();

    if (cacheAllWaveforms)
      dd->enableAllWaveformDiskCache(HDD::joinPath(pWorkingDir, "tmpcache"));
    else
      dd->disableAllWaveformDiskCache();
  }
  catch (exception &e)
  {
    SEISCOMP_ERROR("Cannot load profile %s (%s)", name.c_str(), e.what());
    unload();
    return;
  }

  loaded    = true;
  lastUsage = Core::Time::GMT();

  SEISCOMP_INFO("Profile %s loaded into memory", name.c_str());
}

void RTDD::Profile::unload()
{
  SEISCOMP_INFO("Unloading profile %s", name.c_str());
  dd.reset();
  loaded    = false;
  lastUsage = Core::Time::GMT();
}

void RTDD::Profile::preloadWaveforms()
{
  if (!loaded)
  {
    string msg = Core::stringify(
        "Cannot preload catalog waveforms, profile %s not initialized",
        name.c_str());
    throw runtime_error(msg.c_str());
  }
  dd->preloadWaveforms();
  lastUsage = Core::Time::GMT();
}

void RTDD::Profile::freeResources()
{
  if (!loaded) return;
  dd->unloadWaveforms();
  lastUsage = Core::Time::GMT();
}

std::unique_ptr<HDD::Catalog>
RTDD::Profile::relocateSingleEvent(DataModel::Origin *org)
{
  if (!loaded)
  {
    string msg = Core::stringify(
        "Cannot relocate origin, profile %s not initialized", name.c_str());
    throw runtime_error(msg.c_str());
  }
  lastUsage = Core::Time::GMT();

  DataSource dataSrc(query, cache, eventParameters);

  // we pass the stations information from the background catalog, to avoid
  // wasting time accessing the inventory again for information we already have
  HDD::Catalog orgToRelocate(
      dd->getCatalog().getStations(), map<unsigned, HDD::Catalog::Event>(),
      unordered_multimap<unsigned, HDD::Catalog::Phase>());
  addToCatalog(orgToRelocate, {org}, dataSrc);

  if (org->evaluationMode() == DataModel::MANUAL)
    singleEventClustering.xcorrDetectMissingPhases =
        this->detectMissingPhasesManual;
  else
    singleEventClustering.xcorrDetectMissingPhases =
        this->detectMissingPhasesAuto;

  unique_ptr<HDD::Catalog> rel = dd->relocateSingleEvent(
      orgToRelocate, singleEventClustering, singleEventClustering, solverCfg);

  return rel;
}

std::unique_ptr<HDD::Catalog>
RTDD::Profile::relocateCatalog(const std::string &xcorrFile)
{
  if (!loaded)
  {
    string msg = Core::stringify(
        "Cannot relocate catalog, profile %s not initialized", name.c_str());
    throw runtime_error(msg.c_str());
  }
  lastUsage = Core::Time::GMT();

  HDD::XCorrCache xcorr;
  if (!xcorrFile.empty())
    xcorr = HDD::readXCorrFromFile(dd->getCatalog(), xcorrFile);

  multiEventClustering.xcorrDetectMissingPhases = false;

  unique_ptr<HDD::Catalog> relocatedCat =
      dd->relocateMultiEvents(multiEventClustering, solverCfg, xcorr);

  relocatedCat->writeToFile("reloc-event.csv", "reloc-phase.csv",
                            "reloc-station.csv");

  HDD::writeXCorrToFile(xcorr, dd->getCatalog(), "xcorr.csv");

  SEISCOMP_INFO(
      "Wrote relocated catalog files reloc-event.csv, reloc-phase.csv, "
      "reloc-station.csv and cross-correlation results xcorr.csv");
  return relocatedCat;
}

void RTDD::Profile::evalXCorr(const std::string &xcorrFile)
{
  if (!loaded)
  {
    string msg = Core::stringify(
        "Cannot evalute cross-correlation settings, profile %s not initialized",
        name.c_str());
    throw runtime_error(msg.c_str());
  }
  lastUsage = Core::Time::GMT();

  HDD::XCorrCache xcorr;
  if (!xcorrFile.empty())
    xcorr = HDD::readXCorrFromFile(dd->getCatalog(), xcorrFile);

  dd->evalXCorr(multiEventClustering, printEvalXcorrStats, xcorr);

  HDD::writeXCorrToFile(xcorr, dd->getCatalog(), "xcorr.csv");

  SEISCOMP_INFO("Wrote cross-correlation results xcorr.csv");
}

void RTDD::Profile::dumpWaveforms()
{
  if (!loaded)
  {
    string msg = Core::stringify(
        "Cannot dump catalog waveforms, profile %s not initialized",
        name.c_str());
    throw runtime_error(msg.c_str());
  }
  lastUsage = Core::Time::GMT();
  dd->dumpWaveforms();
}

void RTDD::Profile::dumpClusters()
{
  if (!loaded)
  {
    string msg = Core::stringify(
        "Cannot dump clusters, profile %s not initialized", name.c_str());
    throw runtime_error(msg.c_str());
  }
  lastUsage                   = Core::Time::GMT();
  list<HDD::Catalog> clusters = dd->findClusters(multiEventClustering);
  SEISCOMP_INFO("Found %zu clusters", clusters.size());
  unsigned clusterId = 1;
  for (const HDD::Catalog &cat : clusters)
  {
    SEISCOMP_INFO("Writing cluster %u (%zu events)", clusterId,
                  cat.getEvents().size());
    string prefix = stringify("cluster-%u", clusterId++);
    cat.writeToFile(prefix + "-event.csv", prefix + "-phase.csv",
                    prefix + "-station.csv");
  }
}

// End Profile class

} // namespace Seiscomp
