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

#include "rtdd.h"
#include "csvreader.h"
#include "rtddmsg.h"

#include <seiscomp3/logging/channel.h>
#include <seiscomp3/logging/filerotator.h>

#include <seiscomp3/core/genericrecord.h>
#include <seiscomp3/core/strings.h>
#include <seiscomp3/core/system.h>

#include <seiscomp3/client/inventory.h>
#include <seiscomp3/io/archive/xmlarchive.h>
#include <seiscomp3/io/records/mseedrecord.h>

#include <seiscomp3/datamodel/event.h>
#include <seiscomp3/datamodel/journalentry.h>
#include <seiscomp3/datamodel/magnitude.h>
#include <seiscomp3/datamodel/origin.h>
#include <seiscomp3/datamodel/parameter.h>
#include <seiscomp3/datamodel/parameterset.h>
#include <seiscomp3/datamodel/pick.h>
#include <seiscomp3/datamodel/stationmagnitude.h>
#include <seiscomp3/datamodel/stationmagnitudecontribution.h>
#include <seiscomp3/datamodel/utils.h>

#include <seiscomp3/math/geo.h>

#include <seiscomp3/utils/files.h>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace Seiscomp::Processing;
using namespace Seiscomp::DataModel;
using Seiscomp::Core::stringify;
using PhaseType      = Seiscomp::HDD::Catalog::Phase::Type;
using PhaseSrc       = Seiscomp::HDD::Catalog::Phase::Source;
using CatalogDataSrc = Seiscomp::HDD::ScCatalog::DataSource;

#define NEW_OPT(var, ...) addOption(&var, __VA_ARGS__)
#define NEW_OPT_CLI(var, ...) addOption(&var, nullptr, __VA_ARGS__)

namespace {

using namespace Seiscomp;
using Seiscomp::Core::fromString;

template <class T>
bool configGetTypedList(const Application *app,
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

  bool init(const Application *app, const string &prefix)
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

  bool init(const Application *app, const string &prefix)
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

Core::Time now; // this is tricky, I don't like it

void makeUpper(string &dest, const string &src)
{
  dest = src;
  for (size_t i = 0; i < src.size(); ++i) dest[i] = toupper(src[i]);
}

bool startsWith(const string &haystack,
                const string &needle,
                bool caseSensitive = true)
{
  string _haystack = haystack;
  string _needle   = needle;
  if (!caseSensitive)
  {
    makeUpper(_haystack, haystack);
    makeUpper(_needle, needle);
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

double normalizeAz(double az)
{
  if (az < 0)
    az += 360.0;
  else if (az >= 360.0)
    az -= 360.0;
  return az;
}

double normalizeLon(double lon)
{
  while (lon < -180.0) lon += 360.0;
  while (lon > 180.0) lon -= 360.0;
  return lon;
}

} // unnamed namespace

namespace Seiscomp {

RTDD::Config::Config()
{
  workingDirectory     = "/tmp/rtdd";
  saveProcessingFiles  = false;
  onlyPreferredOrigin  = true;
  allowAutomaticOrigin = true;
  allowManualOrigin    = true;
  profileTimeAlive     = -1;
  cacheWaveforms       = false;
  cacheAllWaveforms    = false;
  debugWaveforms       = false;

  loadProfileWf   = false;
  forceProcessing = false;
  testMode        = false;
  dumpWaveforms   = false;
  fExpiry         = 1.0;

  wakeupInterval = 1; // sec
  logCrontab     = true;
}

RTDD::RTDD(int argc, char **argv) : Application(argc, argv)
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

  NEW_OPT(_config.saveProcessingFiles, "saveProcessingFiles");
  NEW_OPT(_config.onlyPreferredOrigin, "onlyPreferredOrigins");
  NEW_OPT(_config.allowAutomaticOrigin, "automaticOrigins");
  NEW_OPT(_config.allowManualOrigin, "manualOrigins");
  NEW_OPT(_config.activeProfiles, "activeProfiles");

  NEW_OPT(_config.logCrontab, "cron.logging");
  NEW_OPT(_config.delayTimes, "cron.delayTimes");

  NEW_OPT(_config.profileTimeAlive, "performance.profileTimeAlive");
  NEW_OPT(_config.cacheWaveforms, "performance.cacheWaveforms");

  NEW_OPT_CLI(
      _config.relocateCatalog, "Mode", "reloc-catalog",
      "Relocate the catalog passed as argument in multi-event mode. The "
      "input can be a single file (containing seiscomp origin ids) or a file "
      "triplet (station.csv,event.csv,phase.csv). For events stored "
      "in a XML files add the --ep option. Use in combination with --profile",
      true);
  NEW_OPT_CLI(
      _config.originIDs, "Mode", "origin-id,O",
      "Relocate  the origin (or multiple comma-separated origins) in "
      "signle-event mode and send a message. Each origin will be processed "
      "accordingly to the matching profile region unless the --profile option "
      " is used.",
      true);
  NEW_OPT_CLI(
      _config.eventXML, "Mode", "ep",
      "Event parameters XML file for offline processing of contained origins "
      "(implies --test option). Each contained origin will be processed in "
      "signle-event mode unless --reloc-catalog is provided, which enable "
      "multi-event mode.",
      true);
  NEW_OPT_CLI(
      _config.evalXCorr, "Mode", "eval-xcorr",
      "Compute cross-correlation statistics fon the catalog passed as "
      "argument. The input can be a single file (containing seiscomp origin "
      "ids) or a file triplet (station.csv,event.csv,phase.csv). Use in "
      "combination with --profile",
      true);
  NEW_OPT_CLI(_config.loadProfileWf, "Mode", "load-profile-wf",
              "Load catalog waveforms from the configured recordstream and "
              "save them into the profile working directory. Use in "
              "combination with --profile",
              false, true);
  NEW_OPT_CLI(_config.reloadProfileMsg, "Mode", "send-reload-profile-msg",
              "Send a message to any running scrtdd module requesting to "
              "reload a specific profile passed as argument",
              true);
  NEW_OPT_CLI(_config.dumpCatalog, "Catalog", "dump-catalog",
              "Dump the seiscomp event/origin id file passed as argument into "
              "a catalog file triplet (station.csv,event.csv,phase.csv).",
              true);
  NEW_OPT_CLI(_config.mergeCatalogs, "Catalog", "merge-catalogs",
              "Merge in a single catalog all the catalog file triplets "
              "(station1.csv,event1.csv,phase1.csv,station2.csv,event2.csv,"
              "phase2.csv,...) passed as arguments.",
              true);
  NEW_OPT_CLI(_config.forceProfile, "ModeOptions", "profile",
              "To be used in combination with other options: select the "
              "profile configuration to use",
              true);
  NEW_OPT_CLI(_config.fExpiry, "ModeOptions", "expiry,x",
              "Defines the time span in hours after which objects expire.",
              true);
  NEW_OPT_CLI(_config.cacheAllWaveforms, "ModeOptions", "cache-wf-all",
              "All waveforms will be saved to disk cache, even temporarily "
              "ones. Normally only catalog phase waveforms are cached to disk. "
              "This is useful to speed up debugging/testing when the same "
              "origins are repeatedly processed.",
              false, true);
  NEW_OPT_CLI(_config.dumpWaveforms, "ModeOptions", "debug-wf",
              "Enable saving of processed waveforms into the profile working "
              "directory for inspection.",
              false, true);
  NEW_OPT_CLI(_config.testMode, "ModeOptions", "test",
              "Test mode, no messages are sent when relocating a single event",
              false, true);
}

RTDD::~RTDD() {}

void RTDD::createCommandLineDescription()
{
  Application::createCommandLineDescription();
  commandline().addOption("Mode", "dump-config",
                          "Dump the configuration and exit");
  commandline().addOption("ModeOptions", "xmlout",
                          "Enable XML output when combined with "
                          "--reloc-catalog or --oring-id options");
  commandline().addOption<string>(
      "Catalog", "merge-catalogs-keepid",
      "Similar to the --merge-catalogs option but events keep their ids. If "
      "multiple events share the same id, subsequent events will be discarded.",
      nullptr, false);
  commandline().addOption<string>(
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
      nullptr, false);
}

bool RTDD::validateParameters()
{
  Environment *env = Environment::Instance();

  if (!Application::validateParameters()) return false;

  if (commandline().hasOption("merge-catalogs-keepid"))
    _config.mergeCatalogs =
        commandline().option<string>("merge-catalogs-keepid");

  // disable messaging (offline mode) with certain command line options
  if (!_config.eventXML.empty() || !_config.dumpCatalog.empty() ||
      !_config.mergeCatalogs.empty() || !_config.evalXCorr.empty() ||
      !_config.relocateCatalog.empty() || _config.loadProfileWf ||
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
  if (!_config.relocateCatalog.empty() &&
      ::splitString(_config.relocateCatalog, ",").size() ==
          1) // single file containing origin ids
  {
    profileRequireDB = true;
  }

  _config.workingDirectory =
      env->absolutePath(configGetPath("workingDirectory"));

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
      makeUpper(regionType, configGetString(prefix + "regionType"));
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
      prof->singleEventClustering.minEStoIEratio = 2;
      prof->multiEventClustering.minEStoIEratio  = 2;
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
      prof->singleEventClustering.xcorrMaxEvStaDist = 85;
      prof->multiEventClustering.xcorrMaxEvStaDist  = 85;
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
      prof->useTheoreticalAuto =
          configGetBool(prefix + "theoreticalPhaseAutoOrigin");
    }
    catch (...)
    {
      prof->useTheoreticalAuto = true;
    }
    try
    {
      prof->useTheoreticalManual =
          configGetBool(prefix + "theoreticalPhaseManualOrigin");
    }
    catch (...)
    {
      prof->useTheoreticalManual = false;
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
      prof->ddCfg.xcorr[PhaseType::P].maxDelay = 0.350;
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
      prof->ddCfg.xcorr[PhaseType::S].maxDelay = 0.350;
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
      prof->ddCfg.xcorr[PhaseType::S].components = {"T", "Z"};
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
      prof->ddCfg.wfFilter.resampleFreq =
          configGetDouble(prefix + "resampling");
    }
    catch (...)
    {
      prof->ddCfg.wfFilter.resampleFreq = 400;
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
      prof->ddCfg.snr.signalEnd = 0.350;
    }

    prefix = string("profile.") + prof->name + ".solver.";
    try
    {
      prof->ddCfg.ttt.type =
          configGetString(prefix + "travelTimeTable.tableType");
    }
    catch (...)
    {
      prof->ddCfg.ttt.type = "libtau";
    }
    try
    {
      prof->ddCfg.ttt.model =
          configGetString(prefix + "travelTimeTable.tableModel");
    }
    catch (...)
    {
      prof->ddCfg.ttt.model = "iasp91";
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

    prof->ddCfg.recordStreamURL = recordStreamURL();

    // no reason to make those configurable
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
       _config.loadProfileWf))
  {
    SEISCOMP_INFO("Disable database connection");
    setDatabaseEnabled(false, false);
  }

  if (!profilesOK) return false;

  if (commandline().hasOption("dump-config"))
  {
    for (Options::const_iterator it = options().begin(); it != options().end();
         ++it)
    {
      if ((*it)->cfgName)
        cout << (*it)->cfgName;
      else if ((*it)->cliParam)
        cout << "--" << (*it)->cliParam;
      else
        continue;

      cout << ": ";
      (*it)->printStorage(cout);
      cout << endl;
    }

    return false;
  }

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

  // evaluate cross-correlation settings and exit
  if (!_config.evalXCorr.empty())
  {
    HDD::CatalogPtr catalog = getCatalog(_config.evalXCorr);
    ProfilePtr profile      = getProfile(_config.forceProfile);
    if (!catalog || !profile) return false;
    loadProfile(profile, false, catalog);
    profile->evalXCorr();
    return true;
  }

  // load catalog waveforms and exit
  if (_config.loadProfileWf)
  {
    ProfilePtr profile = getProfile(_config.forceProfile);
    if (!profile) return false;
    loadProfile(profile, true);
    return true;
  }

  // dump catalog and exit
  if (!_config.dumpCatalog.empty())
  {
    HDD::ScCatalogPtr cat(new HDD::ScCatalog());
    CatalogDataSrc dataSrc(query(), &_cache, _eventParameters.get());
    if (commandline().hasOption("dump-catalog-options"))
    {
      string options = commandline().option<string>("dump-catalog-options");
      vector<DataModel::OriginPtr> origins =
          fetchOrigins(_config.dumpCatalog, options);
      cat->add(origins, dataSrc);
    }
    else
    {
      cat->add(_config.dumpCatalog, dataSrc);
    }
    cat->writeToFile("event.csv", "phase.csv", "station.csv");
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

    HDD::CatalogPtr outCat = new HDD::Catalog();
    for (size_t i = 0; i < tokens.size(); i += 3)
    {
      SEISCOMP_INFO("Reading and merging %s, %s, %s", tokens[i + 0].c_str(),
                    tokens[i + 1].c_str(), tokens[i + 2].c_str());
      HDD::CatalogPtr cat =
          new HDD::Catalog(tokens[i + 0], tokens[i + 1], tokens[i + 2], true);
      outCat->add(*cat, keepEvId);
    }
    outCat->writeToFile("merged-event.csv", "merged-phase.csv",
                        "merged-station.csv");
    SEISCOMP_INFO(
        "Wrote files merged-event.csv, merged-phase.csv, merged-station.csv");
    return true;
  }

  // relocate full catalog and exit
  if (!_config.relocateCatalog.empty())
  {
    std::unordered_map<unsigned, DataModel::OriginPtr> idmap;
    HDD::CatalogPtr catalog = getCatalog(_config.relocateCatalog, &idmap);
    ProfilePtr profile      = getProfile(_config.forceProfile);
    if (!catalog || !profile) return false;

    // if the input catalog is a list of origin ids, then dump its data
    if (!idmap.empty())
    {
      catalog->writeToFile("event.csv", "phase.csv", "station.csv");
      SEISCOMP_INFO(
          "Wrote input catalog files event.csv, phase.csv, station.csv");
    }

    bool preloadData = profile->multiEventClustering.xcorrMaxEvStaDist != 0 &&
                       profile->multiEventClustering.xcorrMaxInterEvDist != 0;

    loadProfile(profile, preloadData, catalog);

    HDD::CatalogPtr relocatedCat;
    try
    {
      relocatedCat = profile->relocateCatalog();
      relocatedCat->writeToFile("reloc-event.csv", "reloc-phase.csv",
                                "reloc-station.csv");
      SEISCOMP_INFO("Wrote relocated catalog files reloc-event.csv, "
                    "reloc-phase.csv, reloc-station.csv");
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
      DataModel::EventParametersPtr evParam = new DataModel::EventParameters();
      evParam->SetRegistrationEnabled(false); // allow existing publicIDs
      for (const auto &kv : relocatedCat->getEvents())
      {
        HDD::CatalogPtr ev = relocatedCat->extractEvent(kv.second.id, true);
        DataModel::OriginPtr srcOrg;
        try
        {
          srcOrg = idmap.at(kv.second.id); // might not be present
        }
        catch (...)
        {};

        DataModel::OriginPtr newOrg;
        std::vector<DataModel::PickPtr> newOrgPicks;
        convertOrigin(ev, profile, srcOrg.get(), true, false, true, newOrg,
                      newOrgPicks);

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

  // Relocate origins coming from scolv
  RTDDRelocateRequestMessage *reloc_req = RTDDRelocateRequestMessage::Cast(msg);
  if (reloc_req)
  {
    SEISCOMP_INFO("Received relocation request");

    RTDDRelocateResponseMessage reloc_resp;
    ProfilePtr currProfile;
    OriginPtr originToReloc = reloc_req->getOrigin();

    if (originToReloc)
    {
      currProfile = getProfile(originToReloc.get(), reloc_req->getProfile());
      if (!currProfile)
      {
        reloc_resp.setError(
            stringify("No profile available, ignoring origin %s",
                      originToReloc->publicID().c_str()));
      }
    }
    else
    {
      reloc_resp.setError("No origin to relocate has been received");
    }

    // Inform scolv we are going to relocate this origin or not
    reloc_resp.setRequestAccepted(!reloc_resp.hasError());

    if (!connection()->send("SERVICE_PROVIDE", &reloc_resp))
      SEISCOMP_ERROR("Failed sending relocation response");

    if (!reloc_resp.hasError())
    {
      OriginPtr relocatedOrg;
      std::vector<DataModel::PickPtr>
          relocatedOrgPicks; // we cannot return these to scolv
      processOrigin(originToReloc.get(), relocatedOrg, relocatedOrgPicks,
                    currProfile, true, true, true, false);

      if (relocatedOrg)
      {
        reloc_resp.setOrigin(relocatedOrg);
      }
      else
      {
        reloc_resp.setError(stringify("OriginId %s has not been relocated",
                                      originToReloc->publicID().c_str()));
      }

      SEISCOMP_DEBUG("Sending relocation response (%s)",
                     (reloc_resp.hasError() ? reloc_resp.getError()
                                            : "no relocation errors")
                         .c_str());

      if (!connection()->send("SERVICE_PROVIDE", &reloc_resp))
        SEISCOMP_ERROR("Failed sending relocation response");
    }
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
      loadProfile(currProfile, (_config.profileTimeAlive < 0));

    // periodic clean up of profiles
    if (_config.profileTimeAlive >= 0)
    {
      Core::TimeSpan expired = Core::TimeSpan(_config.profileTimeAlive);
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
                                      Core::TimeSpan(_config.delayTimes[i]));

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
    loadProfile(profile, (_config.profileTimeAlive < 0));
  HDD::CatalogPtr relocatedOrg = profile->relocateSingleEvent(org);
  bool includeMagnitude        = org->evaluationMode() == DataModel::MANUAL;
  convertOrigin(relocatedOrg, profile, org, includeMagnitude, true, false,
                newOrg, newOrgPicks);
}

void RTDD::convertOrigin(
    const HDD::CatalogCPtr &relocatedOrg,
    ProfilePtr profile,     // can be nullptr
    DataModel::Origin *org, // can be nullptr
    bool includeMagnitude,
    bool fullMagnitude,
    bool includeExistingPicks,
    DataModel::OriginPtr &newOrg,                 // return value
    std::vector<DataModel::PickPtr> &newOrgPicks) // return value
{
  CatalogDataSrc dataSrc(query(), &_cache, _eventParameters.get());

  // there must be only one event in the catalog, the relocated origin
  const HDD::Catalog::Event &event = relocatedOrg->getEvents().begin()->second;

  newOrg = Origin::Create();

  DataModel::CreationInfo ci;
  ci.setAgencyID(agencyID());
  ci.setAuthor(author());
  ci.setCreationTime(Core::Time::GMT());

  newOrg->setCreationInfo(ci);
  newOrg->setEarthModelID(profile ? profile->earthModelID : "");
  newOrg->setMethodID(profile ? profile->methodID : "RTDD");
  newOrg->setEvaluationMode(EvaluationMode(AUTOMATIC));

  newOrg->setTime(DataModel::TimeQuantity(event.time));

  RealQuantity latitude = DataModel::RealQuantity(event.latitude);
  newOrg->setLatitude(event.latitude);

  RealQuantity longitude =
      DataModel::RealQuantity(normalizeLon(event.longitude));
  newOrg->setLongitude(longitude);

  RealQuantity depth = DataModel::RealQuantity(event.depth);
  newOrg->setDepth(depth);

  if (event.relocInfo.isRelocated)
  {
    DataModel::Comment *comment = new DataModel::Comment();
    comment->setId("scrtddRelocationReport");
    comment->setText(HDD::HypoDD::relocationReport(relocatedOrg));
    newOrg->add(comment);
  }

  auto evPhases = relocatedOrg->getPhases().equal_range(
      event.id); // phases of relocated event
  int usedPhaseCount = 0;
  double meanDist    = 0;
  double minDist     = std::numeric_limits<double>::max();
  double maxDist     = 0;
  vector<double> azi;
  set<string> associatedStations;
  set<string> usedStations;

  // If we know the origin before relocation fetch some information from it
  if (org)
  {
    //
    // store source origin id as comment
    //
    DataModel::Comment *comment = new DataModel::Comment();
    comment->setId("scrtddSourceOrigin");
    comment->setText(org->publicID());
    newOrg->add(comment);

    //
    // Copy magnitude from org if that is Manual
    //
    if (includeMagnitude)
    {
      dataSrc.loadMagnitudes(org, fullMagnitude, fullMagnitude);

      unordered_map<string, string> staMagIdMap;
      if (fullMagnitude)
      {
        for (size_t i = 0; i < org->stationMagnitudeCount(); i++)
        {
          DataModel::StationMagnitude *staMag = org->stationMagnitude(i);
          DataModel::StationMagnitude *newStaMag =
              DataModel::StationMagnitude::Create();
          *newStaMag = *staMag;
          newOrg->add(newStaMag);
          staMagIdMap[staMag->publicID()] = newStaMag->publicID();
        }
      }

      for (size_t i = 0; i < org->magnitudeCount(); i++)
      {
        DataModel::Magnitude *mag    = org->magnitude(i);
        DataModel::Magnitude *newMag = DataModel::Magnitude::Create();
        *newMag                      = *mag;

        if (fullMagnitude)
        {
          for (size_t j = 0; j < mag->stationMagnitudeContributionCount(); j++)
          {
            DataModel::StationMagnitudeContribution *contrib =
                mag->stationMagnitudeContribution(j);
            DataModel::StationMagnitudeContribution *newContrib =
                new DataModel::StationMagnitudeContribution();
            *newContrib = *contrib;
            try
            {
              newContrib->setStationMagnitudeID(
                  staMagIdMap.at(contrib->stationMagnitudeID()));
            }
            catch (...)
            {}
            newMag->add(newContrib);
          }
        }
        newOrg->add(newMag);
      }
    }

    //
    // add all arrivals that were in the original Origin (before relocation)
    //
    for (size_t i = 0; i < org->arrivalCount(); i++)
    {
      DataModel::Arrival *orgArr = org->arrival(i);

      DataModel::Arrival *newArr = new Arrival();
      newArr->setPickID(orgArr->pickID());
      newArr->setPhase(orgArr->phase());
      newArr->setWeight(0.);
      newArr->setTimeUsed(false);

      newOrg->add(newArr);

      DataModel::PickPtr pick = dataSrc.get<DataModel::Pick>(orgArr->pickID());
      if (pick)
      {
        associatedStations.insert(pick->waveformID().networkCode() + "." +
                                  pick->waveformID().stationCode());

        if (includeExistingPicks)
          newOrgPicks.push_back(DataModel::Pick::Cast(pick->clone()));
      }
    }
  }

  // add missing arrivals and fill in all the properties
  for (auto it = evPhases.first; it != evPhases.second; ++it)
  {
    const HDD::Catalog::Phase &phase = it->second;
    bool phaseUsed =
        phase.relocInfo.isRelocated && phase.relocInfo.finalWeight != 0;

    // drop phases discovered via cross-correlation if those phases were not
    // used for the relocations
    if ((phase.procInfo.source == PhaseSrc::THEORETICAL ||
         phase.procInfo.source == PhaseSrc::XCORR) &&
        !phaseUsed)
    {
      continue;
    }

    associatedStations.insert(phase.networkCode + "." + phase.stationCode);

    // check if this phase has been already added
    bool alreadyAdded = false;
    DataModel::Arrival *newArr;

    for (size_t i = 0; i < newOrg->arrivalCount(); i++)
    {
      newArr                  = newOrg->arrival(i);
      DataModel::PickPtr pick = dataSrc.get<DataModel::Pick>(newArr->pickID());

      if (pick && phase.time == pick->time().value() &&
          phase.networkCode == pick->waveformID().networkCode() &&
          phase.stationCode == pick->waveformID().stationCode() &&
          phase.locationCode == pick->waveformID().locationCode() &&
          phase.channelCode == pick->waveformID().channelCode())
      {
        alreadyAdded = true;
        break;
      }
    }

    if (!alreadyAdded)
    {
      // prepare the new pick
      DataModel::PickPtr newPick = Pick::Create();
      newPick->setCreationInfo(ci);
      newPick->setMethodID(profile ? profile->methodID : "RTDD");
      newPick->setEvaluationMode(phase.isManual ? EvaluationMode(MANUAL)
                                                : EvaluationMode(AUTOMATIC));
      DataModel::TimeQuantity pickTime(phase.time);
      pickTime.setLowerUncertainty(phase.lowerUncertainty);
      pickTime.setUpperUncertainty(phase.upperUncertainty);
      newPick->setTime(pickTime);
      newPick->setPhaseHint(DataModel::Phase(phase.type));
      newPick->setWaveformID(
          WaveformStreamID(phase.networkCode, phase.stationCode,
                           phase.locationCode, phase.channelCode, ""));
      newOrgPicks.push_back(newPick);

      // prepare the new arrival
      newArr = new Arrival();
      newArr->setPickID(newPick->publicID());
      newArr->setPhase(phase.type);

      newOrg->add(newArr);
    }

    newArr->setWeight(phase.relocInfo.isRelocated ? phase.relocInfo.finalWeight
                                                  : 0.);
    newArr->setTimeUsed(phaseUsed);
    newArr->setTimeResidual(
        phase.relocInfo.isRelocated ? phase.relocInfo.finalTTResidual : 0.);

    auto search = relocatedOrg->getStations().find(phase.stationId);
    if (search == relocatedOrg->getStations().end())
    {
      SEISCOMP_WARNING("Cannot find station id '%s' referenced by phase '%s'."
                       "Cannot add Arrival to relocated origin",
                       phase.stationId.c_str(), string(phase).c_str());
      continue;
    }
    const HDD::Catalog::Station &station = search->second;

    double distance, az, baz;
    Math::Geo::delazi(event.latitude, event.longitude, station.latitude,
                      station.longitude, &distance, &az, &baz);
    newArr->setAzimuth(normalizeAz(az));
    newArr->setDistance(distance);

    // update stats
    if (newArr->timeUsed())
    {
      usedPhaseCount++;
      meanDist += distance;
      minDist = distance < minDist ? distance : minDist;
      maxDist = distance > maxDist ? distance : maxDist;
      azi.push_back(az);
      usedStations.insert(phase.stationId);
    }
  }

  // finish computing stats
  meanDist /= usedPhaseCount;

  double primaryAz = 360., secondaryAz = 360.;
  if (azi.size() >= 2)
  {
    primaryAz = secondaryAz = 0.;
    sort(azi.begin(), azi.end());
    vector<double>::size_type aziCount = azi.size();
    azi.push_back(azi[0] + 360.);
    azi.push_back(azi[1] + 360.);
    for (vector<double>::size_type i = 0; i < aziCount; i++)
    {
      double gap = azi[i + 1] - azi[i];
      if (gap > primaryAz) primaryAz = gap;
      gap = azi[i + 2] - azi[i];
      if (gap > secondaryAz) secondaryAz = gap;
    }
  }

  // add quality
  DataModel::OriginQuality oq;
  oq.setAssociatedPhaseCount(newOrg->arrivalCount());
  oq.setUsedPhaseCount(usedPhaseCount);
  oq.setAssociatedStationCount(associatedStations.size());
  oq.setUsedStationCount(usedStations.size());
  oq.setStandardError(event.relocInfo.isRelocated ? event.relocInfo.finalRms
                                                  : 0);
  oq.setMedianDistance(meanDist);
  oq.setMinimumDistance(minDist);
  oq.setMaximumDistance(maxDist);
  oq.setAzimuthalGap(primaryAz);
  oq.setSecondaryAzimuthalGap(secondaryAz);
  newOrg->setQuality(oq);
}

HDD::Catalog *
RTDD::getCatalog(const std::string &catalogPath,
                 std::unordered_map<unsigned, DataModel::OriginPtr> *idmap)
{
  std::vector<std::string> tokens = ::splitString(catalogPath, ",");
  try
  {
    if (tokens.size() == 1) // single file containing origin ids
    {
      CatalogDataSrc dataSrc(query(), &_cache, _eventParameters.get());
      HDD::ScCatalog *cat = new HDD::ScCatalog();
      auto _map           = cat->add(tokens[0], dataSrc);
      if (idmap) *idmap = _map;
      return cat;
    }
    else if (tokens.size() == 3) // triplet: station.csv,event.csv,phase.csv
    {
      return new HDD::Catalog(tokens[0], tokens[1], tokens[2], true);
    }
  }
  catch (...)
  {}
  SEISCOMP_ERROR("Cannot load catalog %s", catalogPath.c_str());
  return nullptr;
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
                       bool preloadData,
                       const HDD::CatalogCPtr &alternativeCatalog)
{
  profile->load(query(), &_cache, _eventParameters.get(),
                _config.workingDirectory, _config.saveProcessingFiles,
                _config.cacheWaveforms, _config.cacheAllWaveforms,
                _config.dumpWaveforms, preloadData, alternativeCatalog);
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
    const string &id = row.at("seiscompId");

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
                         bool debugWaveforms,
                         bool preloadData,
                         const HDD::CatalogCPtr &alternativeCatalog)
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
    HDD::CatalogCPtr ddbgc;
    if (alternativeCatalog) // force this catalog
    {
      ddbgc = alternativeCatalog;
    }
    else if (!eventIDFile.empty()) // catalog is a list of origin ids
    {
      CatalogDataSrc dataSrc(query, cache, eventParameters);
      HDD::ScCatalogPtr ddbgc_ = new HDD::ScCatalog();
      ddbgc_->add(eventIDFile, dataSrc);
      ddbgc = ddbgc_;
    }
    else // catalog is extended format station.csv,event.csv,phase.csv
    {
      ddbgc = new HDD::Catalog(stationFile, eventFile, phaFile);
    }

    hypodd = new HDD::HypoDD(ddbgc, ddCfg, pWorkingDir);
    hypodd->setSaveProcessing(saveProcessingFiles);
    hypodd->setUseCatalogWaveformDiskCache(cacheWaveforms);
    hypodd->setWaveformCacheAll(cacheAllWaveforms);
    hypodd->setWaveformDebug(debugWaveforms);

    if (preloadData)
    {
      hypodd->preloadWaveforms();
    }
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

void RTDD::Profile::freeResources()
{
  if (!loaded) return;
  hypodd->unloadTTT();
  hypodd->unloadWaveforms();
  lastUsage = Core::Time::GMT();
}

void RTDD::Profile::unload()
{
  SEISCOMP_INFO("Unloading profile %s", name.c_str());
  hypodd.reset();
  loaded    = false;
  lastUsage = Core::Time::GMT();
}

HDD::CatalogPtr RTDD::Profile::relocateSingleEvent(DataModel::Origin *org)
{
  if (!loaded)
  {
    string msg = Core::stringify(
        "Cannot relocate origin, profile %s not initialized", name.c_str());
    throw runtime_error(msg.c_str());
  }
  lastUsage = Core::Time::GMT();

  CatalogDataSrc dataSrc(query, cache, eventParameters);

  // we pass the stations information from the background catalog, to avoid
  // wasting time accessing the inventory again for information we already have
  HDD::ScCatalogPtr orgToRelocate = new HDD::ScCatalog(
      hypodd->getCatalog()->getStations(), map<unsigned, HDD::Catalog::Event>(),
      unordered_multimap<unsigned, HDD::Catalog::Phase>());
  orgToRelocate->add({org}, dataSrc);

  if (org->evaluationMode() == DataModel::MANUAL)
    hypodd->setUseArtificialPhases(this->useTheoreticalManual);
  else
    hypodd->setUseArtificialPhases(this->useTheoreticalAuto);

  HDD::CatalogPtr rel = hypodd->relocateSingleEvent(
      orgToRelocate, singleEventClustering, singleEventClustering, solverCfg);

  hypodd
      ->unloadTTT(); // free memory and file descriptors (mostly for NLL grids)

  return rel;
}

HDD::CatalogPtr RTDD::Profile::relocateCatalog()
{
  if (!loaded)
  {
    string msg = Core::stringify(
        "Cannot relocate catalog, profile %s not initialized", name.c_str());
    throw runtime_error(msg.c_str());
  }
  lastUsage = Core::Time::GMT();
  hypodd->setUseArtificialPhases(this->useTheoreticalManual);
  return hypodd->relocateMultiEvents(multiEventClustering, solverCfg);
}

void RTDD::Profile::evalXCorr()
{
  if (!loaded)
  {
    string msg = Core::stringify(
        "Cannot evalute cross-correlation settings, profile %s not initialized",
        name.c_str());
    throw runtime_error(msg.c_str());
  }
  lastUsage = Core::Time::GMT();
  hypodd->evalXCorr(multiEventClustering, false);
}

// End Profile class

} // namespace Seiscomp
