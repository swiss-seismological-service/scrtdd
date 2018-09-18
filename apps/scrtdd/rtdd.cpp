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


#include "rtdd.h"
#include "csvreader.h"

#include <seiscomp3/logging/filerotator.h>
#include <seiscomp3/logging/channel.h>

#include <seiscomp3/core/strings.h>
#include <seiscomp3/core/genericrecord.h>
#include <seiscomp3/core/system.h>

#include <seiscomp3/client/inventory.h>
#include <seiscomp3/io/archive/xmlarchive.h>
#include <seiscomp3/io/records/mseedrecord.h>

#include <seiscomp3/datamodel/event.h>
#include <seiscomp3/datamodel/pick.h>
#include <seiscomp3/datamodel/origin.h>
#include <seiscomp3/datamodel/magnitude.h>
#include <seiscomp3/datamodel/utils.h>
#include <seiscomp3/datamodel/parameter.h>
#include <seiscomp3/datamodel/parameterset.h>
#include <seiscomp3/datamodel/journalentry.h>
#include <seiscomp3/datamodel/utils.h>

#include <seiscomp3/math/geo.h>
#include <seiscomp3/math/filter/butterworth.h>

#include <seiscomp3/utils/files.h>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>


using namespace std;
using namespace Seiscomp::Processing;
using namespace Seiscomp::DataModel;


#define JOURNAL_ACTION           "RTDD"
#define JOURNAL_ACTION_COMPLETED "completed"


namespace Seiscomp {

#define NEW_OPT(var, ...) addOption(&var, __VA_ARGS__)
#define NEW_OPT_CLI(var, ...) addOption(&var, NULL, __VA_ARGS__)


namespace {

using Seiscomp::Core::fromString;

// Rectangular region class defining a rectangular region
// by latmin, lonmin, latmax, lonmax.
struct RectangularRegion : public Seiscomp::RTDD::Region
{
	RectangularRegion() {
	}

	bool init(const Application* app, const string &prefix) {
		vector<string> region;
		try { region = app->configGetStrings(prefix + "region"); }
		catch ( ... ) {}

		if ( region.empty() )
			isEmpty = true;
		else {
			isEmpty = false;

			// Parse region
			if ( region.size() != 4 ) {
				SEISCOMP_ERROR("%s: expected 4 values in region definition, got %d",
				               prefix.c_str(), (int)region.size());
				return false;
			}

			if ( !fromString(latMin, region[0]) ||
			     !fromString(lonMin, region[1]) ||
			     !fromString(latMax, region[2]) ||
			     !fromString(lonMax, region[3]) ) {
				SEISCOMP_ERROR("%s: invalid region value(s)", prefix.c_str());
				return false;
			}
		}

		return true;
	}

	bool isInside(double lat, double lon) const {
		if ( isEmpty ) return true;

		double len, dist;

		if ( lat < latMin || lat > latMax ) return false;

		len = lonMax - lonMin;
		if ( len < 0 )
			len += 360.0;

		dist = lon - lonMin;
		if ( dist < 0 )
			dist += 360.0;

		return dist <= len;
	}

	bool isEmpty;
	double latMin, lonMin;
	double latMax, lonMax;
};

// Rectangular region class defining a circular region
// by lat, lon, radius.
struct CircularRegion : public Seiscomp::RTDD::Region
{
	CircularRegion() {
	}

	bool init(const Application* app, const string &prefix) {
		vector<string> region;
		try { region = app->configGetStrings(prefix + "region"); }
		catch ( ... ) {}

		if ( region.empty() )
			isEmpty = true;
		else {
			isEmpty = false;

			// Parse region
			if ( region.size() != 3 ) {
				SEISCOMP_ERROR("%s: expected 3 values in region definition, got %d",
				               prefix.c_str(), (int)region.size());
				return false;
			}

			if ( !fromString(lat, region[0]) ||
			     !fromString(lon, region[1]) ||
			     !fromString(radius, region[2]) ) {
				SEISCOMP_ERROR("%s: invalid region value(s)", prefix.c_str());
				return false;
			}
		}

		return true;
	}

	bool isInside(double lat, double lon) const {
		if ( isEmpty ) return true;

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
	for ( size_t i = 0; i < src.size(); ++i )
		dest[i] = toupper(src[i]);
}


bool startsWith(const string& haystack, const string& needle, bool caseSensitive = true)
{
	string _haystack = haystack;
	string _needle   = needle;
	if ( !caseSensitive )
	{
		makeUpper(_haystack, haystack);
		makeUpper(_needle, needle);
	}
	return _haystack.compare(0, _needle.length(), _needle);
}

}


RTDD::Config::Config()
{
	publicIDPattern = "RTDD.@time/%Y%m%d%H%M%S.%f@.@id@";
	processManualOrigin = true;
	outputPath = "/tmp/rtdd";

	forceProcessing = false;
	testMode = false;
    fExpiry = 1.0;

	wakeupInterval = 10;
	eventMaxIdleTime = 3600;
	logCrontab = true;
	updateDelay = 60;
}



RTDD::RTDD(int argc, char **argv) : Application(argc, argv)
{
	setAutoApplyNotifierEnabled(true);
	setInterpretNotifierEnabled(true);

	setLoadInventoryEnabled(true);
	setLoadConfigModuleEnabled(true);

	setPrimaryMessagingGroup("LOCATION");

	addMessagingSubscription("EVENT");
	addMessagingSubscription("LOCATION");

	_cache.setPopCallback(boost::bind(&RTDD::removedFromCache, this, _1));

	_processingInfoChannel = NULL;
	_processingInfoOutput = NULL;

	NEW_OPT(_config.publicIDPattern, "rtdd.publicIDPattern");
	NEW_OPT(_config.activeProfiles, "rtdd.activeProfiles");
	NEW_OPT(_config.outputPath, "rtdd.outputPath");
	NEW_OPT(_config.processManualOrigin, "rtdd.manualOrigin");

	NEW_OPT(_config.wakeupInterval, "rtdd.cron.wakeupInterval");
	NEW_OPT(_config.eventMaxIdleTime, "rtdd.cron.eventMaxIdleTime");
	NEW_OPT(_config.logCrontab, "rtdd.cron.logging");
	NEW_OPT(_config.updateDelay, "rtdd.cron.updateDelay");
	NEW_OPT(_config.delayTimes, "rtdd.cron.delayTimes");

	NEW_OPT_CLI(_config.testMode, "Mode", "test",
	            "Test mode, no messages are sent", false, true);
	NEW_OPT_CLI(_config.forceProcessing, "Mode", "force",
	            "Force event processing even if a journal entry exists that processing has completed",
	            false, true);
	NEW_OPT_CLI(_config.fExpiry, "Mode", "expiry,x",
	            "Time span in hours after which objects expire", true);
	NEW_OPT_CLI(_config.originID, "Mode", "origin-id,O",
	            "Reprocess the origin(s) and send a message", true);
	NEW_OPT_CLI(_config.eventXML, "Mode", "ep",
	            "Event parameters XML file for offline processing of all contained origins (imply test option)", true);
	NEW_OPT_CLI(_config.forceProfile, "Mode", "profile",
	            "Force this profile to be used", true);
	NEW_OPT_CLI(_config.relocateCatalog, "Mode", "reloc-catalog",
	            "Relocate the catalog of profile passed as argument", true);
	NEW_OPT_CLI(_config.dumpCatalog, "Mode", "dump-catalog",
	            "Dump catalog files from the seiscomp event/origin ids file passed as argument", true);
}



RTDD::~RTDD() {
}


void RTDD::createCommandLineDescription() {
	Application::createCommandLineDescription();
	commandline().addOption("Mode", "dump-config", "Dump the configuration and exit");
}



bool RTDD::validateParameters() 
{
	Environment *env = Environment::Instance();

	if ( !Application::validateParameters() )
		return false;

	// if --ep option is enabled then disable the messaging (offline mode)
	if ( !_config.eventXML.empty())
	{
		setMessagingEnabled(false);
		_config.testMode = true; // we won't send any message
	}

	// If the inventory is provided by an XML file or an event XML
	// is provided, disable the database because we don't need to access it
	if ( !isInventoryDatabaseEnabled() || !_config.eventXML.empty() )
    {
		setDatabaseEnabled(false, false);
    }

	if ( !Util::pathExists(_config.outputPath) ) {
		if ( ! Util::createPath(_config.outputPath) ) {
			SEISCOMP_ERROR("rtdd.outputPath: failed to create path %s",_config.outputPath.c_str());
			return false;
		}
	}

	if ( _config.outputPath.size() > 0 ) {
		if ( _config.outputPath[_config.outputPath.size()-1] != '/' )
			_config.outputPath += '/';
	}

	bool profilesOK = true;

	for ( vector<string>::iterator it = _config.activeProfiles.begin();
	      it != _config.activeProfiles.end(); )
	{

		ProfilePtr prof = new Profile;
		string prefix = string("rtdd.profile.") + *it + ".";

		prof->name = *it;

		try {
			prof->earthModelID = configGetString(prefix + "earthModelID");
		} catch ( ... ) {}
		try {
			prof->methodID = configGetString(prefix + "methodID");
		} catch ( ... ) {}
		if ( ! startsWith(prof->methodID, "RTDD", false) )
		{
			prof->methodID = "RTDD" + prof->methodID;
		}

		string regionType;
		makeUpper(regionType, configGetString(prefix + "regionType"));
		if ( regionType == "RECTANGLE" )
			prof->region = new RectangularRegion;
		else if ( regionType == "CIRCLE" )
			prof->region = new CircularRegion;

		if ( prof->region == NULL ) {
			SEISCOMP_ERROR("rtdd.profile.%s: invalid region type: %s",
			               it->c_str(), regionType.c_str());
			it = _config.activeProfiles.erase(it);
			profilesOK = false;
			continue;
		}

		if ( !prof->region->init(this, prefix) ) {
			SEISCOMP_ERROR("rtdd.profile.%s: invalid region parameters", it->c_str());
			it = _config.activeProfiles.erase(it);
			profilesOK = false;
			continue;
		}

		prefix = string("rtdd.profile.") + *it + ".catalog.";

		prof->eventFile = env->absolutePath(configGetString(prefix + "eventFile"));
		try {
			prof->stationFile = env->absolutePath(configGetString(prefix + "stationFile"));
		} catch ( ... ) {}
		try {
			prof->phaFile = env->absolutePath(configGetString(prefix + "phaFile"));
		} catch ( ... ) {}

		prefix = string("rtdd.profile.") + *it + ".dtct.";
		prof->ddcfg.dtt.minWeight = configGetDouble(prefix + "minWeight");
		prof->ddcfg.dtt.maxESdist = configGetDouble(prefix + "maxESdist");
		prof->ddcfg.dtt.maxIEdist = configGetDouble(prefix + "maxIEdist");
		prof->ddcfg.dtt.minNumNeigh = configGetInt(prefix + "minNumNeigh");
		prof->ddcfg.dtt.maxNumNeigh = configGetInt(prefix + "maxNumNeigh");
		prof->ddcfg.dtt.minDTperEvt = configGetInt(prefix + "minDTperEvt");

		prefix = string("rtdd.profile.") + *it + ".dtcc.";
		prof->ddcfg.xcorr.recordStreamURL = recordStreamURL();
		prof->ddcfg.xcorr.minWeight = configGetDouble(prefix + "minWeight");
		prof->ddcfg.xcorr.maxESdist = configGetDouble(prefix + "maxESdist");
		prof->ddcfg.xcorr.maxIEdist = configGetDouble(prefix + "maxIEdist");
		prof->ddcfg.xcorr.minNumNeigh = configGetInt(prefix + "minNumNeigh");
		prof->ddcfg.xcorr.maxNumNeigh = configGetInt(prefix + "maxNumNeigh");
		prof->ddcfg.xcorr.minDTperEvt = configGetInt(prefix + "minDTperEvt");

		prefix = string("rtdd.profile.") + *it + ".dtcc.crosscorrelation.";
		try {
			prof->ddcfg.xcorr.filterFmin = configGetDouble(prefix + "filterFmin");
		} catch ( ... ) {}
		try {
			prof->ddcfg.xcorr.filterFmax = configGetDouble(prefix + "filterFmax");
		} catch ( ... ) {}
		try {
			prof->ddcfg.xcorr.filterFsamp = configGetDouble(prefix + "filterFsamp");
		} catch ( ... ) {}
		try {
			prof->ddcfg.xcorr.filterOrder = configGetInt(prefix + "filterOrder");
		} catch ( ... ) { prof->ddcfg.xcorr.filterOrder = 3; }
		prof->ddcfg.xcorr.timeBeforePick = configGetDouble(prefix + "timeBeforePick");
		prof->ddcfg.xcorr.timeAfterPick = configGetDouble(prefix + "timeAfterPick");
		prof->ddcfg.xcorr.maxDelay = configGetDouble(prefix + "maxDelay");

		prefix = string("rtdd.profile.") + *it + ".hypodd.";

		try {
			prof->ddcfg.hypodd.exec = env->absolutePath(configGetString(prefix + "execPath"));
		} catch ( ... ) { prof->ddcfg.hypodd.exec = "hypoDD"; }
		prof->ddcfg.hypodd.ctrlFile = env->absolutePath(configGetString(prefix + "controlFile"));

		prefix = string("rtdd.profile.") + *it + ".ph2dt.";

		try {
			prof->ddcfg.ph2dt.exec = env->absolutePath(configGetString(prefix + "execPath"));
		} catch ( ... ) { prof->ddcfg.hypodd.exec = "ph2dt"; }
		try {
			prof->ddcfg.ph2dt.minwght = configGetDouble(prefix + "minwght");
			prof->ddcfg.ph2dt.maxdist = configGetDouble(prefix + "maxdist");
			prof->ddcfg.ph2dt.maxsep = configGetDouble(prefix + "maxsep");
			prof->ddcfg.ph2dt.maxngh = configGetInt(prefix + "maxngh");
			prof->ddcfg.ph2dt.minlnk = configGetInt(prefix + "minlnk");
			prof->ddcfg.ph2dt.minobs = configGetInt(prefix + "minobs");
			prof->ddcfg.ph2dt.maxobs = configGetInt(prefix + "maxobs");
		} catch ( ... ) {}

		_profiles.push_back(prof);

		++it;
	}

	if (!profilesOK) return false;

	if ( commandline().hasOption("dump-config") )
    {
		for ( Options::const_iterator it = options().begin(); it != options().end(); ++it )
        {
			if ( (*it)->cfgName )
				cout << (*it)->cfgName;
			else if ( (*it)->cliParam)
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



bool RTDD::init() {

	if ( !Application::init() )
		return false;

	// Log into processing/info to avoid logging the same information into the global info channel
	_processingInfoChannel = SEISCOMP_DEF_LOGCHANNEL("processing/info", Logging::LL_INFO);
	_processingInfoOutput = new Logging::FileRotatorOutput(Environment::Instance()->logFile("scrtdd-processing-info").c_str(),  60*60*24, 30);

	_processingInfoOutput->subscribe(_processingInfoChannel);

	_inputEvts = addInputObjectLog("event");
	_inputOrgs = addInputObjectLog("origin");
	_outputOrgs = addOutputObjectLog("origin", primaryMessagingGroup());

	_cache.setTimeSpan(Core::TimeSpan(_config.fExpiry*3600.));
	_cache.setDatabaseArchive(query());

	// Enable periodic timer: handleTimeout()
	enableTimer(1);

	// Check each 10 seconds if a new job needs to be started
	_cronCounter = _config.wakeupInterval;

	return true;
}



bool RTDD::run() {

	// dump catalog and exit
	if ( !_config.dumpCatalog.empty() )
	{
		HDD::CatalogPtr cat = new HDD::Catalog(_config.dumpCatalog, query());
		cat->writeToFile("event.csv","phase.csv","station.csv");
		return true;
	}

	// relocate full catalog and exit
	if ( !_config.relocateCatalog.empty() )
	{
		for (list<ProfilePtr>::iterator it = _profiles.begin(); it != _profiles.end(); ++it )
		{
			ProfilePtr profile = *it;
			if ( profile->name == _config.relocateCatalog)
			{
				string workingDir = (boost::filesystem::path(_config.outputPath)/ profile->name).string();
				profile->load(query(), workingDir);
				HDD::CatalogPtr relocatedCat = profile->relocateCatalog();
				profile->unload();
				relocatedCat->writeToFile("event.csv","phase.csv","station.csv");
				break;
			}
		}
		return true;
	}
	
	// relocate passed origin and exit
	if ( !_config.originID.empty() )
	{
		OriginPtr org = _cache.get<Origin>(_config.originID);
		if ( !org ) {
			cerr << "Event " << _config.originID << " not found." << endl;
			return false;
		}

		// Start processing immediately
		_config.delayTimes.clear();
		_config.updateDelay = 0;

		if ( !addProcess(org.get()) )
		    return false;

		return true;
	}

	// relocate xml event and exit
	if ( !_config.eventXML.empty() )
	{
		IO::XMLArchive ar;
		if ( !ar.open(_config.eventXML.c_str()) ) {
			cerr << "Unable to open " << _config.eventXML << endl;
			return false;
		}

		ar >> _eventParameters;
		ar.close();

		if ( !_eventParameters ) {
			cerr << "No event parameters found in " << _config.eventXML << endl;
			return false;
		}

		for(unsigned i = 0; i < _eventParameters->originCount(); i++)
		{
			OriginPtr org = _eventParameters->origin(i);

			// Start processing immediately
			_config.delayTimes.clear();
			_config.updateDelay = 0;

			if ( !addProcess(org.get()) )
			    return false;
		}

		ar.create("-");
		ar.setFormattedOutput(true);
		ar << _eventParameters;
		ar.close();
        return true;
	}

	// real time processing
	return Application::run();
}



void RTDD::done() {
	Application::done();

	// Remove crontab log file if exists
	unlink((Environment::Instance()->logDir() + "/" + name() + ".sched").c_str());

	if ( _processingInfoChannel ) delete _processingInfoChannel;
	if ( _processingInfoOutput )  delete _processingInfoOutput;
}



void RTDD::handleMessage(Core::Message *msg)
{
	Application::handleMessage(msg);

	// Add all events collected by addObject/updateObject
	Todos::iterator it;
	for ( it = _todos.begin(); it != _todos.end(); ++it )
		addProcess(it->get());
	_todos.clear();
}



void RTDD::addObject(const string& parentID, DataModel::Object* object)
{
	updateObject(parentID, object);
}



void RTDD::updateObject(const string &parentID, Object* object)
{
	Origin *origin = Origin::Cast(object);
	if ( origin )
	{
		_cache.feed(origin);
		// process manual origin
		try {
			if (_config.processManualOrigin && origin->evaluationMode() == MANUAL )
			{
				_todos.insert(origin);
				logObject(_inputOrgs, Core::Time::GMT());
			}
		} catch ( ... ) {}
		return;
	}

	Event *event = Event::Cast(object);
	if ( event )
	{
		logObject(_inputEvts, now);

		// store events in _todos so that they will be processed
		if ( !event->registered() )
		{
			EventPtr cached = Event::Find(event->publicID());
			if ( cached )
			{
				_todos.insert(cached.get());
				return;
			}
		}
		_todos.insert(event);
        return;
	}
}



void RTDD::handleTimeout() 
{
	cleanUnusedProfiles();
	runNewJobs();

    //if ( !_config.originID.empty() || !_config.eventXML.empty() )
	//{
	//	quit();
	//}
}



void RTDD::cleanUnusedProfiles() 
{
	// Peridoically clean up profiles unused for some time as they
	// use lots of memory (waveform data)

	Core::TimeSpan expired = Core::TimeSpan(60*60*24); // 1 day

	for (list<ProfilePtr>::iterator it = _profiles.begin(); it != _profiles.end(); ++it )
	{
		ProfilePtr currProfile = *it;
		if ( currProfile->isLoaded() && currProfile->inactiveTime() > expired )
		{
			SEISCOMP_DEBUG("Unload profile %s, inactive for more than %f seconds",
			               currProfile->name.c_str(), expired.length());
			currProfile->unload();
		}
	} 
}



void RTDD::runNewJobs() 
{
	if ( --_cronCounter <= 0 )
	{
		// Reset counter
		_cronCounter = _config.wakeupInterval;

		Crontab::iterator it;
		now = Core::Time::GMT();

		// Update crontab
		for ( it = _crontab.begin(); it != _crontab.end(); )
		{
			Cronjob *job = it->second.get();

			Processes::iterator pit = _processes.find(it->first);
			ProcessPtr proc = (pit == _processes.end()?NULL:pit->second);

			if ( proc == NULL )
			{
				SEISCOMP_WARNING("No processor for cronjob %s", it->first.c_str());
				++it;
				continue;
			}

			// Skip processes where nextRun is not set
			if ( job->runTimes.empty() )
			{
				// Jobs stopped for more than a day now?
				if ( (now - proc->lastRun).seconds() >= _config.eventMaxIdleTime )
				{
					SEISCOMP_DEBUG("Process %s idle time expired, removing",
					               proc->obj->publicID().c_str());
					removeProcess(it, proc.get());
				}
				else
					++it;

				continue;
			}

			Core::Time nextRun = job->runTimes.front();

			// Time of next run in future?
			if ( nextRun > now )
			{
				++it;
				continue;
			}

			// Remove all times in the past
			while ( !job->runTimes.empty() && (job->runTimes.front() <= now) )
				job->runTimes.pop_front();

			// Add eventID to processQueue if not already inserted
			if ( find(_processQueue.begin(), _processQueue.end(), proc) ==
				 _processQueue.end() )
			{
				SEISCOMP_DEBUG("Pushing %s to process queue",
				               proc->obj->publicID().c_str());
				_processQueue.push_back(proc);
			}

			/*
			// No more jobs to start later?
			if ( !job->nextRun.valid() ) {
				_crontab.erase(it++);
				continue;
			}
			*/

			++it;
		}

		// Process event queue
		if ( !_processQueue.empty() )
		{
			ProcessPtr proc = _processQueue.front();
			_processQueue.pop_front();
			startProcess(proc.get());
		}

		// Dump crontab if activated
		if ( _config.logCrontab )
		{
			ofstream of((Environment::Instance()->logDir() + "/" + name() + ".sched").c_str());
			of << "Now: " << now.toString("%F %T") << endl;
			of << "------------------------" << endl;
			of << "[Schedule]" << endl;
			for ( it = _crontab.begin(); it != _crontab.end(); ++it ) {
				if ( !it->second->runTimes.empty() )
					of << it->second->runTimes.front().toString("%F %T") << "\t" << it->first
					   << "\t" << (it->second->runTimes.front()-now).seconds() << endl;
				else
					of << "STOPPED            \t" << it->first << endl;
			}

			// Dump process queue if not empty
			if ( !_processQueue.empty() || _currentProcess ) {
				of << endl << "[Queue]" << endl;

				ProcessQueue::iterator it;
				for ( it = _processQueue.begin(); it != _processQueue.end(); ++it )
					of << "WAITING            \t" << (*it)->obj->publicID() << endl;
				if ( _currentProcess )
					of << "RUNNING            \t" << _currentProcess->obj->publicID() << endl;
			}
		}
	}
}


bool RTDD::addProcess(DataModel::PublicObject* obj)
{
	if (obj == NULL) return false;

	_cache.feed(obj);

	now = Core::Time::GMT();

	// New process?
	ProcessPtr proc;
	Processes::iterator pit = _processes.find(obj->publicID());
	if ( pit == _processes.end() )
	{
		SEISCOMP_DEBUG("Adding process [%s]", obj->publicID().c_str());
		proc = new Process;
		proc->created = now;
		proc->obj = obj;
		_processes[obj->publicID()] = proc;
	}
	else
		proc = pit->second;

	Core::Time nextRun = now + Core::TimeSpan(_config.updateDelay);

	Crontab::iterator it = _crontab.find(obj->publicID());
	if ( it != _crontab.end() )
	{
		// Process currently stopped?
		if ( it->second->runTimes.empty() )
		{
			it->second->runTimes.push_back(nextRun);
			SEISCOMP_DEBUG("Update delay = %ds, next run at %s",
			               _config.updateDelay, nextRun.toString("%FT%T").c_str());
		}
		else
		{
			// Insert next run into queue
			Core::Time first = it->second->runTimes.front();
			if ( (first-nextRun).seconds() >= _config.updateDelay )
				it->second->runTimes.push_front(nextRun);
		}

		return true;
	}

	CronjobPtr job = new Cronjob;

	// Debug to test next run in 20 seconds
	//job->runTimes.push_back(now+Core::TimeSpan(20));

	if ( _config.delayTimes.empty() )
		job->runTimes.push_back(nextRun);
	else
	{
		for ( size_t i = 0; i < _config.delayTimes.size(); ++i )
			job->runTimes.push_back(now + Core::TimeSpan(_config.delayTimes[i]));
	}

	SEISCOMP_DEBUG("%s: adding new cronjob", obj->publicID().c_str());
	_crontab[obj->publicID()] = job;
	handleTimeout();

	return true;
}



bool RTDD::startProcess(Process *proc)
{
	SEISCOMP_DEBUG("Starting process [%s]", proc->obj->publicID().c_str());
	_currentProcess = proc;
	_currentProcess->lastRun = now;

	// if the process contain an origin, process that
	OriginPtr org = Origin::Cast(proc->obj);
	if ( !org )
	{
		// otherwise fetch the preferred origin of the event
		EventPtr obj = Event::Cast(proc->obj);
		OriginPtr org = _cache.get<Origin>(obj->preferredOriginID());
		if ( !org )
		{
			cerr << "Preferred origin " << obj->preferredOriginID() << " not found." << endl;
			_currentProcess = NULL;
			return false;
		}
	}
	process(org.get());
	return true;
}



void RTDD::stopProcess(Process *proc)
{
	Crontab::iterator cit = _crontab.find(proc->obj->publicID());
	if ( cit != _crontab.end() ) cit->second->runTimes.clear();
}



void RTDD::removeProcess(RTDD::Crontab::iterator &it, Process *proc)
{
	// Remove process from process map
	Processes::iterator pit = _processes.find(proc->obj->publicID());
	if ( pit != _processes.end() ) _processes.erase(pit);

	// Remove process from queue
	ProcessQueue::iterator qit = find(_processQueue.begin(), _processQueue.end(), proc);

	if ( qit != _processQueue.end() ) _processQueue.erase(qit);

	// Remove cronjob
	_crontab.erase(it++);
}



void RTDD::process(Origin *origin)
{
	if ( !origin ) return;

	if ( startsWith(origin->methodID(), "RTDD", false) )
	{
		SEISCOMP_INFO("%s: origin was generated by RTDD, skip it",
		              origin->publicID().c_str());
		return;
	}

	if ( isAgencyIDBlocked(objectAgencyID(origin)) )
	{
		SEISCOMP_INFO("%s: origin's agencyID '%s' is blocked",
		              origin->publicID().c_str(), objectAgencyID(origin).c_str());
		return;
	}

	if ( !_config.forceProcessing )
	{
		// Check the origin hasn't been already processed and if it was processed
		// check the processing time is older than origin modification time
		if ( query() )
		{
			DatabaseIterator it;
			JournalEntryPtr entry;
			it = query()->getJournalAction(origin->publicID(), JOURNAL_ACTION);
			while ( (entry = static_cast<JournalEntry*>(*it)) != NULL )
			{
				if ( entry->parameters() == JOURNAL_ACTION_COMPLETED &&
				     entry->created() >= origin->creationInfo().modificationTime() )
				{
					SEISCOMP_INFO("%s: found journal entry \"completely processed\", ignoring origin",
					              origin->publicID().c_str());
					it.close();
					return;
				}
				++it;
			}
			it.close();
			SEISCOMP_DEBUG("No journal entry \"completely processed\" found, go ahead");
		}
	}
	else
		SEISCOMP_DEBUG("Force processing, journal ignored");

	double latitude;
	double longitude;

	try {
		latitude  = origin->latitude().value();
		longitude = origin->longitude().value();
	}
	catch ( ... ) {
		SEISCOMP_WARNING("Ignoring origin %s with unset lat/lon or time",
		                 origin->publicID().c_str());
		return;
	}

	// Find best earth model based on region information and the initial origin
	ProfilePtr currProfile;

	for (list<ProfilePtr>::iterator it = _profiles.begin(); it != _profiles.end(); ++it )
	{
		if ( ! _config.forceProfile.empty() )
		{
			// if user forced a profile, use that
			if ( (*it)->name == _config.forceProfile)
				currProfile = *it;
		}
		else
		{
			// if epicenter is inside the configured region, use it
			if ( (*it)->region->isInside(latitude, longitude) )
				currProfile = *it;
		}
		if (currProfile) break;
	} 

	if ( !currProfile )
	{
		SEISCOMP_ERROR("No profile found for location (lat:%s lon:%s), ignoring origin %s",
		               Core::toString(latitude).c_str(), Core::toString(longitude).c_str(),
		               origin->publicID().c_str());
		return;
	}

	SEISCOMP_DEBUG("Relocating origin '%s'", origin->publicID().c_str());

	OriginPtr newOrg;

	try {
		newOrg = runHypoDD(origin, currProfile);
	}
	catch ( exception &e ) {
		SEISCOMP_ERROR("%s", e.what());
		SEISCOMP_ERROR("Cannot relocate origin %s", origin->publicID().c_str());
		return;
	}

	if ( !newOrg )
	{
		SEISCOMP_ERROR("processing of origin '%s' failed", origin->publicID().c_str());
		return;
	}

	// finished processing, send new origin and update journal
	if (!send(newOrg.get()) )
		SEISCOMP_ERROR("%s: sending of derived origin failed", origin->publicID().c_str());

	SEISCOMP_INFO("Origin %s has been processed, stop process", origin->publicID().c_str());

	stopProcess(_currentProcess.get());

	if ( connection() && !_config.testMode )
	{
		DataModel::Journaling journal;
		JournalEntryPtr entry = new JournalEntry;
		entry->setObjectID(origin->publicID());
		entry->setAction(JOURNAL_ACTION);
		entry->setParameters(JOURNAL_ACTION_COMPLETED);
		entry->setSender(name() + "@" + Core::getHostname());
		entry->setCreated(Core::Time::GMT());
		Notifier::Enable();
		Notifier::Create(journal.publicID(), OP_ADD, entry.get());
		Notifier::Disable();

		Core::MessagePtr msg = Notifier::GetMessage();
		if ( msg ) connection()->send("EVENT", msg.get());
	}

	_currentProcess = NULL;
	handleTimeout();


}



void RTDD::removedFromCache(Seiscomp::DataModel::PublicObject *po) {
    // do nothing
}



bool RTDD::send(Origin *org)
{
	if ( org == NULL ) return false;

	logObject(_outputOrgs, Core::Time::GMT());

	if (!_config.eventXML.empty())
	{
        _eventParameters->add(org);
	}

	if ( _config.testMode ) return true;

	EventParametersPtr ep = new EventParameters;

	bool wasEnabled = Notifier::IsEnabled();
	Notifier::Enable();

	// Insert origin to event parameters
	ep->add(org);

	NotifierMessagePtr msg = Notifier::GetMessage();

	bool result = false;
	if ( connection() )
		result = connection()->send(msg.get());

	Notifier::SetEnabled(wasEnabled);

	return result;
}



OriginPtr RTDD::runHypoDD(Origin *org, ProfilePtr profile)
{
	string workingDir = (boost::filesystem::path(_config.outputPath)/profile->name).string();
	profile->load(query(), workingDir);

	OriginPtr newOrg;

	if ( !_config.publicIDPattern.empty() )
	{
		newOrg = Origin::Create("");
		PublicObject::GenerateId(&*newOrg, _config.publicIDPattern);
	}
	else
		newOrg = Origin::Create();

	DataModel::CreationInfo ci;
	ci.setAgencyID(agencyID());
	ci.setAuthor(author());
	ci.setCreationTime(Core::Time().gmt());

	newOrg->setCreationInfo(ci);
	newOrg->setEarthModelID(profile->earthModelID);
	newOrg->setMethodID(profile->methodID);
	newOrg->setEvaluationMode(EvaluationMode(AUTOMATIC));
	newOrg->setEpicenterFixed(true);

	HDD::CatalogPtr relocatedOrg = profile->relocateSingleEvent(org);
	// there must be only one event in the catalog, the relocated origin
	const HDD::Catalog::Event& event = relocatedOrg->getEvents().begin()->second;

	newOrg->setLatitude(DataModel::RealQuantity(event.latitude));
	newOrg->setLongitude(DataModel::RealQuantity(event.longitude));
	newOrg->setDepth(DataModel::RealQuantity(event.depth));
	newOrg->setTime(DataModel::TimeQuantity(event.time));

	return newOrg;
}



// Profile class

RTDD::Profile::Profile()
{
	loaded = false;
}


void RTDD::Profile::load(DataModel::DatabaseQuery* query, string workingDir)
{
	if ( loaded ) return;
	SEISCOMP_DEBUG("Loading profile %s", name.c_str());

	this->query = query;

	// check if the catalog uses seiscomp event ids or if it contains full information
	HDD::CatalogPtr ddbgc;
	if ( CSV::readWithHeader(eventFile)[0].count("seiscompId") != 0)
		ddbgc = new HDD::Catalog(eventFile, query);
	else
		ddbgc = new HDD::Catalog(stationFile, eventFile, phaFile);

	hypodd = new HDD::HypoDD(ddbgc, ddcfg, workingDir);
	loaded = true;
	lastUsage = Core::Time::GMT();
}


void RTDD::Profile::unload()
{
	hypodd.reset();
	loaded = false;
	lastUsage = Core::Time();
}


HDD::CatalogPtr RTDD::Profile::relocateSingleEvent(Origin *org)
{
	if ( !loaded )
	{
		string msg = Core::stringify("Cannot relocate origin, profile %s not initialized", name);
		throw runtime_error(msg.c_str());
	}
	lastUsage = Core::Time::GMT();

	// FIXME: bad performance since we are not passing the cache but only the query object
	HDD::CatalogPtr orgToRelocate = new HDD::Catalog({org->publicID()}, query);
	return hypodd->relocateSingleEvent(orgToRelocate);
}



HDD::CatalogPtr RTDD::Profile::relocateCatalog()
{
	if ( !loaded )
	{
		string msg = Core::stringify("Cannot relocate catalog, profile %s not initialized", name);
		throw runtime_error(msg.c_str());
	}
	lastUsage = Core::Time::GMT();
	return hypodd->relocateCatalog();
}

// End Profile class

} // Seiscomp

