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


#define SEISCOMP_COMPONENT RTDDLocator
#define EXTERN_MODE

#include "rtddloc.h"
#include "rtddmsg.h"
#include <seiscomp3/logging/log.h>
#include <seiscomp3/core/system.h>

ADD_SC_PLUGIN(
    "Locator implementation using scrtdd (real time double difference)",
    "Luca Scarabello, Swiss Seismological Service ETH Zuerich",
    1, 0, 0
)


using namespace std;
using namespace Seiscomp::Core;
using namespace Seiscomp::DataModel;


namespace Seiscomp {

namespace Seismology {

namespace Plugins {


REGISTER_LOCATOR(RTDDLocator, "RTDD");

RTDDLocator::IDList RTDDLocator::_allowedParameters;

bool RTDDLocator::init(const Config::Config &config)
{
    //
    // allowed parameters    
    //
    if ( _allowedParameters.empty() ) {
        _allowedParameters.push_back("RTDD_HOST");  
        _allowedParameters.push_back("RELOC_TIMEOUT");
    }
    // Set defaults value
    _parameters["RTDD_HOST"] = "localhost";
    _parameters["RELOC_TIMEOUT"] = "300"; // seconds

    //
    // Set up profiles
    //
    _profileNames.clear();

    try { _profileNames = config.getStrings("rtddloc.activeProfiles"); }
    catch ( ... ) { }

    _profileNames.insert(_profileNames.begin(), "automatic");
    _currentProfile = "";

    return true;
}


void RTDDLocator::setProfile(const string &name)
{
    if ( find(_profileNames.begin(), _profileNames.end(), name) ==
         _profileNames.end() )
        return;

    _currentProfile = (name == "automatic" ? "" : name);
}


string RTDDLocator::parameter(const string &name) const
{
    ParameterMap::const_iterator it = _parameters.find(name);
    if ( it != _parameters.end() )
        return it->second;

    return "";
}


bool RTDDLocator::setParameter(const string &name,
                             const string &value)
 {
    ParameterMap::iterator it = _parameters.find(name);
    if ( it == _parameters.end() )
        return false;

    it->second = value;
    return true;
}


Communication::ConnectionPtr RTDDLocator::createConnection()
{
    Communication::ConnectionPtr connection;

    connection = Communication::Connection::Create(_parameters["RTDD_HOST"], "", "SERVICE_REQUEST");
    if ( ! connection )
    {
        throw GeneralException(stringify("Failed to create connection to '%s'", _parameters["RTDD_HOST"].c_str()));
    }

    if ( connection->subscribe("SERVICE_REQUEST") != Status::SEISCOMP_SUCCESS )
    {
        connection.reset();
        throw GeneralException("Failed to subscribe to SERVICE_REQUEST");
    }

    return connection;
}


Origin* RTDDLocator::relocate(const Origin *origin)
{
    if ( origin == NULL ) return NULL;

    SEISCOMP_DEBUG("Relocating origin (%s) with profile '%s'",
                   origin->publicID().c_str(), _currentProfile.c_str());

    //
    // prepare the connection to the messaging
    //
    double msgWaitTimeout = std::stod(_parameters["RELOC_TIMEOUT"]);
    Communication::ConnectionPtr connection = createConnection();

    //
    // send relocation request
    //
    RTDDRelocateRequestMessage msg;
    OriginPtr nonConstOrg = DataModel::Origin::Cast(origin->clone());
    for (size_t i = 0; i < origin->arrivalCount(); i++)
        nonConstOrg->add(DataModel::Arrival::Cast(origin->arrival(i)->clone()));
    msg.setOrigin(nonConstOrg);
    msg.setProfile(_currentProfile);
    connection->send(&msg);

    //
    // wait 1 second for a request accepted response
    //
    string lastError;
    bool requestAccepted = false;
    Core::Time started = Core::Time::GMT();
    while ( (Core::Time::GMT() - started).length() < 1 )
    {
        int error;
        Message *msg = connection->readMessage(false, Communication::Connection::READ_ALL, NULL, &error);
        RTDDRelocateResponseMessage* relocMsg = RTDDRelocateResponseMessage::Cast(msg);
        if ( relocMsg )
        {
            SEISCOMP_DEBUG("Received RTDDRelocateResponseMessage");
            if ( relocMsg->hasError() )
            {
                lastError = relocMsg->getError();
                // do not exit here, maybe there is more than one scrtdd running
            }
            else if ( relocMsg->isRequestAccepted() )
            {
                requestAccepted = true;
                break;
            }
        }
        Core::msleep(0.1);
    }

    //
    // If nobody accepted the request maybe there is no scrtdd listening
    //
    if ( ! requestAccepted )
    {
        if ( lastError.empty() )
            throw LocatorException("No reply from scrtdd. Has scrtdd module been started?");
        else
            throw LocatorException(lastError);
    }

    //
    // wait for a relocation reply
    //
    started = Core::Time::GMT();
    while ( (Core::Time::GMT() - started).length() < msgWaitTimeout)
    {
        int error;
        Message *msg = connection->readMessage(false, Communication::Connection::READ_ALL, NULL, &error);
        RTDDRelocateResponseMessage* relocMsg = RTDDRelocateResponseMessage::Cast(msg);
        if ( relocMsg )
        {
            SEISCOMP_DEBUG("Received RTDDRelocateResponseMessage");
            if ( relocMsg->hasError() )
                throw LocatorException(relocMsg->getError());
            return relocMsg->getOrigin().get();
        }
        Core::msleep(0.1);
    }

    throw LocatorException("No reply from scrtdd. Has scrtdd module been started?");
}


} // namespace Plugins

} // namespace Seismology

} // namespace Seiscomp

