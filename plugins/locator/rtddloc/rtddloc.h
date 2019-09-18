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

#ifndef __RTDD_PLUGIN_RTDDLOC_H__
#define __RTDD_PLUGIN_RTDDLOC_H__

#include <seiscomp3/core/plugin.h>
#include <seiscomp3/seismology/locatorinterface.h>
#include <seiscomp3/communication/connection.h>
#include <string>


namespace Seiscomp {

namespace Seismology {

namespace Plugins {


class RTDDLocator : public Seiscomp::Seismology::LocatorInterface {

    // ----------------------------------------------------------------------
    //  X'truction
    // ----------------------------------------------------------------------
    public:
        //! C'tor
        RTDDLocator() {}

        //! D'tor
        ~RTDDLocator() {}

    // ----------------------------------------------------------------------
    //  Locator interface implementation
    // ----------------------------------------------------------------------
    public:
        //! Initializes the locator.
        virtual bool init(const Config::Config &config);

        //! Returns supported parameters to be changed.
        virtual IDList parameters() const { return _allowedParameters; }

        //! Returns the value of a parameter.
        virtual std::string parameter(const std::string &name) const;

        //! Sets the value of a parameter.
        virtual bool setParameter(const std::string &name,
                                  const std::string &value);
                                  
        //! List available profiles
        virtual IDList profiles() const { return _profileNames; }

        //! specify the profile to be used
        virtual void setProfile(const std::string &name);
        
        //! Returns the implementations capabilities
        virtual int capabilities() const { return NoCapability; }

        //! not supported, only relocation
        virtual DataModel::Origin* locate(PickList& pickList) { return NULL; }
        
        //! not supported, only relocation
        virtual DataModel::Origin* locate(PickList& pickList,
                                          double initLat, double initLon, double initDepth,
                                          const Core::Time &initTime) { return NULL; } 
                                          
        virtual DataModel::Origin* relocate(const DataModel::Origin* origin);

    // ----------------------------------------------------------------------
    //  Private members
    // ----------------------------------------------------------------------
    private:
    
        Communication::ConnectionPtr createConnection();
    
        ParameterMap  _parameters;
        IDList        _profileNames;
        std::string   _currentProfile;
        
        static IDList _allowedParameters;
};


}

}

}

#endif

