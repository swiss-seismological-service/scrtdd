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

#ifndef __RTDD_PLUGIN_RTDDLOC_H__
#define __RTDD_PLUGIN_RTDDLOC_H__

#include <seiscomp/core/plugin.h>
#include <seiscomp/messaging/connection.h>
#include <seiscomp/seismology/locatorinterface.h>
#include <string>

namespace Seiscomp {

namespace Seismology {

namespace Plugins {

class RTDDLocator : public Seiscomp::Seismology::LocatorInterface
{

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
  virtual bool setParameter(const std::string &name, const std::string &value);

  //! List available profiles
  virtual IDList profiles() const { return _profileNames; }

  //! specify the profile to be used
  virtual void setProfile(const std::string &name);

  //! Returns the implementations capabilities
  virtual int capabilities() const { return NoCapability; }

  //! not supported, only relocation
  virtual DataModel::Origin *locate(PickList &pickList) { return NULL; }

  //! not supported, only relocation
  virtual DataModel::Origin *locate(PickList &pickList,
                                    double initLat,
                                    double initLon,
                                    double initDepth,
                                    const Core::Time &initTime)
  {
    return NULL;
  }

  virtual DataModel::Origin *relocate(const DataModel::Origin *origin);

  // ----------------------------------------------------------------------
  //  Private members
  // ----------------------------------------------------------------------
private:
  Client::ConnectionPtr createConnection();

  ParameterMap _parameters;
  IDList _profileNames;
  std::string _currentProfile;

  static IDList _allowedParameters;
};

} // namespace Plugins

} // namespace Seismology

} // namespace Seiscomp

#endif
