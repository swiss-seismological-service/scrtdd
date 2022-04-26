/***************************************************************************
 *   Copyright (C) by ETHZ/SED                                             *
 *                                                                         *
 * This program is free software: you can redistribute it and/or modify    *
 * it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as          *
 * published by the Free Software Foundation, either version 3 of the      *
 * License, or (at your option) any later version.                         *
 *                                                                         *
 * This software is distributed in the hope that it will be useful,        *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 * GNU Affero General Public License for more details.                     *
 *                                                                         *
 *   Developed by Luca Scarabello <luca.scarabello@sed.ethz.ch>            *
 ***************************************************************************/

#ifndef __RTDD_RTDDMSG_H__
#define __RTDD_RTDDMSG_H__

#include <seiscomp/client.h>
#include <seiscomp/core/genericmessage.h>
#include <seiscomp/datamodel/origin.h>
#include <seiscomp/io/archive/xmlarchive.h>

namespace Seiscomp {

DEFINE_SMARTPOINTER(RTDDRelocateRequestMessage);

class SC_SYSTEM_CLIENT_API RTDDRelocateRequestMessage
    : public Seiscomp::Core::Message
{
  DECLARE_SC_CLASS(RTDDRelocateRequestMessage);
  DECLARE_SERIALIZATION;

public:
  RTDDRelocateRequestMessage() : _origin(0), _profile("") {}

  void setOrigin(DataModel::OriginPtr org) { _origin = org; }
  DataModel::OriginPtr getOrigin() const { return _origin; }

  void setProfile(const std::string &name) { _profile = name; }
  std::string getProfile() const { return _profile; }

  //! Implemented interface from Message
  virtual bool empty() const { return false; }

private:
  DataModel::OriginPtr _origin;
  std::string _profile;
};

DEFINE_SMARTPOINTER(RTDDRelocateResponseMessage);

class SC_SYSTEM_CLIENT_API RTDDRelocateResponseMessage
    : public Seiscomp::Core::Message
{
  DECLARE_SC_CLASS(RTDDRelocateResponseMessage);
  DECLARE_SERIALIZATION;

public:
  RTDDRelocateResponseMessage()
      : _relocatedOrigin(0), _error(""), _requestAccepted(false)
  {}

  void setOrigin(DataModel::OriginPtr org) { _relocatedOrigin = org; }
  DataModel::OriginPtr getOrigin() { return _relocatedOrigin; }

  void setError(const std::string err) { _error = err; }
  std::string getError() const { return _error; }
  bool hasError() const { return !_error.empty(); }

  void setRequestAccepted(bool accepted) { _requestAccepted = accepted; }
  bool isRequestAccepted() const { return _requestAccepted; }

  //! Implemented interface from Message
  virtual bool empty() const { return false; }

private:
  DataModel::OriginPtr _relocatedOrigin;
  std::string _error;
  bool _requestAccepted;
};

DEFINE_SMARTPOINTER(RTDDReloadProfileRequestMessage);

class SC_SYSTEM_CLIENT_API RTDDReloadProfileRequestMessage
    : public Seiscomp::Core::Message
{
  DECLARE_SC_CLASS(RTDDReloadProfileRequestMessage)
  DECLARE_SERIALIZATION;

public:
  RTDDReloadProfileRequestMessage() : _profile("") {}

  void setProfile(const std::string &name) { _profile = name; }
  std::string getProfile() const { return _profile; }

  //! Implemented interface from Message
  virtual bool empty() const { return false; }

private:
  std::string _profile;
};

DEFINE_SMARTPOINTER(RTDDReloadProfileResponseMessage);

class SC_SYSTEM_CLIENT_API RTDDReloadProfileResponseMessage
    : public Seiscomp::Core::Message
{
  DECLARE_SC_CLASS(RTDDReloadProfileResponseMessage)
  DECLARE_SERIALIZATION;

public:
  RTDDReloadProfileResponseMessage() : _error("") {}

  void setError(const std::string err) { _error = err; }
  std::string getError() const { return _error; }
  bool hasError() const { return !_error.empty(); }

  //! Implemented interface from Message
  virtual bool empty() const { return false; }

private:
  std::string _error;
};

} // namespace Seiscomp

#endif
