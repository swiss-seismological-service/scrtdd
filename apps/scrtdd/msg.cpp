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

#include "msg.h"

namespace Seiscomp {

void RTDDReloadProfileRequestMessage::serialize(Archive &ar)
{
  Core::Message::serialize(ar);
  if (!ar.success()) return;
  ar &NAMED_OBJECT("profile", _profile);
}

IMPLEMENT_SC_CLASS_DERIVED(RTDDReloadProfileRequestMessage,
                           Message,
                           "rtdd_reload_profile_request_message");

void RTDDReloadProfileResponseMessage::serialize(Archive &ar)
{
  Core::Message::serialize(ar);
  if (!ar.success()) return;
  ar &NAMED_OBJECT("error", _error);
}

IMPLEMENT_SC_CLASS_DERIVED(RTDDReloadProfileResponseMessage,
                           Message,
                           "rtdd_reload_profile_response_message");

} // namespace Seiscomp
