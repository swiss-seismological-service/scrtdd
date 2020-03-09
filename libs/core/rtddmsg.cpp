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
 
#include "rtddmsg.h"

namespace Seiscomp {


void RTDDRelocateRequestMessage::serialize(Archive& ar)
{
    Core::Message::serialize(ar);
    if ( !ar.success() ) return;
    ar & NAMED_OBJECT("origin", _origin);
    ar & NAMED_OBJECT("profile", _profile);
}

IMPLEMENT_SC_CLASS_DERIVED(
    RTDDRelocateRequestMessage, Message, "rtdd_relocate_request_message"
);

void RTDDRelocateResponseMessage::serialize(Archive& ar)
{
    Core::Message::serialize(ar);
    if ( !ar.success() ) return;
    ar & NAMED_OBJECT("relocatedOrigin", _relocatedOrigin);
    ar & NAMED_OBJECT("error", _error);
    ar & NAMED_OBJECT("requestAccepted", _requestAccepted);
}

IMPLEMENT_SC_CLASS_DERIVED(
    RTDDRelocateResponseMessage, Message, "rtdd_relocate_response_message"
);


} // namespace Seiscomp

