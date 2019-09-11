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
 
#include "rtddmsg.h"

namespace Seiscomp {


void RTDDRelocateRequestMessage::serialize(Archive& ar)
{
    Core::Message::serialize(ar);
    if ( !ar.success() ) return;
    ar & NAMED_OBJECT("origin", _origin);
    ar & NAMED_OBJECT("publicID", _publicID);
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
}

IMPLEMENT_SC_CLASS_DERIVED(
    RTDDRelocateResponseMessage, Message, "rtdd_relocate_response_message"
);


} // namespace Seiscomp

