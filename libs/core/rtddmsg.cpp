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

