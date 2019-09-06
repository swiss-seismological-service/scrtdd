#include "rtddmsg.h"

namespace Seiscomp {


void RTDDRelocateRequestMessage::serialize(Archive& ar) {}

IMPLEMENT_SC_CLASS_DERIVED(
    RTDDRelocateRequestMessage, Message, "rtdd_relocate_request_message"
);

void RTDDRelocateResponseMessage::serialize(Archive& ar) {}

IMPLEMENT_SC_CLASS_DERIVED(
    RTDDRelocateResponseMessage, Message, "rtdd_relocate_response_message"
);


} // namespace Seiscomp

