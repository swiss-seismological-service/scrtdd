SET(HDDSC_DIR ${SC3_PACKAGE_SOURCE_DIR}/libs/hdd/adapters)

SET(HDDSC_SOURCES
	${HDDSC_DIR}/scutils.cpp
	${HDDSC_DIR}/scttt.cpp
	${HDDSC_DIR}/sclog.cpp
	${HDDSC_DIR}/scwaveform.cpp
)

SET(HDDSC_HEADERS
	${HDDSC_DIR}/scutils.h
	${HDDSC_DIR}/scttt.h
	${HDDSC_DIR}/sclog.h
	${HDDSC_DIR}/scwaveform.h
)

SC_ADD_LIBRARY(HDDSC hddsc)
#SC_LIB_INSTALL_HEADERS(HDDSC)

