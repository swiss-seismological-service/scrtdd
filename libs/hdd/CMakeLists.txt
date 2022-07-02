SET(ADAPTERS_DIR ${SC3_PACKAGE_SOURCE_DIR}/libs/hdd/adapters)
SET(3RDPARTY_DIR ${SC3_PACKAGE_SOURCE_DIR}/libs/hdd/3rd-party)
SET(LSMR_LSQR_DIR ${3RDPARTY_DIR}/lsqr-lsmr)
SET(GMT_DIR ${3RDPARTY_DIR}/gmt)
INCLUDE_DIRECTORIES(${3RDPARTY_DIR} ${LSMR_LSQR_DIR} ${GMT_DIR})

SET(HDD_SOURCES
	${LSMR_LSQR_DIR}/lsmr.cpp
	${LSMR_LSQR_DIR}/lsqr.cpp
	${GMT_DIR}/map_project.cpp
	${3RDPARTY_DIR}/geo.cpp
	log.cpp
	utils.cpp
	csvreader.cpp
	utctime.cpp
	catalog.cpp
	clustering.cpp
	ttt.cpp
	nllgrid.cpp
	nllttt.cpp
	solver.cpp
	xcorr.cpp
	waveform.cpp
	dd.cpp
	${ADAPTERS_DIR}/scttt.cpp
	${ADAPTERS_DIR}/sclog.cpp
	${ADAPTERS_DIR}/scwaveform.cpp
)

SET(HDD_HEADERS
	${3RDPARTY_DIR}/lrucache.h
	${3RDPARTY_DIR}/geo.h
	ellipsoid.h
	xcorrcache.h
	random.h
	index.h
	trace.h
	timewindow.h
	log.h
	utils.h
	csvreader.h
	utctime.h
	catalog.h
	clustering.h
	ttt.h
	nllgrid.h
	nllttt.h
	solver.h
	xcorr.h
	waveform.h
	dd.h
	${ADAPTERS_DIR}/scttt.h
	${ADAPTERS_DIR}/sclog.h
	${ADAPTERS_DIR}/scconversion.h
)

SC_ADD_LIBRARY(HDD hdd)
SC_LIB_LINK_LIBRARIES_INTERNAL(hdd core)
SC_LIB_INSTALL_HEADERS(HDD)

# Add test directory
IF (SC_GLOBAL_UNITTESTS)
	SUBDIRS(test)
ENDIF (SC_GLOBAL_UNITTESTS)
