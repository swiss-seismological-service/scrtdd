SET(3RDPARTY_DIR ${SC3_PACKAGE_SOURCE_DIR}/libs/hdd/3rd-party)
SET(LSMR_LSQR_DIR ${3RDPARTY_DIR}/lsqr-lsmr)
SET(GMT_DIR ${3RDPARTY_DIR}/gmt)
INCLUDE_DIRECTORIES(${3RDPARTY_DIR} ${LSMR_LSQR_DIR} ${GMT_DIR})

ADD_LIBRARY(hdd SHARED
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
	cvttt.cpp
	solver.cpp
	xcorr.cpp
	waveform.cpp
	dd.cpp
)

INSTALL(TARGETS hdd
	DESTINATION ${SC3_PACKAGE_LIB_DIR}
)

SET(HDD_HEADERS
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
	cvttt.h
	solver.h
	xcorr.h
	waveform.h
	dd.h
)

#INSTALL(
#  FILES ${HDD_HEADERS}
#	DESTINATION ${SC3_PACKAGE_INCLUDE_DIR}/hdd
#)

SET(HDD_HEADERS_3RDPARTY
	${3RDPARTY_DIR}/lrucache.h
	${3RDPARTY_DIR}/geo.h
)

#INSTALL(
#  FILES ${HDD_HEADERS_3RDPARTY}
#	DESTINATION ${SC3_PACKAGE_INCLUDE_DIR}/hdd/3rd-party
#)

SUBDIRS(adapters)

# Add test directory
IF (SC_GLOBAL_UNITTESTS)
	SUBDIRS(test)
ENDIF (SC_GLOBAL_UNITTESTS)
