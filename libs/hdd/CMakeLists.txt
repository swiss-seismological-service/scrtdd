SET(3RDPARTY_DIR ./3rd-party)
SET(LSMR_LSQR_DIR ${3RDPARTY_DIR}/lsqr-lsmr)
SET(GMT_DIR ${3RDPARTY_DIR}/gmt)
INCLUDE_DIRECTORIES(${3RDPARTY_DIR} ${LSMR_LSQR_DIR} ${GMT_DIR} .)

IF (${HDD_SHARED_LIB})
  SET(LIB_TYPE SHARED)
ELSE()
  SET(LIB_TYPE STATIC)
ENDIF()

ADD_LIBRARY(hdd ${LIB_TYPE}
	${LSMR_LSQR_DIR}/lsmr.cpp
	${LSMR_LSQR_DIR}/lsqr.cpp
	${GMT_DIR}/map_project.cpp
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
	LIBRARY DESTINATION ${HDD_LIB_INSTALL_DIR}
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

SET(HDD_HEADERS_3RDPARTY
	${3RDPARTY_DIR}/lrucache.h
)

IF (HDD_USE_BOOST_FS)
	ADD_DEFINITIONS( -DUSE_BOOST_FS )
ENDIF()

IF (HDD_INSTALL_HEADERS)
INSTALL(
  TARGETS ${HDD_HEADERS}
	PUBLIC_HEADER DESTINATION ${HDD_HEADERS_INSTALL_DIR}/hdd
)
INSTALL(
  TARGETS ${HDD_HEADERS_3RDPARTY}
	PUBLIC_HEADER DESTINATION ${HDD_HEADERS_INSTALL_DIR}/hdd/3rd-party
)  
ENDIF()

IF (HDD_TESTS)
	SUBDIRS(test)
ENDIF()

