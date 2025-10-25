SET(LSMR_LSQR_DIR 3rd-party/lsqr-lsmr)
SET(GMT_DIR 3rd-party/gmt)
INCLUDE_DIRECTORIES(${LSMR_LSQR_DIR} ${GMT_DIR})

IF (${HDD_SHARED_LIB})
  SET(LIB_TYPE SHARED)
ELSE()
  SET(LIB_TYPE STATIC)
ENDIF()

ADD_LIBRARY(hdd ${LIB_TYPE}
	${LSMR_LSQR_DIR}/lsmr.cpp
	${LSMR_LSQR_DIR}/lsqr.cpp
	${GMT_DIR}/map_project.cpp
	nll/grid.cpp
	ttt/nllgrid.cpp
	ttt/homogeneous.cpp
	log.cpp
	utils.cpp
	csvreader.cpp
	utctime.cpp
	catalog.cpp
	clustering.cpp
	solver.cpp
	xcorr.cpp
	waveform.cpp
	xcorrcache.cpp
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
	solver.h
	xcorr.h
	waveform.h
	dd.h
)

SET(HDD_HEADERS_TTT
	ttt/nllgrid.h
	ttt/homogeneous.h
)

SET(HDD_HEADERS_NLL
	nll/grid.h
)

SET(HDD_HEADERS_3RDPARTY
	3rd-party/lrucache.h
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
	TARGETS ${HDD_HEADERS_TTT}
	PUBLIC_HEADER DESTINATION ${HDD_HEADERS_INSTALL_DIR}/hdd/ttt
)
INSTALL(
	TARGETS ${HDD_HEADERS_NLL}
	PUBLIC_HEADER DESTINATION ${HDD_HEADERS_INSTALL_DIR}/hdd/nll
)
INSTALL(
	TARGETS ${HDD_HEADERS_3RDPARTY}
	PUBLIC_HEADER DESTINATION ${HDD_HEADERS_INSTALL_DIR}/hdd/3rd-party
)
ENDIF()

IF (HDD_TESTS)
	SUBDIRS(test)
ENDIF()

