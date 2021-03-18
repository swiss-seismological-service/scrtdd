SET(LSMR_LSQR_DIR ${SC3_PACKAGE_SOURCE_DIR}/libs/3rd-party/lsqr-lsmr)
INCLUDE_DIRECTORIES(${LSMR_LSQR_DIR})

SET(HDD_SOURCES
	${LSMR_LSQR_DIR}/lsmr.cpp
	${LSMR_LSQR_DIR}/lsqr.cpp
	utils.cpp
	solver.cpp
	csvreader.cpp
	catalog.cpp
	sccatalog.cpp
	waveform.cpp
	clustering.cpp
	ttt.cpp
	scttt.cpp
	nllttt.cpp
	hypodd.cpp
)

SET(HDD_HEADERS
	ellipsoid.ipp
	xcorrcache.ipp
	utils.h
	clustering.h
	solver.h
	ttt.h
	scttt.h
	nllttt.h
	waveform.h
	csvreader.h
	catalog.h
	sccatalog.h
	hypodd.h
)

SC_ADD_LIBRARY(HDD hdd)
SC_LIB_LINK_LIBRARIES_INTERNAL(hdd core)
SC_LIB_INSTALL_HEADERS(HDD)
