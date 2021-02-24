SET(HDD_SOURCES
	utils.cpp
	lsmr.cpp
	lsqr.cpp
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
	csvreader.h
	catalog.h
	sccatalog.h
	hypodd.h
)

SC_ADD_LIBRARY(HDD hdd)
SC_LIB_LINK_LIBRARIES_INTERNAL(hdd core)
SC_LIB_INSTALL_HEADERS(HDD)
