SET(HDD_SOURCES
	utils.cpp
	lsmr.cpp
	lsqr.cpp
	solver.cpp
	csvreader.cpp
	datasrc.cpp
	catalog.cpp
	waveform.cpp
	clustering.cpp
	ttt.cpp
	hypodd.cpp
)

SET(HDD_HEADERS
	csvreader.h
	datasrc.h
	catalog.h
	hypodd.h
)

SC_ADD_LIBRARY(HDD hdd)
SC_LIB_INSTALL_HEADERS(HDD)
