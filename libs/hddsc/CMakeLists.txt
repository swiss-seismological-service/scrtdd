INCLUDE_DIRECTORIES(.)

SET(HDDSC_SOURCES
	scutils.cpp
	scttt.cpp
	scwaveform.cpp
)

SET(HDDSC_HEADERS
	scutils.h
	scttt.h
	scwaveform.h
)

SC_ADD_LIBRARY(HDDSC hddsc)

IF (HDD_INSTALL_HEADERS)
SC_LIB_INSTALL_HEADERS(HDDSC)
ENDIF()

IF (HDD_TESTS)
	SUBDIRS(test)
ENDIF()

