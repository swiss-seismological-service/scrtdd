SET(RTDD_TARGET scrtdd)

SET(RTDD_SOURCES
		hdd/utils.cpp
		hdd/lsmr.cpp
		hdd/lsqr.cpp
		hdd/solver.cpp
		hdd/csvreader.cpp
		hdd/datasrc.cpp
		hdd/catalog.cpp
		hdd/waveform.cpp
		hdd/clustering.cpp
		hdd/ttt.cpp
		hdd/hypodd.cpp
		app.cpp
		rtdd.cpp
		main.cpp
)

INCLUDE_DIRECTORIES(hdd)
INCLUDE_DIRECTORIES(${RTDDMSG_DIR}) 

SC_ADD_EXECUTABLE(RTDD ${RTDD_TARGET})
SC_LINK_LIBRARIES_INTERNAL(${RTDD_TARGET} client rtddmsg)
SC_INSTALL_INIT(${RTDD_TARGET} ../../../trunk/apps/templates/initd.py)

FILE(GLOB descs "${CMAKE_CURRENT_SOURCE_DIR}/descriptions/*.xml")
INSTALL(FILES ${descs} DESTINATION ${SC3_PACKAGE_APP_DESC_DIR})
