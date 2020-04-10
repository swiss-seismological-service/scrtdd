SET(RTDD_TARGET scrtdd)

SET(RTDD_SOURCES
		lsmr.cpp
		lsqr.cpp
		solver.cpp
		csvreader.cpp
		datasrc.cpp
		catalog.cpp
		wfmngr.cpp
		hypodd.cpp
		app.cpp
		rtdd.cpp
		main.cpp
)

INCLUDE_DIRECTORIES(${RTDDMSG_DIR}) 

SC_ADD_EXECUTABLE(RTDD ${RTDD_TARGET})
SC_LINK_LIBRARIES_INTERNAL(${RTDD_TARGET} client rtddmsg)
SC_INSTALL_INIT(${RTDD_TARGET} ../../../trunk/apps/templates/initd.py)

FILE(GLOB descs "${CMAKE_CURRENT_SOURCE_DIR}/descriptions/*.xml")
INSTALL(FILES ${descs} DESTINATION ${SC3_PACKAGE_APP_DESC_DIR})
