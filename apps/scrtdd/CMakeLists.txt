SET(RTDD_TARGET scrtdd)

SET(
    RTDD_SOURCES
		main.cpp
		app.cpp
		rtdd.cpp
		hypodd.cpp
		csvreader.cpp
)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

SC_ADD_EXECUTABLE(RTDD ${RTDD_TARGET})
SC_LINK_LIBRARIES_INTERNAL(${RTDD_TARGET} client)
SC_INSTALL_INIT(${RTDD_TARGET} ../../../trunk/apps/templates/initd.py)

FILE(GLOB descs "${CMAKE_CURRENT_SOURCE_DIR}/descriptions/*.xml")
INSTALL(FILES ${descs} DESTINATION ${SC3_PACKAGE_APP_DESC_DIR})
