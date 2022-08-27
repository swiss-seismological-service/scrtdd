SET(TESTS
	ellipsoid.cpp
	ttt.cpp
	clustering.cpp
	dd.cpp
	waveform.cpp
)

INCLUDE_DIRECTORIES(${HDD_DIR})

# add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
# add_link_options(-fsanitize=address)

FOREACH(testSrc ${TESTS})
	GET_FILENAME_COMPONENT(testName ${testSrc} NAME_WE)
	SET(testName test_hdd_${testName})
	ADD_EXECUTABLE(${testName} ${testSrc})
	TARGET_LINK_LIBRARIES(${testName} hdd)
	SC_LINK_LIBRARIES_INTERNAL(${testName} unittest hddsc client)
	SC_LINK_LIBRARIES(${testName})

	ADD_TEST(
		NAME ${testName}
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		COMMAND ${testName}
	)

	# Set the data directory to read the test travel time tables
	SET_TESTS_PROPERTIES(${testName} PROPERTIES
		ENVIRONMENT "SEISCOMP_ROOT=${CMAKE_CURRENT_SOURCE_DIR}/data"
	) 
ENDFOREACH(testSrc)
