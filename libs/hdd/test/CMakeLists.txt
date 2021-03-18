SET(TESTS
	waveform.cpp
	ellipsoid.cpp
	clustering.cpp
	ttt.cpp
	solver.cpp
	ddme.cpp
	ddse.cpp
)

INCLUDE_DIRECTORIES(${HDD_DIR})

# add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
# add_link_options(-fsanitize=address)

FOREACH(testSrc ${TESTS})
	GET_FILENAME_COMPONENT(testName ${testSrc} NAME_WE)
	SET(testName test_hdd_${testName})
	ADD_EXECUTABLE(${testName} ${testSrc})
	SC_LINK_LIBRARIES_INTERNAL(${testName} unittest hdd client)
	SC_LINK_LIBRARIES(${testName})

	ADD_TEST(
		NAME ${testName}
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		COMMAND ${testName}
	)
ENDFOREACH(testSrc)
