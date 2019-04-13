cmake_minimum_required(VERSION 2.8.7)

project(ad9361_filter    C)

# if (MSVC)
# 	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../deps/wingetopt/src)
# 	set(GETOPT_C_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../deps/wingetopt/src/getopt.c)
# endif (MSVC)

if (WIN32)
	set(LIBIIO_ORIGINAL_FILENAME )
	set(LIBIIO_RC ${CMAKE_CURRENT_BINARY_DIR}/properties.rc)
	configure_file(../properties.rc.cmakein ${LIBIIO_RC} @ONLY)
endif()
add_executable(ad9361_filter  ad9361_filter.c ${GETOPT_C_FILE} ${LIBIIO_RC})

target_link_libraries(ad9361_filter ad9361)

set(AD9361_APPS_TARGETS ad9361_filter)

# if(PTHREAD_LIBRARIES)
# 	project(iio_adi_xflow_check C)
# 	add_executable(iio_adi_xflow_check iio_adi_xflow_check.c)
# 	target_link_libraries(iio_adi_xflow_check iio ${PTHREAD_LIBRARIES})
# 	set(IIO_TESTS_TARGETS ${IIO_TESTS_TARGETS} iio_adi_xflow_check)
#
# 	target_link_libraries(iio_readdev  ${PTHREAD_LIBRARIES})
# 	target_link_libraries(iio_writedev ${PTHREAD_LIBRARIES})
# endif()

set_target_properties(${IIO_TESTS_TARGETS} PROPERTIES
	C_STANDARD 99
	C_STANDARD_REQUIRED ON
	C_EXTENSIONS OFF
)

# if(NOT SKIP_INSTALL_ALL)
# 	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
# 		install(TARGETS ${IIO_TESTS_TARGETS} RUNTIME DESTINATION /Library/Frameworks/iio.framework/Tools)
# 	else()
# 		install(TARGETS ${IIO_TESTS_TARGETS} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
# 	endif()
# endif()

# Make the test targets available to the main CMakeLists.txt
set(AD9361_APPS_TARGETS ${IIO_TESTS_TARGETS} PARENT_SCOPE)
