cmake_minimum_required(VERSION 2.8.7)

project(ad9361_fm5_phase_sync  C)

add_executable(ad9361_fm5_phase_sync ad9361_fm5_phase_sync.c)
target_link_libraries(ad9361_fm5_phase_sync ad9361 iio m)
set(AD9361_TOOLS_TARGETS ad9361_fm5_phase_sync)

if(NOT SKIP_INSTALL_ALL)
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		install(TARGETS ${AD9361_TOOLS_TARGETS} RUNTIME DESTINATION /Library/Frameworks/ad9361.framework/Tools)
	else()
		install(TARGETS ${AD9361_TOOLS_TARGETS} RUNTIME DESTINATION bin)
	endif()
endif()

# Make the test targets available to the main CMakeLists.txt
set(AD9361_TOOLS_TARGETS ${AD9361_TOOLS_TARGETS} PARENT_SCOPE)
