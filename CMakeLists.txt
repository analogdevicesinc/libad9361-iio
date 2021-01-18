# Copyright (C) 2015 Analog Devices, Inc.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

cmake_minimum_required(VERSION 2.8.7)
project(ad9361 C)

set(LIBAD9361_VERSION_MAJOR 0)
set(LIBAD9361_VERSION_MINOR 2)
set(VERSION ${LIBAD9361_VERSION_MAJOR}.${LIBAD9361_VERSION_MINOR})

if (WIN32)
	string(TIMESTAMP BUILD_YEAR "%Y")
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	option(OSX_PACKAGE "Create a OSX package" ON)
	set(CMAKE_MACOSX_RPATH ON)
	set(SKIP_INSTALL_ALL ${OSX_PACKAGE})
endif()

include(FindGit OPTIONAL)
if (GIT_FOUND)
	execute_process(
		COMMAND ${GIT_EXECUTABLE} rev-parse --show-toplevel
		WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		OUTPUT_VARIABLE LIBAD9361_GIT_REPO
		OUTPUT_STRIP_TRAILING_WHITESPACE
		ERROR_QUIET
	)

if ("${LIBAD9361_GIT_REPO}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
		execute_process(
			COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
			OUTPUT_VARIABLE LIBAD9361_VERSION_GIT
			OUTPUT_STRIP_TRAILING_WHITESPACE
		)
	endif()
endif()

if (NOT LIBAD9361_VERSION_GIT)
	set(LIBAD9361_VERSION_GIT v${VERSION})
endif()

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries")

if (NOT WIN32)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fvisibility=hidden")
endif()

# based on host, set these on macos
if(APPLE)
        #full Single Unix Standard v3 (SUSv3) conformance (the Unix API)
        add_definitions(-D_DARWIN_C_SOURCE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set(CMAKE_REQUIRED_DEFINITIONS "-D_GNU_SOURCE=1")
        add_definitions(-D_GNU_SOURCE=1)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "FreeBSD|DragonFly|OpenBSD|NetBSD")
        set(CMAKE_REQUIRED_DEFINITIONS "-D__BSD_VISIBLE")
        add_definitions(-D__BSD_VISIBLE=1)
endif()
add_definitions(-D_POSIX_C_SOURCE=200809L -D__XSI_VISIBLE=500 -DLIBAD9361_EXPORTS=1)

set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib"
	CACHE PATH "Installation directory for libraries")
mark_as_advanced(INSTALL_LIB_DIR)

find_library(LIBIIO_LIBRARIES iio)
find_path(LIBIIO_INCLUDEDIR iio.h)

set(LIBAD9361_HEADERS ad9361.h)

file(GLOB_RECURSE FD_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/filterdesigner/*.c)
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${LIBIIO_INCLUDEDIR})
add_library(ad9361 ad9361_multichip_sync.c
            ad9361_baseband_auto_rate.c
            ad9361_design_taps.c
            ad9361_calculate_rf_clock_chain.c
            ad9361_fmcomms5_phase_sync.c
            ${FD_SRC_FILES}
            ${LIBAD9361_HEADERS})
enable_testing()
add_subdirectory(test)
set_target_properties(ad9361 PROPERTIES
	VERSION ${VERSION}
	SOVERSION ${LIBAD9361_VERSION_MAJOR}
	FRAMEWORK TRUE
	PUBLIC_HEADER ${LIBAD9361_HEADERS}
	C_STANDARD 99
	C_STANDARD_REQUIRED ON
	C_EXTENSIONS OFF
)
target_link_libraries(ad9361 LINK_PRIVATE ${LIBIIO_LIBRARIES})

if (MSVC)
	set_target_properties(ad9361 PROPERTIES OUTPUT_NAME libad9361)
endif()

configure_file(libad9361-iio.iss.cmakein ${CMAKE_CURRENT_BINARY_DIR}/libad9361-iio.iss @ONLY)

set(LIBAD9361_PC ${CMAKE_CURRENT_BINARY_DIR}/libad9361.pc)
configure_file(libad9361.pc.cmakein ${LIBAD9361_PC} @ONLY)
install(FILES ${LIBAD9361_PC} DESTINATION ${INSTALL_LIB_DIR}/pkgconfig)

if(NOT SKIP_INSTALL_ALL)
	install(TARGETS ad9361
		ARCHIVE DESTINATION lib
		LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
		RUNTIME DESTINATION bin
		FRAMEWORK DESTINATION lib
		PUBLIC_HEADER DESTINATION include)
endif()

set(CMAKE_INSTALL_DOCDIR "" CACHE PATH "documentation root (DATAROOTDIR/doc/${PROJECT_NAME}${LIBAD9361_VERSION_MAJOR}-doc)")
include(GNUInstallDirs)
set(CMAKE_INSTALL_DOCDIR "${CMAKE_INSTALL_DATAROOTDIR}/doc/${PROJECT_NAME}${LIBAD9361_VERSION_MAJOR}-doc")

find_package(Doxygen)
if(DOXYGEN_FOUND)
	option(WITH_DOC "Generate documentation with Doxygen" ON)

	if (WITH_DOC)
		configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
		set(HTML_DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)
		file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/doc DESTINATION ${HTML_DEST_DIR})

		add_custom_command(TARGET ad9361 POST_BUILD
			COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			COMMENT "Generating API documentation with Doxygen" VERBATIM
			)

		if(NOT SKIP_INSTALL_ALL)
			install(DIRECTORY ${HTML_DEST_DIR} DESTINATION ${CMAKE_INSTALL_DOCDIR})
		endif()
	endif()
else()
	message(STATUS "Doxygen not found, API documentation won't be generated")
endif()

add_subdirectory(bindings/matlab)

# Create an installer if compiling for OSX
if(OSX_PACKAGE)
	set(LIBAD9361_PKG ${CMAKE_CURRENT_BINARY_DIR}/libad9361-${VERSION}.pkg)
	set(LIBAD9361_TEMP_PKG ${CMAKE_CURRENT_BINARY_DIR}/libad9361-${VERSION}-temp.pkg)
	set(LIBAD9361_DISTRIBUTION_XML ${CMAKE_CURRENT_BINARY_DIR}/Distribution.xml)
	set(LIBAD9361_FRAMEWORK_DIR ${CMAKE_CURRENT_BINARY_DIR}/ad9361.framework)
	configure_file(Distribution.xml.cmakein ${LIBAD9361_DISTRIBUTION_XML} @ONLY)

	find_program(PKGBUILD_EXECUTABLE
		NAMES pkgbuild
		DOC "OSX Package builder (pkgbuild)")
	mark_as_advanced(PKGBUILD_EXECUTABLE)

	find_program(PRODUCTBUILD_EXECUTABLE
		NAMES productbuild
		DOC "OSX Package builder (productbuild)")
	mark_as_advanced(PRODUCTBUILD_EXECUTABLE)

	add_custom_command(OUTPUT ${LIBAD9361_PKG}
		COMMAND ${PKGBUILD_EXECUTABLE}
			--component ${LIBAD9361_FRAMEWORK_DIR}
			--identifier com.adi.ad9361 --version ${VERSION}
			--install-location /Library/Frameworks ${LIBAD9361_TEMP_PKG}
		COMMAND ${PRODUCTBUILD_EXECUTABLE}
			--distribution ${LIBAD9361_DISTRIBUTION_XML} ${LIBAD9361_PKG}
		COMMAND ${CMAKE_COMMAND} -E remove ${LIBAD9361_TEMP_PKG}
		DEPENDS ad9361 ${LIBAD9361_DISTRIBUTION_XML}
	)

	if (PKGBUILD_EXECUTABLE AND PRODUCTBUILD_EXECUTABLE)
		add_custom_target(libad9361-pkg ALL DEPENDS ${LIBAD9361_PKG})

		install(CODE "execute_process(COMMAND /usr/sbin/installer -pkg ${LIBAD9361_PKG} -target /)")
	else()
		message(WARNING "Missing pkgbuild or productbuild: OSX installer won't be created.")
	endif()
else()
	# Support creating some basic binpkgs via `make package`.
	# Disabled if OSX_PACKAGE is enabled, as tarballs would
	#   end up empty otherwise.
	option(ENABLE_PACKAGING "Create .deb/.rpm or .tar.gz packages via 'make package'" ON)

	if(ENABLE_PACKAGING)
		if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
			include(cmake/DarwinPackaging.cmake)
		endif()
		if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
			include(cmake/LinuxPackaging.cmake)
		endif()
	endif()
endif()
