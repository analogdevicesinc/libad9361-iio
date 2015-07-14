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
set(LIBAD9361_VERSION_MINOR 1)
set(LIBAD9361_VERSION ${LIBAD9361_VERSION_MAJOR}.${LIBAD9361_VERSION_MINOR})

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries")

set(LIBAD9361_HEADERS ad9361.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_library(ad9361 ad9361_multichip_sync.c ${LIBAD9361_HEADERS})
set_target_properties(ad9361 PROPERTIES VERSION ${LIBAD9361_VERSION} SOVERSION ${LIBAD9361_VERSION_MAJOR} FRAMEWORK TRUE PUBLIC_HEADER "${LIBAD9361_HEADERS}")
target_link_libraries(ad9361 LINK_PRIVATE iio)

set(LIBAD9361_PC ${CMAKE_CURRENT_BINARY_DIR}/libad9361.pc)
configure_file(libad9361.pc.cmakein ${LIBAD9361_PC} @ONLY)
install(FILES ${LIBAD9361_PC} DESTINATION lib/pkgconfig)

if(NOT SKIP_INSTALL_ALL)
	install(TARGETS ad9361
		ARCHIVE DESTINATION lib
		LIBRARY DESTINATION lib
		RUNTIME DESTINATION bin
		FRAMEWORK DESTINATION lib
		PUBLIC_HEADER DESTINATION include)
endif()
