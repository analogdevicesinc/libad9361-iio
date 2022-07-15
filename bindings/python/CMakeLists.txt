cmake_minimum_required(VERSION 2.87)
project(libad9361-py NONE)

if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
	if(Python3_EXECUTABLE)
		set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
	endif()
	find_package(PythonInterp 3.6 REQUIRED)

	set(Python3_FOUND ${PYTHONINTERP_FOUND})
	set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
else()
	find_package(Python3 3.6 REQUIRED COMPONENTS Interpreter)
endif()

if (Python3_FOUND)
	message(STATUS "Found Python: Building bindings")
	set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.cmakein)
	set(SETUP_PY ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

	configure_file(${SETUP_PY_IN} ${SETUP_PY})

	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ad9361.py  ${CMAKE_CURRENT_BINARY_DIR}/ad9361.py COPYONLY)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tasks.py  ${CMAKE_CURRENT_BINARY_DIR}/tasks.py COPYONLY)

	add_custom_target(libad9361-py-py ALL DEPENDS ${SETUP_PY} COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} --quiet build)

	if(NOT SKIP_INSTALL_ALL)
		install(CODE "execute_process(WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} install --root=\$ENV{DESTDIR}/ --prefix=${CMAKE_INSTALL_PREFIX})")
	endif()
else()
	message(FATAL_ERROR "Python search failed : Can not build Python bindings")
endif()

if(WITH_DOC)
	find_program(SPHINX_EXECUTABLE
		NAMES sphinx-build
		DOC "Sphinx Documentation Builder (sphinx-doc.org)"
	)
	if (NOT SPHINX_EXECUTABLE)
		message(FATAL_ERROR "Can not build python doc without sphinx-build")
	endif()
	message (STATUS "Building with Python Doc (Sphinx)")
	set(CMAKE_HTML_DEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/html/v${LIBAD9361_VERSION_MAJOR}.${LIBAD9361_VERSION_MINOR}")
	configure_file(
		${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py.in
		${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py)
	add_custom_command(TARGET libad9361-py-py POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E env "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}"
		${SPHINX_EXECUTABLE}
			-b html
			-n
			-c ${CMAKE_CURRENT_SOURCE_DIR}/doc
			-d ${CMAKE_CURRENT_BINARY_DIR}/doctrees
			-w ${CMAKE_BINARY_DIR}/Spx_output_python
			${CMAKE_CURRENT_SOURCE_DIR}/doc
			${CMAKE_HTML_DEST_DIR}/python
		COMMENT "Generating Python binding documentation with Sphinx" VERBATIM
	)
endif()
