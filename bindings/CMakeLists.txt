if (MATLAB_BINDINGS)
	add_subdirectory(matlab)
endif()

if (PYTHON_BINDINGS)
	add_subdirectory(python)
endif()