add_executable(Test filter_designer_test.c)
target_link_libraries(Test ad9361)
add_test(NAME FilterDesignerTest COMMAND Test WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
