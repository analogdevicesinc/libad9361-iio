add_executable(FilterDesignerTest filter_designer_test.c)
target_link_libraries(FilterDesignerTest ad9361)
add_test(NAME FilterDesignerTest COMMAND FilterDesignerTest WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(FMComms5SyncTest fmcomms5_sync_test.c)
target_link_libraries(FMComms5SyncTest ad9361 iio m)
add_test(NAME FMComms5SyncTest COMMAND FMComms5SyncTest WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
