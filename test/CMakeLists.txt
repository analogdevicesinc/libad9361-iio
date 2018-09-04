add_executable(FilterDesignerTest filter_designer_test.c)
target_link_libraries(FilterDesignerTest LINK_PRIVATE ${LIBIIO_LIBRARIES} ad9361)
add_test(NAME FilterDesignerTest
         COMMAND FilterDesignerTest
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(GenerateRatesTest gen_rates_test.c)
target_link_libraries(GenerateRatesTest LINK_PRIVATE ${LIBIIO_LIBRARIES} ad9361)
add_test(NAME GenerateRatesTest
         COMMAND GenerateRatesTest
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(FilterDesignerHardwareTest filter_designer_hw.c)
target_link_libraries(FilterDesignerHardwareTest LINK_PRIVATE ${LIBIIO_LIBRARIES} ad9361)
add_test(NAME FilterDesignerHardwareTest
         COMMAND FilterDesignerHardwareTest
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

find_library(LIBIIO_LIBRARIES iio)
add_executable(FMComms5SyncTest fmcomms5_sync_test.c)
target_link_libraries(FMComms5SyncTest ad9361 ${LIBIIO_LIBRARIES})
add_test(NAME FMComms5SyncTest
         COMMAND FMComms5SyncTest
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
