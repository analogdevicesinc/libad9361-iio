find_library(LIBIIO_LIBRARIES iio)
if(NOT LIBIIO_LIBRARIES)
  message(FATAL_ERROR "libiio library not found")
endif()
find_library(AD9361_LIBRARIES ad9361)
if(NOT AD9361_LIBRARIES)
  message(FATAL_ERROR "libad9361 library not found")
endif()

add_executable(FilterDesignerTest filter_designer_test.c)
target_link_libraries(FilterDesignerTest ad9361 ${LIBIIO_LIBRARIES})
add_test(NAME FilterDesignerTest
         COMMAND FilterDesignerTest
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(GenerateRatesTest gen_rates_test.c)
target_link_libraries(GenerateRatesTest ad9361 ${LIBIIO_LIBRARIES})
add_test(NAME GenerateRatesTest
         COMMAND GenerateRatesTest
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(FilterDesignerHardwareTest filter_designer_hw.c)
target_link_libraries(FilterDesignerHardwareTest ad9361 ${LIBIIO_LIBRARIES})
add_test(NAME FilterDesignerHardwareTest
         COMMAND FilterDesignerHardwareTest
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(FMComms5SyncTest fmcomms5_sync_test.c)
target_link_libraries(FMComms5SyncTest ad9361 ${LIBIIO_LIBRARIES})
add_test(NAME FMComms5SyncTest
         COMMAND FMComms5SyncTest
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
