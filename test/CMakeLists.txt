add_executable(FilterDesignerTest filter_designer_test.c)
target_link_libraries(FilterDesignerTest ad9361)
add_test(NAME FilterDesignerTest
         COMMAND FilterDesignerTest
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(GenerateRatesTest gen_rates_test.c)
target_link_libraries(GenerateRatesTest ad9361)
add_test(NAME GenerateRatesTest
         COMMAND GenerateRatesTest
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

 add_executable(FilterDesignerHardwareTest filter_designer_hw.c)
 target_link_libraries(FilterDesignerHardwareTest ad9361 iio)
 add_test(NAME FilterDesignerHardwareTest
          COMMAND FilterDesignerHardwareTest
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
