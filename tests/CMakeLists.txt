cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -funroll-loops -ffast-math -Wall -O3")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include/)

set(launcher_libs
  vcn_exceptions 
  vcn_statistics_bot
  vcn_math_bot
  dl m
)

message("Configuring test type")
add_library(test
  ${CMAKE_CURRENT_SOURCE_DIR}/test.c
)

message("Configuring Unit-Tests launcher.")
add_executable(launch_unit_tests
  ${CMAKE_CURRENT_SOURCE_DIR}/unit_test_driver.c
  ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.c)
target_link_libraries(launch_unit_tests ${launcher_libs})

message("Configuring Performance-Tests launcher.")
add_executable(launch_performance_tests
  ${CMAKE_CURRENT_SOURCE_DIR}/performance_test_driver.c
  ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.c)
target_link_libraries(launch_performance_tests ${launcher_libs} test)

message("Configuring Complexity-Tests launcher.")
add_executable(launch_complexity_tests
  ${CMAKE_CURRENT_SOURCE_DIR}/complexity_funcs.c
  ${CMAKE_CURRENT_SOURCE_DIR}/complexity_test_driver.c
  ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.c)
target_link_libraries(launch_complexity_tests ${launcher_libs} test)

message("Configuring vcn_test_add()")
add_library(vcn_test_add SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/test.c
  ${CMAKE_CURRENT_SOURCE_DIR}/test_add.c
)
target_link_libraries(vcn_test_add vcn_container_bot test)

message("Configuring Complexity-Tests launcher.")
add_executable(valgrind_test
  ${CMAKE_CURRENT_SOURCE_DIR}/valgrind_template/valgrind_test.c
  )
target_link_libraries(valgrind_test vcn_bots)

add_subdirectory(vcn)