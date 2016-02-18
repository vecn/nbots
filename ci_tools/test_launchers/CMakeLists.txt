cmake_minimum_required (VERSION 2.8.6)
project (Numerical_Bots)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../include/)

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -funroll-loops -ffast-math -Wall -O3 ${OpenMP_C_FLAGS}")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 


set(launcher_libs
  nb_exceptions 
  nb_statistics_bot
  nb_math_bot
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

message("Configuring nb_test_add()")
add_library(nb_test_add SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/test.c
  ${CMAKE_CURRENT_SOURCE_DIR}/test_add.c
)
target_link_libraries(nb_test_add nb_container_bot test)