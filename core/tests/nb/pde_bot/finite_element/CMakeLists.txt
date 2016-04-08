cmake_minimum_required (VERSION 2.8.6)
project (Numerical_Bots)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../ci_tools/test_launchers)

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -Wall -pedantic -funroll-loops -ffast-math -O3")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

message(" ==  FINITE ELEMENT  ==")

message("Configuring Unit-Tests")

message("-- solid_mechanics/static_elasticity2D_UT")
add_library(static_elasticity2D_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/solid_mechanics/static_elasticity2D_UT.c)
target_link_libraries(static_elasticity2D_UT
  nb_test_add
  nb_cfreader_bot
  nb_container_bot
  nb_geometric_bot
  nb_pde_bot
  nb_geometric_bot-cairo #TEMPORAL
  nb_pde_bot-cairo #TEMPORAL
  m)
