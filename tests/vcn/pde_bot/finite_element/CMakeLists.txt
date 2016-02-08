cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../)

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -Wall -pedantic -funroll-loops -ffast-math -O3")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

message(" ==  FINITE ELEMENT  ==")

message("Configuring Unit-Tests")

message("-- solid_mechanics/static_elasticity2D_UT")
add_library(static_elasticity_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/solid_mechanics/static_elasticity2D_UT.c)
target_link_libraries(static_elasticity_UT
  vcn_container_bot
  vcn_geometric_bot
  vcn_pde_bot
  m)
