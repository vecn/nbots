cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -funroll-loops -ffast-math -Wall -pedantic -O3 -fopenmp")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

set(pde_bot
  ${CMAKE_CURRENT_SOURCE_DIR}/material.c
  ${CMAKE_CURRENT_SOURCE_DIR}/boundary_conditions.c
  ${CMAKE_CURRENT_SOURCE_DIR}/finite_element/element.c
  ${CMAKE_CURRENT_SOURCE_DIR}/finite_element/gaussp_to_nodes.c
  ${CMAKE_CURRENT_SOURCE_DIR}/finite_element/elements/trg_linear.c
  ${CMAKE_CURRENT_SOURCE_DIR}/finite_element/solid_mechanics/pipeline.c
  ${CMAKE_CURRENT_SOURCE_DIR}/finite_element/solid_mechanics/static_elasticity2D.c
  ${CMAKE_CURRENT_SOURCE_DIR}/finite_element/solid_mechanics/static_damage2D.c
)

add_library(vcn_pde_bot SHARED ${pde_bot})