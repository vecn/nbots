cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../../include)

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -funroll-loops -ffast-math -Wall -pedantic -O3 -fopenmp")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

add_library(vcn_pde_bot-fem-elements OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/trg_linear.c
)