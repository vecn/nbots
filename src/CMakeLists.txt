cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
  "${CMAKE_CURRENT_SOURCE_DIR}/../CMake/Modules")

message("Finding Cairo (Drawing Library)...")
find_package(Cairo)

message("Finding OpenMP...")
find_package(OpenMP)

set(VCN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/../include
  ${CAIRO_INCLUDE_DIRS}
  )

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -funroll-loops -ffast-math -Wall -O3 ${OpenMP_C_FLAGS}")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

add_subdirectory(vcn)

message("Configuring VCN Math bot")
add_library(vcn_math_bot SHARED
  ${VCN_SOURCE_DIR}/vcn/math_bot.c)
target_link_libraries(vcn_math_bot m)

message("Configuring VCN Statistics bot")
add_library(vcn_statistics_bot SHARED
  ${VCN_SOURCE_DIR}/vcn/statistics_bot.c)
target_link_libraries(vcn_statistics_bot m)

message("Configuring VCN Bots Library (ALL in one)")
add_library(vcn_bots SHARED
  ${VCN_SOURCE_DIR}/vcn/math_bot.c
  ${VCN_SOURCE_DIR}/vcn/statistics_bot.c
  ${VCN_SOURCE_DIR}/vcn/cfreader_cat.c
  ${VCN_SOURCE_DIR}/vcn/fem_bot.c
)
target_link_libraries(vcn_bots
  vcn_container_bot
  vcn_eigen_bot
  vcn_graph_bot
  vcn_geometric_bot
  m 
  )

message("Configuring VCN Bots: Cairo Extension Library")
add_library(vcn_bots-cairo SHARED
  ${VCN_SOURCE_DIR}/vcn/geometric_bot/model/modules2D/exporter_cairo.c
  ${VCN_SOURCE_DIR}/vcn/geometric_bot/mesh/modules2D/exporter_cairo.c
)
target_link_libraries(vcn_bots-cairo ${CAIRO_LIBRARIES})

message("Configuring VCN Exceptions")
add_library(vcn_exceptions SHARED
  ${VCN_SOURCE_DIR}/vcn/exceptions.c
)
target_link_libraries(vcn_exceptions vcn_container_bot)