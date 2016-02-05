cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../)

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -Wall -pedantic -funroll-loops -ffast-math -O3")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

message("==== TESTS OF GEOMETRIC BOT ====")

message("Configuring Unit-Tests")

set(geometric_bot_libs 
  vcn_test_add  vcn_geometric_bot m
)

message("-- point2D_UT")
add_library(point2D_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/point2D_UT.c)
target_link_libraries(point2D_UT ${geometric_bot_libs})

message("-- utils2D_UT")
add_library(utils2D_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/utils2D_UT.c)
target_link_libraries(utils2D_UT ${geometric_bot_libs})

message("-- bins2D_UT")
add_library(bins2D_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/knn/bins2D_UT.c)
target_link_libraries(bins2D_UT
  vcn_container_bot
  ${geometric_bot_libs}
  )

message("-- bins2D_iterator_UT")
add_library(bins2D_iterator_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/knn/bins2D_iterator_UT.c)
target_link_libraries(bins2D_iterator_UT
  vcn_container_bot
  ${geometric_bot_libs}
  )

message("-- dewall_UT")
add_library(dewall_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/dewall_UT.c)
target_link_libraries(dewall_UT ${geometric_bot_libs})

message("-- constrained_delaunay_UT")
add_library(constrained_delaunay_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/constrained_delaunay_UT.c)
target_link_libraries(constrained_delaunay_UT ${geometric_bot_libs})

message("-- mesh2D_UT")
add_library(mesh2D_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/mesh2D_UT.c)
target_link_libraries(mesh2D_UT
  ${geometric_bot_libs}
  vcn_geometric_bot-cairo #TEMPORAL
  )

message("-- image_density_UT")
add_library(image_density_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/modules2D/image_density_UT.c)
target_link_libraries(image_density_UT
  ${geometric_bot_libs}
  vcn_geometric_bot-cairo #TEMPORAL
)