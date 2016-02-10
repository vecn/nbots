cmake_minimum_required (VERSION 2.8.6)
project (Numerical_Bots)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
  "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules")

message("Finding Cairo (Drawing Library)...")
find_package(Cairo)

message("Finding OpenMP...")
find_package(OpenMP)

set(NB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  ${CAIRO_INCLUDE_DIRS}
  )

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -funroll-loops -ffast-math -Wall -O3 ${OpenMP_C_FLAGS}")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

add_subdirectory(nb)

message("Configuring Math bot")
add_library(nb_math_bot SHARED
  ${NB_SOURCE_DIR}/nb/math_bot.c)
target_link_libraries(nb_math_bot m)

message("Configuring Statistics bot")
add_library(nb_statistics_bot SHARED
  ${NB_SOURCE_DIR}/nb/statistics_bot.c)
target_link_libraries(nb_statistics_bot m)

message("Configuring TextFile Bot")
add_library(nb_cfreader_bot SHARED
  ${NB_SOURCE_DIR}/nb/cfreader_cat.c
)

message("Configuring Drawing Bot")
add_library(nb_drawing_bot SHARED
  ${NB_SOURCE_DIR}/nb/visual_cat.c
)

message("Configuring NB Bots Library (ALL in one)")
add_library(nb_bots SHARED
  ${NB_SOURCE_DIR}/nb/math_bot.c
  ${NB_SOURCE_DIR}/nb/statistics_bot.c
)
target_link_libraries(nb_bots LINK_INTERFACE_LIBRARIES
  nb_cfreader_bot
  nb_container_bot
  nb_eigen_bot
  nb_graph_bot
  nb_image_bot
  nb_geometric_bot
  nb_pde_bot
  m 
  )

message("Configuring NB Bots: Cairo Extension Library")
add_library(nb_bots-cairo SHARED
  ${NB_SOURCE_DIR}/nb/geometric_bot/model/modules2D/exporter_cairo.c
  ${NB_SOURCE_DIR}/nb/geometric_bot/mesh/modules2D/exporter_cairo.c
)
target_link_libraries(nb_bots-cairo ${CAIRO_LIBRARIES}
  nb_pde_bot-cairo
  )

message("Configuring NB Exceptions")
add_library(nb_exceptions SHARED
  ${NB_SOURCE_DIR}/nb/exceptions.c
)
target_link_libraries(nb_exceptions nb_container_bot)