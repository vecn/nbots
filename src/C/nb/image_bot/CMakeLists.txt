cmake_minimum_required (VERSION 2.8.6)
project (Numerical_Bots)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -funroll-loops -ffast-math -Wall -pedantic -O3 -fopenmp")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

message("Configuring NB Eigen Bot")

set(image_bot
  ${CMAKE_CURRENT_SOURCE_DIR}/image.c
)

add_library(nb_image_bot SHARED ${image_bot})
target_link_libraries(nb_image_bot m)