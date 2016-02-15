cmake_minimum_required (VERSION 2.8.6)
project (Numerical_Bots)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -funroll-loops -ffast-math -Wall -pedantic -O3")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

message("Configuring NB Container Bot")
add_library(queue_container SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/queue/queue_node.c
  ${CMAKE_CURRENT_SOURCE_DIR}/queue/queue_dst.c
  ${CMAKE_CURRENT_SOURCE_DIR}/queue/queue_iterator.c
)

add_library(stack_container SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/stack/stack_node.c
  ${CMAKE_CURRENT_SOURCE_DIR}/stack/stack_dst.c
  ${CMAKE_CURRENT_SOURCE_DIR}/stack/stack_iterator.c
)

add_library(sorted_container SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/sorted/avl_tree.c
  ${CMAKE_CURRENT_SOURCE_DIR}/sorted/sorted_dst.c
  ${CMAKE_CURRENT_SOURCE_DIR}/sorted/sorted_iterator.c
)

add_library(heap_container SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/heap/heap_tree.c
  ${CMAKE_CURRENT_SOURCE_DIR}/heap/heap_dst.c
  ${CMAKE_CURRENT_SOURCE_DIR}/heap/heap_iterator.c
)

add_library(hash_container SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/hash/hash_dst.c
  ${CMAKE_CURRENT_SOURCE_DIR}/hash/hash_iterator.c
)

set(container_bot
  ${CMAKE_CURRENT_SOURCE_DIR}/array.c
  ${CMAKE_CURRENT_SOURCE_DIR}/container.c
  ${CMAKE_CURRENT_SOURCE_DIR}/iterator.c
)

add_library(nb_container_bot  SHARED ${container_bot})
target_link_libraries(nb_container_bot LINK_INTERFACE_LIBRARIES
  queue_container
  stack_container
  sorted_container
  heap_container
  hash_container
)
  