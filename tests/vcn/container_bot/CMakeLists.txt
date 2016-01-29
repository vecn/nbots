cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../)

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -Wall -pedantic -funroll-loops -ffast-math -O3")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 


message("Configuring Unit-Tests")

message("-- array_UT")
add_library(array_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/array_UT.c)
target_link_libraries(array_UT vcn_test_add vcn_container_bot m)

set(container_UT_libs 
  vcn_test_add 
  vcn_container_bot
  vcn_statistics_bot
)

message("-- container_QUEUE_UT")
add_library(container_QUEUE_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/container_QUEUE_UT.c)
target_link_libraries(container_QUEUE_UT ${container_UT_libs})

message("-- container_STACK_UT")
add_library(container_STACK_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/container_STACK_UT.c)
target_link_libraries(container_STACK_UT ${container_UT_libs})

message("-- container_SORTED_UT")
add_library(container_SORTED_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/container_SORTED_UT.c)
target_link_libraries(container_SORTED_UT ${container_UT_libs})

message("-- container_HEAP_UT")
add_library(container_HEAP_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/container_HEAP_UT.c)
target_link_libraries(container_HEAP_UT ${container_UT_libs})

message("-- container_HASH_UT")
add_library(container_HASH_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/container_HASH_UT.c)
target_link_libraries(container_HASH_UT ${container_UT_libs})

set(iterator_UT_libs 
  vcn_test_add 
  vcn_container_bot
)

message("-- iterator_QUEUE_UT")
add_library(iterator_QUEUE_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/iterator_QUEUE_UT.c)
target_link_libraries(iterator_QUEUE_UT ${iterator_UT_libs})

message("-- iterator_STACK_UT")
add_library(iterator_STACK_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/iterator_STACK_UT.c)
target_link_libraries(iterator_STACK_UT ${iterator_UT_libs})

message("-- iterator_SORTED_UT")
add_library(iterator_SORTED_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/iterator_SORTED_UT.c)
target_link_libraries(iterator_SORTED_UT ${iterator_UT_libs})

message("-- iterator_HEAP_UT")
add_library(iterator_HEAP_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/iterator_HEAP_UT.c)
target_link_libraries(iterator_HEAP_UT ${iterator_UT_libs})

message("-- iterator_HASH_UT")
add_library(iterator_HASH_UT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/iterator_HASH_UT.c)
target_link_libraries(iterator_HASH_UT ${iterator_UT_libs})

message("Configuring Performance-Tests")

message("-- array_PT")
add_library(array_PT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/array_PT.c)
target_link_libraries(array_PT vcn_test_add vcn_container_bot m)


message("Configuring Complexity-Tests")

message("-- array_CT")
add_library(array_CT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/array_CT.c)
target_link_libraries(array_CT vcn_test_add vcn_container_bot m)

set(container_CT_libs 
  vcn_test_add 
  vcn_container_bot
  vcn_statistics_bot
)

message("-- container_QUEUE_CT")
add_library(container_QUEUE_CT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/container_QUEUE_CT.c)
target_link_libraries(container_QUEUE_CT ${container_CT_libs})

message("-- container_STACK_CT")
add_library(container_STACK_CT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/container_STACK_CT.c)
target_link_libraries(container_STACK_CT ${container_CT_libs})

message("-- container_SORTED_CT")
add_library(container_SORTED_CT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/container_SORTED_CT.c)
target_link_libraries(container_SORTED_CT ${container_CT_libs})

message("-- container_HEAP_CT")
add_library(container_HEAP_CT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/container_HEAP_CT.c)
target_link_libraries(container_HEAP_CT ${container_CT_libs})

message("-- container_HASH_CT")
add_library(container_HASH_CT MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/container_HASH_CT.c)
target_link_libraries(container_HASH_CT ${container_CT_libs})