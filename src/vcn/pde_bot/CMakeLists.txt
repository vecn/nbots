cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

add_subdirectory(finite_element)

add_library(vcn_pde_bot SHARED
  $<TARGET_OBJECTS:vcn_pde_bot-fem>
  )
target_link_libraries(vcn_pde_bot m)