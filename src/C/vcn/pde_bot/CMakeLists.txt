cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

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


include_directories(${CAIRO_INCLUDE_DIRS})
add_library(vcn_pde_bot-cairo SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/finite_element/modules/exporter_cairo.c)
target_link_libraries(vcn_pde_bot-cairo
  vcn_pde_bot
  ${CAIRO_LIBRARIES}
  m)