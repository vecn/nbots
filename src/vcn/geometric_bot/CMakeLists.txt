cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -funroll-loops -ffast-math -Wall -pedantic -O3")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

message("Configuring VCN Geometric Bot")

set(knn
  ${CMAKE_CURRENT_SOURCE_DIR}/knn/bins2D.c
  ${CMAKE_CURRENT_SOURCE_DIR}/knn/bins2D_iterator.c
  ${mesh2D}
)

set(model
  ${CMAKE_CURRENT_SOURCE_DIR}/model/model2D_struct.c
  ${CMAKE_CURRENT_SOURCE_DIR}/model/model2D.c
  ${CMAKE_CURRENT_SOURCE_DIR}/model/modules2D/verifier.c
  ${CMAKE_CURRENT_SOURCE_DIR}/model/modules2D/exporter_asymptote.c
  ${CMAKE_CURRENT_SOURCE_DIR}/model/modules2D/regularizer.c
  ${CMAKE_CURRENT_SOURCE_DIR}/model/modules2D/simplifier.c
  ${CMAKE_CURRENT_SOURCE_DIR}/model/modules2D/vtx.c
  ${CMAKE_CURRENT_SOURCE_DIR}/model/modules2D/edge.c
  ${CMAKE_CURRENT_SOURCE_DIR}/model/modules2D/blender.c
)

set(mesh2D
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/mesh2D_structs.c
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/mesh2D.c
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/ruppert/ruppert.c
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/modules2D/image_density.c
#  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/modules2D/graph_generator.c
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/modules2D/pruner.c
#  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/modules2D/extra_fem.c
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/elements2D/triangles.c
#  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/elements2D/polygons.c
#  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/elements2D/disks.c
)

set(mesh
  ${CMAKE_CURRENT_SOURCE_DIR}/point2D.c
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/dewall.c
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/constrained_delaunay.c
  ${mesh2D}
)

set(geometric_bot
  ${CMAKE_CURRENT_SOURCE_DIR}/utils2D.c
  ${knn}
  ${model}
  ${mesh}
)

add_library(vcn_geometric_bot SHARED ${geometric_bot})
target_link_libraries(vcn_geometric_bot
  vcn_eigen_bot
  vcn_graph_bot
  vcn_image_bot
  m)


set(geometric_bot-cairo
  ${CMAKE_CURRENT_SOURCE_DIR}/model/modules2D/exporter_cairo.c
  ${CMAKE_CURRENT_SOURCE_DIR}/mesh/modules2D/exporter_cairo.c
)
include_directories(${CAIRO_INCLUDE_DIRS})
add_library(vcn_geometric_bot-cairo SHARED ${geometric_bot-cairo})
target_link_libraries(vcn_geometric_bot-cairo
  vcn_geometric_bot
  ${CAIRO_LIBRARIES}
  m)