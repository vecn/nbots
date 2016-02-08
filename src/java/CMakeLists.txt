cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

find_package(Java)
include(UseJava)

set(CMAKE_JNI_TARGET TRUE)

message("Configuring (Java) VCNGeometricBot")

set(geometricBot
  ${CMAKE_CURRENT_SOURCE_DIR}/vcn/geometricBot/Box.java
  ${CMAKE_CURRENT_SOURCE_DIR}/vcn/geometricBot/Mesh.java
  ${CMAKE_CURRENT_SOURCE_DIR}/vcn/geometricBot/Model.java
  ${CMAKE_CURRENT_SOURCE_DIR}/vcn/geometricBot/GeometricBot.java
  ${CMAKE_CURRENT_SOURCE_DIR}/vcn/geometricBot/GeometricBotDraw.java
  )

add_jar(GeometricBot SOURCES ${geometricBot})

set(pdeBot
  ${CMAKE_CURRENT_SOURCE_DIR}/vcn/pdeBot/Material.java
  ${CMAKE_CURRENT_SOURCE_DIR}/vcn/pdeBot/BoundaryConditions.java
  ${CMAKE_CURRENT_SOURCE_DIR}/vcn/pdeBot/finiteElement/MeshResults.java
  ${CMAKE_CURRENT_SOURCE_DIR}/vcn/pdeBot/finiteElement/solidMechanics/StaticElasticity2D.java
  )

add_jar(PdeBot SOURCES ${pdeBot})