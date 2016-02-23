cmake_minimum_required (VERSION 2.8.6)
project (Numerical_Bots)

find_package(Java)
include(UseJava)

set(CMAKE_JNI_TARGET TRUE)

message("Configuring (Java) NBGeometricBot")

set(geometricBot
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/geometricBot/Box.java
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/geometricBot/Mesh.java
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/geometricBot/Model.java
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/geometricBot/ModelStatus.java
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/geometricBot/GeometricBot.java
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/geometricBot/GeometricBotDraw.java
  )

add_jar(GeometricBot SOURCES ${geometricBot})

set(pdeBot
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/pdeBot/Material.java
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/pdeBot/BoundaryConditions.java
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/pdeBot/BCType.java
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/pdeBot/finiteElement/MeshResults.java
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/pdeBot/finiteElement/MeshResultsDraw.java
  ${CMAKE_CURRENT_SOURCE_DIR}/nb/pdeBot/finiteElement/solidMechanics/StaticElasticity2D.java
  )

add_jar(PdeBot SOURCES ${pdeBot})