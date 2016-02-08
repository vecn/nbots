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

set(LIB_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_jar(VCNGeometricBot
  SOURCES ${geometricBot})
install_jar(VCNGeometricBot ${LIB_INSTALL_DIR}/VCNGeometricBot)