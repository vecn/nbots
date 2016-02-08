cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

find_package(Java)
include(UseJava)

message("Configuring Java Unit-Tests")

message("-- TestGeometricBot")
add_jar(TestGeometricBot
  SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/TestGeometricBot.java
  INCLUDE_JARS GeometricBot
  ENTRY_POINT TestGeometricBot)

message("-- TestPdeBot")
add_jar(TestPdeBot
  SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/TestPdeBot.java
  INCLUDE_JARS PdeBot GeometricBot
  ENTRY_POINT TestPdeBot)