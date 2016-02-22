cmake_minimum_required (VERSION 2.8.6)
project (Numerical_Bots)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
  "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules")

set(JDK_INCLUDE 
  /opt/java/jdk1.8.0_51/include/
  /opt/java/jdk1.8.0_51/include/linux/
  /usr/local/java/jdk1.8.0_51/include/
  /usr/local/java/jdk1.8.0_51/include/linux/
)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/../../include
  ${JDK_INCLUDE}
  )

if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_C_FLAGS 
    "${CMAKE_C_FLAGS} -std=c11 -funroll-loops -ffast-math -Wall -O3 -pedantic")
  message("Additional default options were set for gcc.")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 

message("Configuring Java Native Interface (JNI)")

message("Configuring JNI - Geometric Bot")
add_library(nb_geometric_bot-jni SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometric_bot/jMesh.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometric_bot/jModel.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometric_bot/JNI_jModel.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometric_bot/JNI_jGeometricBot.c
  )
target_link_libraries(nb_geometric_bot-jni
  nb_geometric_bot
  m)

message("Configuring JNI - PDE Bot")
add_library(nb_pde_bot-jni SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometric_bot/jMesh.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometric_bot/jModel.c
  ${CMAKE_CURRENT_SOURCE_DIR}/src/pde_bot/JNI_jPdeBot.c
  )
target_link_libraries(nb_pde_bot-jni
  nb_geometric_bot
  nb_pde_bot
  m)