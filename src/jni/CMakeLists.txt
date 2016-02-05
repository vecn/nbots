cmake_minimum_required (VERSION 2.8.6)
project (VCN_Bots)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} 
  "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules")

set(JDK_INCLUDE 
  /opt/java/jdk1.8.0_51/include/
  /opt/java/jdk1.8.0_51/include/linux/
  /usr/local/java/jdk1.8.0_51/include/
  /usr/local/java/jdk1.8.0_51/include/linux/
)

include_directories(
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
add_library(vcn_geometric_bot-jni SHARED
  ${CMAKE_CURRENT_SOURCE_DIR}/geometric_bot/geometric_bot.c
  )
target_link_libraries(vcn_geometric_bot-jni
  vcn_geometric_bot
  m)