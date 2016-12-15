cmake_minimum_required(VERSION 3.0)
project(NBots)

# Build Library
file(GLOB_RECURSE NB_SOURCES "src/*.c")

set(NB_EXCLUDE_IMPORTED_LIBS false)
if (NB_EXCLUDE_IMPORTED_LIBS)
  message("== NB: Exclude imported libs")
  add_definitions(-DNB_EXCLUDE_IMPORTED_LIBS)
  set(NB_IMPORTED_DIRS "/imported_libs/")
  foreach(NB_SOURCE ${NB_SOURCES})
    string(FIND ${NB_SOURCE} ${NB_IMPORTED_DIRS} NB_EXCLUDE_DIR_FOUND)
    if (NOT ${NB_EXCLUDE_DIR_FOUND} EQUAL -1)
      list(REMOVE_ITEM NB_SOURCES ${NB_SOURCE})
    endif ()
  endforeach(NB_SOURCE)
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(include)

add_library(objlib OBJECT ${NB_SOURCES})
target_compile_options(objlib PUBLIC -std=c99 -Wall)
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1) # -fPIC
add_library(nbots SHARED $<TARGET_OBJECTS:objlib>)
add_library(nbots_static STATIC $<TARGET_OBJECTS:objlib>)
set_target_properties(nbots_static PROPERTIES OUTPUT_NAME nbots)

#file(GLOB_RECURSE sources_test utest/*.c)

#file(COPY ${data} DESTINATION res)

#include(CTest)

#install(TARGETS nbots DESTINATION /home/victor/Desktop)
#install(DIRECTORY res DESTINATION /home/victor/Desktop)

#set(CPACK_PACKAGE_NAME "nbots")
#set(CPACK_PACKAGE_VERSION "0.3.0")
#set(CPACK_MONOLITHIC_INSTALL 1)
#include(CPack)