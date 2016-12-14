cmake_minimum_required(VERSION 3.0)
project(NBots)

file(GLOB_RECURSE sources src/*.c)
#file(GLOB_RECURSE sources_test src/test/*.cpp)

list(REMOVE_ITEM sources src/nb/graph_bot/imported_libs/*) # Temporal

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(include)
add_library(nbots ${sources})
target_compile_options(nbots PUBLIC -std=c99 -Wall)
file(COPY ${data} DESTINATION res)

install(TARGETS nbots DESTINATION /home/victor/Desktop)
install(DIRECTORY res DESTINATION /home/victor/Desktop)

set(CPACK_PACKAGE_NAME "nbots")
set(CPACK_PACKAGE_VERSION "0.3.0")
set(CPACK_MONOLITHIC_INSTALL 1)
include(CPack)