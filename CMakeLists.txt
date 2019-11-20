cmake_minimum_required(VERSION 2.6)
project(nsga2)

set(CMAKE_CXX_STANDARD 11)

include_directories ("${PROJECT_SOURCE_DIR}")
include_directories ("${PROJECT_SOURCE_DIR}/Header")
include_directories ("${PROJECT_SOURCE_DIR}/Source")

if(MSVC)
    file(GLOB_RECURSE Header *.h *.hpp)
    file(GLOB_RECURSE Source *.c *.cpp)

    source_group("Header" FILES ${Header})
    source_group("Source" FILES ${Source})
    add_executable(nsga2 ${PROJECT_SOURCE_DIR} ${Header} ${Source})
elseif(UNIX)
    add_executable(nsga2 main.cpp)
endif()