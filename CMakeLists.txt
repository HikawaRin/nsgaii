cmake_minimum_required(VERSION 2.6)
project(nsga2)

set(CMAKE_CXX_STANDARD 11)

include_directories ("${PROJECT_SOURCE_DIR}")

if(MSVC)
    set(CMAKE_EXE_LINKER_FLAGS  "/FORCE")
    file(GLOB_RECURSE Header Header/*.h Header/*.hpp)
    file(GLOB_RECURSE Source Source/*.c Source/*.cpp)

    add_executable(nsga2 "main.cpp" ${Header} ${Source})

    source_group("Header" FILES ${Header})
    source_group("Source" FILES ${Source})
elseif(UNIX)
    add_executable(nsga2 main.cpp)
endif()
