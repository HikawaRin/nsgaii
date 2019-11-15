cmake_minimum_required(VERSION 2.6)
project(nsga2)

set(CMAKE_CXX_STANDARD 11)

include_directories(Header)
include_directories(Source)
include_directories(Input)
include_directories(Include)

AUX_SOURCE_DIRECTORY(./ DIR_SRCS)
add_executable(nsga2 ${DIR_SRCS})
