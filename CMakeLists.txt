cmake_minimum_required(VERSION 3.20)
project(comp_programming_fun CXX)

# no build type means someone forgot -- pick the fast one, not the slow default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/flags.cmake)

# header-only common/ lives at the root -- every problem includes it the same
# way, "common/bench.hpp", no relative-path gymnastics.
include_directories(${CMAKE_SOURCE_DIR})

# `ctest` runs every problem's test binary from one command.
enable_testing()

# problems register themselves. drop a new NNNN-slug/ with its own CMakeLists
# and it joins the build -- no editing this file to add a problem.
file(GLOB problem_dirs ${CMAKE_SOURCE_DIR}/problems/*)
foreach(dir ${problem_dirs})
    if(EXISTS ${dir}/CMakeLists.txt)
        add_subdirectory(${dir})
    endif()
endforeach()
