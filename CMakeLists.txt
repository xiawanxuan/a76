cmake_minimum_required(VERSION 3.16)
project(RoadbedFreezeThawSim
    VERSION 1.0.0
    DESCRIPTION "寒区道路工程季节冻融循环耦合仿真"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/W4 /permissive- /EHsc)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
else()
    add_compile_options(-Wall -Wextra -Wpedantic -O2)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

option(BUILD_TESTS "Build unit tests" ON)
option(USE_BLOCK_LOADING "Enable block-wise memory loading for large meshes" ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost 1.70 REQUIRED COMPONENTS program_options filesystem system serialization)
find_package(Boost 1.70 QUIET OPTIONAL_COMPONENTS unit_test_framework)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_subdirectory(src)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

install(DIRECTORY examples/ DESTINATION examples)
