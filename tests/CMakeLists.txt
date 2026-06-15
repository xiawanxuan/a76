cmake_minimum_required(VERSION 3.16)

find_package(Boost 1.70 REQUIRED COMPONENTS unit_test_framework)

add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories(${CMAKE_SOURCE_DIR}/include)

set(TEST_FILES
    test_mesh2d.cpp
    test_thermohydro.cpp
    test_assembler.cpp
    test_solver.cpp
    test_vtkwriter.cpp
    test_config.cpp
)

set(TEST_DATA_DIR ${CMAKE_SOURCE_DIR}/examples)
add_definitions(-DTEST_DATA_DIR="${TEST_DATA_DIR}")

foreach(TEST_FILE ${TEST_FILES})
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_FILE})
    target_link_libraries(${TEST_NAME}
        PRIVATE
            roadbed_core
            Boost::unit_test_framework
            Boost::filesystem
            Boost::system
    )
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endforeach()
