set(CORE_SOURCES
    Mesh2D.cpp
    ThermoHydro.cpp
    Assembler.cpp
    NewtonRaphsonSolver.cpp
    VTKWriter.cpp
    Config.cpp
    BlockLoader.cpp
)

set(CORE_HEADERS
    ${CMAKE_SOURCE_DIR}/include/Mesh2D.h
    ${CMAKE_SOURCE_DIR}/include/ThermoHydro.h
    ${CMAKE_SOURCE_DIR}/include/Assembler.h
    ${CMAKE_SOURCE_DIR}/include/NewtonRaphsonSolver.h
    ${CMAKE_SOURCE_DIR}/include/VTKWriter.h
    ${CMAKE_SOURCE_DIR}/include/Config.h
    ${CMAKE_SOURCE_DIR}/include/BlockLoader.h
    ${CMAKE_SOURCE_DIR}/include/Types.h
)

add_library(roadbed_core STATIC ${CORE_SOURCES} ${CORE_HEADERS})
target_link_libraries(roadbed_core
    PUBLIC
        Eigen3::Eigen
        Boost::boost
        Boost::filesystem
        Boost::system
)

target_include_directories(roadbed_core PUBLIC ${CMAKE_SOURCE_DIR}/include)

if(USE_BLOCK_LOADING)
    target_compile_definitions(roadbed_core PUBLIC USE_BLOCK_LOADING)
endif()

add_executable(RoadbedSim main.cpp)
target_link_libraries(RoadbedSim
    PRIVATE
        roadbed_core
        Boost::program_options
)
