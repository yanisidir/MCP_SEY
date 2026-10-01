cmake_minimum_required(VERSION 3.16)
project(MCP_SEY LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Geant4 REQUIRED ui_all vis_all)

add_executable(MCP_SEY
    MCP_SEY.cc
    src/ActionInitialization.cc
    src/DetectorConstruction.cc
    src/EventAction.cc
    src/FurmanPiviModel.cc
    src/FurmanPiviPhysics.cc
    src/FurmanPiviProcess.cc
    src/GaussianSource.cc
    src/PrimaryGeneratorAction.cc
    src/RunAction.cc
    src/SteppingAction.cc
    src/TrackingAction.cc
    src/WuModel.cc)

target_include_directories(MCP_SEY PRIVATE include)
target_compile_features(MCP_SEY PRIVATE cxx_std_17)
target_link_libraries(MCP_SEY PRIVATE ${Geant4_LIBRARIES})

# Macros disponibles a cote de l'executable ; recopiees a chaque edition.
file(GLOB MCP_SEY_MACROS ${CMAKE_CURRENT_SOURCE_DIR}/macros/*.mac)
foreach(macro ${MCP_SEY_MACROS})
    get_filename_component(name ${macro} NAME)
    configure_file(${macro} ${CMAKE_CURRENT_BINARY_DIR}/macros/${name} COPYONLY)
endforeach()
