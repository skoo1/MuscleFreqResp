cmake_minimum_required(VERSION 3.20)
project(muscle LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenSim REQUIRED)

add_executable(muscle_app
    src/main.cpp
    src/MMM_active.cpp
    src/MMM_passive.cpp
    src/TMM_active.cpp
    src/TMM_passive.cpp
    src/Thelen2003MuscleWrapper.cpp
    src/Thelen2003MuscleHelper.cpp
    src/Millard2012EquilibriumMuscleHelper.cpp
)

target_include_directories(muscle_app
    PRIVATE
    ${OpenSim_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(muscle_app
    PRIVATE
    ${OpenSim_LIBRARIES}
)
