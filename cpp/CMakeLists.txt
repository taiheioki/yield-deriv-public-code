cmake_minimum_required(VERSION 3.30)

project(rcmc
    VERSION 1.0
    DESCRIPTION "Code for yield derivative"
    LANGUAGES CXX
)

cmake_policy(SET CMP0167 NEW)
find_package(Boost CONFIG REQUIRED COMPONENTS program_options timer)
find_package(Eigen3 REQUIRED NO_MODULE)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED gmp gmpxx IMPORTED_TARGET)

add_library(rcmc STATIC src/chemistry.cpp src/io.cpp src/network.cpp src/rate_constant_matrix.cpp src/rcmc.cpp)

target_include_directories(rcmc PUBLIC include)
target_compile_features(rcmc PUBLIC cxx_std_17)
target_compile_options(rcmc PUBLIC -Wall -Wextra)
target_link_libraries(rcmc PUBLIC Boost::program_options Boost::timer Eigen3::Eigen PkgConfig::GMP)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(deriv bin/deriv.cpp)
target_link_libraries(deriv rcmc)

add_executable(pop bin/pop.cpp)
target_link_libraries(pop rcmc)
