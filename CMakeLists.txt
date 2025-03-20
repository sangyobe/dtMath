#
# The root cmake script for dtMath, the C++ linear algebra library.
#
# date      2024. 07. 08.
#
# Copyrignt RoboticsLab ART All rights reserved.
#
# #######################################################

cmake_minimum_required(VERSION 3.13)

project(dtMath
    LANGUAGES CXX
    VERSION 1.5.0
    DESCRIPTION "A math library."
    HOMEPAGE_URL "https://gitlabee.hmg-corp.io/rlab/art/ctrlpart/project/arch/03-art-framework/dtmath.git"
)
set(PROJECT_NAMESPACE ${PROJECT_NAME})

cmake_policy(VERSION 3.13)
cmake_policy(SET CMP0077 NEW)

# Add CMake modules
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Check if ARTF_INSTALL_DIR is set
file(TO_CMAKE_PATH $ENV{ARTF_INSTALL_DIR} ARTF_INSTALL_DIR)

if("${ARTF_INSTALL_DIR}" STREQUAL "")
    message(FATAL_ERROR "ARTF_INSTALL_DIR is not set.")
else()
    message(STATUS "ARTF_INSTALL_DIR: $ENV{ARTF_INSTALL_DIR}")
endif()

set(CMAKE_PREFIX_PATH "${ARTF_INSTALL_DIR};${CMAKE_PREFIX_PATH}")
set(CMAKE_INSTALL_PREFIX "${ARTF_INSTALL_DIR}")

# --------------------------------------------------------
# Set build options
# --------------------------------------------------------
# Set all project options
# include(ProjectOptions)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build options
option(BUILD_TEST "Build tests for ${PROJECT_NAME}" OFF)
option(BUILD_USER "Build user application" OFF)

# --------------------------------------------------------
# Fecth all dependent library
# --------------------------------------------------------
# include(ProjectDependencies)

# --------------------------------------------------------
# Build
# --------------------------------------------------------
# Tests to build
if(BUILD_TEST)
    enable_testing()
    add_subdirectory(test)
endif()

# User application
if(BUILD_USER)
    add_subdirectory(user)
endif()

# --------------------------------------------------------
# Install
# --------------------------------------------------------
install(DIRECTORY ./include/dtMath
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)

# --------------------------------------------------------
# Finalize configuration and print options
# --------------------------------------------------------
message(STATUS "BUILD_TEST                          : ${BUILD_TEST}")
message(STATUS "BUILD_USER                          : ${BUILD_USER}")