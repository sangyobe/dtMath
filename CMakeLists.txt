#
# The root cmake script for dtMath, the C++ linear algebra library.
#
# date      2024. 07. 08.
#
# Copyrignt RoboticsLab ART All rights reserved. 
# 
########################################################

cmake_minimum_required(VERSION 3.13)

project(dtMath
    LANGUAGES CXX
	VERSION 1.2.0
	DESCRIPTION "A math library."
	HOMEPAGE_URL "https://gitlabee.hmg-corp.io/rlab/art/ctrlpart/project/arch/03-art-framework/dtmath.git"
)
set(PROJECT_NAMESPACE ${PROJECT_NAME})

cmake_policy(VERSION 3.13)
cmake_policy(SET CMP0077 NEW)

# Add CMake modules
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")


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