#
#
# This file is a sample RobotCtrl project embedding mujoco simulation.
# Check whether the model opens up on mujoco.
# author    Joonhee Jo
# email     allusivejune@gmail.com
# date      2023. 02. 20.
#
# 
########################################################

cmake_minimum_required(VERSION 3.16...3.25)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
# Default to GLVND if available.
set(CMAKE_POLICY_DEFAULT_CMP0072 NEW)

# This line has to appear before 'PROJECT' in order to be able to disable incremental linking
set(MSVC_INCREMENTAL_DEFAULT ON)

project(
	user
	VERSION 1.0.0
	DESCRIPTION "${PROJECT_NAME} binaries for a math library."
	HOMEPAGE_URL "https://cody-escm.autoever.com/rlab/art/ctrlpart/personal/joonhee-jo/dt-math-cmake.git"
)
set(PROJECT_NAMESPACE ${PROJECT_NAME})
set(PARENT_PROJECT_NAME ${PROJECT_NAME})

# Add CMake modules
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")


# Set all project options
include(ProjectOptions)

# Set all build options
option(TEST_PROJECT "Build tests for ${PROJECT_NAME}" ON)


# Fecth all dependant library
include(ProjectDependancies)


file(GLOB_RECURSE TARGET_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/user/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/user/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/user/*.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/user/*.ipp
)

add_library(${PROJECT_NAME} STATIC ${TARGET_SRC})

set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")

target_include_directories(${PROJECT_NAME} 
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
		$<INSTALL_INTERFACE:include>
	PRIVATE user
)

target_link_libraries(${PROJECT_NAME} PRIVATE osqp)

add_library(${PROJECT_NAMESPACE}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})


##########################################################
# 5. Tests to build 
# This part can verify this example builds. For your code,
# you'll probably want tests too

if(TEST_PROJECT)
	message(STATUS "Build test: ${BUILD_TEST}")
	add_subdirectory(test)
endif()

