

if(NOT TEST_PROJECT)
	return()
endif()

cmake_minimum_required(VERSION 3.16)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
# Default to GLVND if available.
set(CMAKE_POLICY_DEFAULT_CMP0072 NEW)

# This line has to appear before 'PROJECT' in order to be able to disable incremental linking
set(MSVC_INCREMENTAL_DEFAULT ON)

project(
    test
    VERSION 1.0.0
    DESCRIPTION "test main file for main integraion"
    HOMEPAGE_URL ""
)



file(GLOB_RECURSE EXECUTE_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/user/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/user/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/user/*.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/user/*.ipp
)

add_executable(${PROJECT_NAME} 
    ${EXECUTE_SRC}
)

target_link_libraries(${PROJECT_NAME}
    ${PARENT_PROJECT_NAME}
)

target_include_directories(${PROJECT_NAME} PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/../user
)