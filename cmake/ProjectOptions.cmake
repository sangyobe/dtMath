#
# This file is an example of cmakelist for compile/build.
# This file contains the project options.
# author    Joonhee Jo
# date      2023. 02. 10.
# 

# C++ version to use
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Configuration/build type
option(BUILD_SHARED_LIBS "Build shared libraries." ON)

get_property(MultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(MultiConfig)
	if(NOT CMAKE_CONFIGURATION_TYPES)
		set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING
		"Choose the type of builds, options are: Debug Release RelWithDebInfo MinSizeRel. (default: Release;Debug)"
		FORCE)
	endif()
	message(STATUS "Configuration types: ${CMAKE_CONFIGURATION_TYPES}")
else()
	if(NOT CMAKE_BUILD_TYPE)
		set(CMAKE_BUILD_TYPE "Release" CACHE STRING
		"Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel. (default: Release)"
		FORCE)
	endif()
	message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()

# Layout build dir like install dir
include(GNUInstallDirs)
option(BUILD_SHARED_LIBS "Build shared libraries." ON)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_INSTALL_RPATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# for multi-config build system (e.g. Xcode, Ninja Multi-Config)
foreach(OutputConfig IN LISTS CMAKE_CONFIGURATION_TYPES)
	string(TOUPPER ${OutputConfig} OUTPUTCONFIG)
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OutputConfig}/${CMAKE_INSTALL_LIBDIR})
	set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OutputConfig}/${CMAKE_INSTALL_LIBDIR})
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OutputConfig}/${CMAKE_INSTALL_BINDIR})
endforeach()


# RPath : path search for executables or libraries.
option(ENABLE_RPATH "Enable RPath support when installing." ON)
mark_as_advanced(ENABLE_RPATH)

# external libraries by FindOrFetch macro will be downloaded into the following path.
set(EXTERNAL_LIBRARY_FETCH_PATH ${CMAKE_CURRENT_SOURCE_DIR}/extlib)
