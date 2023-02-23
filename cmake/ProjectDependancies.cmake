#
# This file is an example of cmakelist for compile/build.
# This file contains the project options.
# author    Joonhee Jo
# date      2023. 02. 02.
# 

# set the library version compitible to the project
set(PROJECT_DEP_VERSION_osqp
	v0.6.2
	CACHE STRING "Version of `osqp` to be fetched."
)

set(PROJECT_DEP_VERSION_qdldl
	v0.1.6
	CACHE STRING "Version of `osqp` to be fetched."
)

set(PROJECT_DEP_VERSION_gtest
	58d77fa8070e8cec2dc1ed015d66b454c8d78850 # release-1.12.1
	CACHE STRING "Version of `gtest` to be fetched."
)

# set the variable as advanced
mark_as_advanced(PROJECT_DEP_VERSION_osqp)
mark_as_advanced(PROJECT_DEP_VERSION_qdldl)
mark_as_advanced(PROJECT_DEP_VERSION_gtest)

include(FetchContent)
include(FindOrFetch)

# We force all the dependencies to be compiled as static libraries.
# TODO(fraromano) Revisit this choice when adding support for install.
set(BUILD_SHARED_LIBS_OLD ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS
	OFF
	CACHE INTERNAL "Build SHARED libraries"
)

# fetch osqp library
findorfetch(
	USE_SYSTEM_PACKAGE
	OFF
	PACKAGE_NAME
	osqp
	LIBRARY_NAME
	osqp
	GIT_REPO
	https://github.com/osqp/osqp.git
	GIT_TAG
	${PROJECT_DEP_VERSION_osqp}
	TARGETS
	osqp
	EXCLUDE_FROM_ALL
)

# fetch qdldl library
findorfetch(
	USE_SYSTEM_PACKAGE
	OFF
	PACKAGE_NAME
	qdldl
	LIBRARY_NAME
	qdldl
	GIT_REPO
	https://github.com/qdldl/qdldl.git
	GIT_TAG
	${PROJECT_DEP_VERSION_qdldl}
	TARGETS
	qdldl
	EXCLUDE_FROM_ALL
)

if(TEST_PROJECT)
	findorfetch(
		USE_SYSTEM_PACKAGE
		OFF
		PACKAGE_NAME
		gtest
		LIBRARY_NAME
		googletest
		GIT_REPO
		https://github.com/google/googletest.git
		GIT_TAG
		${PROJECT_DEP_VERSION_gtest}
		TARGETS
		gtest
		gmock
		gtest_main
		EXCLUDE_FROM_ALL
	)
endif()

# Reset BUILD_SHARED_LIBS to its previous value
set(BUILD_SHARED_LIBS
	${BUILD_SHARED_LIBS_OLD}
	CACHE BOOL "Build project as a shared library" FORCE
)
unset(BUILD_SHARED_LIBS_OLD)
