# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/robotics/Downloads/_JhJo/1_works/00_Joonhee Jo/dt-math-cmake/extlib/osqp/osqp-src"
  "/home/robotics/Downloads/_JhJo/1_works/00_Joonhee Jo/dt-math-cmake/extlib/osqp/osqp-build"
  "/home/robotics/Downloads/_JhJo/1_works/00_Joonhee Jo/dt-math-cmake/extlib/osqp/osqp-subbuild/osqp-populate-prefix"
  "/home/robotics/Downloads/_JhJo/1_works/00_Joonhee Jo/dt-math-cmake/extlib/osqp/osqp-subbuild/osqp-populate-prefix/tmp"
  "/home/robotics/Downloads/_JhJo/1_works/00_Joonhee Jo/dt-math-cmake/extlib/osqp/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp"
  "/home/robotics/Downloads/_JhJo/1_works/00_Joonhee Jo/dt-math-cmake/extlib/osqp/osqp-subbuild/osqp-populate-prefix/src"
  "/home/robotics/Downloads/_JhJo/1_works/00_Joonhee Jo/dt-math-cmake/extlib/osqp/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/robotics/Downloads/_JhJo/1_works/00_Joonhee Jo/dt-math-cmake/extlib/osqp/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/robotics/Downloads/_JhJo/1_works/00_Joonhee Jo/dt-math-cmake/extlib/osqp/osqp-subbuild/osqp-populate-prefix/src/osqp-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
