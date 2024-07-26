
# Introduction

dtMath is a linear algebra library. It is simple and fast. It is header-only C++ library.

# Installation

dtMath is header-only library. Get source from gitlab repository and install it using the following CMake command.
Before installing dtMath, ARTF_INSTALL_DIR environment variable must be set.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target install
```

# Build user application

To build and test user application, set BUILD_USER cmake option ON.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_USER=ON
cmake --build build --target install
```

Run user application and see the result.
```
./build/user/apptmath
```
Output:
```
  1.000000   2.000000   3.000000 
  4.000000   5.000000   6.000000 
  7.000000   8.000000   9.000000 
  1.000000   4.000000   7.000000 
  2.000000   5.000000   8.000000 
  3.000000   6.000000   9.000000 
```


## Structure

```
.
├── cmake
│    ├── FindOrFetch.cmake
│    ├── ProjectDependancies.cmake
│    ├── ProjectOptions.cmake
│    └── TargetAddRpath.cmake
├── CMakeLists.txt
├── README.md
├── include
│    └── dtMath
│        ├── dtCARE.h
│        ├── dtCARE.tpp
│        ├── dtCommaInit.h
│        ├── dtCommaInit.tpp
│        ├── dtCscMatrix.h
│        ├── dtCscMatrix.tpp
│        ├── dtDefine.h
│        ├── dtLDLT.h
│        ├── dtLDLT.tpp
│        ├── dtLLT.h
│        ├── dtLLT.tpp
│        ├── dtLowerTriangular.h
│        ├── dtLowerTriangular.tpp
│        ├── dtMath.h
│        ├── dtMatrix3.h
│        ├── dtMatrix3.tpp
│        ├── dtMatrix.h
│        ├── dtMatrix.tpp
│        ├── dtNoPivLU.h
│        ├── dtNoPivLU.tpp
│        ├── dtPartialPivLU.h
│        ├── dtPartialPivLU.tpp
│        ├── dtQR.h
│        ├── dtQR.tpp
│        ├── dtQuadProg.h
│        ├── dtQuadProg.tpp
│        ├── dtQuaternion.h
│        ├── dtQuaternion.tpp
│        ├── dtRotation.h
│        ├── dtRotation.tpp
│        ├── dtSVD.h
│        ├── dtSVD.tpp
│        ├── dtTransform.h
│        ├── dtTransform.tpp
│        ├── dtUpperTriangular.h
│        ├── dtUpperTriangular.tpp
│        ├── dtVector3.h
│        ├── dtVector3.tpp
│        ├── dtVector4.h
│        ├── dtVector4.tpp
│        ├── dtVector6.h
│        ├── dtVector6.tpp
│        ├── dtVector.h
│        └── dtVector.tpp
├── test
│   ├── CMakeLists.txt
│   └── user
│        ├── dhTerm.cpp
│        ├── dhTerm.h
│        ├── main.cpp
│        ├── testCompareClass.cpp
│        ├── testCompareClass.h
│        ├── testContiAlgebraicRiccatiEq.cpp
│        ├── testContiAlgebraicRiccatiEq.h
│        ├── testCscMatrix.cpp
│        ├── testCscMatrix.h
│        ├── testDecomposition.cpp
│        ├── testDecomposition.h
│        ├── testEulerAngles.cpp
│        ├── testEulerAngles.h
│        ├── testGnuPlot.cpp
│        ├── testGnuPlot.h
│        ├── testMatrix3.cpp
│        ├── testMatrix3.h
│        ├── testMatrix.cpp
│        ├── testMatrix.h
│        ├── testMatrixInverse.cpp
│        ├── testMatrixInverse.h
│        ├── testMPC.cpp
│        ├── testMPC.h
│        ├── testOriErr.cpp
│        ├── testOriErr.h
│        ├── testPrint.cpp
│        ├── testPrint.h
│        ├── testQuadProg.cpp
│        ├── testQuadProg.h
│        ├── testQuaternion.cpp
│        ├── testQuaternion.h
│        ├── testRotation.cpp
│        ├── testRotation.h
│        ├── testTransform.cpp
│        ├── testTransform.h
│        ├── testVector3.cpp
│        ├── testVector3.h
│        ├── testVector4.cpp
│        ├── testVector4.h
│        ├── testVector6.cpp
│        ├── testVector6.h
│        ├── testVector.cpp
│        └── testVector.h
└── user
     ├── CMakeLists.txt
     └── main.cpp

```
