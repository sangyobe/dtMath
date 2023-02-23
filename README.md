


# Introduction

This project is a converted project of dtMath_VS to use cmake. You can clone or fetch this library.



## Structure

```
.
├── cmake
│   ├── FindOrFetch.cmake
│   ├── ProjectDependancies.cmake
│   ├── ProjectOptions.cmake
│   └── TargetAddRpath.cmake
├── CMakeLists.txt
├── doc
│   ├── notes on quaternions_Simo Sarkka.pdf
│   ├── quaternion computation_Neil Dantam.pdf
│   └── quaternion kinematics for the error-state kalman filter_Joan Sola.pdf
├── extlib
│   ├── googletest
│   ├── osqp
│   ├── osqp.dll
│   ├── osqp.lib
│   ├── qdldl.dll
│   └── qdldl.lib
├── README.md
├── test
│   ├── CMakeLists.txt
│   └── user
│       ├── dhTerm.cpp
│       ├── dhTerm.h
│       ├── main.cpp
│       ├── testCompareClass.cpp
│       ├── testCompareClass.h
│       ├── testContiAlgebraicRiccatiEq.cpp
│       ├── testContiAlgebraicRiccatiEq.h
│       ├── testCscMatrix.cpp
│       ├── testCscMatrix.h
│       ├── testDecomposition.cpp
│       ├── testDecomposition.h
│       ├── testEulerAngles.cpp
│       ├── testEulerAngles.h
│       ├── testGnuPlot.cpp
│       ├── testGnuPlot.h
│       ├── testMatrix3.cpp
│       ├── testMatrix3.h
│       ├── testMatrix.cpp
│       ├── testMatrix.h
│       ├── testMatrixInverse.cpp
│       ├── testMatrixInverse.h
│       ├── testMPC.cpp
│       ├── testMPC.h
│       ├── testOriErr.cpp
│       ├── testOriErr.h
│       ├── testPrint.cpp
│       ├── testPrint.h
│       ├── testQuadProg.cpp
│       ├── testQuadProg.h
│       ├── testQuaternion.cpp
│       ├── testQuaternion.h
│       ├── testRotation.cpp
│       ├── testRotation.h
│       ├── testTransform.cpp
│       ├── testTransform.h
│       ├── testVector3.cpp
│       ├── testVector3.h
│       ├── testVector4.cpp
│       ├── testVector4.h
│       ├── testVector6.cpp
│       ├── testVector6.h
│       ├── testVector.cpp
│       └── testVector.h
└── user
    ├── dtMath
    │   ├── dtCARE.h
    │   ├── dtCARE.tpp
    │   ├── dtCommaInit.h
    │   ├── dtCommaInit.tpp
    │   ├── dtCscMatrix.h
    │   ├── dtCscMatrix.tpp
    │   ├── dtDefine.h
    │   ├── dtLDLT.h
    │   ├── dtLDLT.tpp
    │   ├── dtLLT.h
    │   ├── dtLLT.tpp
    │   ├── dtLowerTriangular.h
    │   ├── dtLowerTriangular.tpp
    │   ├── dtMath.h
    │   ├── dtMatrix3.h
    │   ├── dtMatrix3.tpp
    │   ├── dtMatrix.h
    │   ├── dtMatrix.tpp
    │   ├── dtNoPivLU.h
    │   ├── dtNoPivLU.tpp
    │   ├── dtPartialPivLU.h
    │   ├── dtPartialPivLU.tpp
    │   ├── dtQR.h
    │   ├── dtQR.tpp
    │   ├── dtQuadProg.h
    │   ├── dtQuadProg.tpp
    │   ├── dtQuaternion.h
    │   ├── dtQuaternion.tpp
    │   ├── dtRotation.h
    │   ├── dtRotation.tpp
    │   ├── dtSVD.h
    │   ├── dtSVD.tpp
    │   ├── dtTransform.h
    │   ├── dtTransform.tpp
    │   ├── dtUpperTriangular.h
    │   ├── dtUpperTriangular.tpp
    │   ├── dtVector3.h
    │   ├── dtVector3.tpp
    │   ├── dtVector4.h
    │   ├── dtVector4.tpp
    │   ├── dtVector6.h
    │   ├── dtVector6.tpp
    │   ├── dtVector.h
    │   └── dtVector.tpp
    ├── robotics
    │   ├── djBezier.h
    │   ├── djBezierOri.h
    │   ├── djBezierOri.tpp
    │   ├── djBezier.tpp
    │   ├── dtBezier.h
    │   ├── dtBezier.tpp
    │   ├── dtGaitCPG.h
    │   ├── dtGaitCPG.tpp
    │   ├── dtGaitPlanner.h
    │   ├── dtGaitPlanner.tpp
    │   ├── dtStanceBezier.h
    │   ├── dtStanceBezier.tpp
    │   ├── dtSwingBezier.h
    │   └── dtSwingBezier.tpp
    └── utility
        ├── dhTimeCheck.cpp
        ├── dhTimeCheck.h
        ├── dtGnuPlot.h
        └── dtGnuPlot.tpp
```
