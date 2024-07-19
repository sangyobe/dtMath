/*!
\file       dtMath.h
\brief      dtMath
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\date       2020. 10. 21
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#pragma once

#include "dtDefine.h"

#include "dtCommaInit.h"
#include "dtMathMem.h"
#include "dtStorage.h"

#include "dtMatrix3.h"
#include "dtVector3.h"
#include "dtVector4.h"
#include "dtVector6.h"

#include "dtQuaternion.h"
#include "dtRotation.h"
#include "dtTransform.h"

#include "dtCscMatrix.h"
#include "dtMatrix.h"
#include "dtVector.h"

#include "dtFullPivLU.h"    // LU Decomposition with complete pivoting, Doolittle form
#include "dtLDLT.h"         // Choleski(LDLT) Decomposition
#include "dtLLT.h"          // Choleski(LLT) Decomposition
#include "dtLowerMatrix.h"  // Lower triangular matrix solution
#include "dtNoPivLU.h"      // LU Decomposition without pivoting, Doolittle form
#include "dtPartialPivLU.h" // LU Decompostion using partial pivoting and Gaussian-elimination, Doolittle form
#include "dtQR.h"           // QR Decomposition without pivoting, using Householder refection method
#include "dtSVD.h"          // Singular Value Decomposition
#include "dtUpperMatrix.h"  // Upper triangular matrix solution
////#include "dtColPivQR.h" // To do!

#include "dtCscAMD.h" // Approximate Minimum Degree (AMD) for sparse matrix
#include "dtCscLLT.h" // Cholesky (LLT) Decomposition for sparse matrices

#include "dtCARE.h" // Continuous Algebraic Riccati Equation Solver
#include "dtQuadProg.h"
