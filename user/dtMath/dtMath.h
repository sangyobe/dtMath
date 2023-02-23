/*!
\file       dtMath.h
\brief      dtMath
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#pragma once

#include "dtDefine.h"

#include "dtCommaInit.h"

#include "dtVector3.h"
#include "dtVector4.h"
#include "dtVector6.h"
#include "dtMatrix3.h"

#include "dtQuaternion.h"
#include "dtRotation.h"
#include "dtTransform.h"

#include "dtVector.h"
#include "dtMatrix.h"
#include "dtCscMatrix.h"

#include "dtLowerTriangular.h"// Lower triangular matrix solution
#include "dtUpperTriangular.h"// Upper triangular matrix solution
#include "dtNoPivLU.h"        // LU Decomposition without pivoting, Doolittle form
#include "dtPartialPivLU.h"   // LU Decompostion using partial pivoting and Gaussian-elimination, Doolittle form
#include "dtLLT.h"            // Choleski(LLT) Decomposition
#include "dtLDLT.h"           // Choleski(LDLT) Decomposition
#include "dtQR.h"             // QR Decomposition without pivoting, using Householder refection method
#include "dtSVD.h"            // Singular Value Decomposition
////#include "dtFullPivLU.h" // To do!
////#include "dtColPivQR.h" // To do!

#include "dtQuadProg.h"
#include "dtCARE.h"           // Continuous Algebraic Riccati Equation Solver
