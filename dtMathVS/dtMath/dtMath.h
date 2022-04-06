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

#include "dtCommaInit.hpp"

#include "dtVector3.hpp"
#include "dtVector4.hpp"
#include "dtVector6.hpp"
#include "dtMatrix3.hpp"

#include "dtQuaternion.hpp"
#include "dtRotation.hpp"
#include "dtTransform.hpp"

#include "dtVector.hpp"
#include "dtMatrix.hpp"
#include "dtCscMatrix.hpp"

#include "dtLowerTriangular.hpp"// Lower triangular matrix solution
#include "dtUpperTriangular.hpp"// Upper triangular matrix solution
#include "dtNoPivLU.hpp"        // LU Decomposition without pivoting, Doolittle form
#include "dtPartialPivLU.hpp"   // LU Decompostion using partial pivoting and Gaussian-elimination, Doolittle form
#include "dtLLT.hpp"            // Choleski(LLT) Decomposition
#include "dtLDLT.hpp"           // Choleski(LDLT) Decomposition
#include "dtQR.hpp"             // QR Decomposition without pivoting, using Householder refection method
#include "dtSVD.hpp"            // Singular Value Decomposition
////#include "dtFullPivLU.hpp" // To do!
////#include "dtColPivQR.hpp" // To do!

#include "dtQuadProg.hpp"
#include "dtCARE.hpp"           // Continuous Algebraic Riccati Equation Solver
