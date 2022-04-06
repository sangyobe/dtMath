/*!
\file       dtContiRiccatiEq.hpp
\brief      dtMath, The continuous-time algebraic Riccati Equation Class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
\see        https://github.com/RobotLocomotion/drake/tree/master/math
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#pragma once

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtMatrix;

// Computes the unique stabilizing solution X to the continuous-time algebraic
// Riccati equation:
// 
// equation: X A + A'X - X B R^{-1} B' X + Q = 0
// where: A(n x n), B(n x m), Q(n x n), R(m x m)
//
// @throws std::runtime_error if R is not positive definite.
//
// Based on the Matrix Sign Function method outlined in this paper:
// http://www.engr.iupui.edu/~skoskie/ECE684/Riccati_algorithms.pdf
//

template <uint16_t m_dimN, uint16_t m_dimM, typename m_type = float>
class CdtCARE
{
private:
    //CdtMatrix<2 * m_dimN, 2 * m_dimN, m_type> m_mH;
    CdtMatrix<2 * m_dimN, 2 * m_dimN, m_type> m_mZ;
    CdtMatrix<2 * m_dimN, 2 * m_dimN, m_type> m_mZold;

    CdtMatrix<m_dimN, m_dimN, m_type> m_mW11;
    CdtMatrix<m_dimN, m_dimN, m_type> m_mW12;
    CdtMatrix<m_dimN, m_dimN, m_type> m_mW21;
    CdtMatrix<m_dimN, m_dimN, m_type> m_mW22;

    CdtMatrix<2 * m_dimN, m_dimN> m_mLhs;
    CdtMatrix<2 * m_dimN, m_dimN> m_mRhs;
    CdtMatrix<m_dimN, m_dimN> m_mEye;

    uint16_t m_maxIteration;
    m_type m_tolerance;
    const m_type m_power = (m_type)(-1) / (2 * m_dimN);

public:
    CdtCARE();
    ~CdtCARE() {}
    void SetSolveCriteria(m_type tolerance, uint16_t maxIteration);
    CdtMatrix<m_dimN, m_dimN, m_type> Solve(
        const CdtMatrix<m_dimN, m_dimN, m_type> &A,
        const CdtMatrix<m_dimN, m_dimM, m_type> &B,
        const CdtMatrix<m_dimN, m_dimN, m_type> &Q,
        const CdtMatrix<m_dimM, m_dimM, m_type> &R);

    int8_t Solve(
        const CdtMatrix<m_dimN, m_dimN, m_type> &A,
        const CdtMatrix<m_dimN, m_dimM, m_type> &B,
        const CdtMatrix<m_dimN, m_dimN, m_type> &Q,
        const CdtMatrix<m_dimM, m_dimM, m_type> &R,
        CdtMatrix<m_dimN, m_dimN, m_type> &X);
};

#include "dtCARE.ipp"