/*!
\file       dtCARE.h
\brief      dtMath, The Continuous-time Algebraic Riccati Equation Class, cmake version
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\see        https://github.com/RobotLocomotion/drake/tree/master/math
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCARE_H_
#define DTMATH_DTCARE_H_

#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
#include <stdint.h>
#include <stdio.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;

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

template <uint16_t t_dimN, uint16_t t_dimM, typename t_type = float>
class CARE
{
private:
    // Matrix<2 * t_dimN, 2 * t_dimN, t_type> m_mH;
    Matrix<2 * t_dimN, 2 * t_dimN, t_type> m_mZ;
    Matrix<2 * t_dimN, 2 * t_dimN, t_type> m_mZold;

    Matrix<t_dimN, t_dimN, t_type> m_mW11;
    Matrix<t_dimN, t_dimN, t_type> m_mW12;
    Matrix<t_dimN, t_dimN, t_type> m_mW21;
    Matrix<t_dimN, t_dimN, t_type> m_mW22;

    Matrix<2 * t_dimN, t_dimN, t_type> m_mLhs;
    Matrix<2 * t_dimN, t_dimN, t_type> m_mRhs;
    Matrix<t_dimN, t_dimN, t_type> m_mEye;

    uint16_t m_maxIteration;
    t_type m_tolerance;
    const t_type m_power = (t_type)(-1) / (2 * t_dimN);

public:
    CARE();
    ~CARE() {}
    void SetSolveCriteria(t_type tolerance, uint16_t maxIteration);
    Matrix<t_dimN, t_dimN, t_type> Solve(
        const Matrix<t_dimN, t_dimN, t_type> &A,
        const Matrix<t_dimN, t_dimM, t_type> &B,
        const Matrix<t_dimN, t_dimN, t_type> &Q,
        const Matrix<t_dimM, t_dimM, t_type> &R);

    int8_t Solve(
        const Matrix<t_dimN, t_dimN, t_type> &A,
        const Matrix<t_dimN, t_dimM, t_type> &B,
        const Matrix<t_dimN, t_dimN, t_type> &Q,
        const Matrix<t_dimM, t_dimM, t_type> &R,
        Matrix<t_dimN, t_dimN, t_type> &X);
};

} // namespace Math
} // namespace dt

#include "dtCARE.tpp"

#endif // DTMATH_DTCARE_H_