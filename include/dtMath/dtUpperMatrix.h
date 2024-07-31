/*!
\file       dtUpperMatrix.h
\brief      dtMath, Upper triangular matrix solver class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTUPPER_MATRIX_H_
#define DTMATH_DTUPPER_MATRIX_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
#include <stdint.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

namespace dt
{
namespace Math
{

template <uint16_t t_row, typename t_type> class Vector;
template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;

template <uint16_t t_row, uint16_t t_col, typename t_type = float>
class UpperMatrix
{
public:
    UpperMatrix() {}

    // for General Upper Triangular Matrix
    static int8_t Solve(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, Vector<t_col, t_type> &x);
    static Vector<t_col, t_type> Solve(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, int8_t *isOk = nullptr);
    static int8_t Inverse(Matrix<t_row, t_col, t_type> U, Matrix<t_row, t_col, t_type> &invU);
    static Matrix<t_row, t_col, t_type> Inverse(Matrix<t_row, t_col, t_type> U, int8_t *isOk = nullptr);

    // For Partial or Full Dynamic Memory Allocation

    static int8_t Solve(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, Vector<0, t_type> &x);
    static int8_t Solve(Matrix<t_row, t_col, t_type> &U, Vector<0, t_type> &b, Vector<0, t_type> &x);
    static int8_t Solve(Matrix<0, 0, t_type> &U, Vector<0, t_type> &b, Vector<0, t_type> &x);

    // for Unit Upper Triangular Matrix
    static int8_t SolveUnit(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, Vector<t_col, t_type> &x);
    static Vector<t_col, t_type> SolveUnit(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, int8_t *isOk = nullptr);
    static int8_t InverseUnit(Matrix<t_row, t_col, t_type> U, Matrix<t_row, t_col, t_type> &invU);
    static Matrix<t_row, t_col, t_type> InverseUnit(Matrix<t_row, t_col, t_type> U, int8_t *isOk = nullptr);

    // For Partial or Full Dynamic Memory Allocation

    static int8_t SolveUnit(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, Vector<0, t_type> &x);
    static int8_t SolveUnit(Matrix<t_row, t_col, t_type> &U, Vector<0, t_type> &b, Vector<0, t_type> &x);
    static int8_t SolveUnit(Matrix<0, 0, t_type> &U, Vector<0, t_type> &b, Vector<0, t_type> &x);
};

} // namespace Math
} // namespace dt

#include "dtUpperMatrix.tpp"
#include "dtUpperMatrix0.h"

#endif // DTMATH_DTUPPER_MATRIX_H_
