/*!
\file       dtLowerMatrix0.h
\brief      dtMath, Lower triangular matrix solver class with dynamic memory
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 05. 21
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLOWER_MATRIX0_H_
#define DTMATH_DTLOWER_MATRIX0_H_

#include "dtDefine.h"
#include "dtLowerMatrix.h"

#if defined(_WIN32) || defined(__linux__)
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

template <typename t_type>
class LowerMatrix<0, 0, t_type>
{
public:
    LowerMatrix() {}

    /* For General Lower Triangular Matrix */
    // Solve for 'x' in Lx = b with dynamic memory vector as return
    template <uint16_t row, uint16_t col>
    static Vector<0, t_type> Solve(Matrix<row, col, t_type> &L, Vector<row, t_type> &b, int8_t *isOk = nullptr);
    template <uint16_t row, uint16_t col>
    static Vector<0, t_type> Solve(Matrix<row, col, t_type> &L, Vector<0, t_type> &b, int8_t *isOk = nullptr);
    static Vector<0, t_type> Solve(Matrix<0, 0, t_type> &L, Vector<0, t_type> &b, int8_t *isOk = nullptr);

    // Inverse involving dynamic memory matrix class
    static int8_t Inverse(Matrix<0, 0, t_type> L, Matrix<0, 0, t_type> &invL);
    static Matrix<0, 0, t_type> Inverse(Matrix<0, 0, t_type> L, int8_t *isOk = nullptr);

    /* For Lower Triangular Matrix with Unit Diagonal */
    // Solve for 'x' in Lx = b with dynamic memory vector as return
    template <uint16_t row, uint16_t col>
    static Vector<0, t_type> SolveUnit(Matrix<row, col, t_type> &L, Vector<row, t_type> &b, int8_t *isOk = nullptr);
    template <uint16_t row, uint16_t col>
    static Vector<0, t_type> SolveUnit(Matrix<row, col, t_type> &L, Vector<0, t_type> &b, int8_t *isOk = nullptr);
    static Vector<0, t_type> SolveUnit(Matrix<0, 0, t_type> &L, Vector<0, t_type> &b, int8_t *isOk = nullptr);

    // Inverse involving dynamic memory matrix class
    static int8_t InverseUnit(Matrix<0, 0, t_type> L, Matrix<0, 0, t_type> &invL);
    static Matrix<0, 0, t_type> InverseUnit(Matrix<0, 0, t_type> L, int8_t *isOk = nullptr);
};

} // namespace Math
} // namespace dt

#include "dtLowerMatrix0.tpp"

#endif // DTMATH_DTLOWER_MATRIX0_H_
