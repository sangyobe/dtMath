/*!
\file       dtNoPivLU.h
\brief      dtMath, LU Decomposition without pivoting(Doolittle form) class, A = LU
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTNO_PIV_LU_H_
#define DTMATH_DTNO_PIV_LU_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
#include <stdint.h>
#include <string.h>
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
template <typename t_type, uint16_t t_row, uint16_t t_col> class Matrix3;

template <uint16_t t_row, uint16_t t_col, typename t_type = float>
class NoPivLU
{
private:
    t_type m_elem[t_row * t_col];
    t_type m_inv[t_row * t_col];
    int8_t m_isOk;
    int8_t Compute(); // Compute Lower/Upper Triangular Matrix, Doolittle form

public:
    NoPivLU();
    NoPivLU(const t_type *element, const size_t n_byte);
    NoPivLU(const Matrix<t_row, t_col, t_type> &m);
    NoPivLU(const Matrix3<t_type, t_row, t_col> &m);

    int8_t Compute(const t_type *element, const size_t n_byte); // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const Matrix<t_row, t_col, t_type> &m);      // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const Matrix3<t_type, t_row, t_col> &m);     // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t IsOk() { return m_isOk; }

    Matrix<t_row, t_col, t_type> GetMatrix() const;  // return matrix A including L and U matrix
    Matrix<t_row, t_col, t_type> GetMatrixL() const; // return Lower Triangular matrix
    Matrix<t_row, t_col, t_type> GetMatrixU() const; // return Upper Triangular matrix

    int8_t Solve(const Vector<t_row, t_type> &b, Vector<t_col, t_type> &x);              // Solve x = (LU)^-1 * b
    Vector<t_col, t_type> Solve(const Vector<t_row, t_type> &b, int8_t *isOk = nullptr); // Solve x = (LU)^-1 * b

    int8_t Inverse(Matrix<t_row, t_col, t_type> &inv);            // Inverse matrix of LU matrix
    int8_t Inverse(Matrix3<t_type, t_row, t_col> &inv);           // Inverse matrix of LU matrix
    Matrix<t_row, t_col, t_type> Inverse(int8_t *isOk = nullptr); // Inverse matrix of LU matrix

    int8_t InverseArray(t_type *inv);             // Inverse array of LU matrix
    t_type *InverseArray(int8_t *isOk = nullptr); // Inverse array of LU matrix
};

} // namespace Math
} // namespace dt

#include "dtNoPivLU.tpp"
#include "dtNoPivLU0.h"

#endif // DTMATH_DTNO_PIV_LU_H_
