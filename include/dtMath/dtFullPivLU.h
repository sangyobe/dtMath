/*!
\file       dtFullPivLU.hpp
\brief      dtMath, LU Decomposition with complete pivoting(Doolittle form) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\date       2023. 3. 21
\version    1.0.0
\warning    Do Not delete this comment for document history!
*/

#ifndef DTMATH_DTFULL_PIV_LU_H_
#define DTMATH_DTFULL_PIV_LU_H_

#pragma once

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
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
class FullPivLU
{
private:
    t_type m_elem[t_row * t_col];
    // t_type m_inv[t_row * t_col];

    // row pivot vector, col pivot vector, permutation vector
    int m_pivot[t_row], m_pivotCol[t_col], m_permCol[t_row];
    int8_t m_isOk;
    int8_t Compute(); // Compute Lower/Upper Triangular Matrix, Doolittle form

public:
    FullPivLU();                                           // Default Constructor
    FullPivLU(const t_type *element, const size_t n_byte); // Constructor with Array Input and Array Size
    FullPivLU(const Matrix<t_row, t_col, t_type> &m);      // Constructor with Matrix Argument
    FullPivLU(const Matrix3<t_type, t_row, t_col> &m);     // Constructor with Matrix3 Argument

    int8_t Compute(const t_type *element, const size_t n_byte); // Compute with Array Input and Array Size
    int8_t Compute(const Matrix<t_row, t_col, t_type> &m);      // Compute with Matrix Argument
    int8_t Compute(const Matrix3<t_type, t_row, t_col> &m);     // Compute with Matrix3 Argument

    int8_t IsOk() { return m_isOk; }

    Matrix<t_row, t_col, t_type> GetMatrix() const;  // return matrix A including L/U matrix
    Matrix<t_row, t_col, t_type> GetMatrixL() const; // return Lower Triangular matrix
    Matrix<t_row, t_col, t_type> GetMatrixU() const; // return Upper Triangular matrix
    Matrix<t_row, t_col, t_type> GetMatrixP() const; // return Permutation matrix P
    Matrix<t_row, t_col, t_type> GetMatrixQ() const; // return Permutation matrix Q

    int8_t Solve(const Vector<t_row, t_type> &b, Vector<t_col, t_type> &x);              // Solve x = (LU)^-1 * b
    Vector<t_col, t_type> Solve(const Vector<t_row, t_type> &b, int8_t *isOk = nullptr); // Solve x = (LU)^-1 * b

    int8_t Inverse(Matrix<t_row, t_col, t_type> &inv);            // Inverse matrix of LU matrix
    int8_t Inverse(Matrix3<t_type, t_row, t_col> &inv);           // Inverse matrix of LU matrix
    Matrix<t_row, t_col, t_type> Inverse(int8_t *isOk = nullptr); // Inverse matrix of LU matrix
};

} // namespace Math
} // namespace dt

#include "dtFullPivLU.tpp"
#include "dtFullPivLU0.h"

#endif // DTMATH_DTFULL_PIV_LU_H_
