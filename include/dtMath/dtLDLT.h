/*!
\file       dtLDLT.h
\brief      dtMath, Cholesky decomposition(L*D*L^T form) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLDLT_H_
#define DTMATH_DTLDLT_H_

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
class LDLT
{
private:
    t_type m_elem[t_row * t_col];
    int8_t m_isOk;

public:
    LDLT();
    LDLT(const t_type *element, const size_t n_byte);
    LDLT(const Matrix<t_row, t_col, t_type> &m);
    LDLT(const Matrix3<t_type, t_row, t_col> &m);

    int8_t Compute();                                           // Compute Cholesky Decomposition, L*D*L^T form
    int8_t Compute(const t_type *element, const size_t n_byte); // Compute Cholesky Decomposition, L*D*L^T form
    int8_t Compute(const Matrix<t_row, t_col, t_type> &m);      // Compute Cholesky Decomposition, L*D*L^T form
    int8_t Compute(const Matrix3<t_type, t_row, t_col> &m);     // Compute Cholesky Decomposition, L*D*L^T form
    int8_t IsOk() { return m_isOk; }

    Matrix<t_row, t_col, t_type> GetMatrix() const;  // return matrix A including L/D/U matrix
    Matrix<t_row, t_col, t_type> GetMatrixL() const; // return Lower Triangular matrix
    Matrix<t_row, t_col, t_type> GetMatrixD() const; // return Diagonal matrix
    Matrix<t_row, t_col, t_type> GetMatrixU() const; // return Upper Triangular matrix

    template <uint16_t col>
    int8_t Solve(const Matrix<t_row, col, t_type> &b, Matrix<t_col, col, t_type> &x); // Solve x = (LDU)^-1 * b
    int8_t Solve(const Vector<t_row, t_type> &b, Vector<t_col, t_type> &x);           // Solve x = (LDU)^-1 * b
    template <uint16_t col>
    Matrix<t_col, col, t_type> Solve(const Matrix<t_row, col, t_type> &b, int8_t *isOk = nullptr); // Solve x = (LDU)^-1 * b
    Vector<t_col, t_type> Solve(const Vector<t_row, t_type> &b, int8_t *isOk = nullptr);           // Solve x = (LDU)^-1 * b

    int8_t Inverse(Matrix<t_row, t_col, t_type> &inv);            // Inverse matrix of LDU matrix
    int8_t Inverse(Matrix3<t_type, t_row, t_col> &inv);           // Inverse matrix of LDU matrix
    Matrix<t_row, t_col, t_type> Inverse(int8_t *isOk = nullptr); // Inverse matrix of LDU matrix

    int8_t InverseArray(t_type *inv);             // Inverse array of LDU matrix
    t_type *InverseArray(int8_t *isOk = nullptr); // Inverse array of LDU matrix
};

} // namespace Math
} // namespace dt

#include "dtLDLT.tpp"
#include "dtLDLT0.h"

#endif // DTMATH_DTLDLT_H_
