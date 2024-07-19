/*!
\file       dtHouseholderQR.h
\brief      dtMath, QR Decomposition without pivoting(Householder method) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTQR_H_
#define DTMATH_DTQR_H_

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

template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Matrix3;

template <uint16_t t_row, uint16_t t_col, typename t_type = float>
class QR
{
private:
    t_type m_elem[t_row * t_col];
    t_type m_R[t_row * t_col]; // Hn * ... * H2 * H1 * A
    t_type m_Q[t_row * t_row]; // Q = H1 * H2 * ... * Hn
    int8_t m_isOk;

public:
    QR();
    QR(const t_type *element, const size_t n_byte);
    QR(const Matrix<t_row, t_col, t_type> &m);
    QR(const Matrix3<t_type, t_row, t_col> &m);

    int8_t Compute();
    int8_t Compute(const t_type *element, const size_t n_byte); // Compute Q & R Matrix, Using Householder reflections
    int8_t Compute(const Matrix<t_row, t_col, t_type> &m);      // Compute Q & R Matrix, Using Householder reflections
    int8_t Compute(const Matrix3<t_type, t_row, t_col> &m);     // Compute Q & R Matrix, Using Householder reflections
    int8_t IsOk() { return m_isOk; }

    Matrix<t_row, t_row, t_type> GetMatrixQ() const; // return matrix Q, orthogonal matrix
    Matrix<t_row, t_col, t_type> GetMatrixR() const; // return matrix R, Upper Triangular matrix
};

} // namespace Math
} // namespace dt

#include "dtQR.tpp"

#endif // DTMATH_DTQR_H_
