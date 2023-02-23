/*!
\file       dtHouseholderQR.h
\brief      dtMath, QR Decomposition without pivoting(Householder method) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2021. 04. 19
\version    1.0.0
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

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class CdtQR
{
private:
    m_type m_elem[m_row * m_col];
    m_type m_R[m_row*m_col] = { 0, }; // Hn * ... * H2 * H1 * A
    m_type m_Q[m_row*m_row] = { 0, }; // Q = H1 * H2 * ... * Hn
    int8_t m_isOk;

public:
    CdtQR();
    CdtQR(const m_type *element, const size_t n_byte);
    CdtQR(const CdtMatrix<m_row, m_col, m_type> &m);
    CdtQR(const CdtMatrix3<m_type, m_row, m_col> &m);

    int8_t Compute();
    int8_t Compute(const m_type *element, const size_t n_byte); // Compute Q & R Matrix, Using Householder reflections
    int8_t Compute(const CdtMatrix<m_row, m_col, m_type> &m);   // Compute Q & R Matrix, Using Householder reflections
    int8_t Compute(const CdtMatrix3<m_type, m_row, m_col> &m);  // Compute Q & R Matrix, Using Householder reflections
    int8_t IsOk() { return m_isOk; }

    CdtMatrix<m_row, m_row, m_type> GetMatrixQ() const;         // return matrix Q, orthogonal matrix
    CdtMatrix<m_row, m_col, m_type> GetMatrixR() const;         // return matrix R, Upper Triangular matrix
};

#include "dtQR.tpp"

#endif // DTMATH_DTQR_H_
