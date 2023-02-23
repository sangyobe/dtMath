/*!
\file       dtNoPivLU.h
\brief      dtMath, LU Decomposition without pivoting(Doolittle form) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTNO_PIV_LU_H_
#define DTMATH_DTNO_PIV_LU_H_

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
class CdtNoPivLU
{
private:
    m_type m_elem[m_row * m_col];
    m_type m_inv[m_row * m_col];
    int8_t m_isOk;

public:
    CdtNoPivLU();
    CdtNoPivLU(const m_type *element, const size_t n_byte);
    CdtNoPivLU(const CdtMatrix<m_row, m_col, m_type> &m);
    CdtNoPivLU(const CdtMatrix3<m_type, m_row, m_col> &m);

    int8_t Compute();                                           // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const m_type *element, const size_t n_byte); // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const CdtMatrix<m_row, m_col, m_type> &m);   // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const CdtMatrix3<m_type, m_row, m_col> &m);  // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t IsOk() { return m_isOk; }

    CdtMatrix<m_row, m_col, m_type> GetMatrix() const;          // return matrix A including L and U matrix
    CdtMatrix<m_row, m_col, m_type> GetMatrixL() const;         // return Lower Triangular matrix
    CdtMatrix<m_row, m_col, m_type> GetMatrixU() const;         // return Upper Triangular matrix

    int8_t Solve(const CdtVector<m_row, m_type> &b, CdtVector<m_col, m_type> &x);             // Solve x = (LU)^-1 * b
    CdtVector<m_col, m_type> Solve(const CdtVector<m_row, m_type> &b, int8_t *isOk = nullptr);// Solve x = (LU)^-1 * b

    int8_t Inverse(CdtMatrix<m_row, m_col, m_type> &inv);           // Inverse matrix of LU matrix
    int8_t Inverse(CdtMatrix3<m_type, m_row, m_col> &inv);          // Inverse matrix of LU matrix
    CdtMatrix<m_row, m_col, m_type> Inverse(int8_t *isOk = nullptr);// Inverse matrix of LU matrix

    int8_t InverseArray(m_type *inv);                               // Inverse array of LU matrix
    m_type* InverseArray(int8_t *isOk = nullptr);                   // Inverse array of LU matrix
};

#include "dtNoPivLU.tpp"

#endif // DTMATH_DTNO_PIV_LU_H_
