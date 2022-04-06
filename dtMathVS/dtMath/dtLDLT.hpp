/*!
\file       dtLDLT.hpp
\brief      dtMath, Cholesky decomposition(L*D*L^T form) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

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

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class CdtLDLT
{
private:
    m_type m_elem[m_row * m_col];
    int8_t m_isOk;

public:
    CdtLDLT();
    CdtLDLT(const m_type* element, const size_t n_byte);
    CdtLDLT(const CdtMatrix<m_row, m_col, m_type> &m);
    CdtLDLT(const CdtMatrix3<m_type, m_row, m_col> &m);

    int8_t Compute();                                           // Compute Cholesky Decomposition, L*D*L^T form
    int8_t Compute(const m_type* element, const size_t n_byte); // Compute Cholesky Decomposition, L*D*L^T form
    int8_t Compute(const CdtMatrix<m_row, m_col, m_type> &m);   // Compute Cholesky Decomposition, L*D*L^T form
    int8_t Compute(const CdtMatrix3<m_type, m_row, m_col> &m);  // Compute Cholesky Decomposition, L*D*L^T form
    int8_t IsOk() { return m_isOk; }

    CdtMatrix<m_row, m_col, m_type> GetMatrix() const;  // return matrix A including L/D/U matrix
    CdtMatrix<m_row, m_col, m_type> GetMatrixL() const; // return Lower Triangular matrix
    CdtMatrix<m_row, m_col, m_type> GetMatrixD() const; // return Diagonal matrix
    CdtMatrix<m_row, m_col, m_type> GetMatrixU() const; // return Upper Triangular matrix

    template <uint16_t col>
    int8_t Solve(const CdtMatrix<m_row, col, m_type> &b, CdtMatrix<m_col, col, m_type> &x);   // Solve x = (LDU)^-1 * b
    int8_t Solve(const CdtVector<m_row, m_type> &b, CdtVector<m_col, m_type> &x);             // Solve x = (LDU)^-1 * b
    template <uint16_t col>
    CdtMatrix<m_col, col, m_type> Solve(const CdtMatrix<m_row, col, m_type> &b, int8_t *isOk = nullptr);  // Solve x = (LDU)^-1 * b
    CdtVector<m_col, m_type> Solve(const CdtVector<m_row, m_type> &b, int8_t *isOk = nullptr);            // Solve x = (LDU)^-1 * b

    int8_t Inverse(CdtMatrix<m_row, m_col, m_type> &inv);           // Inverse matrix of LDU matrix
    int8_t Inverse(CdtMatrix3<m_type, m_row, m_col> &inv);          // Inverse matrix of LDU matrix
    CdtMatrix<m_row, m_col, m_type> Inverse(int8_t *isOk = nullptr);// Inverse matrix of LDU matrix

    int8_t InverseArray(m_type *inv);               // Inverse array of LDU matrix
    m_type* InverseArray(int8_t *isOk = nullptr);   // Inverse array of LDU matrix
};

#include "dtLDLT.ipp"
