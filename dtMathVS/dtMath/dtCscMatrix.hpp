/*!
\file       dtCscMatrix.hpp
\brief      dtMath, Compressed Sparse Column Matrix (m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2022. 03. 28
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#pragma once

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

template <uint16_t m_row, typename m_type> class CdtVector;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtMatrix;

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class CdtCscMatrix
{
private:
    int m_elemNum = 0;
    m_type m_elem[m_row * m_col];
    int m_rowIdx[m_row * m_col];    // row indices
    int m_colPtr[m_col + 1];        // column index pointer
    CdtCscMatrix(const m_type *element, const int elemNum, const int *rowIdx, const int *colPtr);

public:
    CdtCscMatrix();
    CdtCscMatrix(const CdtMatrix<m_row, m_col, m_type>& m);
    CdtCscMatrix(const CdtCscMatrix& m);
    ~CdtCscMatrix() {}

    const m_type* const GetDataAddr() const;
    const int* const GetRowIdx() const;
    const int* const GetColPtr() const;
    uint16_t GetRowSize() const { return m_row; }   // size of row
    uint16_t GetColSize() const { return m_col; }   // size of colum
    CdtVector<m_col, m_type> GetRowVec(const uint16_t idxRow) const;
    CdtVector<m_row, m_type> GetColVec(const uint16_t idxCol) const;
    int8_t GetRowVec(const uint16_t idxRow, CdtVector<m_col, m_type>& v) const;
    int8_t GetColVec(const uint16_t idxCol, CdtVector<m_row, m_type>& v) const;
    CdtMatrix<m_row, m_col, m_type> GetDenseMat() const;
    CdtCscMatrix<m_col, m_row, m_type> Transpose() const;

    m_type GetNorm() const;
    m_type GetSqNorm() const;

    /* Assignment operators */
    CdtCscMatrix& operator  =(const CdtCscMatrix& m);                          // matrix  = matrix
    CdtCscMatrix& operator *=(const m_type s);                                 // matrix *= scalar
    CdtCscMatrix& operator /=(const m_type s);                                 // matrix /= scalar

    /* Arithmetic operators */
    CdtCscMatrix operator *(const m_type s) const;                              // matrix * scalar
    CdtCscMatrix operator /(const m_type s) const;                              // matrix / scalar

    CdtVector<m_row, m_type> operator *(const CdtVector<m_col, m_type>& v) const;   // matrix * vector
    CdtVector<m_col, m_type> TposeVec(const CdtVector<m_row, m_type>& v) const;     // matrix^T * vector

    void Print(const char endChar = 0);

    /* Friend classes */
    template <uint16_t row, uint16_t col, typename type> friend class CdtCscMatrix;
};

#include "dtCscMatrix.ipp"
