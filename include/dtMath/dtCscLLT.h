/*!
\file      dtCscLLT.h
\brief      dtMath, Cholesky decomposition(L*L^T form) for Sparse Matrix Class
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date      2024. 06. 11
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCSC_LLT_H_
#define DTMATH_DTCSC_LLT_H_

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

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class CscLLT
{
private:
    /* Input Parameters for Sparse Cholesky Factorization */
    int m_elemNum = 0;            // no. of non zero elements
    m_type m_elem[m_row * m_col]; // non zero elements of the sparse matrix
    int m_rowIdx[m_row * m_col];  // row indices
    int m_colPtr[m_col + 1];      // column index pointer
    int m_parent[m_col];          // etree for the LLT Pattern

    /* Output Parameters for Sparse Cholesky Factorization */
    int elemNum_s = 0; // no. of non zero elements
    m_type elem_s[m_row * m_col] = {
        0,
    }; // non zero elements of the sparse matrix
    int rowIdx_s[m_row * m_col] = {
        0,
    }; // row indices
    int colPtr_s[m_col + 1] = {
        0,
    }; // column pointers

    /* Initial Variables for AMD Implementation */
    CscAMD<m_row, m_col, m_type> Order;
    int AMDFlag = 0;
    int *P;
    int *Pinv;

    /* etree and reach functions for Sparse Cholesky Factorization */
    void etree();                                                    // Computing the elimination tree
    int ereach(int k, int *s, int *w, int *colPtr, m_type *x) const; // Computing the reach of the elimination tree

    int8_t Compute();

public:
    CscLLT(int order = 0);                                                                  // Default constructor with AMD or natural ordering; Order = 1: AMD Ordering ; Order = 0 or otherwise: natural ordering
    CscLLT(const CscMatrix<m_row, m_col, m_type> &m, int order = 0);                        // Constructor with sparse matrix and order selection as argument
    CscLLT(const m_type *element, const int elemNum, const int *rowIdx, const int *colPtr); // Constructor with sparse matrix parameters (natural ordering)

    /* Compute Functions for Cholesky Factorization - Compute methods gives the U matrix */

    int8_t Compute(const m_type *element, const int elemNum, const int *rowIdx, const int *colPtr, int flag = 0);
    int8_t Compute(const CscMatrix<m_row, m_col, m_type> &m, int flag = 0);

    /* Getting L and LT matrices in sparse form */
    CscMatrix<m_row, m_col, m_type> GetMatrixL() const; // Return L sparse matrix
    CscMatrix<m_col, m_row, m_type> GetMatrixU() const; // Return U sparse matrix

    /* Solve functions using dense and sparse matrices */
    int8_t Solve(const Vector<m_row, m_type> &b, Vector<m_col, m_type> &x);
    Vector<m_row, m_type> Solve(const Vector<m_row, m_type> &b) const;
};

} // namespace Math
} // namespace dt

#include "dtCscLLT.tpp"

#endif // DTMATH_DTCSC_LLT_H_
