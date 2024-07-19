/*!
\file       dtPartialPivLU.h
\brief      dtMath, LU Decomposition with partial pivoting(Doolittle form) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTPARTIAL_PIV_LU0_TPP_
#define DTMATH_DTPARTIAL_PIV_LU0_TPP_

#include "dtMathMem.h"
#include "dtPartialPivLU0.h"

#include <assert.h>

namespace dt
{
namespace Math
{

template <typename t_type>
inline PartialPivLU<0, 0, t_type>::PartialPivLU()
    : m_row(0), m_col(0), m_size(0), m_elem(nullptr), m_inv(nullptr), m_pivot(nullptr), m_isOk(0)
{
    assert(m_row == m_col && "PartialPivLU is only for square (and moreover invertible) matrices");
}

template <typename t_type>
inline PartialPivLU<0, 0, t_type>::PartialPivLU(const uint16_t row, const uint16_t col)
    : m_row(row), m_col(col), m_size(row * col)
{
    assert(row == col && "PartialPivLU is only for square (and moreover invertible) matrices");

    m_elem = MemAlloc<t_type>(m_size);
    m_inv = MemAlloc<t_type>(m_size);
    m_pivot = MemAllocZeroInit<int>(m_row);
}

template <typename t_type>
inline PartialPivLU<0, 0, t_type>::PartialPivLU(const uint16_t row, const uint16_t col, const t_type *element)
    : m_row(row), m_col(col), m_size(row * col)
{
    assert(row == col && "PartialPivLU is only for square (and moreover invertible) matrices");

    m_elem = MemAlloc<t_type>(m_size);
    m_inv = MemAlloc<t_type>(m_size);
    m_pivot = MemAllocZeroInit<int>(m_row);

    memcpy(m_elem, element, sizeof(t_type) * m_size);
    Compute();
}

template <typename t_type>
inline PartialPivLU<0, 0, t_type>::PartialPivLU(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_row == m.m_col && "PartialPivLU is only for square (and moreover invertible) matrices");

    m_row = m.m_row;
    m_col = m.m_col;
    m_size = m_row * m_col;

    m_elem = MemAlloc<t_type>(m_size);
    m_inv = MemAlloc<t_type>(m_size);
    m_pivot = MemAllocZeroInit<int>(m_row);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    Compute();
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline PartialPivLU<0, 0, t_type>::PartialPivLU(const Matrix<row, col, t_type> &m)
    : m_row(row), m_col(col), m_size(row * col)
{
    static_assert(row == col, "PartialPivLU is only for square (and moreover invertible) matrices");

    m_elem = MemAlloc<t_type>(m_size);
    m_inv = MemAlloc<t_type>(m_size);
    m_pivot = MemAllocZeroInit<int>(m_row);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    Compute();
}

template <typename t_type>
inline PartialPivLU<0, 0, t_type>::PartialPivLU(const Matrix3<t_type, 3, 3> &m)
    : m_row(3), m_col(3), m_size(9)
{
    m_elem = MemAlloc<t_type>(9);
    m_inv = MemAlloc<t_type>(9);
    m_pivot = MemAllocZeroInit<int>(3);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * 9);
    Compute();
}

template <typename t_type>
inline PartialPivLU<0, 0, t_type>::~PartialPivLU()
{
    if (m_elem)
    {
        MemFree(m_elem);
        MemFree(m_inv);
        MemFree(m_pivot);
        m_elem = nullptr;
        m_inv = nullptr;
        m_pivot = nullptr;
    }
}

template <typename t_type>
inline int8_t PartialPivLU<0, 0, t_type>::Compute()
{
    uint16_t i, j, k;
    t_type *pMi, *pMk, *p_pivotRow = nullptr;
    t_type max, absElem;
    t_type pivotRow[m_col];

    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        /* Pivoting */
        // find the pivot row
        m_pivot[i] = i;
        max = std::abs(*(pMi + i));
        for (k = i + 1, pMk = pMi + m_col; k < m_row; k++, pMk += m_col)
        {
            if (max < (absElem = std::abs(*(pMk + i))))
            {
                max = absElem;
                m_pivot[i] = k;
                p_pivotRow = pMk; // pMk is pivot row
            }
        }

        // interchange the two rows.
        if (m_pivot[i] != i)
        {
            memcpy(pivotRow, p_pivotRow, sizeof(t_type) * m_col);
            memcpy(p_pivotRow, pMi, sizeof(t_type) * m_col);
            memcpy(pMi, pivotRow, sizeof(t_type) * m_col);
        }

        // matrix is singular, return error
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            m_isOk = 0;
            return -1;
        }

        /* LU Decompostion using Gaussian Elimination */
        for (k = i + 1, pMk = pMi + m_col; k < m_row; pMk += m_col, k++)
        {
            // find the lower triangular matrix elements for column i.
            *(pMk + i) /= *(pMi + i);

            // update the upper triangular matrix for remaining matrix
            for (j = i + 1; j < m_col; j++)
                *(pMk + j) -= *(pMk + i) * *(pMi + j);
        }
    }

    m_isOk = 1;
    return 0;
}

template <typename t_type>
inline int8_t PartialPivLU<0, 0, t_type>::Compute(const uint16_t row, const uint16_t col, const t_type *element)
{
    assert(row == col && "PartialPivLU is only for square (and moreover invertible) matrices");
    assert(row != 0 && "row and col must be non zero");

    if (m_row != row || m_col != col)
    {
        if (m_elem) MemFree(m_elem);
        if (m_inv) MemFree(m_inv);
        if (m_pivot) MemFree(m_pivot);

        m_row = row;
        m_col = col;
        m_size = row * col;
        m_elem = MemAlloc<t_type>(m_size);
        m_inv = MemAlloc<t_type>(m_size);
        m_pivot = MemAlloc<int>(m_row);
    }

    memcpy(m_elem, element, sizeof(t_type) * m_size);
    memset(m_pivot, 0, sizeof(int) * m_row);

    return Compute();
}

template <typename t_type>
inline int8_t PartialPivLU<0, 0, t_type>::Compute(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_row == m.m_col && "PartialPivLU is only for square (and moreover invertible) matrices");
    assert(m.m_row != 0 && "row and col must be non zero");

    if (m_row != m.m_row || m_col != m.m_col)
    {
        if (m_elem) MemFree(m_elem);
        if (m_inv) MemFree(m_inv);
        if (m_pivot) MemFree(m_pivot);

        m_row = m.m_row;
        m_col = m.m_col;
        m_size = m_row * m_col;
        m_elem = MemAlloc<t_type>(m_size);
        m_inv = MemAlloc<t_type>(m_size);
        m_pivot = MemAlloc<int>(m_row);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    memset(m_pivot, 0, sizeof(int) * m_row);

    return Compute();
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t PartialPivLU<0, 0, t_type>::Compute(const Matrix<row, col, t_type> &m)
{
    static_assert(row == col, "PartialPivLU is only for square (and moreover invertible) matrices");

    if (m_row != row || m_col != col)
    {
        if (m_elem) MemFree(m_elem);
        if (m_inv) MemFree(m_inv);
        if (m_pivot) MemFree(m_pivot);

        m_row = row;
        m_col = col;
        m_size = row * col;
        m_elem = MemAlloc<t_type>(m_size);
        m_inv = MemAlloc<t_type>(m_size);
        m_pivot = MemAlloc<int>(m_row);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    memset(m_pivot, 0, sizeof(int) * m_row);

    return Compute();
}

template <typename t_type>
inline int8_t PartialPivLU<0, 0, t_type>::Compute(const Matrix3<t_type, 3, 3> &m)
{
    if (m_row != 3 || m_col != 3)
    {
        if (m_elem) MemFree(m_elem);
        if (m_inv) MemFree(m_inv);
        if (m_pivot) MemFree(m_pivot);

        m_row = 3;
        m_col = 3;
        m_size = 9;
        m_elem = MemAlloc<t_type>(9);
        m_inv = MemAlloc<t_type>(9);
        m_pivot = MemAlloc<int>(3);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * 9);
    memset(m_pivot, 0, sizeof(int) * 3);

    return Compute();
}

template <typename t_type>
inline t_type PartialPivLU<0, 0, t_type>::Determinant()
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    uint16_t offset = m_row + 1;
    t_type det = 1;

    for (uint16_t i = 0; i < m_row; i++)
        det *= m_elem[i * offset];

    return det;
}

template <typename t_type>
inline Matrix<0, 0, t_type> PartialPivLU<0, 0, t_type>::GetMatrix() const
{
    return Matrix<0, 0, t_type>(m_row, m_col, m_elem);
}

template <typename t_type>
inline Matrix<0, 0, t_type> PartialPivLU<0, 0, t_type>::GetMatrixL() const
{
    uint16_t i, j;
    t_type L[m_size];
    memset(L, 0, sizeof(L));

    /* Set diagonal elements as 1 */
    for (i = 0; i < m_row; i++)
    {
        L[i * (m_col + 1)] = 1;
    }

    /* Update remaining matrix from m_elem to L*/
    for (i = 1; i < m_row; i++)
        for (j = 0; j < i; j++)
            L[i * m_col + j] = m_elem[i * m_col + j];

    return Matrix<0, 0, t_type>(m_row, m_col, L);
}

template <typename t_type>
inline Matrix<0, 0, t_type> PartialPivLU<0, 0, t_type>::GetMatrixU() const
{
    uint16_t i, j;
    t_type U[m_size];
    memset(U, 0, sizeof(U));

    for (i = 0; i < m_row; i++)
        for (j = i; j < m_col; j++)
            U[i * m_col + j] = m_elem[i * m_col + j];

    return Matrix<0, 0, t_type>(m_row, m_col, U);
}

template <typename t_type>
inline Matrix<0, 0, t_type> PartialPivLU<0, 0, t_type>::GetMatrixP() const
{
    uint16_t i;
    t_type P[m_size];
    t_type row[m_col];
    memset(P, 0, sizeof(P));

    for (i = 0; i < m_row; i++)
        P[i * (m_col + 1)] = 1;

    for (i = 0; i < m_row; i++)
    {
        if (m_pivot[i] != i)
        {
            memcpy(row, &P[i * m_col], sizeof(t_type) * m_col);
            memcpy(&P[i * m_col], &P[m_pivot[i] * m_col], sizeof(t_type) * m_col);
            memcpy(&P[m_pivot[i] * m_col], row, sizeof(t_type) * m_col);
        }
    }

    return Matrix<0, 0, t_type>(m_row, m_col, P);
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t PartialPivLU<0, 0, t_type>::Solve(const Vector<row, t_type> &b, Vector<col, t_type> &x)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");
    assert(m_col == col && "Col demensions do not match");

    if (!m_isOk)
    {
        assert(false && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, k;
    t_type *pMi;
    t_type tmp;
    t_type vb[row];

    memcpy(vb, b.m_elem, sizeof(vb));

    /* Solve Ly = b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < row; pMi += col, i++)
    {
        if (m_pivot[i] != i)
        {
            tmp = vb[i];
            vb[i] = vb[m_pivot[i]];
            vb[m_pivot[i]] = tmp;
        }

        x.m_elem[i] = vb[i];
        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Ux = y */
    // interchange the row of vector b along the original position
    // Solve the upper triangular (backward substitution)
    for (i = row - 1, pMi = m_elem + (row - 1) * col; i >= 0; i--, pMi -= col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <typename t_type>
inline int8_t PartialPivLU<0, 0, t_type>::Solve(const Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == b.m_row && "Check dimensions");
    assert(m_col == x.m_row && "Check dimensions");
    assert(b.m_elem != nullptr && "Memory has not been allocated");
    assert(x.m_elem != nullptr && "Memory has not been allocated");

    if (!m_isOk)
    {
        assert(false && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, k;
    t_type *pMi;
    t_type tmp;
    t_type vb[m_row];

    memcpy(vb, b.m_elem, sizeof(vb));

    /* Solve Ly = b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        if (m_pivot[i] != i)
        {
            tmp = vb[i];
            vb[i] = vb[m_pivot[i]];
            vb[m_pivot[i]] = tmp;
        }

        x.m_elem[i] = vb[i];
        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Ux = y */
    // interchange the row of vector b along the original position
    // Solve the upper triangular (backward substitution)
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> PartialPivLU<0, 0, t_type>::Solve(const Vector<row, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    Vector<0, t_type> x(m_col);

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return x;
    }

    int i, k;
    t_type *pMi;
    t_type tmp;
    t_type vb[row];

    memcpy(vb, b.m_elem, sizeof(vb));

    /* Solve Ly =b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < row; pMi += m_col, i++)
    {
        if (m_pivot[i] != i)
        {
            tmp = vb[i];
            vb[i] = vb[m_pivot[i]];
            vb[m_pivot[i]] = tmp;
        }

        x.m_elem[i] = vb[i];
        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Ux = y */
    // interchange the row of vector b along the original position
    // Solve the upper triangular (backward substitution)
    for (i = row - 1, pMi = m_elem + (row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x.m_elem[i] /= *(pMi + i);
    }

    if (isOk) *isOk = 1;
    return x;
}

template <typename t_type>
inline Vector<0, t_type> PartialPivLU<0, 0, t_type>::Solve(const Vector<0, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == b.m_row && "Row dimensions do not matched");
    assert(b.m_elem != nullptr && "Memory has not been allocated");

    Vector<0, t_type> x(m_col);

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return x;
    }

    int i, k;
    t_type *pMi;
    t_type tmp;
    t_type vb[m_row];

    memcpy(vb, b.m_elem, sizeof(vb));

    /* Solve Ly =b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        if (m_pivot[i] != i)
        {
            tmp = vb[i];
            vb[i] = vb[m_pivot[i]];
            vb[m_pivot[i]] = tmp;
        }

        x.m_elem[i] = vb[i];
        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Ux = y */
    // interchange the row of vector b along the original position
    // Solve the upper triangular (backward substitution)
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x.m_elem[i] /= *(pMi + i);
    }

    if (isOk) *isOk = 1;
    return x;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t PartialPivLU<0, 0, t_type>::Inverse(Matrix<row, col, t_type> &inv)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");
    assert(m_col == col && "Col demensions do not match");

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    uint16_t colIdx[col]; // column index for interchange the current col with the pivot col
    uint16_t tmpColIdx;
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[row * col]{0};
    t_type invU[row * col]{0};

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < row; i++)
    {
        invL[i * (col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + col;
    p_invMi = invL + col;
    for (i = 1; i < row; i++, p_Mi += col, p_invMi += col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + col;
            for (k = j + 1; k < i; k++, p_invMk += col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < row; k++, p_Mi += (col + 1), p_invMk += (col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + col * (row - 2);
    p_invMi = invU + col * (row - 2);
    for (i = row - 2; i >= 0; i--, p_Mi -= col, p_invMi -= col)
    {
        for (j = col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + col;
            for (k = i + 1; k <= j; k++, p_invMk += col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            if (m_pivot[j] != j)
            {
                tmpColIdx = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tmpColIdx;
            }

            for (k = 0; k < col; k++)
                inv.m_elem[i * col + colIdx[j]] += invU[i * col + k] * invL[k * col + j];

            colIdx[j] = j;
        }
    }

    return 0;
}

template <typename t_type>
inline int8_t PartialPivLU<0, 0, t_type>::Inverse(Matrix<0, 0, t_type> &inv)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == inv.m_row && "Row dimensions do not matched");
    assert(m_col == inv.m_col && "Col demensions do not match");

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    uint16_t colIdx[m_col]; // column index for interchange the current col with the pivot col
    uint16_t tmpColIdx;
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[m_size];
    t_type invU[m_size];

    memset(invL, 0, sizeof(invL));
    memset(invU, 0, sizeof(invU));

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < m_row; i++)
    {
        invL[i * (m_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + m_col;
    p_invMi = invL + m_col;
    for (i = 1; i < m_row; i++, p_Mi += m_col, p_invMi += m_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += m_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + m_col;
            for (k = j + 1; k < i; k++, p_invMk += m_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < m_row; k++, p_Mi += (m_col + 1), p_invMk += (m_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + m_col * (m_row - 2);
    p_invMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, p_Mi -= m_col, p_invMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + m_col;
            for (k = i + 1; k <= j; k++, p_invMk += m_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < m_row; i++)
    {
        for (j = 0; j < m_col; j++)
        {
            if (m_pivot[j] != j)
            {
                tmpColIdx = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tmpColIdx;
            }

            for (k = 0; k < m_col; k++)
                inv.m_elem[i * m_col + colIdx[j]] += invU[i * m_col + k] * invL[k * m_col + j];

            colIdx[j] = j;
        }
    }

    return 0;
}

template <typename t_type>
inline int8_t PartialPivLU<0, 0, t_type>::Inverse(Matrix3<t_type, 3, 3> &inv)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(m_col == 3 && "Col demensions do not match");

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    uint16_t colIdx[3]; // column index for interchange the current col with the pivot col
    uint16_t tmpColIdx;
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[9]{0};
    t_type invU[9]{0};

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < 3; i++)
    {
        invL[i * (3 + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + 3;
    p_invMi = invL + 3;
    for (i = 1; i < 3; i++, p_Mi += 3, p_invMi += 3)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += 3)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + 3;
            for (k = j + 1; k < i; k++, p_invMk += 3)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < 3; k++, p_Mi += 4, p_invMk += 4)
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + 3;
    p_invMi = invU + 3;
    for (i = 1; i >= 0; i--, p_Mi -= 3, p_invMi -= 3)
    {
        for (j = 2; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + 3;
            for (k = i + 1; k <= j; k++, p_invMk += 3)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            if (m_pivot[j] != j)
            {
                tmpColIdx = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tmpColIdx;
            }

            for (k = 0; k < 3; k++)
                inv.m_elem[i * 3 + colIdx[j]] += invU[i * 3 + k] * invL[k * 3 + j];

            colIdx[j] = j;
        }
    }

    return 0;
}

template <typename t_type>
inline Matrix<0, 0, t_type> PartialPivLU<0, 0, t_type>::Inverse(int8_t *isOk)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return Matrix<0, 0, t_type>();
    }

    int i, j, k;
    uint16_t colIdx[m_col]; // column index for interchange the current col with the pivot col
    uint16_t tmpColIdx;
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[m_size];
    t_type invU[m_size];

    memset(invL, 0, sizeof(invL));
    memset(invU, 0, sizeof(invU));

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < m_row; i++)
    {
        invL[i * (m_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + m_col;
    p_invMi = invL + m_col;
    for (i = 1; i < m_row; i++, p_Mi += m_col, p_invMi += m_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += m_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + m_col;
            for (k = j + 1; k < i; k++, p_invMk += m_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < m_row; k++, p_Mi += (m_col + 1), p_invMk += (m_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + m_col * (m_row - 2);
    p_invMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, p_Mi -= m_col, p_invMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + m_col;
            for (k = i + 1; k <= j; k++, p_invMk += m_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    memset(m_inv, 0, sizeof(t_type) * m_size);
    for (i = 0; i < m_row; i++)
    {
        for (j = 0; j < m_col; j++)
        {
            if (m_pivot[j] != j)
            {
                tmpColIdx = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tmpColIdx;
            }

            for (k = 0; k < m_col; k++)
                m_inv[i * m_col + colIdx[j]] += invU[i * m_col + k] * invL[k * m_col + j];

            colIdx[j] = j;
        }
    }

    if (isOk) *isOk = 1;
    return Matrix<0, 0, t_type>(m_row, m_col, m_inv);
}

template <typename t_type>
inline int8_t PartialPivLU<0, 0, t_type>::InverseArray(t_type *inv)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    uint16_t colIdx[m_col]; // column index for interchange the current col with the pivot col
    uint16_t tmpColIdx;
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[m_size];
    t_type invU[m_size];

    memset(invL, 0, sizeof(invL));
    memset(invU, 0, sizeof(invU));

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < m_row; i++)
    {
        invL[i * (m_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + m_col;
    p_invMi = invL + m_col;
    for (i = 1; i < m_row; i++, p_Mi += m_col, p_invMi += m_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += m_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + m_col;
            for (k = j + 1; k < i; k++, p_invMk += m_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < m_row; k++, p_Mi += (m_col + 1), p_invMk += (m_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + m_col * (m_row - 2);
    p_invMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, p_Mi -= m_col, p_invMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + m_col;
            for (k = i + 1; k <= j; k++, p_invMk += m_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < m_row; i++)
    {
        for (j = 0; j < m_col; j++)
        {
            if (m_pivot[j] != j)
            {
                tmpColIdx = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tmpColIdx;
            }

            for (k = 0; k < m_col; k++)
                inv[i * m_col + colIdx[j]] += invU[i * m_col + k] * invL[k * m_col + j];

            colIdx[j] = j;
        }
    }

    return 0;
}

template <typename t_type>
inline t_type *PartialPivLU<0, 0, t_type>::InverseArray(int8_t *isOk)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return nullptr;
    }

    int i, j, k;
    uint16_t colIdx[m_col]; // column index for interchange the current col with the pivot col
    uint16_t tmpColIdx;
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[m_size];
    t_type invU[m_size];

    memset(invL, 0, sizeof(invL));
    memset(invU, 0, sizeof(invU));
    memset(m_inv, 0, sizeof(t_type) * m_size);

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < m_row; i++)
    {
        invL[i * (m_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + m_col;
    p_invMi = invL + m_col;
    for (i = 1; i < m_row; i++, p_Mi += m_col, p_invMi += m_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += m_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + m_col;
            for (k = j + 1; k < i; k++, p_invMk += m_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < m_row; k++, p_Mi += (m_col + 1), p_invMk += (m_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + m_col * (m_row - 2);
    p_invMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, p_Mi -= m_col, p_invMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + m_col;
            for (k = i + 1; k <= j; k++, p_invMk += m_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < m_row; i++)
    {
        for (j = 0; j < m_col; j++)
        {
            if (m_pivot[j] != j)
            {
                tmpColIdx = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tmpColIdx;
            }

            for (k = 0; k < m_col; k++)
                m_inv[i * m_col + colIdx[j]] += invU[i * m_col + k] * invL[k * m_col + j];

            colIdx[j] = j;
        }
    }

    if (isOk) *isOk = 1;
    return m_inv;
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTPARTIAL_PIV_LU0_TPP_