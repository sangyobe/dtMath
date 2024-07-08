/*!
\file       dtPartialPivLU.h
\brief      dtMath, LU Decomposition with partial pivoting(Doolittle form) class, PA = LU
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTPARTIAL_PIV_LU_TPP_
#define DTMATH_DTPARTIAL_PIV_LU_TPP_

#include "dtPartialPivLU.h"

#include <cassert>

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline PartialPivLU<t_row, t_col, t_type>::PartialPivLU() : m_elem(), m_inv(), m_pivot(), m_isOk(0)
{
    static_assert(t_row == t_col, "PartialPivLU is only for square (and moreover invertible) matrices");
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline PartialPivLU<t_row, t_col, t_type>::PartialPivLU(const t_type *element, const size_t n_byte) : m_inv(), m_pivot()
{
    static_assert(t_row == t_col, "PartialPivLU is only for square (and moreover invertible) matrices");

    if ((sizeof(t_type) * t_row * t_col) != n_byte)
        m_isOk = 0;
    else
    {
        memset(m_elem, element, n_byte);
        Compute();
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline PartialPivLU<t_row, t_col, t_type>::PartialPivLU(const Matrix<t_row, t_col, t_type> &m) : m_inv(), m_pivot()
{
    static_assert(t_row == t_col, "PartialPivLU is only for square (and moreover invertible) matrices");

    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline PartialPivLU<t_row, t_col, t_type>::PartialPivLU(const Matrix3<t_type, t_row, t_col> &m) : m_inv(), m_pivot()
{
    static_assert(t_row == t_col, "PartialPivLU is only for square (and moreover invertible) matrices");

    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t PartialPivLU<t_row, t_col, t_type>::Compute()
{
    uint16_t i, j, k;
    t_type *pMi, *pMk, *p_pivotRow = nullptr;
    t_type max, absElem;
    t_type pivotRow[t_col];

    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
    {
        /* Pivoting */
        // find the pivot row
        m_pivot[i] = i;
        max = std::abs(*(pMi + i));
        for (k = i + 1, pMk = pMi + t_col; k < t_row; k++, pMk += t_col)
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
            memcpy(pivotRow, p_pivotRow, sizeof(t_type) * t_col);
            memcpy(p_pivotRow, pMi, sizeof(t_type) * t_col);
            memcpy(pMi, pivotRow, sizeof(t_type) * t_col);
        }

        // matrix is singular, return error
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            m_isOk = 0;
            return -1;
        }

        /* LU Decompostion using Gaussian Elimination */
        for (k = i + 1, pMk = pMi + t_col; k < t_row; pMk += t_col, k++)
        {
            // find the lower triangular matrix elements for column i.
            *(pMk + i) /= *(pMi + i);

            // update the upper triangular matrix for remaining matrix
            for (j = i + 1; j < t_col; j++)
                *(pMk + j) -= *(pMk + i) * *(pMi + j);
        }
    }

    m_isOk = 1;
    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t PartialPivLU<t_row, t_col, t_type>::Compute(const t_type *element, const size_t n_byte)
{
    if ((sizeof(t_type) * t_row * t_col) != n_byte)
    {
        m_isOk = 0;
        return -1;
    }

    memcpy(m_elem, element, n_byte);
    memset(m_pivot, 0, sizeof(m_pivot));
    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t PartialPivLU<t_row, t_col, t_type>::Compute(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    memset(m_pivot, 0, sizeof(m_pivot));
    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t PartialPivLU<t_row, t_col, t_type>::Compute(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    memset(m_pivot, 0, sizeof(m_pivot));
    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type PartialPivLU<t_row, t_col, t_type>::Determinant()
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    uint16_t offset = t_row + 1;
    t_type det = 1;

    for (uint16_t i = 0; i < t_row; i++)
        det *= m_elem[i * offset];

    return det;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> PartialPivLU<t_row, t_col, t_type>::GetMatrix() const
{
    return Matrix<t_row, t_col, t_type>(m_elem);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> PartialPivLU<t_row, t_col, t_type>::GetMatrixL() const
{
    uint16_t i, j;
    t_type L[t_row * t_col]{0};

    /* Set diagonal elements as 1 */
    for (i = 0; i < t_row; i++)
    {
        L[i * (t_col + 1)] = 1;
    }

    /* Update remaining matrix from m_elem to L*/
    for (i = 1; i < t_row; i++)
        for (j = 0; j < i; j++)
            L[i * t_col + j] = m_elem[i * t_col + j];

    return Matrix<t_row, t_col, t_type>(L);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> PartialPivLU<t_row, t_col, t_type>::GetMatrixU() const
{
    uint16_t i, j;
    t_type U[t_row * t_col]{0};

    for (i = 0; i < t_row; i++)
        for (j = i; j < t_col; j++)
            U[i * t_col + j] = m_elem[i * t_col + j];

    return Matrix<t_row, t_col, t_type>(U);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> PartialPivLU<t_row, t_col, t_type>::GetMatrixP() const
{
    uint16_t i;
    t_type P[t_row * t_col]{0};
    t_type row[t_col]{0};

    for (i = 0; i < t_row; i++)
        P[i * (t_col + 1)] = 1;

    for (i = 0; i < t_row; i++)
    {
        if (m_pivot[i] != i)
        {
            memcpy(row, &P[i * t_col], sizeof(t_type) * t_col);
            memcpy(&P[i * t_col], &P[m_pivot[i] * t_col], sizeof(t_type) * t_col);
            memcpy(&P[m_pivot[i] * t_col], row, sizeof(t_type) * t_col);
        }
    }

    return Matrix<t_row, t_col, t_type>(P);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t PartialPivLU<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, Vector<t_col, t_type> &x)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    if (!m_isOk)
    {
        assert(false && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, k;
    t_type *pMi;
    t_type tmp;
    t_type vb[t_row];

    memcpy(vb, b.m_elem, sizeof(vb));

    /* Solve Ly = b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
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
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < t_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t PartialPivLU<t_row, t_col, t_type>::Solve(const Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    assert(b.m_elem != nullptr && "Memory has not been allocated");
    assert(x.m_elem != nullptr && "Memory has not been allocated");
    assert(b.m_row == t_row && "Check dimensions");
    assert(x.m_row == t_col && "Check dimensions");

    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    if (!m_isOk)
    {
        assert(false && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, k;
    t_type *pMi;
    t_type tmp;
    t_type vb[t_row];

    memcpy(vb, b.m_elem, sizeof(vb));

    /* Solve Ly = b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
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
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < t_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> PartialPivLU<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return Vector<t_col, t_type>();
    }

    int i, k;
    t_type *pMi;
    t_type x[t_col]{0};
    t_type tmp;
    t_type vb[t_row];

    memcpy(vb, b.m_elem, sizeof(vb));

    /* Solve Ly =b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
    {
        if (m_pivot[i] != i)
        {
            tmp = vb[i];
            vb[i] = vb[m_pivot[i]];
            vb[m_pivot[i]] = tmp;
        }

        x[i] = vb[i];
        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];
    }

    /* Solve Ux = y */
    // interchange the row of vector b along the original position
    // Solve the upper triangular (backward substitution)
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < t_col; k++)
            x[i] -= *(pMi + k) * x[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x[i] /= *(pMi + i);
    }

    if (isOk) *isOk = 1;
    return Vector<t_col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> PartialPivLU<t_row, t_col, t_type>::Solve(const Vector<0, t_type> &b, int8_t *isOk)
{
    assert(b.m_elem != nullptr && "Memory has not been allocated");
    assert(b.m_row == t_row && "Check dimensions");

    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return Vector<t_col, t_type>();
    }

    int i, k;
    t_type *pMi;
    t_type x[t_col]{0};
    t_type tmp;
    t_type vb[t_row];

    memcpy(vb, b.m_elem, sizeof(vb));

    /* Solve Ly =b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
    {
        if (m_pivot[i] != i)
        {
            tmp = vb[i];
            vb[i] = vb[m_pivot[i]];
            vb[m_pivot[i]] = tmp;
        }

        x[i] = vb[i];
        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];
    }

    /* Solve Ux = y */
    // interchange the row of vector b along the original position
    // Solve the upper triangular (backward substitution)
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < t_col; k++)
            x[i] -= *(pMi + k) * x[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x[i] /= *(pMi + i);
    }

    if (isOk) *isOk = 1;
    return Vector<t_col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t PartialPivLU<t_row, t_col, t_type>::Inverse(Matrix<t_row, t_col, t_type> &inv)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    uint16_t colIdx[t_col]; // column index for interchange the current col with the pivot col
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};
    int tempCol;

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < t_row; i++)
    {
        invL[i * (t_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + t_col;
    p_invMi = invL + t_col;
    for (i = 1; i < t_row; i++, p_Mi += t_col, p_invMi += t_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += t_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + t_col;
            for (k = j + 1; k < i; k++, p_invMk += t_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < t_row; k++, p_Mi += (t_col + 1), p_invMk += (t_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + t_col * (t_row - 2);
    p_invMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, p_Mi -= t_col, p_invMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + t_col;
            for (k = i + 1; k <= j; k++, p_invMk += t_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < t_row; i++)
    {
        for (j = 0; j < t_col; j++)
        {
            if (m_pivot[j] != j)
            {
                tempCol = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tempCol;
            }

            for (k = 0; k < t_col; k++)
                inv.m_elem[i * t_col + colIdx[j]] += invU[i * t_col + k] * invL[k * t_col + j];

            colIdx[j] = j;
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t PartialPivLU<t_row, t_col, t_type>::Inverse(Matrix3<t_type, t_row, t_col> &inv)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    uint16_t colIdx[t_col]; // column index for interchange the current col with the pivot col
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};
    int tempCol;

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < t_row; i++)
    {
        invL[i * (t_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + t_col;
    p_invMi = invL + t_col;
    for (i = 1; i < t_row; i++, p_Mi += t_col, p_invMi += t_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += t_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + t_col;
            for (k = j + 1; k < i; k++, p_invMk += t_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < t_row; k++, p_Mi += (t_col + 1), p_invMk += (t_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + t_col * (t_row - 2);
    p_invMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, p_Mi -= t_col, p_invMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + t_col;
            for (k = i + 1; k <= j; k++, p_invMk += t_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < t_row; i++)
    {
        for (j = 0; j < t_col; j++)
        {
            if (m_pivot[j] != j)
            {
                tempCol = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tempCol;
            }

            for (k = 0; k < t_col; k++)
                inv.m_elem[i * t_col + colIdx[j]] += invU[i * t_col + k] * invL[k * t_col + j];

            colIdx[j] = j;
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> PartialPivLU<t_row, t_col, t_type>::Inverse(int8_t *isOk)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return Matrix<t_row, t_col, t_type>();
    }

    int i, j, k;
    uint16_t colIdx[t_col]; // column index for interchange the current col with the pivot col
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};
    int tempCol;

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < t_row; i++)
    {
        invL[i * (t_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + t_col;
    p_invMi = invL + t_col;
    for (i = 1; i < t_row; i++, p_Mi += t_col, p_invMi += t_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += t_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + t_col;
            for (k = j + 1; k < i; k++, p_invMk += t_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < t_row; k++, p_Mi += (t_col + 1), p_invMk += (t_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + t_col * (t_row - 2);
    p_invMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, p_Mi -= t_col, p_invMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + t_col;
            for (k = i + 1; k <= j; k++, p_invMk += t_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    memset(m_inv, 0, sizeof(t_type) * t_row * t_col);
    for (i = 0; i < t_row; i++)
    {
        for (j = 0; j < t_col; j++)
        {
            if (m_pivot[j] != j)
            {
                tempCol = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tempCol;
            }

            for (k = 0; k < t_col; k++)
                m_inv[i * t_col + colIdx[j]] += invU[i * t_col + k] * invL[k * t_col + j];

            colIdx[j] = j;
        }
    }

    if (isOk) *isOk = 1;
    return Matrix<t_row, t_col, t_type>(m_inv);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t PartialPivLU<t_row, t_col, t_type>::InverseArray(t_type *inv)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    uint16_t colIdx[t_col]; // column index for interchange the current col with the pivot col
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};
    int tempCol;

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < t_row; i++)
    {
        invL[i * (t_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + t_col;
    p_invMi = invL + t_col;
    for (i = 1; i < t_row; i++, p_Mi += t_col, p_invMi += t_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += t_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + t_col;
            for (k = j + 1; k < i; k++, p_invMk += t_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < t_row; k++, p_Mi += (t_col + 1), p_invMk += (t_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + t_col * (t_row - 2);
    p_invMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, p_Mi -= t_col, p_invMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + t_col;
            for (k = i + 1; k <= j; k++, p_invMk += t_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < t_row; i++)
    {
        for (j = 0; j < t_col; j++)
        {
            if (m_pivot[j] != j)
            {
                tempCol = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tempCol;
            }

            for (k = 0; k < t_col; k++)
                inv[i * t_col + colIdx[j]] += invU[i * t_col + k] * invL[k * t_col + j];

            colIdx[j] = j;
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type *PartialPivLU<t_row, t_col, t_type>::InverseArray(int8_t *isOk)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return nullptr;
    }

    int i, j, k;
    uint16_t colIdx[t_col]; // column index for interchange the current col with the pivot col
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};
    int tempCol;
    memset(m_inv, 0, sizeof(t_type) * t_row * t_col);

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < t_row; i++)
    {
        invL[i * (t_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + t_col;
    p_invMi = invL + t_col;
    for (i = 1; i < t_row; i++, p_Mi += t_col, p_invMi += t_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += t_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + t_col;
            for (k = j + 1; k < i; k++, p_invMk += t_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < t_row; k++, p_Mi += (t_col + 1), p_invMk += (t_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + t_col * (t_row - 2);
    p_invMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, p_Mi -= t_col, p_invMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + t_col;
            for (k = i + 1; k <= j; k++, p_invMk += t_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < t_row; i++)
    {
        for (j = 0; j < t_col; j++)
        {
            if (m_pivot[j] != j)
            {
                tempCol = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tempCol;
            }

            for (k = 0; k < t_col; k++)
                m_inv[i * t_col + colIdx[j]] += invU[i * t_col + k] * invL[k * t_col + j];

            colIdx[j] = j;
        }
    }

    if (isOk) *isOk = 1;
    return m_inv;
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTPARTIAL_PIV_LU_TPP_