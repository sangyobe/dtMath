/*!
\file       dtFullPivLU.tpp
\brief      dtMath, LU Decomposition with partial pivoting(Doolittle form) class for Dynamic Allocation
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 05. 14
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTFULL_PIV_LU0_TPP_
#define DTMATH_DTFULL_PIV_LU0_TPP_

#include "dtFullPivLU0.h"
#include "dtMathMem.h"

#include <assert.h>

namespace dt
{
namespace Math
{

/** Initializes the FullPivLU object (Default Constructor)
 */
template <typename t_type>
inline FullPivLU<0, 0, t_type>::FullPivLU()
    : m_row(0), m_col(0), m_size(0), m_elem(nullptr), m_inv(nullptr), m_pivot(nullptr), m_pivotCol(nullptr), m_permCol(nullptr), m_isOk(0)
{
    assert(m_row == m_col && "FullPivLU is only for square (and moreover invertible) matrices");
}

/** Initializes the FullPivLU object
 *
 * Function Arguments:
 * row       -> no. of rows of the matrix
 * col       -> no. of columns of the matrix
 * element   -> vector (array) containing elements of the matrix
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline FullPivLU<0, 0, t_type>::FullPivLU(const uint16_t row, const uint16_t col, const t_type *element)
    : m_row(row), m_col(col), m_size(row * col)
{
    assert(row == col && "FullPivLU is only for square (and moreover invertible) matrices");

    m_elem = MemAlloc<t_type>(m_size);
    m_inv = MemAlloc<t_type>(m_size);
    m_pivot = MemAllocZeroInit<int>(m_row);
    m_pivotCol = MemAllocZeroInit<int>(m_row);
    m_permCol = MemAllocZeroInit<int>(m_row);

    memcpy(m_elem, element, sizeof(t_type) * m_size);
    Compute();
}

/** Initializes the FullPivLU object
 *
 * Function Arguments:
 *
 * m   ->  dynamic memory Matrix object
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline FullPivLU<0, 0, t_type>::FullPivLU(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_row == m.m_col && "FullPivLU is only for square (and moreover invertible) matrices");

    m_row = m.m_row;
    m_col = m.m_col;
    m_size = m_row * m_col;

    m_elem = MemAlloc<t_type>(m_size);
    m_inv = MemAlloc<t_type>(m_size);
    m_pivot = MemAllocZeroInit<int>(m_row);
    m_pivotCol = MemAllocZeroInit<int>(m_row);
    m_permCol = MemAllocZeroInit<int>(m_row);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    Compute();
}

/** Initializes the FullPivLU object
 *
 * Function Arguments:
 *
 * m   ->  Matrix object
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
template <uint16_t row, uint16_t col>
inline FullPivLU<0, 0, t_type>::FullPivLU(const Matrix<row, col, t_type> &m)
    : m_row(row), m_col(col), m_size(row * col)
{
    static_assert(row == col, "FullPivLU is only for square (and moreover invertible) matrices");

    m_elem = MemAlloc<t_type>(m_size);
    m_inv = MemAlloc<t_type>(m_size);
    m_pivot = MemAllocZeroInit<int>(m_row);
    m_pivotCol = MemAllocZeroInit<int>(m_row);
    m_permCol = MemAllocZeroInit<int>(m_row);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    Compute();
}

/** Initializes the FullPivLU object
 *
 * Function Arguments:
 *
 * m   ->  Matrix3 object
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline FullPivLU<0, 0, t_type>::FullPivLU(const Matrix3<t_type, 3, 3> &m)
    : m_row(3), m_col(3), m_size(9)
{
    m_elem = MemAlloc<t_type>(9);
    m_inv = MemAlloc<t_type>(9);
    m_pivot = MemAllocZeroInit<int>(3);
    m_pivotCol = MemAllocZeroInit<int>(3);
    m_permCol = MemAllocZeroInit<int>(3);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * 9);
    Compute();
}

/** Deconstructor for the FullPivLU object
 */
template <typename t_type>
inline FullPivLU<0, 0, t_type>::~FullPivLU()
{
    if (m_elem)
    {
        MemFree(m_elem);
        MemFree(m_inv);
        MemFree(m_pivot);
        MemFree(m_pivotCol);
        MemFree(m_permCol);
        m_elem = nullptr;
        m_inv = nullptr;
        m_pivot = nullptr;
        m_pivotCol = nullptr;
        m_permCol = nullptr;
    }
}

/** Compute the PAQ = LU factorization with full pivoting
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t FullPivLU<0, 0, t_type>::Compute()
{
    int i, j, k;
    t_type *pMi, *pMk, *p_pivotRow = nullptr;
    t_type max_row, max_col, absElem;
    t_type pivotRow[m_row];
    t_type temp_clm;
    int tempPerm;

    for (i = 0; i < m_row; i++) // filling the Q permutation vector as identity
        m_permCol[i] = i;

    for (i = 0, pMi = m_elem; i < m_row; pMi += m_row, i++)
    {
        // Pivoting /
        // find the pivot row
        m_pivot[i] = i;
        max_row = std::abs(*(pMi + i));

        for (j = i, pMk = pMi; j < m_row; j++, pMk += m_row)
        {
            for (k = i; k < m_row; k++) // Loop for checking maximum value in each element of row j
            {
                if (max_row < (absElem = std::abs(*(pMk + k))))
                {
                    max_row = absElem;
                    m_pivot[i] = j;
                    p_pivotRow = pMk; // pMk is pivot row
                }
            }
        }

        if (m_pivot[i] != i)
        {
            memcpy(pivotRow, p_pivotRow, sizeof(t_type) * m_col);
            memcpy(p_pivotRow, pMi, sizeof(t_type) * m_col);
            memcpy(pMi, pivotRow, sizeof(t_type) * m_col);
        }

        m_pivotCol[i] = i;                                    // Stores the column to be swapped temporarily as the ith column
        max_col = std::abs(*(pMi + i));                       // Temporarily stores the value of diagonal "i x i" as maximum value
        for (k = i, pMk = pMi + i; k < m_col - 1; k++, pMk++) // Loop for finding for highest value in the "i"th row starting from i x i value
        {
            if (max_col < (absElem = std::abs(*(pMk + 1)))) // Check for maximum value for each element in row "i" (Rows already arranged in previous step)
            {
                max_col = absElem; // Replace maximum value in "max_col"
                m_pivotCol[i] = k + 1;
            }
        }

        if (m_pivotCol[i] != i)
        {

            for (k = 0, pMk = m_elem; k < m_col; k++, pMk += m_col)
            {

                temp_clm = *(pMk + i);
                *(pMk + i) = *(pMk + m_pivotCol[i]); // These three lines swap the "k"th column with the column with the next maximum value
                *(pMk + m_pivotCol[i]) = temp_clm;
            }
            tempPerm = m_permCol[i];
            m_permCol[i] = m_permCol[m_pivotCol[i]]; // These three lines swap the "i"th row of the permutation vector
            m_permCol[m_pivotCol[i]] = tempPerm;
        }

        // matrix is singular, return error
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            m_isOk = 0;
            assert(false && "The matrix is singular");
            return -1;
        }

        // LU Decompostion using Gaussian Elimination //
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

/** Compute the LU factorization with full pivoting
 *
 * Function Arguments:
 *
 * row       -> rows of the matrix
 * col       -> columns of the matrix
 * element   -> vector (array) containing elements of the matrix
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t FullPivLU<0, 0, t_type>::Compute(const uint16_t row, const uint16_t col, const t_type *element)
{
    assert(row == col && "FullPivLU is only for square (and moreover invertible) matrices");
    assert(row != 0 && "row and col must be non zero");

    if (m_row != row || m_col != col)
    {
        if (m_elem) MemFree(m_elem);
        if (m_inv) MemFree(m_inv);
        if (m_pivot) MemFree(m_pivot);
        if (m_pivotCol) MemFree(m_pivotCol);
        if (m_permCol) MemFree(m_permCol);

        m_row = row;
        m_col = col;
        m_size = row * col;
        m_elem = MemAlloc<t_type>(m_size);
        m_inv = MemAlloc<t_type>(m_size);
        m_pivot = MemAlloc<int>(m_row);
        m_pivotCol = MemAlloc<int>(m_row);
        m_permCol = MemAlloc<int>(m_row);
    }

    memcpy(m_elem, element, sizeof(t_type) * m_size);
    memset(m_pivot, 0, sizeof(int) * m_row);
    memset(m_pivotCol, 0, sizeof(int) * m_row);
    memset(m_permCol, 0, sizeof(int) * m_row);

    return Compute();
}

/** Compute the LU factorization with full pivoting
 *
 * Function Arguments:
 *
 * m   ->  dynamic memory Matrix object
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t FullPivLU<0, 0, t_type>::Compute(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_row == m.m_col && "FullPivLU is only for square (and moreover invertible) matrices");
    assert(m.m_row != 0 && "row and col must be non zero");

    if (m_row != m.m_row || m_col != m.m_col)
    {
        if (m_elem) MemFree(m_elem);
        if (m_inv) MemFree(m_inv);
        if (m_pivot) MemFree(m_pivot);
        if (m_pivotCol) MemFree(m_pivotCol);
        if (m_permCol) MemFree(m_permCol);

        m_row = m.m_row;
        m_col = m.m_col;
        m_size = m_row * m_col;
        m_elem = MemAlloc<t_type>(m_size);
        m_inv = MemAlloc<t_type>(m_size);
        m_pivot = MemAlloc<int>(m_row);
        m_pivotCol = MemAlloc<int>(m_row);
        m_permCol = MemAlloc<int>(m_row);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    memset(m_pivot, 0, sizeof(int) * m_row);
    memset(m_pivotCol, 0, sizeof(int) * m_row);
    memset(m_permCol, 0, sizeof(int) * m_row);

    return Compute();
}

/** Compute the LU factorization with full pivoting
 *
 * Function Arguments:
 *
 * m   ->  Matrix object
 *
 * \return    -> 0 on success, -1 on failure
 */

template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t FullPivLU<0, 0, t_type>::Compute(const Matrix<row, col, t_type> &m)
{
    static_assert(row == col, "FullPivLU is only for square (and moreover invertible) matrices");

    if (m_row != row || m_col != col)
    {
        if (m_elem) MemFree(m_elem);
        if (m_inv) MemFree(m_inv);
        if (m_pivot) MemFree(m_pivot);
        if (m_pivotCol) MemFree(m_pivotCol);
        if (m_permCol) MemFree(m_permCol);

        m_row = row;
        m_col = col;
        m_size = row * col;
        m_elem = MemAlloc<t_type>(m_size);
        m_inv = MemAlloc<t_type>(m_size);
        m_pivot = MemAlloc<int>(m_row);
        m_pivotCol = MemAlloc<int>(m_row);
        m_permCol = MemAlloc<int>(m_row);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    memset(m_pivot, 0, sizeof(int) * m_row);
    memset(m_pivotCol, 0, sizeof(int) * m_row);
    memset(m_permCol, 0, sizeof(int) * m_row);

    return Compute();
}

/** Compute the LU factorization with full pivoting
 *
 * Function Arguments:
 *
 * m   ->  Matrix3 object
 *
 * \return    -> 0 on success, -1 on failure
 */

template <typename t_type>
inline int8_t FullPivLU<0, 0, t_type>::Compute(const Matrix3<t_type, 3, 3> &m)
{
    if (m_row != 3 || m_col != 3)
    {

        if (m_elem) MemFree(m_elem);
        if (m_inv) MemFree(m_inv);
        if (m_pivot) MemFree(m_pivot);
        if (m_pivotCol) MemFree(m_pivotCol);
        if (m_permCol) MemFree(m_permCol);

        m_row = 3;
        m_col = 3;
        m_size = 9;
        m_elem = MemAlloc<t_type>(9);
        m_inv = MemAlloc<t_type>(9);
        m_pivot = MemAlloc<int>(3);
        m_pivotCol = MemAlloc<int>(3);
        m_permCol = MemAlloc<int>(3);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * 9);
    memset(m_pivot, 0, sizeof(int) * 3);
    memset(m_pivotCol, 0, sizeof(int) * 3);
    memset(m_permCol, 0, sizeof(int) * 3);

    return Compute();
}

/** Get the complete matrix after factorization with with full pivoting
 *
 * Get the complete solved matrix after the PAQ = LU factorization
 *
 * \return    -> LU matrix elements in a full matrix in dynamic Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> FullPivLU<0, 0, t_type>::GetMatrix() const
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return Matrix<0, 0, t_type>();
    }
    return Matrix<0, 0, t_type>(m_row, m_col, m_elem);
}

/** Get the L matrix after factorization with full pivoting
 *
 * Get the lower triangular matrix after the PAQ = LU factorization with full pivoting
 *
 * \return    -> L matrix in the PAQ = LU equation in dynamic Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> FullPivLU<0, 0, t_type>::GetMatrixL() const
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return Matrix<0, 0, t_type>();
    }
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

/** Get the U matrix after factorization with full pivoting
 *
 * Get the upper triangular matrix after the PAQ = LU factorization with full pivoting
 *
 * \return    -> U matrix in the PAQ = LU equation in dynamic Matrix class
 */

template <typename t_type>
inline Matrix<0, 0, t_type> FullPivLU<0, 0, t_type>::GetMatrixU() const
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return Matrix<0, 0, t_type>();
    }
    uint16_t i, j;
    t_type U[m_size];
    memset(U, 0, sizeof(U));

    for (i = 0; i < m_row; i++)
        for (j = i; j < m_col; j++)
            U[i * m_col + j] = m_elem[i * m_col + j];

    return Matrix<0, 0, t_type>(m_row, m_col, U);
}

/** Get the P matrix
 *
 * Get the row permutation matrix 'P' in the PAQ = LU equation
 *
 * \return    -> P matrix in the PAQ = LU equation in dynamic Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> FullPivLU<0, 0, t_type>::GetMatrixP() const
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return Matrix<0, 0, t_type>();
    }

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

/** Get the Q matrix
 *
 * Get the column permutation matrix 'Q' in the PAQ = LU equation
 *
 * \return    -> Q matrix in the PAQ = LU equation in dynamic Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> FullPivLU<0, 0, t_type>::GetMatrixQ() const
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return Matrix<0, 0, t_type>();
    }

    uint16_t i, j;
    t_type Q[m_size];
    t_type temp_clm;
    t_type *pMk;

    memset(Q, 0, sizeof(Q));

    for (i = 0, pMk = Q; i < m_row; pMk += m_row, i++)
        *(pMk + i) = 1;

    for (i = 0; i < m_row; i++)
    {
        if (m_pivotCol[i] != i)
        {
            for (j = 0, pMk = Q; j < m_col; j++, pMk += m_col)
            {

                temp_clm = *(pMk + i);
                *(pMk + i) = *(pMk + m_pivotCol[i]); // These three lines swap the "j"th column with the m_pivotCol[i]th column
                *(pMk + m_pivotCol[i]) = temp_clm;
            }
        }
    }

    return Matrix<0, 0, t_type>(m_row, m_col, Q);
}

/** Perform the P*A*Q*x = L*U*x = b operation
 *
 * Function Arguments:
 *
 * b ->  The 'b' vector in the L*U*x = b equation in Vector class
 * x ->  The 'x' vector in the L*U*x = b equation in Vector class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t FullPivLU<0, 0, t_type>::Solve(const Vector<row, t_type> &b, Vector<col, t_type> &x)
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
    t_type x_out[col]{0};

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

        x_out[i] = vb[i];
        for (k = 0; k < i; k++)
            x_out[i] -= *(pMi + k) * x_out[k];
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
            x_out[i] -= *(pMi + k) * x_out[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x_out[i] /= *(pMi + i);
        x.m_elem[m_permCol[i]] = x_out[i];
    }

    return 0;
}

/** Perform the P*A*Q*x = L*U*x = b operation

* Function Arguments:
*
* b ->  The 'b' vector in the L*U*x = b equation in dynamic Vector class
* x ->  The 'x' vector in the L*U*x = b equation in dynamic Vector class
*
* \return    -> 0 on success, -1 on failure
*/
template <typename t_type>
inline int8_t FullPivLU<0, 0, t_type>::Solve(const Vector<0, t_type> &b, Vector<0, t_type> &x)
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
    t_type x_out[m_col];

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

        x_out[i] = vb[i];
        for (k = 0; k < i; k++)
            x_out[i] -= *(pMi + k) * x_out[k];
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
            x_out[i] -= *(pMi + k) * x_out[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x_out[i] /= *(pMi + i);
        x.m_elem[m_permCol[i]] = x_out[i];
    }

    return 0;
}

/** Perform the P*A*Q*x = L*U*x = b operation
 *
 * Function Arguments:
 *
 * b     ->  The 'b' vector in the L*U*x = b equation in Vector class
 * isOk  ->  the Ok flag to indicate that the process completed properly
 *
 * \return    -> The 'x' vector in the L*U*x = b in dynamic Vector class
 */
template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> FullPivLU<0, 0, t_type>::Solve(const Vector<row, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    Vector<0, t_type> x_out(m_col);

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return x_out;
    }

    int i, k;
    t_type *pMi;
    t_type x[m_col];
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

        x[i] = vb[i];
        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];
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
            x[i] -= *(pMi + k) * x[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x[i] /= *(pMi + i);
        x_out.m_elem[m_permCol[i]] = x[i]; //  x_out = Q * x
    }

    if (isOk) *isOk = 1;
    return x_out;
}

/** Perform the P*A*Q*x = L*U*x = b operation
 *
 * Function Arguments:
 *
 * b     ->  The 'b' vector in the L*U*x = b equation in dynamic Vector class
 * isOk  ->  the Ok flag to indicate that the process completed properly
 *
 * \return    -> The 'x' vector in the L*U*x = b in dynamic Vector class
 */
template <typename t_type>
inline Vector<0, t_type> FullPivLU<0, 0, t_type>::Solve(const Vector<0, t_type> &b, int8_t *isOk)
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

    Vector<0, t_type> x_out(m_col);

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return x_out;
    }

    int i, k;
    t_type *pMi;
    t_type tmp;
    t_type x[m_col];
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

        x[i] = vb[i];
        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];
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
            x[i] -= *(pMi + k) * x[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x[i] /= *(pMi + i);
        x_out.m_elem[m_permCol[i]] = x[i]; //  x_out = Q * x
    }

    if (isOk) *isOk = 1;
    return x_out;
}

/** Inverse of the A matrix
 *
 * Function Arguments:
 *
 * inv    -> the output inverse matrix in Matrix class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t FullPivLU<0, 0, t_type>::Inverse(Matrix<row, col, t_type> &inv)
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
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[row * col]{0};
    t_type invU[row * col]{0};
    t_type QinvPAQ[row * col]{0};
    t_type *p_pivotRow;
    int tempCol;

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
                tempCol = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tempCol;
            }

            for (k = 0; k < col; k++)
                inv.m_elem[i * col + colIdx[j]] += invU[i * col + k] * invL[k * col + j];

            colIdx[j] = j;
        }
    }
    for (i = 0, p_Mi = inv.m_elem; i < row; i++, p_Mi += col)
    {
        p_pivotRow = QinvPAQ + (m_permCol[i] * col);
        memcpy(p_pivotRow, p_Mi, sizeof(t_type) * col);

        // printf("m_pivot[%d] = %d, colIdx[%d] = %d \n", i, m_pivot[i], i, colIdx[i]);
    }
    memcpy(inv.m_elem, QinvPAQ, sizeof(t_type) * col * row);

    return 0;
}

/** Inverse of the A matrix
 *
 * Function Arguments:
 *
 * inv    -> the output inverse matrix in dynamic memory Matrix class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t FullPivLU<0, 0, t_type>::Inverse(Matrix<0, 0, t_type> &inv)
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
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[m_size];
    t_type invU[m_size];
    t_type QinvPAQ[m_size];
    t_type *p_pivotRow;
    int tempCol;

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
                tempCol = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tempCol;
            }

            for (k = 0; k < m_col; k++)
                inv.m_elem[i * m_col + colIdx[j]] += invU[i * m_col + k] * invL[k * m_col + j];

            colIdx[j] = j;
        }
    }
    for (i = 0, p_Mi = inv.m_elem; i < m_row; i++, p_Mi += m_col)
    {
        p_pivotRow = QinvPAQ + (m_permCol[i] * m_col);
        memcpy(p_pivotRow, p_Mi, sizeof(t_type) * m_col);

        // printf("m_pivot[%d] = %d, colIdx[%d] = %d \n", i, m_pivot[i], i, colIdx[i]);
    }
    memcpy(inv.m_elem, QinvPAQ, sizeof(t_type) * m_col * m_row);

    return 0;
}

/** Inverse of the A matrix
 *
 * Function Arguments:
 *
 * inv    -> the output inverse matrix in Matrix3 class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t FullPivLU<0, 0, t_type>::Inverse(Matrix3<t_type, 3, 3> &inv)
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
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[9]{0};
    t_type invU[9]{0};
    t_type QinvPAQ[9]{0};
    t_type *p_pivotRow;
    int tempCol;

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
                tempCol = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tempCol;
            }

            for (k = 0; k < 3; k++)
                inv.m_elem[i * 3 + colIdx[j]] += invU[i * 3 + k] * invL[k * 3 + j];

            colIdx[j] = j;
        }
    }

    for (i = 0, p_Mi = inv.m_elem; i < 3; i++, p_Mi += 3)
    {
        p_pivotRow = QinvPAQ + (m_permCol[i] * 3);
        memcpy(p_pivotRow, p_Mi, sizeof(t_type) * 3);

        // printf("m_pivot[%d] = %d, colIdx[%d] = %d \n", i, m_pivot[i], i, colIdx[i]);
    }

    memcpy(inv.m_elem, QinvPAQ, sizeof(t_type) * 9);

    return 0;
}

/** Inverse of the A matrix
 *
 * Function Arguments:
 *
 * isOk  ->  the Ok flag to indicate that the process completed properly
 *
 * \return    -> the inverse matrix in dynamic memory Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> FullPivLU<0, 0, t_type>::Inverse(int8_t *isOk)
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
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[m_size];
    t_type invU[m_size];
    t_type QinvPAQ[m_size];
    t_type m_inv[m_size];
    t_type *p_pivotRow;
    int tempCol;

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
                tempCol = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tempCol;
            }

            for (k = 0; k < m_col; k++)
                m_inv[i * m_col + colIdx[j]] += invU[i * m_col + k] * invL[k * m_col + j];

            colIdx[j] = j;
        }
    }

    for (i = 0, p_Mi = m_inv; i < m_row; i++, p_Mi += m_col)
    {
        p_pivotRow = QinvPAQ + (m_permCol[i] * m_col);
        memcpy(p_pivotRow, p_Mi, sizeof(t_type) * m_col);

        // printf("m_pivot[%d] = %d, colIdx[%d] = %d \n", i, m_pivot[i], i, colIdx[i]);
    }

    memcpy(m_inv, QinvPAQ, sizeof(t_type) * m_col * m_row);

    if (isOk) *isOk = 1;
    return Matrix<0, 0, t_type>(m_row, m_col, m_inv);
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTFULL_PIV_LU0_TPP_
