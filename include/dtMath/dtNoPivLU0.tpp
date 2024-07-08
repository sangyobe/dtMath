/*!
\file       dtNoPivLU.tpp
\brief      dtMath, Cholesky decomposition(L*L^T form) Class with Dynamic Memory Allocation
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 05. 17
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTNO_PIV_LU0_TPP_
#define DTMATH_DTNO_PIV_LU0_TPP_
#include "dtMathMem.h"
#include "dtNoPivLU.h"

#include <assert.h>

namespace dt
{
namespace Math
{

/** Initializes the NoPivLU object (Default Constructor)
 */
template <typename t_type>
inline NoPivLU<0, 0, t_type>::NoPivLU()
    : m_row(0), m_col(0), m_size(0), m_elem(nullptr), m_isOk(0)
{
    assert(m_row == m_col && "NoPivLU is only for square (and moreover invertible) matrices");
}

/** Initializes the NoPivLU object
 *
 * Function Arguments:
 *
 * row     -> no. of rows of the matrix
 * col     -> no. of cols of the matrix
 * element -> array containing elements of the matrix
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline NoPivLU<0, 0, t_type>::NoPivLU(const uint16_t row, const uint16_t col, const t_type *element)
    : m_row(row), m_col(col), m_size(row * col)
{
    assert(row == col && "NoPivLU is only for square (and moreover invertible) matrices");

    m_elem = MemAlloc<t_type>(m_size);

    memcpy(m_elem, element, sizeof(t_type) * m_size);
    Compute();
}

/** Initializes the NoPivLU object
 *
 * Function Arguments:
 *
 * m  -> matrix in the dynamic memory Matrix class form
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline NoPivLU<0, 0, t_type>::NoPivLU(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_row == m.m_col && "NoPivLU is only for square (and moreover invertible) matrices");

    m_row = m.m_row;
    m_col = m.m_col;
    m_size = m_row * m_col;

    m_elem = MemAlloc<t_type>(m_size);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    Compute();
}

/** Initializes the NoPivLU object

* Function Arguments:
*
* m  -> matrix in the Matrix class form
*
* \return    -> 0 on success, -1 on failure
*/
template <typename t_type>
template <uint16_t row, uint16_t col>
inline NoPivLU<0, 0, t_type>::NoPivLU(const Matrix<row, col, t_type> &m)
    : m_row(row), m_col(col), m_size(row * col)
{
    static_assert(row == col, "NoPivLU is only for square (and moreover invertible) matrices");

    m_elem = MemAlloc<t_type>(m_size);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    Compute();
}

/** Initializes the NoPivLU object
 *
 * Function Arguments:
 *
 * m  -> matrix in the Matrix3 class form
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline NoPivLU<0, 0, t_type>::NoPivLU(const Matrix3<t_type, 3, 3> &m)
    : m_row(3), m_col(3), m_size(9)
{
    m_elem = MemAlloc<t_type>(9);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * 9);
    Compute();
}

/** Deconstructor for the NoPivLU object
 */
template <typename t_type>
inline NoPivLU<0, 0, t_type>::~NoPivLU()
{
    if (m_elem)
    {
        MemFree(m_elem);
        m_elem = nullptr;
    }
}

/** Compute the L and U Triangular Matrix
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t NoPivLU<0, 0, t_type>::Compute()
{
    int x, i, j, k;
    t_type *pMx, *pMi, *pMk;

    for (x = 0, pMx = m_elem; x < m_row; pMx += m_col, x++)
    {
        /* find the uppper triangular matrix element for row x */
        for (j = x; j < m_col; j++)
        {
            // U(x,j) = A(x,j) - sum of (L(x,k) * U(k,j))
            for (k = 0, pMk = m_elem; k < x; pMk += m_col, k++)
                *(pMx + j) -= *(pMx + k) * *(pMk + j); // in-place
        }

        if (std::abs(*(pMx + x)) <= std::numeric_limits<t_type>::epsilon())
        {
            m_isOk = 0;
            assert(false && "The Matrix is Singular");
            return -1; // matrix is singular, if Diagonal is 0
        }

        /* find the lower triangular matrix element for col x */
        for (i = x + 1, pMi = pMx + m_row; i < m_row; pMi += m_row, i++)
        {
            // L(i,x) = [A(i,x) - sum of (L(i,k) * U(k,x))] / Uxx
            for (k = 0, pMk = m_elem; k < x; pMk += m_row, k++)
                *(pMi + x) -= *(pMi + k) * *(pMk + x);
            *(pMi + x) /= *(pMx + x);
        }
    }

    m_isOk = 1;
    return 0;
}

/** Compute the L and U Triangular Matrix
 *
 * Function Arguments:
 *
 * m  -> matrix in the Matrix class form
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t NoPivLU<0, 0, t_type>::Compute(const Matrix<row, col, t_type> &m)
{
    static_assert(row == col, "NoPivLU should be square (and moreover invertible) matrices");

    if (m_row != row || m_col != col)
    {
        if (m_elem) MemFree(m_elem);

        m_row = row;
        m_col = col;
        m_size = row * col;
        m_elem = MemAlloc<t_type>(m_size);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);

    return Compute();
}

/** Compute the L and U Triangular Matrix
 *
 * Function Arguments:
 *
 * row     -> no. of rows of the matrix
 * col     -> no. of cols of the matrix
 * element -> array containing elements of the matrix
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t NoPivLU<0, 0, t_type>::Compute(const uint16_t row, const uint16_t col, const t_type *element)
{
    assert(row == col && "NoPivLU should be square (and moreover invertible) matrices");
    assert(row != 0 && "row and col must be non zero");

    if (m_row != row || m_col != col)
    {
        if (m_elem) MemFree(m_elem);

        m_row = row;
        m_col = col;
        m_size = row * col;
        m_elem = MemAlloc<t_type>(m_size);
    }

    memcpy(m_elem, element, sizeof(t_type) * m_size);

    return Compute();
}

/** Compute the L and U Triangular Matrix
 *
 * Function Arguments:
 *
 * m  -> matrix in the dynamic memory Matrix class form
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t NoPivLU<0, 0, t_type>::Compute(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_row == m.m_col && "NoPivLU should be square (and moreover invertible) matrices");
    assert(m.m_row != 0 && "row and col must be non zero");

    if (m_row != m.m_row || m_col != m.m_col)
    {
        if (m_elem) MemFree(m_elem);

        m_row = m.m_row;
        m_col = m.m_col;
        m_size = m_row * m_col;
        m_elem = MemAlloc<t_type>(m_size);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);

    return Compute();
}

/** Compute the L and U Triangular Matrix
 *
 * Function Arguments:
 *
 * m  -> matrix in the Matrix3 class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t NoPivLU<0, 0, t_type>::Compute(const Matrix3<t_type, 3, 3> &m)
{

    if (m_row != 3 || m_col != 3)
    {
        if (m_elem) MemFree(m_elem);

        m_row = 3;
        m_col = 3;
        m_size = 9;
        m_elem = MemAlloc<t_type>(9);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * 9);
    return Compute();
}

/** Get the complete matrix after factorization with with no pivoting
 *
 * Get the complete solved matrix after the A = LU factorization
 *
 * \return    -> LU matrix elements in a full matrix in dynamic Matrix class
 */

template <typename t_type>
inline Matrix<0, 0, t_type> NoPivLU<0, 0, t_type>::GetMatrix() const
{
    return Matrix<0, 0, t_type>(3, 3, m_elem);
}

/** Get the L matrix after factorization with full pivoting
 *
 * Get the lower triangular matrix after the A = LU factorization with full pivoting
 *
 * \return    -> L matrix in the A = LU equation in dynamic Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> NoPivLU<0, 0, t_type>::GetMatrixL() const
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

/** Get the U matrix after factorization with full pivoting
 *
 * Get the upper triangular matrix after the A = LU factorization with full pivoting
 *
 * \return    -> U matrix in the A = LU equation in dynamic Matrix class
 */

template <typename t_type>
inline Matrix<0, 0, t_type> NoPivLU<0, 0, t_type>::GetMatrixU() const
{
    uint16_t i, j;
    t_type U[m_size];
    memset(U, 0, sizeof(U));

    for (i = 0; i < m_row; i++)
        for (j = i; j < m_col; j++)
            U[i * m_col + j] = m_elem[i * m_col + j];

    return Matrix<0, 0, t_type>(m_row, m_col, U);
}

/** Solve the Ax = LUx = b problem
 *
 * Function Arguments:
 *
 * b  -> the b vector in from the Vector class
 * x  -> the x vector in from the Vector class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t NoPivLU<0, 0, t_type>::Solve(const Vector<row, t_type> &b, Vector<col, t_type> &x)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix
    //       U is a upper triangular matrix
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

    /* Solve Ly = b */
    // Solve the lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < row; pMi += col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Ux = y */
    // Solve the upper triangular (backward substitution), LT = U
    for (i = row - 1, pMi = m_elem + (row - 1) * col; i >= 0; i--, pMi -= col)
    {
        for (k = i + 1; k < col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

/** Solve the Ax = LUx = b problem
 *
 * Function Arguments:
 *
 * b  -> the b vector in from the dynamic memory Vector class
 * x  -> the x vector in from the dynamic memory Vector class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t NoPivLU<0, 0, t_type>::Solve(const Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix
    //       U is a upper triangular matrix
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

    /* Solve Ly = b */
    // Solve the lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Ux = y */
    // Solve the upper triangular (backward substitution), LT = U
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

/** Solve the Ax = LUx = b problem
 *
 * Function Arguments:
 *
 * b    -> the b vector in from the Vector class
 * isOk -> pointer to indicate if method is correct (default value = nullptr)
 *
 * \return    -> the x vector in from the dynamic Vector class
 */
template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> NoPivLU<0, 0, t_type>::Solve(const Vector<row, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix
    //       U is a upper triangular matrix
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

    /* Solve Ly = b */
    // Solve the lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < row; pMi += m_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Ux = y */
    // Solve the upper triangular (backward substitution), LT = U
    for (i = row - 1, pMi = m_elem + (row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    if (isOk) *isOk = 1;
    return x;
}

/** Solve the Ax = LUx = b problem
 *
 * Function Arguments:
 *
 * b    -> the b vector in from the dynamic memory Vector class
 * isOk -> pointer to indicate if method is correct (default value = nullptr)
 *
 * \return    -> the x vector in from the dynamic Vector class
 */
template <typename t_type>
inline Vector<0, t_type> NoPivLU<0, 0, t_type>::Solve(const Vector<0, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix
    //       U is a upper triangular matrix
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

    /* Solve Ly = b */
    // Solve the lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Ux = y */
    // Solve the upper triangular (backward substitution), LT = U
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    if (isOk) *isOk = 1;
    return x;
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
inline int8_t NoPivLU<0, 0, t_type>::Inverse(Matrix<row, col, t_type> &inv)
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
    t_type *pMi;
    t_type *pInvMi, *pInvMj, *pInvMk;
    t_type sum;
    t_type invL[row * col]{0};
    t_type invU[row * col]{0};

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    for (i = 0; i < row; i++)
        invL[i * (col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + col;
    pInvMi = invL + col;
    for (i = 1; i < row; i++, pMi += col, pInvMi += col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + col;
            for (k = j + 1; k < i; k++, pInvMk += col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < row; k++, pMi += (col + 1), pInvMk += (col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + col * (row - 2);
    pInvMi = invU + col * (row - 2);
    for (i = row - 2; i >= 0; i--, pMi -= col, pInvMi -= col)
    {
        for (j = col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + col;
            for (k = i + 1; k <= j; k++, pInvMk += col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < row; i++)
        for (j = 0; j < col; j++)
            for (k = 0; k < col; k++)
                inv.m_elem[i * col + j] += invU[i * col + k] * invL[k * col + j];

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
inline int8_t NoPivLU<0, 0, t_type>::Inverse(Matrix<0, 0, t_type> &inv)
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
    t_type *pMi;
    t_type *pInvMi, *pInvMj, *pInvMk;
    t_type sum;
    t_type invL[m_size];
    t_type invU[m_size];

    memset(invL, 0, sizeof(invL));
    memset(invU, 0, sizeof(invU));

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    for (i = 0; i < m_row; i++)
        invL[i * (m_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + m_col;
    pInvMi = invL + m_col;
    for (i = 1; i < m_row; i++, pMi += m_col, pInvMi += m_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += m_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + m_col;
            for (k = j + 1; k < i; k++, pInvMk += m_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < m_row; k++, pMi += (m_col + 1), pInvMk += (m_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + m_col * (m_row - 2);
    pInvMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, pMi -= m_col, pInvMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + m_col;
            for (k = i + 1; k <= j; k++, pInvMk += m_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < m_row; i++)
        for (j = 0; j < m_col; j++)
            for (k = 0; k < m_col; k++)
                inv.m_elem[i * m_col + j] += invU[i * m_col + k] * invL[k * m_col + j];

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
inline int8_t NoPivLU<0, 0, t_type>::Inverse(Matrix3<t_type, 3, 3> &inv)
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
    t_type *pMi;
    t_type *pInvMi, *pInvMj, *pInvMk;
    t_type sum;
    t_type invL[9]{0};
    t_type invU[9]{0};

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    for (i = 0; i < 3; i++)
        invL[i * 4] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + 3;
    pInvMi = invL + 3;
    for (i = 1; i < 3; i++, pMi += 3, pInvMi += 3)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += 3)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + 3;
            for (k = j + 1; k < i; k++, pInvMk += 3)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < 3; k++, pMi += (4), pInvMk += (4))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + 3;
    pInvMi = invU + 3;
    for (i = 1; i >= 0; i--, pMi -= 3, pInvMi -= 3)
    {
        for (j = 2; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + 3;
            for (k = i + 1; k <= j; k++, pInvMk += 3)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++)
                inv.m_elem[i * 3 + j] += invU[i * 3 + k] * invL[k * 3 + j];

    return 0;
}

/** Inverse of the A matrix
 *
 * Function Arguments:
 *
 * isOk    -> the output inverse matrix in Matrix
 *
 * \return    -> the output inverse matrix in dynamic memory Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> NoPivLU<0, 0, t_type>::Inverse(int8_t *isOk)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return Matrix<0, 0, t_type>();
    }

    int i, j, k;
    t_type *pMi;
    t_type *pInvMi, *pInvMj, *pInvMk;
    t_type sum;
    t_type invL[m_size];
    t_type invU[m_size];
    t_type inv[m_size];

    memset(invL, 0, sizeof(t_type) * m_size);
    memset(invU, 0, sizeof(t_type) * m_size);
    memset(inv, 0, sizeof(t_type) * m_size);

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    for (i = 0; i < m_row; i++)
        invL[i * (m_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + m_col;
    pInvMi = invL + m_col;
    for (i = 1; i < m_row; i++, pMi += m_col, pInvMi += m_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += m_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + m_col;
            for (k = j + 1; k < i; k++, pInvMk += m_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < m_row; k++, pMi += (m_col + 1), pInvMk += (m_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + m_col * (m_row - 2);
    pInvMi = invU + m_col * (m_row - 2);
    for (i = m_row - 2; i >= 0; i--, pMi -= m_col, pInvMi -= m_col)
    {
        for (j = m_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + m_col;
            for (k = i + 1; k <= j; k++, pInvMk += m_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < m_row; i++)
        for (j = 0; j < m_col; j++)
            for (k = 0; k < m_col; k++)
                inv[i * m_col + j] += invU[i * m_col + k] * invL[k * m_col + j];

    return Matrix<0, 0, t_type>(m_row, m_col, inv);
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTNO_PIV_LU0_TPP_
