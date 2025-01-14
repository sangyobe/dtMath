/*!
\file       dtLDLT.tpp
\brief      dtMath, Cholesky decomposition (L*D*L^T form) Class with Dynamic Memory Allocation
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 05. 16
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLDLT0_TPP_
#define DTMATH_DTLDLT0_TPP_

#include "dtLDLT0.h"
#include "dtMathMem.h"

#include <assert.h>

namespace dt
{
namespace Math
{

/** Initializes the LDLT object (Default Constructor)
 */

template <typename t_type>
inline LDLT<0, 0, t_type>::LDLT()
    : m_row(0), m_col(0), m_size(0), m_elem(nullptr), m_isOk(0)
{
    assert(m_row == m_col && "LDLT is only for square (and moreover invertible) matrices");
}

/** Initializes the LDLT object
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
inline LDLT<0, 0, t_type>::LDLT(const uint16_t row, const uint16_t col, const t_type *element)
    : m_row(row), m_col(col), m_size(row * col)
{
    assert(row == col && "LDLT is only for square (and moreover invertible) matrices");

    m_elem = MemAlloc<t_type>(m_size);

    memcpy(m_elem, element, sizeof(t_type) * m_size);
    Compute();
}

/** Initializes the LDLT object
 *
 * Function Arguments:
 *
 * m  -> matrix in the dynamic memory Matrix class form
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline LDLT<0, 0, t_type>::LDLT(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_row == m.m_col && "LDLT is only for square (and moreover invertible) matrices");

    m_row = m.m_row;
    m_col = m.m_col;
    m_size = m_row * m_col;

    m_elem = MemAlloc<t_type>(m_size);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    Compute();
}

/** Initializes the LDLT object

* Function Arguments:
*
* m  -> matrix in the Matrix class form
*
* \return    -> 0 on success, -1 on failure
*/
template <typename t_type>
template <uint16_t row, uint16_t col>
inline LDLT<0, 0, t_type>::LDLT(const Matrix<row, col, t_type> &m)
    : m_row(row), m_col(col), m_size(row * col)
{
    static_assert(row == col, "LDLT is only for square (and moreover invertible) matrices");

    m_elem = MemAlloc<t_type>(m_size);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
    Compute();
}

/** Initializes the LDLT object
 *
 * Function Arguments:
 *
 * m  -> matrix in the Matrix3 class form
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline LDLT<0, 0, t_type>::LDLT(const Matrix3<t_type, 3, 3> &m)
    : m_row(3), m_col(3), m_size(9)
{
    m_elem = MemAlloc<t_type>(9);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * 9);
    Compute();
}

/** Deconstructor for the LDLT object
 */
template <typename t_type>
inline LDLT<0, 0, t_type>::~LDLT()
{
    if (m_elem)
    {
        MemFree(m_elem);
        m_elem = nullptr;
    }
}

/** Compute the L, LT Triangular Matrix, and D Diagonal Matrix
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t LDLT<0, 0, t_type>::Compute()
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type Lik;

    for (i = 1, pMi = m_elem + m_col; i < m_row; pMi += m_col, i++)
    {
        // Calculate the part of L*D Matrix from Lij*Djj
        // L*D matrix are used for Calculating Dii
        // here *(pMi + j) = Lij*Djj
        for (j = 0, pMj = m_elem; j < i; j++, pMj += m_row)
            for (k = 0; k < j; k++)
                *(pMi + j) -= *(pMi + k) * *(pMj + k);

        // Calculate the diagonal element Dii
        // Calculate the Mij from Lik and Store the Mjk
        for (k = 0, pMk = m_elem; k < i; pMk += m_col, k++)
        {
            Lik = *(pMi + k) / *(pMk + k);  // *(pMi + k) = Lik * Dkk
            *(pMi + i) -= *(pMi + k) * Lik; // Dii = Aii - sum(Lik * Dkk * Lik)
            *(pMi + k) = Lik;               // Lik
            *(pMk + i) = Lik;               // Lki
        }

        // If diagonal element is not positive, return the error,
        // the matrix is not positive definite symmetric.
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            m_isOk = 0;
            assert(false && "The Matrix is not Symmetric Positive Definite");
            return -1;
        }
    }

    m_isOk = 1;
    return 0;
}

/** Compute the L, LT Triangular Matrix, and D Diagonal Matrix
 *
 * Function Arguments:
 *
 * m  -> matrix in the Matrix class form
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t LDLT<0, 0, t_type>::Compute(const Matrix<row, col, t_type> &m)
{
    static_assert(row == col, "LDLT is only for square (and moreover invertible) matrices");

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

/** Compute the L, LT Triangular Matrix, and D Diagonal Matrix
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
inline int8_t LDLT<0, 0, t_type>::Compute(const uint16_t row, const uint16_t col, const t_type *element)
{
    assert(row == col && "LDLT is only for square (and moreover invertible) matrices");
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

/** Compute the L, LT Triangular Matrix, and D Diagonal Matrix
 *
 * Function Arguments:
 *
 * m  -> matrix in the dynamic memory Matrix class form
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t LDLT<0, 0, t_type>::Compute(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_row == m.m_col && "LDLT is only for square (and moreover invertible) matrices");
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

/** Compute the L, LT Triangular Matrix, and D Diagonal Matrix
 *
 * Function Arguments:
 *
 * m  -> matrix in the Matrix3 class form
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t LDLT<0, 0, t_type>::Compute(const Matrix3<t_type, 3, 3> &m)
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

/** Return the full matrix with L, LT Triangular Matrix, and D Diagonal Matrix
 *
 * \return    -> the full matrix with L, D, and LT
 */
template <typename t_type>
inline Matrix<0, 0, t_type> LDLT<0, 0, t_type>::GetMatrix() const
{
    return Matrix<0, 0, t_type>(m_row, m_col, m_elem);
}

/** Return the Lower Triangular 'L' Matrix
 *
 * \return    -> lower triangular 'L' matrix in dynamic memory Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> LDLT<0, 0, t_type>::GetMatrixL() const
{
    uint16_t i, j;
    t_type L[m_size];
    memset(L, 0, sizeof(L));

    /* Set diagonal elements as 1 */
    for (i = 0; i < m_row; i++)
        L[i * m_col + i] = 1;

    for (i = 1; i < m_row; i++)
        for (j = 0; j < i; j++)
            L[i * m_col + j] = m_elem[i * m_col + j];

    return Matrix<0, 0, t_type>(m_row, m_col, L);
}

/** Return the Diagonal Triangular 'D' Matrix
 *
 * \return    -> diagonal 'D' matrix in dynamic memory Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> LDLT<0, 0, t_type>::GetMatrixD() const
{
    int i;
    t_type D[m_size];
    memset(D, 0, sizeof(D));

    for (i = 0; i < m_row; i++)
        D[i * m_col + i] = m_elem[i * m_col + i];

    return Matrix<0, 0, t_type>(m_row, m_col, D);
}

/** Return the Upper Triangular 'U' Matrix
 *
 * \return    -> upper triangular 'U' matrix in dynamic memory Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> LDLT<0, 0, t_type>::GetMatrixU() const
{
    uint16_t i, j;
    t_type U[m_size];
    memset(U, 0, sizeof(U));

    for (i = 0; i < m_row; i++)
        U[i * m_col + i] = 1;

    for (i = 1; i < m_row; i++)
        for (j = 0; j < i; j++)
            U[i + m_col * j] = m_elem[i * m_col + j];

    return Matrix<0, 0, t_type>(m_row, m_col, U);
}

/** Solve the Ax = LDUx = b problem
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
inline int8_t LDLT<0, 0, t_type>::Solve(const Vector<row, t_type> &b, Vector<col, t_type> &x)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
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

    /* Solve Lz = b */
    // Solve the unit lower triangular (forward substitution), here x is z
    for (i = 0, pMi = m_elem; i < row; pMi += col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Dy = z */
    // Solve the diagonal matrix, here x is y
    for (i = 0, pMi = m_elem; i < row; i++, pMi += col)
    {
        x.m_elem[i] /= *(pMi + i);
    }

    /* Solve Ux = y */
    // Solve the unit upper triangular (backward substitution)
    for (i = row - 1, pMi = m_elem + (row - 1) * col; i >= 0; i--, pMi -= col)
    {
        for (k = i + 1; k < col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    return 0;
}

/** Solve the Ax = LDUx = b problem
 *
 * Function Arguments:
 *
 * b  -> the b vector in from the dynamic memory Vector class
 * x  -> the x vector in from the dynamic memory Vector class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t LDLT<0, 0, t_type>::Solve(const Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
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

    /* Solve Lz = b */
    // Solve the unit lower triangular (forward substitution), here x is z
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Dy = z */
    // Solve the diagonal matrix, here x is y
    for (i = 0, pMi = m_elem; i < m_row; i++, pMi += m_col)
    {
        x.m_elem[i] /= *(pMi + i);
    }

    /* Solve Ux = y */
    // Solve the unit upper triangular (backward substitution)
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    return 0;
}

/** Solve the Ax = LDUx = b problem
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
inline Vector<0, t_type> LDLT<0, 0, t_type>::Solve(const Vector<row, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
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

    /* Solve Lz = b */
    // Solve the unit lower triangular (forward substitution), here x is z
    for (i = 0, pMi = m_elem; i < row; pMi += m_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Dy = z */
    // Solve the diagonal matrix, here x is y
    for (i = 0, pMi = m_elem; i < row; i++, pMi += m_col)
    {
        x.m_elem[i] /= *(pMi + i);
    }

    /* Solve Ux = y */
    // Solve the unit upper triangular (backward substitution)
    for (i = row - 1, pMi = m_elem + (row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    if (isOk) *isOk = 1;
    return x;
}

/** Solve the Ax = LDUx = b problem
 *
 * Function Arguments:
 *
 * b    -> the b vector in from the dynamic memory Vector class
 * isOk -> pointer to indicate if method is correct (default value = nullptr)
 *
 * \return    -> the x vector in from the dynamic Vector class
 */
template <typename t_type>
inline Vector<0, t_type> LDLT<0, 0, t_type>::Solve(const Vector<0, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
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

    /* Solve Lz = b */
    // Solve the unit lower triangular (forward substitution), here x is z
    for (i = 0, pMi = m_elem; i < m_row; pMi += m_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Dy = z */
    // Solve the diagonal matrix, here x is y
    for (i = 0, pMi = m_elem; i < m_row; i++, pMi += m_col)
    {
        x.m_elem[i] /= *(pMi + i);
    }

    /* Solve Ux = y */
    // Solve the unit upper triangular (backward substitution)
    for (i = m_row - 1, pMi = m_elem + (m_row - 1) * m_col; i >= 0; i--, pMi -= m_col)
    {
        for (k = i + 1; k < m_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
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
inline int8_t LDLT<0, 0, t_type>::Inverse(Matrix<row, col, t_type> &inv)
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
    t_type *pMi, *pMj, *pMk;

    memcpy(inv.m_elem, m_elem, sizeof(t_type) * row * col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv.m_elem + col; i < row; i++, pMi += col)
    {
        for (j = 0, pMj = inv.m_elem; j < i; pMj += col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + col; k < i; k++, pMk += col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv.m_elem; j < col; j++, pMj += row)
    {
        for (i = j, pMi = pMj; i < row; pMi += col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + col; k < row; k++, pMk += col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

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
inline int8_t LDLT<0, 0, t_type>::Inverse(Matrix<0, 0, t_type> &inv)
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
    t_type *pMi, *pMj, *pMk;

    memcpy(inv.m_elem, m_elem, sizeof(t_type) * m_row * m_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv.m_elem + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv.m_elem; j < i; pMj += m_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + m_col; k < i; k++, pMk += m_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv.m_elem; j < m_col; j++, pMj += m_row)
    {
        for (i = j, pMi = pMj; i < m_row; pMi += m_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + m_col; k < m_row; k++, pMk += m_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

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
inline int8_t LDLT<0, 0, t_type>::Inverse(Matrix3<t_type, 3, 3> &inv)
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
    t_type *pMi, *pMj, *pMk;
    t_type sum;

    memcpy(inv.m_elem, m_elem, sizeof(t_type) * 9);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv.m_elem + 3; i < 3; i++, pMi += 3)
    {
        for (j = 0, pMj = inv.m_elem; j < i; pMj += 3, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + 3; k < i; k++, pMk += 3)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv.m_elem; j < 3; j++, pMj += m_row)
    {
        for (i = j, pMi = pMj; i < 3; pMi += 3, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + 3; k < 3; k++, pMk += 3)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

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
inline Matrix<0, 0, t_type> LDLT<0, 0, t_type>::Inverse(int8_t *isOk)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return Matrix<0, 0, t_type>();
    }

    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type sum;
    t_type inv[m_size];

    memcpy(inv, m_elem, sizeof(t_type) * m_row * m_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv + m_col; i < m_row; i++, pMi += m_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += m_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + m_col; k < i; k++, pMk += m_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv; j < m_col; j++, pMj += m_row)
    {
        for (i = j, pMi = pMj; i < m_row; pMi += m_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + m_col; k < m_row; k++, pMk += m_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

    return Matrix<0, 0, t_type>(m_row, m_col, inv);
}

} // namespace Math

} // namespace dt

#endif // DTMATH_DTLDLT0_TPP_