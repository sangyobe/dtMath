/*!
\file       dtSVD0.tpp
\brief      dtMath, Singular Value Decomposition solver class for dynamic memory allocation
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 05. 22
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTSVD0_TPP_
#define DTMATH_DTSVD0_TPP_

#include "dtSVD.h"

namespace dt
{
namespace Math
{

/** Initializes the SVD object (Default Constructor)
 */
template <typename t_type>
inline SVD<0, 0, t_type>::SVD()
    : m_row(0), m_col(0), m_size(0), m_A(nullptr), m_U(nullptr), m_S(nullptr), m_V(nullptr),
      m_superDiagonal(nullptr), m_inv(nullptr)
{
    m_isOk = 0;
}

/** Initializes the SVD object
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
inline SVD<0, 0, t_type>::SVD(const uint16_t row, const uint16_t col, const t_type *element)
    : m_row(row), m_col(col), m_size(row * col)
{
    m_A = MemAlloc<t_type>(m_size);
    m_U = MemAlloc<t_type>((m_row >= m_col) ? m_size : m_col * m_col);
    m_S = MemAlloc<t_type>((m_row >= m_col) ? m_col : m_row);
    m_V = MemAlloc<t_type>((m_row >= m_col) ? m_col * m_col : m_size);
    m_superDiagonal = MemAlloc<t_type>((m_row >= m_col) ? m_row : m_col);
    m_inv = MemAlloc<t_type>(m_size);

    if (m_row >= m_col)
    {
        memcpy(m_A, element, sizeof(t_type) * m_size);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn0();
        if (GivensReductionToDiagonalForm_Mxn0())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_Mxn0();
        m_isOk = 1;
    }
    else
    { // Transposed A
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < m_row; ++irow)
        {
            for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + m_row * icol] = element[irow * m_col + icol];
                m_A[irow + m_row * (icol + 1)] = element[irow * m_col + icol + 1];
                m_A[irow + m_row * (icol + 2)] = element[irow * m_col + icol + 2];
                m_A[irow + m_row * (icol + 3)] = element[irow * m_col + icol + 3];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + m_row * icol] = element[irow * m_col + icol];
            }
        }
        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN0();
        if (GivensReductionToDiagonalForm_mxN0())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_mxN0();
        m_isOk = 1;
    }
}

/** Initializes the SVD object

* Function Arguments:
*
* m  -> matrix in the Matrix class form
*
* \return    -> 0 on success, -1 on failure
*/
template <typename t_type>
template <uint16_t row, uint16_t col>
inline SVD<0, 0, t_type>::SVD(const Matrix<row, col, t_type> &m)
    : m_row(row), m_col(col), m_size(row * col)
{
    m_A = MemAlloc<t_type>(m_size);
    m_U = MemAlloc<t_type>((m_row >= m_col) ? m_size : m_col * m_col);
    m_S = MemAlloc<t_type>((m_row >= m_col) ? m_col : m_row);
    m_V = MemAlloc<t_type>((m_row >= m_col) ? m_col * m_col : m_size);
    m_superDiagonal = MemAlloc<t_type>((m_row >= m_col) ? m_row : m_col);
    m_inv = MemAlloc<t_type>(m_size);

    if (m_row >= m_col)
    {
        memcpy(m_A, m.m_elem, sizeof(t_type) * m_size);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn0();
        if (GivensReductionToDiagonalForm_Mxn0())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_Mxn0();
        m_isOk = 1;
    }
    else
    { // Transposed A
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < m_row; ++irow)
        {
            for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
                m_A[irow + m_row * (icol + 1)] = m.m_elem[irow * m_col + icol + 1];
                m_A[irow + m_row * (icol + 2)] = m.m_elem[irow * m_col + icol + 2];
                m_A[irow + m_row * (icol + 3)] = m.m_elem[irow * m_col + icol + 3];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN0();
        if (GivensReductionToDiagonalForm_mxN0())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_mxN0();
        m_isOk = 1;
    }
}

/** Initializes the SVD object
 *
 * Function Arguments:
 *
 * m  -> matrix in the Matrix3 class form
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline SVD<0, 0, t_type>::SVD(const Matrix3<t_type, 3, 3> &m)
    : m_row(3), m_col(3), m_size(9)
{
    m_A = MemAlloc<t_type>(9);
    m_U = MemAlloc<t_type>(9);
    m_S = MemAlloc<t_type>(3);
    m_V = MemAlloc<t_type>(9);
    m_superDiagonal = MemAlloc<t_type>(3);
    m_inv = MemAlloc<t_type>(9);

    memcpy(m_A, m.m_elem, sizeof(t_type) * 9);

    // Compute SVD
    HouseholdersReductionToBidiagonalForm_Mxn0();
    if (GivensReductionToDiagonalForm_Mxn0())
    {
        m_isOk = 0;
        return;
    }
    SortByDecreasingSingularValues_Mxn0();
    m_isOk = 1;
}

/** Initializes the SVD object
 *
 * Function Arguments:
 *
 * m  -> matrix in the dynamic memory Matrix class form
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline SVD<0, 0, t_type>::SVD(const Matrix<0, 0, t_type> &m)
{
    m_row = m.m_row;
    m_col = m.m_col;
    m_size = m_row * m_col;

    m_A = MemAlloc<t_type>(m_size);
    m_U = MemAlloc<t_type>((m_row >= m_col) ? m_size : m_col * m_col);
    m_S = MemAlloc<t_type>((m_row >= m_col) ? m_col : m_row);
    m_V = MemAlloc<t_type>((m_row >= m_col) ? m_col * m_col : m_size);
    m_superDiagonal = MemAlloc<t_type>((m_row >= m_col) ? m_row : m_col);
    m_inv = MemAlloc<t_type>(m_size);

    if (m_row >= m_col)
    {
        memcpy(m_A, m.m_elem, sizeof(t_type) * m_size);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn0();
        if (GivensReductionToDiagonalForm_Mxn0())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_Mxn0();
        m_isOk = 1;
    }
    else
    { // Transposed A
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < m_row; ++irow)
        {
            for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
                m_A[irow + m_row * (icol + 1)] = m.m_elem[irow * m_col + icol + 1];
                m_A[irow + m_row * (icol + 2)] = m.m_elem[irow * m_col + icol + 2];
                m_A[irow + m_row * (icol + 3)] = m.m_elem[irow * m_col + icol + 3];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN0();
        if (GivensReductionToDiagonalForm_mxN0())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_mxN0();
        m_isOk = 1;
    }
}

/** Deconstructor for the SVD object
 */
template <typename t_type>
inline SVD<0, 0, t_type>::~SVD()
{
    if (m_A)
    {
        MemFree(m_A);
        MemFree(m_U);
        MemFree(m_S);
        MemFree(m_V);
        MemFree(m_superDiagonal);
        MemFree(m_inv);

        m_A = nullptr;
        m_U = nullptr;
        m_S = nullptr;
        m_V = nullptr;
        m_superDiagonal = nullptr;
        m_inv = nullptr;
    }
}

/** Compute the SVD Factorization
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
inline int8_t SVD<0, 0, t_type>::Compute(const uint16_t row, const uint16_t col, const t_type *element)
{
    assert(row != 0 && "row and col must be non zero");
    assert(col != 0 && "row and col must be non zero");

    if (m_row != row || m_col != col)
    {
        if (m_A) MemFree(m_A);
        if (m_U) MemFree(m_U);
        if (m_S) MemFree(m_S);
        if (m_V) MemFree(m_V);
        if (m_superDiagonal) MemFree(m_superDiagonal);
        if (m_inv) MemFree(m_inv);

        m_row = row;
        m_col = col;
        m_size = row * col;
        m_A = MemAlloc<t_type>(m_size);
        m_U = MemAlloc<t_type>((m_row >= m_col) ? m_size : m_col * m_col);
        m_S = MemAlloc<t_type>((m_row >= m_col) ? m_col : m_row);
        m_V = MemAlloc<t_type>((m_row >= m_col) ? m_col * m_col : m_size);
        m_superDiagonal = MemAlloc<t_type>((m_row >= m_col) ? m_row : m_col);
        m_inv = MemAlloc<t_type>(m_size);
    }

    memcpy(m_A, element, sizeof(t_type) * m_size);
    memset(m_U, 0, sizeof(t_type) * ((m_row >= m_col) ? m_size : m_col * m_col));
    memset(m_S, 0, sizeof(t_type) * ((m_row >= m_col) ? m_col : m_row));
    memset(m_V, 0, sizeof(t_type) * ((m_row >= m_col) ? m_col * m_col : m_size));
    memset(m_superDiagonal, 0, sizeof(t_type) * ((m_row >= m_col) ? m_row : m_col));
    memset(m_inv, 0, sizeof(t_type) * m_size);

    if (m_row >= m_col)
    {
        memcpy(m_A, element, sizeof(t_type) * m_size);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn0();
        if (GivensReductionToDiagonalForm_Mxn0())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_Mxn0();
        m_isOk = 1;
    }
    else
    { // Transposed A
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < m_row; ++irow)
        {
            for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + m_row * icol] = element[irow * m_col + icol];
                m_A[irow + m_row * (icol + 1)] = element[irow * m_col + icol + 1];
                m_A[irow + m_row * (icol + 2)] = element[irow * m_col + icol + 2];
                m_A[irow + m_row * (icol + 3)] = element[irow * m_col + icol + 3];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + m_row * icol] = element[irow * m_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN0();
        if (GivensReductionToDiagonalForm_mxN0())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_mxN0();
        m_isOk = 1;
    }

    return 0;
}

/** Compute the SVD factorization
 *
 * Function Arguments:
 *
 * m   ->  Matrix object
 *
 * \return    -> 0 on success, -1 on failure
 */

template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t SVD<0, 0, t_type>::Compute(const Matrix<row, col, t_type> &m)
{
    assert(row != 0 && "row and col must be non zero");
    assert(col != 0 && "row and col must be non zero");

    if (m_row != row || m_col != col)
    {
        if (m_A) MemFree(m_A);
        if (m_U) MemFree(m_U);
        if (m_S) MemFree(m_S);
        if (m_V) MemFree(m_V);
        if (m_superDiagonal) MemFree(m_superDiagonal);
        if (m_inv) MemFree(m_inv);

        m_row = row;
        m_col = col;
        m_size = m_row * m_col;
        m_A = MemAlloc<t_type>(m_size);
        m_U = MemAlloc<t_type>((m_row >= m_col) ? m_size : m_col * m_col);
        m_S = MemAlloc<t_type>((m_row >= m_col) ? m_col : m_row);
        m_V = MemAlloc<t_type>((m_row >= m_col) ? m_col * m_col : m_size);
        m_superDiagonal = MemAlloc<t_type>((m_row >= m_col) ? m_row : m_col);
        m_inv = MemAlloc<t_type>(m_size);
    }

    memcpy(m_A, m.m_elem, sizeof(t_type) * m_size);
    memset(m_U, 0, sizeof(t_type) * ((m_row >= m_col) ? m_size : m_col * m_col));
    memset(m_S, 0, sizeof(t_type) * ((m_row >= m_col) ? m_col : m_row));
    memset(m_V, 0, sizeof(t_type) * ((m_row >= m_col) ? m_col * m_col : m_size));
    memset(m_superDiagonal, 0, sizeof(t_type) * ((m_row >= m_col) ? m_row : m_col));
    memset(m_inv, 0, sizeof(t_type) * m_size);

    if (m_row >= m_col)
    {
        memcpy(m_A, m.m_elem, sizeof(t_type) * m_size);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn0();
        if (GivensReductionToDiagonalForm_Mxn0())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_Mxn0();
        m_isOk = 1;
    }
    else
    { // Transposed A
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < m_row; ++irow)
        {
            for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
                m_A[irow + m_row * (icol + 1)] = m.m_elem[irow * m_col + icol + 1];
                m_A[irow + m_row * (icol + 2)] = m.m_elem[irow * m_col + icol + 2];
                m_A[irow + m_row * (icol + 3)] = m.m_elem[irow * m_col + icol + 3];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN0();
        if (GivensReductionToDiagonalForm_mxN0())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_mxN0();
        m_isOk = 1;
    }

    return 0;
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
inline int8_t SVD<0, 0, t_type>::Compute(const Matrix3<t_type, 3, 3> &m)
{
    if (m_row != 3 || m_col != 3)
    {
        if (m_A) MemFree(m_A);
        if (m_U) MemFree(m_U);
        if (m_S) MemFree(m_S);
        if (m_V) MemFree(m_V);
        if (m_superDiagonal) MemFree(m_superDiagonal);
        if (m_inv) MemFree(m_inv);

        m_row = 3;
        m_col = 3;
        m_size = 9;
        m_A = MemAlloc<t_type>(9);
        m_U = MemAlloc<t_type>(9);
        m_S = MemAlloc<t_type>(3);
        m_V = MemAlloc<t_type>(9);
        m_superDiagonal = MemAlloc<t_type>(3);
        m_inv = MemAlloc<t_type>(9);
    }

    memcpy(m_A, m.m_elem, sizeof(t_type) * 9);
    memset(m_U, 0, sizeof(t_type) * 9);
    memset(m_S, 0, sizeof(t_type) * 3);
    memset(m_V, 0, sizeof(t_type) * 9);
    memset(m_superDiagonal, 0, sizeof(t_type) * 3);
    memset(m_inv, 0, sizeof(t_type) * 9);

    if (m_row >= m_col)
    {
        memcpy(m_A, m.m_elem, sizeof(t_type) * m_size);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn0();
        if (GivensReductionToDiagonalForm_Mxn0())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_Mxn0();
        m_isOk = 1;
    }
    else
    { // Transposed A
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < m_row; ++irow)
        {
            for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
                m_A[irow + m_row * (icol + 1)] = m.m_elem[irow * m_col + icol + 1];
                m_A[irow + m_row * (icol + 2)] = m.m_elem[irow * m_col + icol + 2];
                m_A[irow + m_row * (icol + 3)] = m.m_elem[irow * m_col + icol + 3];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN0();
        if (GivensReductionToDiagonalForm_mxN0())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_mxN0();
        m_isOk = 1;
    }

    return 0;
}

/** Compute the SVD factorization
 *
 * Function Arguments:
 *
 * m   ->  dynamic memory Matrix object
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t SVD<0, 0, t_type>::Compute(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_row != 0 && "row and col must be non zero");
    assert(m.m_col != 0 && "row and col must be non zero");

    if (m_row != m.m_row || m_col != m.m_col)
    {
        if (m_A) MemFree(m_A);
        if (m_U) MemFree(m_U);
        if (m_S) MemFree(m_S);
        if (m_V) MemFree(m_V);
        if (m_superDiagonal) MemFree(m_superDiagonal);
        if (m_inv) MemFree(m_inv);

        m_row = m.m_row;
        m_col = m.m_col;
        m_size = m_row * m_col;
        m_A = MemAlloc<t_type>(m_size);
        m_U = MemAlloc<t_type>((m_row >= m_col) ? m_size : m_col * m_col);
        m_S = MemAlloc<t_type>((m_row >= m_col) ? m_col : m_row);
        m_V = MemAlloc<t_type>((m_row >= m_col) ? m_col * m_col : m_size);
        m_superDiagonal = MemAlloc<t_type>((m_row >= m_col) ? m_row : m_col);
        m_inv = MemAlloc<t_type>(m_size);
    }

    memcpy(m_A, m.m_elem, sizeof(t_type) * m_size);
    memset(m_U, 0, sizeof(t_type) * ((m_row >= m_col) ? m_size : m_col * m_col));
    memset(m_S, 0, sizeof(t_type) * ((m_row >= m_col) ? m_col : m_row));
    memset(m_V, 0, sizeof(t_type) * ((m_row >= m_col) ? m_col * m_col : m_size));
    memset(m_superDiagonal, 0, sizeof(t_type) * ((m_row >= m_col) ? m_row : m_col));
    memset(m_inv, 0, sizeof(t_type) * m_size);

    if (m_row >= m_col)
    {
        memcpy(m_A, m.m_elem, sizeof(t_type) * m_size);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn0();
        if (GivensReductionToDiagonalForm_Mxn0())
        {
            m_isOk = 0;
            return 0;
        }
        SortByDecreasingSingularValues_Mxn0();
        m_isOk = 1;
    }
    else
    { // Transposed A
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < m_row; ++irow)
        {
            for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
                m_A[irow + m_row * (icol + 1)] = m.m_elem[irow * m_col + icol + 1];
                m_A[irow + m_row * (icol + 2)] = m.m_elem[irow * m_col + icol + 2];
                m_A[irow + m_row * (icol + 3)] = m.m_elem[irow * m_col + icol + 3];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + m_row * icol] = m.m_elem[irow * m_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN0();
        if (GivensReductionToDiagonalForm_mxN0())
        {
            m_isOk = 0;
            return 0;
        }
        SortByDecreasingSingularValues_mxN0();
        m_isOk = 1;
    }

    return 0;
}

/** Get the U matrix after factorization
 *
 * Get the U matrix after the A = USV^T factorization
 *
 * \return    -> U matrix in the A = USV^T equation in dynamic Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> SVD<0, 0, t_type>::GetMatrixU() const
{
    if (m_row > m_col)
    {
        t_type U[m_row * m_row];
        memset(U, 0, sizeof(U));

        for (int i = 0; i < m_row; i++)
            memcpy(&U[i * m_row], &m_U[i * m_col], sizeof(t_type) * m_col);

        return Matrix<0, 0, t_type>(m_row, m_row, U);
    }
    else
    {
        return Matrix<0, 0, t_type>(m_row, m_row, m_U);
    }
}

/** Get the S matrix after factorization
 *
 * Get the S matrix after the A = USV^T factorization
 *
 * \return    -> S matrix in the A = USV^T equation in dynamic Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> SVD<0, 0, t_type>::GetMatrixS() const
{
    t_type S[m_row * m_col];
    memset(S, 0, sizeof(S));

    for (int i = 0; i < m_row; i++)
    {
        S[i * m_col + i] = m_S[i];
        if (i == m_col) break;
    }

    return Matrix<0, 0, t_type>(m_row, m_col, S);
}

/** Get the V matrix after factorization
 *
 * Get the V matrix after the A = USV^T factorization
 *
 * \return    -> V matrix in the A = USV^T equation in dynamic Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> SVD<0, 0, t_type>::GetMatrixV() const
{
    if (m_row < m_col)
    {
        t_type V[m_col * m_col];
        memset(V, 0, sizeof(V));

        for (int i = 0; i < m_col; i++)
            memcpy(&V[i * m_col], &m_V[i * m_row], sizeof(t_type) * m_row);

        return Matrix<0, 0, t_type>(m_col, m_col, V);
    }
    else
    {
        return Matrix<0, 0, t_type>(m_col, m_col, m_V);
    }
}

/** Perform the x = (USVT)^-1 * b solve operation
 *
 * Function Arguments:
 *
 * b ->  The 'b' matrix viz a set of vectors in the x = (USVT)^-1 * b equation in dynamic Matrix class
 * x ->  The 'x' matrix viz a set of vectors in the x = (USVT)^-1 * b equation in dynamic Matrix class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t SVD<0, 0, t_type>::Solve(const Matrix<0, 0, t_type> &b, Matrix<0, 0, t_type> &x, t_type tolerance)
{
    assert(b.m_col == x.m_col && "Rows of b and x do not match");
    assert(m_row == b.m_row && "Check size of matrix b");
    assert(m_col == x.m_row && "Check size of matrix x");

    int i, j, k, c;
    t_type *pU, *pV;
    t_type sum;
    int col = x.m_col;

    if (!m_isOk)
    {
        assert(false && "SVD factorization not computed");
        return -1;
    }

    if (m_row >= m_col)
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_col;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < m_col; i++, pV += m_col)
            {
                x.m_elem[i * col + c] = 0;
                for (j = 0; j < m_col; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_col)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x.m_elem[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_row;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < m_col; i++, pV += m_row)
            {
                x.m_elem[i * col + c] = 0;
                for (j = 0; j < m_row; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_row)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x.m_elem[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }

    return 0;
}

/** Perform the x = (USVT)^-1 * b solve operation
 *
 * Function Arguments:
 *
 * b ->  The 'b' vector in the x = (USVT)^-1 * b equation in dynamic Matrix class
 * x ->  The 'x' vector in the x = (USVT)^-1 * b equation in dynamic Matrix class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t SVD<0, 0, t_type>::Solve(const Vector<0, t_type> &b, Vector<0, t_type> &x, t_type tolerance)
{
    assert(m_row == b.m_row && "Check size of matrix b");
    assert(m_col == x.m_row && "Check size of matrix x");

    int i, j, k;
    t_type *pU, *pV;
    t_type sum;

    if (!m_isOk)
    {
        assert(false && "SVD factorization not computed");
        return -1;
    }

    if (m_row >= m_col)
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_col;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < m_col; i++, pV += m_col)
        {
            x.m_elem[i] = 0;
            for (j = 0; j < m_col; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_col)
                        sum += *(pU + j) * b.m_elem[k];

                    x.m_elem[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_row;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < m_col; i++, pV += m_row)
        {
            x.m_elem[i] = 0;
            for (j = 0; j < m_row; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_row)
                        sum += *(pU + j) * b.m_elem[k];

                    x.m_elem[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }

    return 0;
}

/** Perform the x = (USVT)^-1 * b solve operation
 *
 * Function Arguments:
 *
 * b ->  The 'b' matrix viz a set of vectors in the x = (USVT)^-1 * b equation in dynamic Matrix class
 *
 * \return    -> The 'x' matrix viz a set of vectors in the x = (USVT)^-1 * b equation in dynamic Matrix class
 */
template <typename t_type>
inline Matrix<0, 0, t_type> SVD<0, 0, t_type>::Solve(const Matrix<0, 0, t_type> &b, int8_t *isOk, t_type tolerance)
{
    assert(m_row == b.m_row && "Check size of matrix b");

    if (isOk)
        *isOk = 1;

    int i, j, k, c;
    t_type *pU, *pV;
    t_type sum;
    int col = b.m_col;
    Matrix<0, 0, t_type> x(m_col, col);
    // t_type x[m_col * col] = {
    //     0,
    // };
    if (!m_isOk && isOk)
    {
        *isOk = 0;
        assert(false && "SVD factorization not computed");
        return Matrix<0, 0, t_type>(m_col, col);
    }

    if (m_row >= m_col)
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_col;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < m_col; i++, pV += m_col)
            {
                x.m_elem[i * col + c] = 0;
                for (j = 0; j < m_col; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_col)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x.m_elem[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_row;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < m_col; i++, pV += m_row)
            {
                x.m_elem[i * col + c] = 0;
                for (j = 0; j < m_row; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < m_row; k++, pU += m_row)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x.m_elem[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }

    return x;
}

/** Perform the x = (USVT)^-1 * b solve operation
 *
 * Function Arguments:
 *
 * b ->  The 'b' vector in the x = (USVT)^-1 * b equation in Matrix class
 *
 * \return    -> The 'x' vector in the x = (USVT)^-1 * b equation in dynamic Matrix class
 */
template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> SVD<0, 0, t_type>::Solve(const Vector<row, t_type> &b, int8_t *isOk, t_type tolerance)
{
    assert(m_row == row && "Check size of matrix b");

    int i, j, k;
    t_type *pU, *pV;
    t_type sum;
    Vector<0, t_type> x(m_col);

    if (!m_isOk)
    {
        assert(false && "SVD factorization not computed");
        return Vector<0, t_type>();
    }

    if (row >= m_col)
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_col;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < m_col; i++, pV += m_col)
        {
            x.m_elem[i] = 0;
            for (j = 0; j < m_col; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < row; k++, pU += m_col)
                        sum += *(pU + j) * b.m_elem[k];

                    x.m_elem[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)row;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < m_col; i++, pV += row)
        {
            x.m_elem[i] = 0;
            for (j = 0; j < row; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < row; k++, pU += row)
                        sum += *(pU + j) * b.m_elem[k];

                    x.m_elem[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }

    return x;
}

/** Perform the x = (USVT)^-1 * b solve operation
 *
 * Function Arguments:
 *
 * b ->  The 'b' vector in the x = (USVT)^-1 * b equation in dynamic Matrix class
 *
 * \return    -> The 'x' vector in the x = (USVT)^-1 * b equation in dynamic Matrix class
 */
template <typename t_type>
Vector<0, t_type> SVD<0, 0, t_type>::Solve(const Vector<0, t_type> &b, int8_t *isOk, t_type tolerance)
{
    assert(m_row == b.m_row && "Check size of matrix b");

    int row = b.m_row;
    int i, j, k;
    t_type *pU, *pV;
    t_type sum;
    Vector<0, t_type> x(m_col);

    if (!m_isOk)
    {
        assert(false && "SVD factorization not computed");
        return Vector<0, t_type>();
    }

    if (row >= m_col)
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_col;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < m_col; i++, pV += m_col)
        {
            x.m_elem[i] = 0;
            for (j = 0; j < m_col; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < row; k++, pU += m_col)
                        sum += *(pU + j) * b.m_elem[k];

                    x.m_elem[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)row;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < m_col; i++, pV += row)
        {
            x.m_elem[i] = 0;
            for (j = 0; j < row; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < row; k++, pU += row)
                        sum += *(pU + j) * b.m_elem[k];

                    x.m_elem[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }

    return x;
}

/** Pseudoinverse of the A matrix
 *
 * Function Arguments:
 *
 * inv    -> the output inverse matrix in dynamic Matrix class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <typename t_type>
inline int8_t SVD<0, 0, t_type>::Inverse(Matrix<0, 0, t_type> &inv, t_type tolerance)
{
    assert(m_row == inv.m_col && "Check size of matrix");
    assert(m_col == inv.m_row && "Check size of matrix");

    int i, j, k;
    t_type *pU, *pV, *pInvA;
    t_type tol;
    inv.m_row = m_col;
    inv.m_col = m_row;

    if (!m_isOk)
    {
        assert(false && "SVD factorization not computed");
        return -1;
    }

    if (m_row >= m_col)
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_col;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv.m_elem; i < m_col; i++, pV += m_col)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_col; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }
    else
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_row;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv.m_elem; i < m_col; i++, pV += m_row)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_row; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }

    return 0;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t SVD<0, 0, t_type>::Inverse(Matrix<row, col, t_type> &inv, t_type tolerance) // Inverse matrix of USVT matrix
{
    assert(m_row == col && "Check size of matrix");
    assert(m_col == row && "Check size of matrix");

    int i, j, k;
    t_type *pU, *pV, *pInvA;
    t_type tol;

    if (!m_isOk)
    {
        assert(false && "SVD factorization not computed");
        return -1;
    }

    if (m_row >= m_col)
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_col;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv.m_elem; i < m_col; i++, pV += m_col)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_col; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }
    else
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_row;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv.m_elem; i < m_col; i++, pV += m_row)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_row; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }

    return 0;
}

template <typename t_type>
inline Matrix<0, 0, t_type> SVD<0, 0, t_type>::Inverse(int8_t *isOk, t_type tolerance)
{
    int i, j, k;
    t_type *pU, *pV, *pInvA;
    t_type tol;

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        assert(false && "SVD factorization not computed");
        return Matrix<0, 0, t_type>(m_col, m_row);
    }

    memset(m_inv, 0, sizeof(t_type) * m_row * m_col);

    if (m_row >= m_col)
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_col;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = m_inv; i < m_col; i++, pV += m_col)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_col; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }
    else
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)m_row;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = m_inv; i < m_col; i++, pV += m_row)
        {
            for (j = 0, pU = m_U; j < m_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < m_row; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }

    return Matrix<0, 0, t_type>(m_col, m_row, m_inv);
}

/**** Supporting Codes for Computing SVD ****/
template <typename t_type>
inline void SVD<0, 0, t_type>::HouseholdersReductionToBidiagonalForm_Mxn0()
{
    int i, j, k, ip1;
    t_type s, s2, si, scale;
    t_type *pU, *pUi, *pV, *pVi;
    t_type halfNormSquared;

    // Copy A to U
    memcpy(m_U, m_A, sizeof(t_type) * m_row * m_col);

    m_S[0] = 0;
    s = 0;
    scale = 0;

    for (i = 0, pUi = m_U, ip1 = 1; i < m_col; pUi += m_col, i++, ip1++)
    {
        m_superDiagonal[i] = scale * s;

        // Perform Householder transform on columns.
        // Calculate the normed squared of the i-th column vector starting at row i.
        for (j = i, pU = pUi, scale = 0; j < m_row; j++, pU += m_col)
            scale += std::abs(*(pU + i));

        if (scale > std::numeric_limits<t_type>::epsilon())
        {
            for (j = i, pU = pUi, s2 = 0; j < m_row; j++, pU += m_col)
            {
                *(pU + i) /= scale;
                s2 += *(pU + i) * *(pU + i);
            }

            // Chose sign of s which maximizes the norm
            s = (*(pUi + i) < 0) ? std::sqrt(s2) : -std::sqrt(s2);

            // Calculate -2/u'u
            halfNormSquared = *(pUi + i) * s - s2;

            // Transform remaining columns by the Householder transform.
            *(pUi + i) -= s;

            for (j = ip1; j < m_col; j++)
            {
                for (k = i, si = 0, pU = pUi; k < m_row; k++, pU += m_col)
                    si += *(pU + i) * *(pU + j);

                si /= halfNormSquared;

                for (k = i, pU = pUi; k < m_row; k++, pU += m_col)
                    *(pU + j) += si * *(pU + i);
            }
        }

        for (j = i, pU = pUi; j < m_row; j++, pU += m_col)
            *(pU + i) *= scale;

        m_S[i] = s * scale;

        // Perform Householder transform on rows.
        // Calculate the normed squared of the i-th row vector starting at column i.
        s = 0;
        scale = 0;
        if (i >= m_row || i == (m_col - 1))
            continue;

        for (j = ip1; j < m_col; j++)
            scale += std::abs(*(pUi + j));

        if (scale > std::numeric_limits<t_type>::epsilon())
        {
            for (j = ip1, s2 = 0; j < m_col; j++)
            {
                *(pUi + j) /= scale;
                s2 += *(pUi + j) * *(pUi + j);
            }

            // Chose sign of s which maximizes the norm
            s = (*(pUi + ip1) < 0) ? std::sqrt(s2) : -std::sqrt(s2);

            // Calculate -2/u'u
            halfNormSquared = *(pUi + ip1) * s - s2;

            // Transform the rows by the Householder transform.
            *(pUi + ip1) -= s;

            for (k = ip1; k < m_col; k++)
                m_superDiagonal[k] = *(pUi + k) / halfNormSquared;

            if (i < (m_row - 1))
            {
                for (j = ip1, pU = pUi + m_col; j < m_row; j++, pU += m_col)
                {
                    for (k = ip1, si = 0; k < m_col; k++)
                        si += *(pUi + k) * *(pU + k);

                    for (k = ip1; k < m_col; k++)
                        *(pU + k) += si * m_superDiagonal[k];
                }
            }
            for (k = ip1; k < m_col; k++)
                *(pUi + k) *= scale;
        }
    }

    /* Update V */
    pUi = m_U + m_col * (m_col - 2);
    pVi = m_V + m_col * (m_col - 1);
    *(pVi + m_col - 1) = 1;
    s = m_superDiagonal[m_col - 1];
    pVi -= m_col;

    for (i = m_col - 2, ip1 = m_col - 1; i >= 0; i--, pUi -= m_col, pVi -= m_col, ip1--)
    {
        if (std::abs(s) > std::numeric_limits<t_type>::epsilon())
        {
            pV = pVi + m_col;

            for (j = ip1; j < m_col; j++, pV += m_col)
                *(pV + i) = (*(pUi + j) / *(pUi + ip1)) / s;

            for (j = ip1; j < m_col; j++)
            {
                si = 0;

                for (k = ip1, pV = pVi + m_col; k < m_col; k++, pV += m_col)
                    si += *(pUi + k) * *(pV + j);

                for (k = ip1, pV = pVi + m_col; k < m_col; k++, pV += m_col)
                    *(pV + j) += si * *(pV + i);
            }
        }

        pV = pVi + m_col;
        for (j = ip1; j < m_col; j++, pV += m_col)
        {
            *(pVi + j) = 0;
            *(pV + i) = 0;
        }

        *(pVi + i) = 1;
        s = m_superDiagonal[i];
    }

    /* Update U */
    pUi = m_U + m_col * (m_col - 1);
    for (i = m_col - 1, ip1 = m_col; i >= 0; ip1 = i, i--, pUi -= m_col)
    {
        s = m_S[i];

        for (j = ip1; j < m_col; j++)
            *(pUi + j) = 0;

        if (std::abs(s) > std::numeric_limits<t_type>::epsilon())
        {
            for (j = ip1; j < m_col; j++)
            {
                si = 0;
                pU = pUi + m_col;

                for (k = ip1; k < m_row; k++, pU += m_col)
                    si += *(pU + i) * *(pU + j);

                si = (si / *(pUi + i)) / s;
                for (k = i, pU = pUi; k < m_row; k++, pU += m_col)
                    *(pU + j) += si * *(pU + i);
            }

            for (j = i, pU = pUi; j < m_row; j++, pU += m_col)
                *(pU + i) /= s;
        }
        else
        {
            for (j = i, pU = pUi; j < m_row; j++, pU += m_col)
                *(pU + i) = 0;
        }

        *(pUi + i) += 1;
    }
}

template <typename t_type>
inline int8_t SVD<0, 0, t_type>::GivensReductionToDiagonalForm_Mxn0()
{
    t_type epsilon;
    t_type c, s;
    t_type f, g, h;
    t_type x, y, z;
    t_type *pU, *pV;
    int i, j, k, m;
    int rotation_test;
    int iteration_count;

    for (i = 0, x = 0; i < m_col; i++)
    {
        y = std::abs(m_S[i]) + std::abs(m_superDiagonal[i]);
        if (x < y)
            x = y;
    }

    epsilon = x * std::numeric_limits<t_type>::epsilon();

    for (k = m_col - 1; k >= 0; k--)
    {
        iteration_count = 0;
        while (1)
        {
            rotation_test = 1;
            for (m = k; m >= 0; m--)
            {
                if (std::abs(m_superDiagonal[m]) <= epsilon)
                {
                    rotation_test = 0;
                    break;
                }

                if (std::abs(m_S[m - 1]) <= epsilon)
                    break;
            }

            if (rotation_test)
            {
                c = 0;
                s = 1;
                for (i = m; i <= k; i++)
                {
                    f = s * m_superDiagonal[i];
                    m_superDiagonal[i] *= c;

                    if (std::abs(f) <= epsilon)
                        break;

                    g = m_S[i];
                    h = std::sqrt(f * f + g * g);
                    m_S[i] = h;
                    c = g / h;
                    s = -f / h;

                    for (j = 0, pU = m_U; j < m_row; j++, pU += m_col)
                    {
                        y = *(pU + m - 1);
                        z = *(pU + i);
                        *(pU + m - 1) = y * c + z * s;
                        *(pU + i) = -y * s + z * c;
                    }
                }
            }

            z = m_S[k];

            if (m == k)
            {
                if (z < 0)
                {
                    m_S[k] = -z;
                    for (j = 0, pV = m_V; j < m_col; j++, pV += m_col)
                        *(pV + k) = -*(pV + k);
                }
                break;
            }
            else
            {
                if (iteration_count >= MAX_ITERATION_COUNT)
                    return -1;

                iteration_count++;
                x = m_S[m];
                y = m_S[k - 1];
                g = m_superDiagonal[k - 1];
                h = m_superDiagonal[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
                g = std::sqrt(f * f + 1);

                if (f < 0)
                    g = -g;

                f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
                // Next QR Transformtion
                c = 1;
                s = 1;

                for (i = m + 1; i <= k; i++)
                {
                    g = m_superDiagonal[i];
                    y = m_S[i];
                    h = s * g;
                    g *= c;
                    z = std::sqrt(f * f + h * h);
                    m_superDiagonal[i - 1] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = -x * s + g * c;
                    h = y * s;
                    y *= c;

                    for (j = 0, pV = m_V; j < m_col; j++, pV += m_col)
                    {
                        x = *(pV + i - 1);
                        z = *(pV + i);
                        *(pV + i - 1) = x * c + z * s;
                        *(pV + i) = -x * s + z * c;
                    }

                    z = std::sqrt(f * f + h * h);
                    m_S[i - 1] = z;

                    if (std::abs(z) > std::numeric_limits<t_type>::epsilon())
                    {
                        c = f / z;
                        s = h / z;
                    }

                    f = c * g + s * y;
                    x = -s * g + c * y;

                    for (j = 0, pU = m_U; j < m_row; j++, pU += m_col)
                    {
                        y = *(pU + i - 1);
                        z = *(pU + i);
                        *(pU + i - 1) = c * y + s * z;
                        *(pU + i) = -s * y + c * z;
                    }
                }
                m_superDiagonal[m] = 0;
                m_superDiagonal[k] = f;
                m_S[k] = x;
            }
        }
    }

    return 0;
}

template <typename t_type>
inline void SVD<0, 0, t_type>::SortByDecreasingSingularValues_Mxn0()
{
    int i, j, maxIdx;
    t_type temp;
    t_type *pM1, *pM2;

    for (i = 0; i < m_col - 1; i++)
    {
        maxIdx = i;

        for (j = i + 1; j < m_col; j++)
        {
            if (m_S[j] > m_S[maxIdx])
                maxIdx = j;
        }

        if (maxIdx == i)
            continue;

        temp = m_S[i];
        m_S[i] = m_S[maxIdx];
        m_S[maxIdx] = temp;
        pM1 = m_U + maxIdx;
        pM2 = m_U + i;

        for (j = 0; j < m_row; j++, pM1 += m_col, pM2 += m_col)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }

        pM1 = m_V + maxIdx;
        pM2 = m_V + i;

        for (j = 0; j < m_col; j++, pM1 += m_col, pM2 += m_col)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }
    }
}

template <typename t_type>
inline void SVD<0, 0, t_type>::HouseholdersReductionToBidiagonalForm_mxN0()
{
    int i, j, k, ip1;
    t_type s, s2, si, scale;
    t_type *pV, *pVi, *pU, *pUi;
    t_type halfNormSquared;

    // Copy A to V
    memcpy(m_V, m_A, sizeof(t_type) * m_col * m_row);

    m_S[0] = 0;
    s = 0;
    scale = 0;

    for (i = 0, pVi = m_V, ip1 = 1; i < m_row; pVi += m_row, i++, ip1++)
    {
        m_superDiagonal[i] = scale * s;

        // Perform Householder transform on columns.
        // Calculate the normed squared of the i-th column vector starting at row i.
        for (j = i, pV = pVi, scale = 0; j < m_col; j++, pV += m_row)
            scale += std::abs(*(pV + i));

        if (scale > std::numeric_limits<t_type>::epsilon())
        {
            for (j = i, pV = pVi, s2 = 0; j < m_col; j++, pV += m_row)
            {
                *(pV + i) /= scale;
                s2 += *(pV + i) * *(pV + i);
            }

            // Chose sign of s which maximizes the norm
            s = (*(pVi + i) < 0) ? std::sqrt(s2) : -std::sqrt(s2);

            // Calculate -2/v'v
            halfNormSquared = *(pVi + i) * s - s2;

            // Transform remaining columns by the Householder transform.
            *(pVi + i) -= s;

            for (j = ip1; j < m_row; j++)
            {
                for (k = i, si = 0, pV = pVi; k < m_col; k++, pV += m_row)
                    si += *(pV + i) * *(pV + j);

                si /= halfNormSquared;

                for (k = i, pV = pVi; k < m_col; k++, pV += m_row)
                    *(pV + j) += si * *(pV + i);
            }
        }
        for (j = i, pV = pVi; j < m_col; j++, pV += m_row)
            *(pV + i) *= scale;
        m_S[i] = s * scale;

        // Perform Householder transform on rows.
        // Calculate the normed squared of the i-th row vector starting at column i.
        s = 0;
        scale = 0;
        if (i >= m_col || i == (m_row - 1))
            continue;

        for (j = ip1; j < m_row; j++)
            scale += std::abs(*(pVi + j));

        if (scale > std::numeric_limits<t_type>::epsilon())
        {
            for (j = ip1, s2 = 0; j < m_row; j++)
            {
                *(pVi + j) /= scale;
                s2 += *(pVi + j) * *(pVi + j);
            }

            // Chose sign of s which maximizes the norm
            s = (*(pVi + ip1) < 0) ? std::sqrt(s2) : -std::sqrt(s2);

            // Calculate -2/v'v
            halfNormSquared = *(pVi + ip1) * s - s2;

            // Transform the rows by the Householder transform.
            *(pVi + ip1) -= s;

            for (k = ip1; k < m_row; k++)
                m_superDiagonal[k] = *(pVi + k) / halfNormSquared;

            if (i < (m_col - 1))
            {
                for (j = ip1, pV = pVi + m_row; j < m_col; j++, pV += m_row)
                {
                    for (k = ip1, si = 0; k < m_row; k++)
                        si += *(pVi + k) * *(pV + k);

                    for (k = ip1; k < m_row; k++)
                        *(pV + k) += si * m_superDiagonal[k];
                }
            }
            for (k = ip1; k < m_row; k++)
                *(pVi + k) *= scale;
        }
    }

    /* Update U */
    pVi = m_V + m_row * (m_row - 2);
    pUi = m_U + m_row * (m_row - 1);
    *(pUi + m_row - 1) = 1;
    s = m_superDiagonal[m_row - 1];
    pUi -= m_row;

    for (i = m_row - 2, ip1 = m_row - 1; i >= 0; i--, pVi -= m_row, pUi -= m_row, ip1--)
    {
        if (std::abs(s) > std::numeric_limits<t_type>::epsilon())
        {
            pU = pUi + m_row;

            for (j = ip1; j < m_row; j++, pU += m_row)
                *(pU + i) = (*(pVi + j) / *(pVi + ip1)) / s;

            for (j = ip1; j < m_row; j++)
            {
                si = 0;

                for (k = ip1, pU = pUi + m_row; k < m_row; k++, pU += m_row)
                    si += *(pVi + k) * *(pU + j);

                for (k = ip1, pU = pUi + m_row; k < m_row; k++, pU += m_row)
                    *(pU + j) += si * *(pU + i);
            }
        }

        pU = pUi + m_row;
        for (j = ip1; j < m_row; j++, pU += m_row)
        {
            *(pUi + j) = 0;
            *(pU + i) = 0;
        }

        *(pUi + i) = 1;
        s = m_superDiagonal[i];
    }

    /* Update V */
    pVi = m_V + m_row * (m_row - 1);
    for (i = m_row - 1, ip1 = m_row; i >= 0; ip1 = i, i--, pVi -= m_row)
    {
        s = m_S[i];

        for (j = ip1; j < m_row; j++)
            *(pVi + j) = 0;

        if (std::abs(s) > std::numeric_limits<t_type>::epsilon())
        {
            for (j = ip1; j < m_row; j++)
            {
                si = 0;
                pV = pVi + m_row;

                for (k = ip1; k < m_col; k++, pV += m_row)
                    si += *(pV + i) * *(pV + j);

                si = (si / *(pVi + i)) / s;
                for (k = i, pV = pVi; k < m_col; k++, pV += m_row)
                    *(pV + j) += si * *(pV + i);
            }

            for (j = i, pV = pVi; j < m_col; j++, pV += m_row)
                *(pV + i) /= s;
        }
        else
        {
            for (j = i, pV = pVi; j < m_col; j++, pV += m_row)
                *(pV + i) = 0;
        }

        *(pVi + i) += 1;
    }
}

template <typename t_type>
inline int8_t SVD<0, 0, t_type>::GivensReductionToDiagonalForm_mxN0()
{
    t_type epsilon;
    t_type c, s;
    t_type f, g, h;
    t_type x, y, z;
    t_type *pV, *pU;
    int i, j, k, m;
    int rotation_test;
    int iteration_count;

    for (i = 0, x = 0; i < m_row; i++)
    {
        y = std::abs(m_S[i]) + std::abs(m_superDiagonal[i]);
        if (x < y)
            x = y;
    }

    epsilon = x * std::numeric_limits<t_type>::epsilon();

    for (k = m_row - 1; k >= 0; k--)
    {
        iteration_count = 0;
        while (1)
        {
            rotation_test = 1;
            for (m = k; m >= 0; m--)
            {
                if (std::abs(m_superDiagonal[m]) <= epsilon)
                {
                    rotation_test = 0;
                    break;
                }

                if (std::abs(m_S[m - 1]) <= epsilon)
                    break;
            }

            if (rotation_test)
            {
                c = 0;
                s = 1;
                for (i = m; i <= k; i++)
                {
                    f = s * m_superDiagonal[i];
                    m_superDiagonal[i] *= c;

                    if (std::abs(f) <= epsilon)
                        break;

                    g = m_S[i];
                    h = std::sqrt(f * f + g * g);
                    m_S[i] = h;
                    c = g / h;
                    s = -f / h;

                    for (j = 0, pV = m_V; j < m_col; j++, pV += m_row)
                    {
                        y = *(pV + m - 1);
                        z = *(pV + i);
                        *(pV + m - 1) = y * c + z * s;
                        *(pV + i) = -y * s + z * c;
                    }
                }
            }

            z = m_S[k];

            if (m == k)
            {
                if (z < 0)
                {
                    m_S[k] = -z;
                    for (j = 0, pU = m_U; j < m_row; j++, pU += m_row)
                        *(pU + k) = -*(pU + k);
                }
                break;
            }
            else
            {
                if (iteration_count >= MAX_ITERATION_COUNT)
                    return -1;

                iteration_count++;
                x = m_S[m];
                y = m_S[k - 1];
                g = m_superDiagonal[k - 1];
                h = m_superDiagonal[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
                g = std::sqrt(f * f + 1);

                if (f < 0)
                    g = -g;

                f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
                // Next QR Transformtion
                c = 1;
                s = 1;

                for (i = m + 1; i <= k; i++)
                {
                    g = m_superDiagonal[i];
                    y = m_S[i];
                    h = s * g;
                    g *= c;
                    z = std::sqrt(f * f + h * h);
                    m_superDiagonal[i - 1] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = -x * s + g * c;
                    h = y * s;
                    y *= c;

                    for (j = 0, pU = m_U; j < m_row; j++, pU += m_row)
                    {
                        x = *(pU + i - 1);
                        z = *(pU + i);
                        *(pU + i - 1) = x * c + z * s;
                        *(pU + i) = -x * s + z * c;
                    }

                    z = std::sqrt(f * f + h * h);
                    m_S[i - 1] = z;

                    if (std::abs(z) > std::numeric_limits<t_type>::epsilon())
                    {
                        c = f / z;
                        s = h / z;
                    }

                    f = c * g + s * y;
                    x = -s * g + c * y;

                    for (j = 0, pV = m_V; j < m_col; j++, pV += m_row)
                    {
                        y = *(pV + i - 1);
                        z = *(pV + i);
                        *(pV + i - 1) = c * y + s * z;
                        *(pV + i) = -s * y + c * z;
                    }
                }
                m_superDiagonal[m] = 0;
                m_superDiagonal[k] = f;
                m_S[k] = x;
            }
        }
    }

    return 0;
}

template <typename t_type>
inline void SVD<0, 0, t_type>::SortByDecreasingSingularValues_mxN0()
{
    int i, j, maxIdx;
    t_type temp;
    t_type *pM1, *pM2;

    for (i = 0; i < m_row - 1; i++)
    {
        maxIdx = i;

        for (j = i + 1; j < m_row; j++)
        {
            if (m_S[j] > m_S[maxIdx])
                maxIdx = j;
        }

        if (maxIdx == i)
            continue;

        temp = m_S[i];
        m_S[i] = m_S[maxIdx];
        m_S[maxIdx] = temp;
        pM1 = m_V + maxIdx;
        pM2 = m_V + i;

        for (j = 0; j < m_col; j++, pM1 += m_row, pM2 += m_row)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }

        pM1 = m_U + maxIdx;
        pM2 = m_U + i;

        for (j = 0; j < m_row; j++, pM1 += m_row, pM2 += m_row)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }
    }
}

} // namespace Math

} // namespace dt

#endif