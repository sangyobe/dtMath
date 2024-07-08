/*!
\file       dtLowerMatrix0.tpp
\brief      dtMath, Lower triangular matrix solver class with dynamic memory
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 05. 21
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLOWER_MATRIX0_TPP_
#define DTMATH_DTLOWER_MATRIX0_TPP_

#include "dtLowerMatrix.h"

namespace dt
{
namespace Math
{

/** Solve the Lx = b problem for x vector
 *
 * Function Arguments:
 *
 * L    -> lower triangular Matrx
 * b    -> the right hand side 'b' vector
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Vector<0, t_type> LowerMatrix<0, 0, t_type>::Solve(Matrix<row, col, t_type> &L, Vector<row, t_type> &b, int8_t *isOk)
{
    assert(L.m_elem != nullptr && "Memory has not been allocated to matrix 'L'");
    assert(b.m_elem != nullptr && "Memory has not been allocated to vector 'b'");

    Vector<0, t_type> x(row);
    if (row != col)
    {
        assert(false && "Check row and column dimensions of 'L' and 'b' do not match");
        if (isOk) *isOk = 0;
        return x;
    }

    t_type *pMi = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a lower triangular matrix.
    for (int i = 0; i < row; i++, pMi += col)
    {
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            if (isOk)
                *isOk = 0; // singular;
            assert(false && "The 'L' matrix is singular");
            return x;
        }

        x.m_elem[i] = b.m_elem[i];

        for (int k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    if (isOk)
        *isOk = 1;

    return x;
}

/** Solve the Lx = b problem for x vector
 *
 * Function Arguments:
 *
 * L    -> lower triangular Matrx
 * b    -> the right hand side 'b' vector in dynamic memory Vector class
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Vector<0, t_type> LowerMatrix<0, 0, t_type>::Solve(Matrix<row, col, t_type> &L, Vector<0, t_type> &b, int8_t *isOk)
{
    assert(L.m_elem != nullptr && "Memory has not been allocated to matrix 'L'");
    assert(col == b.m_row && "'L' Column dimensions and 'b' row dimensions do not matched");
    assert(b.m_elem != nullptr && "Memory has not been allocated to vector 'b'");

    Vector<0, t_type> x(row);
    if (row != col)
    {
        assert(false && "Row and column dimensions of 'L' do not match");
        if (isOk) *isOk = 0;
        return x;
    }

    t_type *pMi = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a lower triangular matrix.
    for (int i = 0; i < row; i++, pMi += col)
    {
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            if (isOk)
                *isOk = 0; // singular;
            assert(false && "The 'L' matrix is singular");
            return x;
        }

        x.m_elem[i] = b.m_elem[i];

        for (int k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    if (isOk)
        *isOk = 1;

    return x;
}

/** Solve the Lx = b problem for x vector
 *
 * Function Arguments:
 *
 * L    -> lower triangular Matrx in dynamic Matrix class
 * b    -> the right hand side 'b' vector in dynamic memory Vector class
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
inline Vector<0, t_type> LowerMatrix<0, 0, t_type>::Solve(Matrix<0, 0, t_type> &L, Vector<0, t_type> &b, int8_t *isOk)
{
    assert(L.m_elem != nullptr && "Memory has not been allocated to matrix 'L'");
    assert(L.m_col == b.m_row && "'L' Column dimensions and 'b' row dimensions do not matched");
    assert(b.m_elem != nullptr && "Memory has not been allocated to vector 'b'");

    int row = L.m_row;
    int col = L.m_col;

    Vector<0, t_type> x(row);
    if (row != L.m_col)
    {
        assert(false && "Row and column dimensions of 'L' do not match");
        if (isOk) *isOk = 0;
        return x;
    }

    t_type *pMi = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a lower triangular matrix.
    for (int i = 0; i < row; i++, pMi += col)
    {
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            if (isOk)
                *isOk = 0; // singular;
            assert(false && "The 'L' matrix is singular");
            return x;
        }

        x.m_elem[i] = b.m_elem[i];

        for (int k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    if (isOk)
        *isOk = 1;

    return x;
}

/** Solve the inverse of 'L' matrix
 *
 * Function Arguments:
 *
 * L    -> lower triangular Matrx in dynamic Matrix class
 * invL -> the matrix to store the inverse of U matrix in dynamic Matrix class
 *
 * \return    -> 0 on success; -1 on failure
 */

template <typename t_type>
int8_t LowerMatrix<0, 0, t_type>::Inverse(Matrix<0, 0, t_type> L, Matrix<0, 0, t_type> &invL)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type sum;

    int row = L.m_row;
    int col = L.m_col;

    assert(row == col && "Row and column dimensions of 'L' do not match");

    /* Invert the diagonal elements */
    for (k = 0, pMk = L.m_elem; k < row; pMk += (col + 1), k++)
    {
        if (*pMk == 0) // To do dhl
        {
            assert(false && "The 'L' matrix is singular");
            return -1;
        }
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix L, for row i */
    for (i = 1, pMi = L.m_elem + col; i < row; i++, pMi += col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    invL = L;

    return 0;
}

/** Solve the inverse of 'L' matrix
 *
 * Function Arguments:
 *
 * L    -> lower triangular Matrx in dynamic Matrix class
 * isOK -> the pointer flag for success
 *
 * \return    -> the inverse of L matrix in dynamic Matrix class
 */

template <typename t_type>
inline Matrix<0, 0, t_type> LowerMatrix<0, 0, t_type>::Inverse(Matrix<0, 0, t_type> L, int8_t *isOk)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type sum;

    int row = L.m_row;
    int col = L.m_col;

    if (row != col)
    {
        if (isOk)
            *isOk = 0;
        assert(false && "Row and column dimensions of 'L' do not match");
        return L;
    }

    /* Invert the diagonal elements */
    for (k = 0, pMk = L.m_elem; k < row; pMk += (col + 1), k++)
    {
        if (*pMk == 0) // To do dhl
        {
            if (isOk)
                *isOk = 0;
            assert(false && "The 'L' matrix is singular");
            return L;
        }
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix L, for row i */
    for (i = 1, pMi = L.m_elem + col; i < row; i++, pMi += col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    if (isOk)
        *isOk = 1;

    return L;
}

/** Solve the Lx = b problem for x vector with L having a unit diagonal
 *
 * Function Arguments:
 *
 * L    -> lower triangular Matrx in with unit diagonal
 * b    -> the right hand side 'b' vector
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Vector<0, t_type> LowerMatrix<0, 0, t_type>::SolveUnit(Matrix<row, col, t_type> &L, Vector<row, t_type> &b, int8_t *isOk)
{
    assert(L.m_elem != nullptr && "Memory has not been allocated to matrix 'L'");
    assert(b.m_elem != nullptr && "Memory has not been allocated to vector 'b'");

    Vector<0, t_type> x(row);
    if (row != col)
    {
        assert(false && "Check row and column dimensions of 'L' and 'b' do not match");
        if (isOk) *isOk = 0;
        return x;
    }

    int i, k;
    t_type *pL = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a unit lower triangular matrix.
    x.m_elem[0] = b.m_elem[0];
    for (k = 1, pL += col; k < row; pL += col, k++)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = 0, x.m_elem[k] = b.m_elem[k]; i < k; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pL + i);
    }

    return x;
}

/** Solve the Lx = b problem for x vector with L having a unit diagonal
 *
 * Function Arguments:
 *
 * L    -> lower triangular Matrx in with unit diagonal
 * b    -> the right hand side 'b' vector in dynamic memory Vector class
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Vector<0, t_type> LowerMatrix<0, 0, t_type>::SolveUnit(Matrix<row, col, t_type> &L, Vector<0, t_type> &b, int8_t *isOk)
{
    assert(L.m_elem != nullptr && "Memory has not been allocated to matrix 'L'");
    assert(col == b.m_row && "'L' Column dimensions and 'b' row dimensions do not matched");
    assert(b.m_elem != nullptr && "Memory has not been allocated to vector 'b'");

    Vector<0, t_type> x(row);
    if (row != col)
    {
        assert(false && "Row and column dimensions of 'L' do not match");
        if (isOk) *isOk = 0;
        return x;
    }

    int i, k;
    t_type *pL = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a unit lower triangular matrix.
    x.m_elem[0] = b.m_elem[0];
    for (k = 1, pL += col; k < row; pL += col, k++)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = 0, x.m_elem[k] = b.m_elem[k]; i < k; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pL + i);
    }
    return x;
}

/** Solve the Lx = b problem for x vector with L having a unit diagonal
 *
 * Function Arguments:
 *
 * L    -> lower triangular Matrx in dynamic Matrix class with unit diagonal
 * b    -> the right hand side 'b' vector in dynamic memory Vector class
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
inline Vector<0, t_type> LowerMatrix<0, 0, t_type>::SolveUnit(Matrix<0, 0, t_type> &L, Vector<0, t_type> &b, int8_t *isOk)
{
    assert(L.m_elem != nullptr && "Memory has not been allocated to matrix 'L'");
    assert(L.m_col == b.m_row && "'L' Column dimensions and 'b' row dimensions do not matched");
    assert(b.m_elem != nullptr && "Memory has not been allocated to vector 'b'");

    int row = L.m_row;
    int col = L.m_col;

    Vector<0, t_type> x(row);
    if (row != L.m_col)
    {
        assert(false && "Row and column dimensions of 'L' do not match");
        if (isOk) *isOk = 0;
        return x;
    }
    int i, k;
    t_type *pL = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a unit lower triangular matrix.
    x.m_elem[0] = b.m_elem[0];
    for (k = 1, pL += col; k < row; pL += col, k++)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = 0, x.m_elem[k] = b.m_elem[k]; i < k; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pL + i);
    }
    return x;
}

/** Solve the inverse of 'L' matrix having unit diagonal
 *
 * Function Arguments:
 *
 * L    -> lower triangular Matrx in dynamic Matrix class with unit diagonal
 * invL -> the matrix to store the inverse of L matrix in dynamic Matrix class
 *
 * \return    -> 0 on success
 */

template <typename t_type>
inline int8_t LowerMatrix<0, 0, t_type>::InverseUnit(Matrix<0, 0, t_type> L, Matrix<0, 0, t_type> &invL)
{
    int row = L.m_row;
    int col = L.m_col;

    assert(row == col && "Row and column dimensions of 'L' do not match");

    int i, j, k;
    t_type *pMi, *pMj, *pMk;

    /* Invert the subdiagonal part of the matrix L, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = L.m_elem + col; i < row; i++, pMi += col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + col; k < i; k++, pMk += col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    invL = L;
    return 0;
}

/** Solve the inverse of 'L' matrix with unit diagonal
 *
 * Function Arguments:
 *
 * L    -> lower triangular Matrx with unit diagonal in dynamic Matrix class
 * isOK -> the pointer flag for success
 *
 * \return    -> the inverse of L matrix in dynamic Matrix class
 */

template <typename t_type>
inline Matrix<0, 0, t_type> LowerMatrix<0, 0, t_type>::InverseUnit(Matrix<0, 0, t_type> L, int8_t *isOk)
{
    int row = L.m_row;
    int col = L.m_col;

    if (row != col)
    {
        if (isOk)
            *isOk = 0;
        assert(false && "Row and column dimensions of 'L' do not match");
        return L;
    }

    int i, j, k;
    t_type *pMi, *pMj, *pMk;

    /* Invert the subdiagonal part of the matrix L, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = L.m_elem + col; i < row; i++, pMi += col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + col; k < i; k++, pMk += col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    return L;
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTLOWER_MATRIX0_TPP_