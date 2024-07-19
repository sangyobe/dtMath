/*!
\file       dtUpperMatrix0.tpp
\brief      dtMath, Upper triangular matrix solver class with dynamic memory
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 05. 21
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTUPPER_MATRIX0_TPP_
#define DTMATH_DTUPPER_MATRIX0_TPP_

#include "dtUpperMatrix0.h"

namespace dt
{
namespace Math
{

/** Solve the Ux = b problem for x vector
 *
 * Function Arguments:
 *
 * U    -> upper triangular Matrx
 * b    -> the right hand side 'b' vector
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Vector<0, t_type> UpperMatrix<0, 0, t_type>::Solve(Matrix<row, col, t_type> &U, Vector<row, t_type> &b, int8_t *isOk)
{

    Vector<0, t_type> x(row);
    if (row != col)
    {
        assert(false && "Check row and column dimensions of 'U' and 'b' do not match");
        if (isOk) *isOk = 0;
        return x;
    }

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    for (k = col - 1, pU += row * (col - 1); k >= 0; pU -= col, k--)
    {
        if (std::abs(*(pU + k)) <= std::numeric_limits<t_type>::epsilon())
        {
            assert(false && "The 'U' matrix is singular");
            return -1; // The matrix U is singular
        }

        x.m_elem[k] = b.m_elem[k];

        for (i = k + 1; i < col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);

        x.m_elem[k] /= *(pU + k);
    }

    if (isOk)
        *isOk = 1;

    return x;
}

/** Solve the Ux = b problem for x vector
 *
 * Function Arguments:
 *
 * U    -> upper triangular Matrx
 * b    -> the right hand side 'b' vector in dynamic memory Vector class
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Vector<0, t_type> UpperMatrix<0, 0, t_type>::Solve(Matrix<row, col, t_type> &U, Vector<0, t_type> &b, int8_t *isOk)
{
    assert(col == b.m_row && "'L' Column dimensions and 'b' row dimensions do not matched");
    assert(b.m_elem != nullptr && "Memory has not been allocated to vector 'b'");

    Vector<0, t_type> x(row);
    if (row != col)
    {
        assert(false && "Row and column dimensions of 'U' do not match");
        if (isOk) *isOk = 0;
        return x;
    }

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    for (k = col - 1, pU += row * (col - 1); k >= 0; pU -= col, k--)
    {
        if (std::abs(*(pU + k)) <= std::numeric_limits<t_type>::epsilon())
        {
            assert(false && "The 'U' matrix is singular");
            return -1; // The matrix U is singular
        }

        x.m_elem[k] = b.m_elem[k];

        for (i = k + 1; i < col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);

        x.m_elem[k] /= *(pU + k);
    }

    if (isOk)
        *isOk = 1;

    return x;
}

/** Solve the Ux = b problem for x vector
 *
 * Function Arguments:
 *
 * U    -> upper triangular Matrx in dynamic Matrix class
 * b    -> the right hand side 'b' vector in dynamic memory Vector class
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
inline Vector<0, t_type> UpperMatrix<0, 0, t_type>::Solve(Matrix<0, 0, t_type> &U, Vector<0, t_type> &b, int8_t *isOk)
{
    assert(U.m_elem != nullptr && "Memory has not been allocated to matrix 'U'");
    assert(U.m_col == b.m_row && "'U' Column dimensions and 'b' row dimensions do not matched");
    assert(b.m_elem != nullptr && "Memory has not been allocated to vector 'b'");

    int row = U.m_row;
    int col = U.m_col;

    Vector<0, t_type> x(row);
    if (row != U.m_col)
    {
        assert(false && "Row and column dimensions of 'U' do not match");
        if (isOk) *isOk = 0;
        return x;
    }

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    for (k = col - 1, pU += row * (col - 1); k >= 0; pU -= col, k--)
    {
        if (std::abs(*(pU + k)) <= std::numeric_limits<t_type>::epsilon())
        {
            assert(false && "The 'U' matrix is singular");
            return -1; // The matrix U is singular
        }

        x.m_elem[k] = b.m_elem[k];

        for (i = k + 1; i < col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);

        x.m_elem[k] /= *(pU + k);
    }

    if (isOk)
        *isOk = 1;
    return x;
}

/** Solve the inverse of 'U' matrix
 *
 * Function Arguments:
 *
 * U    -> upper triangular Matrx in dynamic Matrix class
 * invU -> the matrix to store the inverse of U matrix in dynamic Matrix class
 *
 * \return    -> 0 on success; -1 on failure
 */

template <typename t_type>
int8_t UpperMatrix<0, 0, t_type>::Inverse(Matrix<0, 0, t_type> U, Matrix<0, 0, t_type> &invU)
{
    int i, j, k;
    t_type *pMi, *pMk;
    t_type sum;

    int row = U.m_row;
    int col = U.m_col;

    assert(row == col && "Row and column dimensions of 'U' do not match");

    /* Invert the diagonal elements */
    for (k = 0, pMk = U.m_elem; k < row; pMk += (col + 1), k++)
    {
        if (std::abs(*pMk) <= std::numeric_limits<t_type>::epsilon())
        {
            assert(false && "The 'U' matrix is singular");
            return -1;
        }
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix U, for row i */
    for (i = row - 2, pMi = U.m_elem + col * (row - 2); i >= 0; pMi -= col, i--)
    {
        for (j = col - 1; j > i; j--)
        {
            sum = 0;
            for (k = i + 1, pMk = pMi + col; k <= j; pMk += col, k++)
            {
                sum += *(pMi + k) * *(pMk + j);
            }
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    invU = U;

    return 0;
}

/** Solve the inverse of 'U' matrix
 *
 * Function Arguments:
 *
 * U    -> upper triangular Matrx in dynamic Matrix class
 * isOK -> the pointer flag for success
 *
 * \return    -> the inverse of U matrix in dynamic Matrix class
 */

template <typename t_type>
inline Matrix<0, 0, t_type> UpperMatrix<0, 0, t_type>::Inverse(Matrix<0, 0, t_type> U, int8_t *isOk)
{
    int i, j, k;
    t_type *pMi, *pMk;
    t_type sum;

    int row = U.m_row;
    int col = U.m_col;

    if (row != col)
    {
        if (isOk)
            *isOk = 0;
        assert(false && "Row and column dimensions of 'L' do not match");
        return U;
    }

    /* Invert the diagonal elements */
    for (k = 0, pMk = U.m_elem; k < row; pMk += (col + 1), k++)
    {
        if (std::abs(*pMk) <= std::numeric_limits<t_type>::epsilon())
        {
            if (isOk)
                *isOk = 0;
            return U;
        }
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix U, for row i */
    for (i = row - 2, pMi = U.m_elem + col * (row - 2); i >= 0; pMi -= col, i--)
    {
        for (j = col - 1; j > i; j--)
        {
            sum = 0;
            for (k = i + 1, pMk = pMi + col; k <= j; pMk += col, k++)
            {
                sum += *(pMi + k) * *(pMk + j);
            }
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    if (isOk)
        *isOk = 1;

    return U;
}

/** Solve the Ux = b problem for x vector with U having a unit diagonal
 *
 * Function Arguments:
 *
 * U    -> upper triangular Matrx in with unit diagonal
 * b    -> the right hand side 'b' vector
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Vector<0, t_type> UpperMatrix<0, 0, t_type>::SolveUnit(Matrix<row, col, t_type> &U, Vector<row, t_type> &b, int8_t *isOk)
{

    Vector<0, t_type> x(row);
    if (row != col)
    {
        assert(false && "Check row and column dimensions of 'U' and 'b' do not match");
        if (isOk) *isOk = 0;
        return x;
    }

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a unit upper triangular matrix.
    x.m_elem[col - 1] = b.m_elem[row - 1];
    for (k = col - 2, pU += row * (col - 2); k >= 0; pU -= col, k--)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = k + 1; i < col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);
    }

    return x;
}

/** Solve the Ux = b problem for x vector with U having a unit diagonal
 *
 * Function Arguments:
 *
 * U    -> upper triangular Matrx in with unit diagonal
 * b    -> the right hand side 'b' vector in dynamic memory Vector class
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Vector<0, t_type> UpperMatrix<0, 0, t_type>::SolveUnit(Matrix<row, col, t_type> &U, Vector<0, t_type> &b, int8_t *isOk)
{
    assert(col == b.m_row && "'L' Column dimensions and 'b' row dimensions do not matched");
    assert(b.m_elem != nullptr && "Memory has not been allocated to vector 'b'");

    Vector<0, t_type> x(row);
    if (row != col)
    {
        assert(false && "Row and column dimensions of 'U' do not match");
        if (isOk) *isOk = 0;
        return x;
    }

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a unit upper triangular matrix.
    x.m_elem[col - 1] = b.m_elem[row - 1];
    for (k = col - 2, pU += row * (col - 2); k >= 0; pU -= col, k--)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = k + 1; i < col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);
    }
    return x;
}

/** Solve the Ux = b problem for x vector with U having a unit diagonal
 *
 * Function Arguments:
 *
 * U    -> upper triangular Matrx in dynamic Matrix class with unit diagonal
 * b    -> the right hand side 'b' vector in dynamic memory Vector class
 * isOK -> the pointer flag for success
 *
 * \return    -> the output x vector in dynamic memory Vector class
 */

template <typename t_type>
inline Vector<0, t_type> UpperMatrix<0, 0, t_type>::SolveUnit(Matrix<0, 0, t_type> &U, Vector<0, t_type> &b, int8_t *isOk)
{
    assert(U.m_elem != nullptr && "Memory has not been allocated to matrix 'L'");
    assert(U.m_col == b.m_row && "'L' Column dimensions and 'b' row dimensions do not matched");
    assert(b.m_elem != nullptr && "Memory has not been allocated to vector 'b'");

    int row = U.m_row;
    int col = U.m_col;

    Vector<0, t_type> x(row);
    if (row != U.m_col)
    {
        assert(false && "Row and column dimensions of 'U' do not match");
        if (isOk) *isOk = 0;
        return x;
    }

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a unit upper triangular matrix.
    x.m_elem[col - 1] = b.m_elem[row - 1];
    for (k = col - 2, pU += row * (col - 2); k >= 0; pU -= col, k--)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = k + 1; i < col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);
    }
    return x;
}

/** Solve the inverse of 'U' matrix having unit diagonal
 *
 * Function Arguments:
 *
 * U    -> upper triangular Matrx in dynamic Matrix class with unit diagonal
 * invU -> the matrix to store the inverse of U matrix in dynamic Matrix class
 *
 * \return    -> 0 on success
 */

template <typename t_type>
inline int8_t UpperMatrix<0, 0, t_type>::InverseUnit(Matrix<0, 0, t_type> U, Matrix<0, 0, t_type> &invU)
{
    int row = U.m_row;
    int col = U.m_col;

    assert(row == col && "Row and column dimensions of 'U' do not match");

    int i, j, k;
    t_type *pMi, *pMk;

    /* Invert the subdiagonal part of the matrix U, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = row - 2, pMi = U.m_elem + col * (row - 2); i >= 0; pMi -= col, i--)
    {
        for (j = col - 1; j > i; j--)
        {
            *(pMi + j) = -*(pMi + j);
            for (k = i + 1, pMk = pMi + col; k < j; pMk += col, k++)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    invU = U;

    return 0;
}

/** Solve the inverse of 'U' matrix with unit diagonal
 *
 * Function Arguments:
 *
 * U    -> upper triangular Matrx with unit diagonal in dynamic Matrix class
 * isOK -> the pointer flag for success
 *
 * \return    -> the inverse of U matrix in dynamic Matrix class
 */

template <typename t_type>
inline Matrix<0, 0, t_type> UpperMatrix<0, 0, t_type>::InverseUnit(Matrix<0, 0, t_type> U, int8_t *isOk)
{
    int row = U.m_row;
    int col = U.m_col;

    if (row != col)
    {
        if (isOk)
            *isOk = 0;
        assert(false && "Row and column dimensions of 'U' do not match");
        return U;
    }

    int i, j, k;
    t_type *pMi, *pMk;

    /* Invert the subdiagonal part of the matrix U, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = row - 2, pMi = U.m_elem + col * (row - 2); i >= 0; pMi -= col, i--)
    {
        for (j = col - 1; j > i; j--)
        {
            *(pMi + j) = -*(pMi + j);
            for (k = i + 1, pMk = pMi + col; k < j; pMk += col, k++)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    return U;
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTUPPER_MATRIX0_TPP_