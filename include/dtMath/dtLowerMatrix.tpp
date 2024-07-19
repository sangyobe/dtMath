/*!
\file       dtLowerMatrix.h
\brief      dtMath, Lower triangular matrix solver class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLOWER_MATRIX_TPP_
#define DTMATH_DTLOWER_MATRIX_TPP_

#include "dtLowerMatrix.h"

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LowerMatrix<t_row, t_col, t_type>::Solve(Matrix<t_row, t_col, t_type> &L, Vector<t_row, t_type> &b, Vector<t_col, t_type> &x)
{
    if (t_row != t_col)
        return -1;

    t_type *pMi = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a lower triangular matrix.
    for (int i = 0; i < t_row; i++, pMi += t_col)
    {
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
            return -1; // singular;

        x.m_elem[i] = b.m_elem[i];

        for (int k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> LowerMatrix<t_row, t_col, t_type>::Solve(Matrix<t_row, t_col, t_type> &L, Vector<t_row, t_type> &b, int8_t *isOk)
{
    if (t_row != t_col)
    {
        if (isOk)
            *isOk = 0;
        return Vector<t_col, t_type>();
    }

    t_type x[t_col] = {
        0,
    };
    t_type *pMi = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a lower triangular matrix.
    for (int i = 0; i < t_row; i++, pMi += t_col)
    {
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            if (isOk)
                *isOk = 0; // singular;
            return Vector<t_col, t_type>();
        }

        x[i] = b.m_elem[i];

        for (int k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];

        x[i] /= *(pMi + i);
    }

    if (isOk)
        *isOk = 1;

    return Vector<t_col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LowerMatrix<t_row, t_col, t_type>::Inverse(Matrix<t_row, t_col, t_type> L, Matrix<t_row, t_col, t_type> &invL)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type sum;

    if (t_row != t_col)
        return -1;

    /* Invert the diagonal elements */
    for (k = 0, pMk = L.m_elem; k < t_row; pMk += (t_col + 1), k++)
    {
        if (std::abs(*pMk) <= std::numeric_limits<t_type>::epsilon())
            return -1;
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix L, for row i */
    for (i = 1, pMi = L.m_elem + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += t_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += t_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    invL = L;

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> LowerMatrix<t_row, t_col, t_type>::Inverse(Matrix<t_row, t_col, t_type> L, int8_t *isOk)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type sum;

    if (t_row != t_col)
    {
        if (isOk)
            *isOk = 0;
        return Matrix<t_row, t_col, t_type>();
    }

    /* Invert the diagonal elements */
    for (k = 0, pMk = L.m_elem; k < t_row; pMk += (t_col + 1), k++)
    {
        if (*pMk == 0) // To do dhl
        {
            if (isOk)
                *isOk = 0;
            return Matrix<t_row, t_col, t_type>();
        }
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix L, for row i */
    for (i = 1, pMi = L.m_elem + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += t_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += t_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    if (isOk)
        *isOk = 1;

    return L;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LowerMatrix<t_row, t_col, t_type>::Solve(Matrix<t_row, t_col, t_type> &L, Vector<t_row, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == t_col && "Matrix row and column dimensions do not matched");
    assert(t_col == x.m_row && "Vector row and matrix column dimensions do not matched");

    t_type *pMi = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a lower triangular matrix.
    for (int i = 0; i < t_row; i++, pMi += t_col)
    {
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            assert(false && "The 'L' matrix is singular");
            return -1;
        }

        x.m_elem[i] = b.m_elem[i];

        for (int k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LowerMatrix<t_row, t_col, t_type>::Solve(Matrix<t_row, t_col, t_type> &L, Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == t_col && "Matrix row and column dimensions do not matched");
    assert(t_col == b.m_row && "'b' vector row and matrix column dimensions do not matched");
    assert(t_col == x.m_row && "x vector row and matrix column dimensions do not matched");

    t_type *pMi = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a lower triangular matrix.
    for (int i = 0; i < t_row; i++, pMi += t_col)
    {
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            assert(false && "The 'L' matrix is singular");
            return -1;
        }

        x.m_elem[i] = b.m_elem[i];

        for (int k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LowerMatrix<t_row, t_col, t_type>::Solve(Matrix<0, 0, t_type> &L, Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == L.m_col && "Template values must match matrix dimensions");
    assert(t_col == L.m_row && "Template values must match matrix dimensions");
    assert(L.m_row == L.m_col && "Matrix row and column dimensions do not matched");
    assert(L.m_col == b.m_row && "'b' vector row and matrix column dimensions do not matched");
    assert(b.m_row == x.m_row && "x vector row and matrix column dimensions do not matched");

    t_type *pMi = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a lower triangular matrix.
    for (int i = 0; i < t_row; i++, pMi += t_col)
    {
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            assert(false && "The 'L' matrix is singular");
            return -1;
        }

        x.m_elem[i] = b.m_elem[i];

        for (int k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LowerMatrix<t_row, t_col, t_type>::SolveUnit(Matrix<t_row, t_col, t_type> &L, Vector<t_row, t_type> &b, Vector<t_col, t_type> &x)
{
    if (t_row != t_col)
        return -1;

    int i, k;
    t_type *pL = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a unit lower triangular matrix.
    x.m_elem[0] = b.m_elem[0];
    for (k = 1, pL += t_col; k < t_row; pL += t_col, k++)
    {
        x.m_elem[i] = b.m_elem[i];
        for (i = 0, x.m_elem[k] = b.m_elem[k]; i < k; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pL + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> LowerMatrix<t_row, t_col, t_type>::SolveUnit(Matrix<t_row, t_col, t_type> &L, Vector<t_row, t_type> &b, int8_t *isOk)
{
    if (t_row != t_col)
    {
        if (isOk)
            *isOk = 0;
        return Vector<t_col, t_type>();
    }

    int i, k;
    t_type x[t_col] = {
        0,
    };
    t_type *pL = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a unit lower triangular matrix.
    x[0] = b.m_elem[0];
    for (k = 1, pL += t_col; k < t_row; pL += t_col, k++)
    {
        x[i] = b.m_elem[i];
        for (i = 0, x[k] = b.m_elem[k]; i < k; i++)
            x[k] -= x[i] * *(pL + i);
    }

    return Vector<t_col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LowerMatrix<t_row, t_col, t_type>::InverseUnit(Matrix<t_row, t_col, t_type> L, Matrix<t_row, t_col, t_type> &invL)
{
    if (t_row != t_col)
        return -1;

    int i, j, k;
    t_type *pMi, *pMj, *pMk;

    /* Invert the subdiagonal part of the matrix L, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = L.m_elem + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += t_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + t_col; k < i; k++, pMk += t_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    invL = L;

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> LowerMatrix<t_row, t_col, t_type>::InverseUnit(Matrix<t_row, t_col, t_type> L, int8_t *isOk)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;

    /* Invert the subdiagonal part of the matrix L, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = L.m_elem + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = L.m_elem; j < i; pMj += t_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + t_col; k < i; k++, pMk += t_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    if (isOk)
        *isOk = 1;

    return L;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
int8_t LowerMatrix<t_row, t_col, t_type>::SolveUnit(Matrix<t_row, t_col, t_type> &L, Vector<t_row, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == t_col && "'L' Matrix's row and column dimensions do not matched");
    assert(t_col == x.m_row && "'x' vector's row and 'L' matrix's column dimensions do not matched");

    int i, k;
    t_type *pL = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a unit lower triangular matrix.
    x.m_elem[0] = b.m_elem[0];
    for (k = 1, pL += t_col; k < t_row; pL += t_col, k++)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = 0, x.m_elem[k] = b.m_elem[k]; i < k; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pL + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
int8_t LowerMatrix<t_row, t_col, t_type>::SolveUnit(Matrix<t_row, t_col, t_type> &L, Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == t_col && "Matrix row and column dimensions do not matched");
    assert(t_col == b.m_row && "'b' vector row and matrix column dimensions do not matched");
    assert(t_col == x.m_row && "x vector row and matrix column dimensions do not matched");

    int i, k;
    t_type *pL = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a unit lower triangular matrix.
    x.m_elem[0] = b.m_elem[0];
    for (k = 1, pL += t_col; k < t_row; pL += t_col, k++)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = 0, x.m_elem[k] = b.m_elem[k]; i < k; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pL + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
int8_t LowerMatrix<t_row, t_col, t_type>::SolveUnit(Matrix<0, 0, t_type> &L, Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == L.m_col && "Template values must match matrix dimensions");
    assert(t_col == L.m_row && "Template values must match matrix dimensions");
    assert(L.m_row == L.m_col && "Matrix row and column dimensions do not matched");
    assert(L.m_col == b.m_row && "'b' vector row and matrix column dimensions do not matched");
    assert(b.m_row == x.m_row && "x vector row and matrix column dimensions do not matched");

    int i, k;
    t_type *pL = L.m_elem;

    /* Solve Lx = b */
    // by forward-substitution, where L is a unit lower triangular matrix.
    x.m_elem[0] = b.m_elem[0];
    for (k = 1, pL += t_col; k < t_row; pL += t_col, k++)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = 0, x.m_elem[k] = b.m_elem[k]; i < k; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pL + i);
    }

    return 0;
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTLOWER_MATRIX_TPP_
