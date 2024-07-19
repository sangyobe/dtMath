/*!
\file       dtUpperMatrix.tpp
\brief      dtMath, Upper triangular matrix solver class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTUPPER_MATRIX_TPP_
#define DTMATH_DTUPPER_MATRIX_TPP_

#include "dtUpperMatrix.h"

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t UpperMatrix<t_row, t_col, t_type>::Solve(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, Vector<t_col, t_type> &x)
{
    if (t_row != t_col)
        return -1;

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    for (k = t_col - 1, pU += t_row * (t_col - 1); k >= 0; pU -= t_col, k--)
    {
        if (std::abs(*(pU + k)) <= std::numeric_limits<t_type>::epsilon())
            return -1; // The matrix U is singular

        x.m_elem[k] = b.m_elem[k];

        for (i = k + 1; i < t_col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);

        x.m_elem[k] /= *(pU + k);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> UpperMatrix<t_row, t_col, t_type>::Solve(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, int8_t *isOk)
{
    if (t_row != t_col)
    {
        if (isOk)
            *isOk = 0;
        return Vector<t_col, t_type>();
    }

    int k, i;
    t_type x[t_col] = {
        0,
    };
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    for (k = t_col - 1, pU += t_row * (t_col - 1); k >= 0; pU -= t_col, k--)
    {
        if (std::abs(*(pU + k)) <= std::numeric_limits<t_type>::epsilon())
        {
            // The matrix U is singular
            if (isOk)
                *isOk = 0;
            return Vector<t_col, t_type>();
        }

        x[k] = b.m_elem[k];

        for (i = k + 1; i < t_col; i++)
            x[k] -= x[i] * *(pU + i);

        x[k] /= *(pU + k);
    }

    if (isOk)
        *isOk = 1;

    return Vector<t_col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t UpperMatrix<t_row, t_col, t_type>::Inverse(Matrix<t_row, t_col, t_type> U, Matrix<t_row, t_col, t_type> &invU)
{
    int i, j, k;
    t_type *pMi, *pMk;
    t_type sum;

    if (t_row != t_col)
        return -1;

    /* Invert the diagonal elements */
    for (k = 0, pMk = U.m_elem; k < t_row; pMk += (t_col + 1), k++)
    {
        if (std::abs(*pMk) <= std::numeric_limits<t_type>::epsilon())
            return -1;
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix U, for row i */
    for (i = t_row - 2, pMi = U.m_elem + t_col * (t_row - 2); i >= 0; pMi -= t_col, i--)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            for (k = i + 1, pMk = pMi + t_col; k <= j; pMk += t_col, k++)
            {
                sum += *(pMi + k) * *(pMk + j);
            }
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    invU = U;

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> UpperMatrix<t_row, t_col, t_type>::Inverse(Matrix<t_row, t_col, t_type> U, int8_t *isOk)
{
    int i, j, k;
    t_type *pMi, *pMk;
    t_type sum;

    if (t_row != t_col)
    {
        if (isOk)
            *isOk = 0;
        return Matrix<t_row, t_col, t_type>();
    }

    /* Invert the diagonal elements */
    for (k = 0, pMk = U.m_elem; k < t_row; pMk += (t_col + 1), k++)
    {
        if (std::abs(*pMk) <= std::numeric_limits<t_type>::epsilon())
        {
            if (isOk)
                *isOk = 0;
            return Matrix<t_row, t_col, t_type>();
        }
        else
            *pMk = 1 / *pMk;
    }

    /* Invert the remaining matrix L, for row i */
    for (i = t_row - 2, pMi = U.m_elem + t_col * (t_row - 2); i >= 0; pMi -= t_col, i--)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            for (k = i + 1, pMk = pMi + t_col; k <= j; pMk += t_col, k++)
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

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t UpperMatrix<t_row, t_col, t_type>::Solve(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == t_col && "Matrix row and column dimensions do not matched");
    assert(t_col == x.m_row && "Vector row and matrix column dimensions do not matched");

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    for (k = t_col - 1, pU += t_row * (t_col - 1); k >= 0; pU -= t_col, k--)
    {
        if (std::abs(*(pU + k)) <= std::numeric_limits<t_type>::epsilon())
        {
            assert(false && "The 'U' matrix is singular");
            return -1; // The matrix U is singular
        }

        x.m_elem[k] = b.m_elem[k];

        for (i = k + 1; i < t_col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);

        x.m_elem[k] /= *(pU + k);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t UpperMatrix<t_row, t_col, t_type>::Solve(Matrix<t_row, t_col, t_type> &U, Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == t_col && "Matrix row and column dimensions do not matched");
    assert(t_col == b.m_row && "'b' vector row and matrix column dimensions do not matched");
    assert(t_col == x.m_row && "x vector row and matrix column dimensions do not matched");

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    for (k = t_col - 1, pU += t_row * (t_col - 1); k >= 0; pU -= t_col, k--)
    {
        if (std::abs(*(pU + k)) <= std::numeric_limits<t_type>::epsilon())
        {
            assert(false && "The 'U' matrix is singular");
            return -1; // The matrix U is singular
        }

        x.m_elem[k] = b.m_elem[k];

        for (i = k + 1; i < t_col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);

        x.m_elem[k] /= *(pU + k);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t UpperMatrix<t_row, t_col, t_type>::Solve(Matrix<0, 0, t_type> &U, Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == U.m_col && "Template values must match matrix dimensions");
    assert(t_col == U.m_row && "Template values must match matrix dimensions");
    assert(U.m_row == U.m_col && "Matrix row and column dimensions do not matched");
    assert(U.m_col == b.m_row && "'b' vector row and matrix column dimensions do not matched");
    assert(b.m_row == x.m_row && "x vector row and matrix column dimensions do not matched");

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    for (k = t_col - 1, pU += t_row * (t_col - 1); k >= 0; pU -= t_col, k--)
    {
        if (std::abs(*(pU + k)) <= std::numeric_limits<t_type>::epsilon())
        {
            assert(false && "The 'U' matrix is singular");
            return -1; // The matrix U is singular
        }

        x.m_elem[k] = b.m_elem[k];

        for (i = k + 1; i < t_col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);

        x.m_elem[k] /= *(pU + k);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t UpperMatrix<t_row, t_col, t_type>::SolveUnit(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, Vector<t_col, t_type> &x)
{
    if (t_row != t_col)
        return -1;

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a unit upper triangular matrix.
    x.m_elem[t_col - 1] = b.m_elem[t_row - 1];
    for (k = t_col - 2, pU += t_row * (t_col - 2); k >= 0; pU -= t_col, k--)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = k + 1; i < t_col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> UpperMatrix<t_row, t_col, t_type>::SolveUnit(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, int8_t *isOk)
{
    if (t_row != t_col)
    {
        if (isOk)
            *isOk = 0;
        return Vector<t_col, t_type>();
    }

    int k, i;
    t_type x[t_col] = {
        0,
    };
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a unit upper triangular matrix.
    x[t_col - 1] = b.m_elem[t_row - 1];
    for (k = t_col - 2, pU += t_row * (t_col - 2); k >= 0; pU -= t_col, k--)
    {
        x[k] = b.m_elem[k];
        for (i = k + 1; i < t_col; i++)
            x[k] -= x[i] * *(pU + i);
    }

    if (isOk)
        *isOk = 1;

    return Vector<t_col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t UpperMatrix<t_row, t_col, t_type>::InverseUnit(Matrix<t_row, t_col, t_type> U, Matrix<t_row, t_col, t_type> &invU)
{
    if (t_row != t_col)
        return -1;

    int i, j, k;
    t_type *pMi, *pMk;

    /* Invert the subdiagonal part of the matrix U, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = t_row - 2, pMi = U.m_elem + t_col * (t_row - 2); i >= 0; pMi -= t_col, i--)
    {
        for (j = t_col - 1; j > i; j--)
        {
            *(pMi + j) = -*(pMi + j);
            for (k = i + 1, pMk = pMi + t_col; k < j; pMk += t_col, k++)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    invU = U;

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> UpperMatrix<t_row, t_col, t_type>::InverseUnit(Matrix<t_row, t_col, t_type> U, int8_t *isOk)
{
    if (t_row != t_col)
    {
        if (isOk)
            *isOk = 0;
        return Matrix<t_row, t_col, t_type>();
    }

    int i, j, k;
    t_type *pMi, *pMk;

    /* Invert the subdiagonal part of the matrix U, for row i */
    // where the diagonal elements are assumed to be 1
    for (i = t_row - 2, pMi = U.m_elem + t_col * (t_row - 2); i >= 0; pMi -= t_col, i--)
    {
        for (j = t_col - 1; j > i; j--)
        {
            *(pMi + j) = -*(pMi + j);
            for (k = i + 1, pMk = pMi + t_col; k < j; pMk += t_col, k++)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    if (isOk)
        *isOk = 1;

    return U;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t UpperMatrix<t_row, t_col, t_type>::SolveUnit(Matrix<t_row, t_col, t_type> &U, Vector<t_row, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == t_col && "Matrix row and column dimensions do not matched");
    assert(t_col == x.m_row && "Vector row and matrix column dimensions do not matched");

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a unit upper triangular matrix.
    x.m_elem[t_col - 1] = b.m_elem[t_row - 1];
    for (k = t_col - 2, pU += t_row * (t_col - 2); k >= 0; pU -= t_col, k--)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = k + 1; i < t_col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t UpperMatrix<t_row, t_col, t_type>::SolveUnit(Matrix<t_row, t_col, t_type> &U, Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == t_col && "Matrix row and column dimensions do not matched");
    assert(t_col == b.m_row && "'b' vector row and matrix column dimensions do not matched");
    assert(t_col == x.m_row && "x vector row and matrix column dimensions do not matched");

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    x.m_elem[t_col - 1] = b.m_elem[t_row - 1];
    for (k = t_col - 2, pU += t_row * (t_col - 2); k >= 0; pU -= t_col, k--)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = k + 1; i < t_col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t UpperMatrix<t_row, t_col, t_type>::SolveUnit(Matrix<0, 0, t_type> &U, Vector<0, t_type> &b, Vector<0, t_type> &x)
{
    assert(t_row == U.m_col && "Template values must match matrix dimensions");
    assert(t_col == U.m_row && "Template values must match matrix dimensions");
    assert(U.m_row == U.m_col && "Matrix row and column dimensions do not matched");
    assert(U.m_col == b.m_row && "'b' vector row and matrix column dimensions do not matched");
    assert(b.m_row == x.m_row && "x vector row and matrix column dimensions do not matched");

    int k, i;
    t_type *pU = U.m_elem;

    /* Solve Ux = b */
    // by backward-substitution, where U is a upper triangular matrix.
    x.m_elem[t_col - 1] = b.m_elem[t_row - 1];
    for (k = t_col - 2, pU += t_row * (t_col - 2); k >= 0; pU -= t_col, k--)
    {
        x.m_elem[k] = b.m_elem[k];
        for (i = k + 1; i < t_col; i++)
            x.m_elem[k] -= x.m_elem[i] * *(pU + i);
    }

    return 0;
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTUPPER_MATRIX_TPP_