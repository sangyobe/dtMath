/*!
\file       dtLLT.h
\brief      dtMath, Cholesky decomposition(L*L^T form) Class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLLT_TPP_
#define DTMATH_DTLLT_TPP_

#include "dtLLT.h"

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline LLT<t_row, t_col, t_type>::LLT()
{
    memset(m_elem, 0, sizeof(t_type) * t_row * t_col);
    m_isOk = 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline LLT<t_row, t_col, t_type>::LLT(const t_type *element, const size_t n_byte)
{
    if ((sizeof(t_type) * t_row * t_col) != n_byte)
        m_isOk = 0;
    else
    {
        memset(m_elem, element, n_byte);
        Compute();
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline LLT<t_row, t_col, t_type>::LLT(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline LLT<t_row, t_col, t_type>::LLT(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LLT<t_row, t_col, t_type>::Compute()
{
    if (t_row != t_col)
    {
        m_isOk = 0;
        return -1;
    }

    int i, j, k;
    t_type *pMi, *pMj, *pMjj, *pMjk;

    for (j = 0, pMj = m_elem; j < t_col; pMj += t_row, j++)
    {
        /* Calculate the Diagonal element in colum j */
        pMjj = pMj + j;
        for (k = 0, pMjk = pMj; k < j; pMjk += 1, k++)
            *pMjj -= *pMjk * *pMjk;

        // If diagonal element is not positive, return the error,
        // the matrix is not positive definite symmetric.
        if (*pMjj <= std::numeric_limits<t_type>::epsilon())
        {
            m_isOk = 0;
            return -1;
        }

        *pMjj = std::sqrt(*pMjj);

        /* Calculate the lower triangular matrix for colum j */
        for (i = j + 1, pMi = pMj + t_col; i < t_row; pMi += t_col, i++)
        {
            for (k = 0; k < j; k++)
                *(pMi + j) -= *(pMi + k) * *(pMj + k);

            *(pMi + j) /= *pMjj;     // Lower Triangular Matrix, in-place
            *(pMj + i) = *(pMi + j); // Upper Triangular Matrix, in-place
        }
    }

    m_isOk = 1;
    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LLT<t_row, t_col, t_type>::Compute(const t_type *element, const size_t n_byte)
{
    if ((sizeof(t_type) * t_row * t_col) != n_byte)
    {
        m_isOk = 0;
        return -1;
    }

    memset(m_elem, element, n_byte);
    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LLT<t_row, t_col, t_type>::Compute(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LLT<t_row, t_col, t_type>::Compute(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> LLT<t_row, t_col, t_type>::GetMatrix() const
{
    return Matrix<t_row, t_col, t_type>(m_elem);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> LLT<t_row, t_col, t_type>::GetMatrixL() const
{
    int i, j;
    t_type L[t_row * t_col]{0};

    for (i = 0; i < t_row; i++)
        for (j = 0; j <= i; j++)
            L[i * t_col + j] = m_elem[i * t_col + j];

    return Matrix<t_row, t_col, t_type>(L);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> LLT<t_row, t_col, t_type>::GetMatrixU() const
{
    int i, j;
    t_type U[t_row * t_col]{0};

    for (i = 0; i < t_row; i++)
        for (j = 0; j < i; j++)
            U[i + t_col * j] = m_elem[i * t_col + j];

    return Matrix<t_row, t_col, t_type>(U);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t col>
inline int8_t LLT<t_row, t_col, t_type>::Solve(const Matrix<t_row, col, t_type> &b, Matrix<t_col, col, t_type> &x)
{
    // Solve, Ax = LLTx = LUx = b
    // where L is a lower triangular matrix
    //       LT is a transposed lower triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k, j;
    t_type *pMi;

    if (!m_isOk)
        return -1;

    for (j = 0; j < col; j++)
    {
        /* Solve Ly = b */
        // Solve the lower triangular matrix for y (forward substitution), here x is y
        for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
        {
            x.m_elem[i * col + j] = b.m_elem[i * col + j];

            for (k = 0; k < i; k++)
                x.m_elem[i * col + j] -= *(pMi + k) * x.m_elem[k * col + j];

            x.m_elem[i * col + j] /= *(pMi + i);
        }

        /* Solve LTx = y */
        // Solve the upper triangular (backward substitution), LT = U
        for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
        {
            for (k = i + 1; k < t_col; k++)
                x.m_elem[i * col + j] -= *(pMi + k) * x.m_elem[k * col + j];

            x.m_elem[i * col + j] /= *(pMi + i);
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LLT<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, Vector<t_col, t_type> &x)
{
    // Solve, Ax = LLTx = LUx = b
    // where L is a lower triangular matrix
    //       LT is a transposed lower triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    t_type *pMi;

    if (!m_isOk)
        return -1;

    /* Solve Ly = b */
    // Solve the lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    /* Solve LTx = y */
    // Solve the upper triangular (backward substitution), LT = U
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        for (k = i + 1; k < t_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t col>
inline Matrix<t_col, col, t_type> LLT<t_row, t_col, t_type>::Solve(const Matrix<t_row, col, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LLTx = LUx = b
    // where L is a lower triangular matrix
    //       LT is a transposed lower triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k, j;
    t_type *pMi;
    t_type x[t_col * col]{0};

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return Matrix<t_col, col, t_type>(x);
    }

    for (j = 0; j < col; j++)
    {
        /* Solve Ly = b */
        // Solve the lower triangular matrix for y (forward substitution), here x is y
        for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
        {
            x[i * col + j] = b.m_elem[i * col + j];

            for (k = 0; k < i; k++)
                x[i * col + j] -= *(pMi + k) * x[k * col + j];

            x[i * col + j] /= *(pMi + i);
        }

        /* Solve LTx = y */
        // Solve the upper triangular (backward substitution), LT = U
        for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
        {
            for (k = i + 1; k < t_col; k++)
                x[i * col + j] -= *(pMi + k) * x[k * col + j];

            x[i * col + j] /= *(pMi + i);
        }
    }

    return Matrix<t_col, col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> LLT<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LLTx = LUx = b
    // where L is a lower triangular matrix
    //       LT is a transposed lower triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    t_type *pMi;
    t_type x[t_col]{0};

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return Vector<t_col, t_type>(x);
    }

    /* Solve Ly = b */
    // Solve the lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
    {
        x[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];

        x[i] /= *(pMi + i);
    }

    /* Solve LTx = y */
    // Solve the upper triangular (backward substitution), LT = U
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        for (k = i + 1; k < t_col; k++)
            x[i] -= *(pMi + k) * x[k];

        x[i] /= *(pMi + i);
    }

    return Vector<t_col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LLT<t_row, t_col, t_type>::Inverse(Matrix<t_row, t_col, t_type> &inv)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type sum;

    if (!m_isOk)
        return -1;

    memcpy(inv.m_elem, m_elem, sizeof(t_type) * t_row * t_col);

    /* Calculate the inverse of L */
    // Invert the diagonal elements of the lower triangular matrix L.
    for (k = 0, pMk = inv.m_elem; k < t_row; pMk += (t_col + 1), k++)
    {
        *pMk = 1 / *pMk;
    }

    // Invert the remaining lower triangular matrix L row by row.
    for (i = 1, pMi = inv.m_elem + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv.m_elem; j < i; pMj += t_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += t_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    /* Calculate the inverse of LLT, inv(LT) * inv(L) */
    // inv(LLT) is also positive definite symmetric matrix
    for (i = 0, pMi = inv.m_elem; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv.m_elem; j <= i; j++, pMj += t_row)
        {
            sum = 0;
            for (k = i, pMk = pMi; k < t_col; k++, pMk += t_col)
                sum += *(pMk + i) * *(pMk + j);

            *(pMi + j) = sum; // upper parts
            *(pMj + i) = sum; // lower parts
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LLT<t_row, t_col, t_type>::Inverse(Matrix3<t_type, t_row, t_col> &inv)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type sum;

    if (!m_isOk)
        return -1;

    memcpy(inv.m_elem, m_elem, sizeof(t_type) * t_row * t_col);

    /* Calculate the inverse of L */
    // Invert the diagonal elements of the lower triangular matrix L.
    for (k = 0, pMk = inv.m_elem; k < t_row; pMk += (t_col + 1), k++)
    {
        *pMk = 1 / *pMk;
    }

    // Invert the remaining lower triangular matrix L row by row.
    for (i = 1, pMi = inv.m_elem + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv.m_elem; j < i; pMj += t_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += t_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    /* Calculate the inverse of LLT, inv(LT) * inv(L) */
    // inv(LLT) is also positive definite symmetric matrix
    for (i = 0, pMi = inv.m_elem; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv.m_elem; j <= i; j++, pMj += t_row)
        {
            sum = 0;
            for (k = i, pMk = pMi; k < t_col; k++, pMk += t_col)
                sum += *(pMk + i) * *(pMk + j);

            *(pMi + j) = sum; // upper parts
            *(pMj + i) = sum; // lower parts
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> LLT<t_row, t_col, t_type>::Inverse(int8_t *isOk)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type sum;
    t_type inv[t_row * t_col]{0};

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return Matrix<t_row, t_col, t_type>();
    }

    memcpy(inv, m_elem, sizeof(t_type) * t_row * t_col);

    /* Calculate the inverse of L */
    // Invert the diagonal elements of the lower triangular matrix L.
    for (k = 0, pMk = inv; k < t_row; pMk += (t_col + 1), k++)
    {
        *pMk = 1 / *pMk;
    }

    // Invert the remaining lower triangular matrix L row by row.
    for (i = 1, pMi = inv + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += t_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += t_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    /* Calculate the inverse of LLT, inv(LT) * inv(L) */
    // inv(LLT) is also positive definite symmetric matrix
    for (i = 0, pMi = inv; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv; j <= i; j++, pMj += t_row)
        {
            sum = 0;
            for (k = i, pMk = pMi; k < t_col; k++, pMk += t_col)
                sum += *(pMk + i) * *(pMk + j);

            *(pMi + j) = sum; // upper parts
            *(pMj + i) = sum; // lower parts
        }
    }

    return Matrix<t_row, t_col, t_type>(inv);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LLT<t_row, t_col, t_type>::InverseArray(t_type *inv)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type sum;

    if (!m_isOk)
        return -1;

    memcpy(inv, m_elem, sizeof(t_type) * t_row * t_col);

    /* Calculate the inverse of L */
    // Invert the diagonal elements of the lower triangular matrix L.
    for (k = 0, pMk = inv; k < t_row; pMk += (t_col + 1), k++)
    {
        *pMk = 1 / *pMk;
    }

    // Invert the remaining lower triangular matrix L row by row.
    for (i = 1, pMi = inv + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += t_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += t_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    /* Calculate the inverse of LLT, inv(LT) * inv(L) */
    // inv(LLT) is also positive definite symmetric matrix
    for (i = 0, pMi = inv; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv; j <= i; j++, pMj += t_row)
        {
            sum = 0;
            for (k = i, pMk = pMi; k < t_col; k++, pMk += t_col)
                sum += *(pMk + i) * *(pMk + j);

            *(pMi + j) = sum; // upper parts
            *(pMj + i) = sum; // lower parts
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type *LLT<t_row, t_col, t_type>::InverseArray(int8_t *isOk)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type sum;
    t_type inv[t_row * t_col]{0};

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return inv;
    }

    memcpy(inv, m_elem, sizeof(t_type) * t_row * t_col);

    /* Calculate the inverse of L */
    // Invert the diagonal elements of the lower triangular matrix L.
    for (k = 0, pMk = inv; k < t_row; pMk += (t_col + 1), k++)
    {
        *pMk = 1 / *pMk;
    }

    // Invert the remaining lower triangular matrix L row by row.
    for (i = 1, pMi = inv + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += t_col, j++)
        {
            sum = 0;
            for (k = j, pMk = pMj; k < i; k++, pMk += t_col)
                sum += *(pMi + k) * *(pMk + j);
            *(pMi + j) = -*(pMi + i) * sum;
        }
    }

    /* Calculate the inverse of LLT, inv(LT) * inv(L) */
    // inv(LLT) is also positive definite symmetric matrix
    for (i = 0, pMi = inv; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv; j <= i; j++, pMj += t_row)
        {
            sum = 0;
            for (k = i, pMk = pMi; k < t_col; k++, pMk += t_col)
                sum += *(pMk + i) * *(pMk + j);

            *(pMi + j) = sum; // upper parts
            *(pMj + i) = sum; // lower parts
        }
    }

    return inv;
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTLLT_TPP_
