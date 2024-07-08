/*!
\file       dtLDLT.h
\brief      dtMath, Cholesky decomposition(L*D*L^T form) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTLDLT_TPP_
#define DTMATH_DTLDLT_TPP_

#include "dtLDLT.h"

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline LDLT<t_row, t_col, t_type>::LDLT()
{
    memset(m_elem, 0, sizeof(t_type) * t_row * t_col);
    m_isOk = 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline LDLT<t_row, t_col, t_type>::LDLT(const t_type *element, const size_t n_byte)
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
inline LDLT<t_row, t_col, t_type>::LDLT(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline LDLT<t_row, t_col, t_type>::LDLT(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LDLT<t_row, t_col, t_type>::Compute()
{
    if (t_row != t_col)
    {
        m_isOk = 0;
        return -1;
    }

    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type Lik;

    for (i = 1, pMi = m_elem + t_col; i < t_row; pMi += t_col, i++)
    {
        // Calculate the part of L*D Matrix from Lij*Djj
        // L*D matrix are used for Calculating Dii
        // here *(pMi + j) = Lij*Djj
        for (j = 0, pMj = m_elem; j < i; j++, pMj += t_row)
            for (k = 0; k < j; k++)
                *(pMi + j) -= *(pMi + k) * *(pMj + k);

        // Calculate the diagonal element Dii
        // Calculate the Mij from Lik and Store the Mjk
        for (k = 0, pMk = m_elem; k < i; pMk += t_col, k++)
        {
            Lik = *(pMi + k) / *(pMk + k);  // *(pMi + k) = Lik * Dkk
            *(pMi + i) -= *(pMi + k) * Lik; // Dii = Aii - sum(Lik * Dkk * Lik)
            *(pMi + k) = Lik;               // Lik
            *(pMk + i) = Lik;               // Lki
        }

        // If diagonal element is not positive, return the error,
        // the matrix is not positive definite symmetric.
        if (*(pMi + i) <= std::numeric_limits<t_type>::epsilon())
        {
            m_isOk = 0;
            return -1;
        }
    }

    m_isOk = 1;
    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LDLT<t_row, t_col, t_type>::Compute(const t_type *element, const size_t n_byte)
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
inline int8_t LDLT<t_row, t_col, t_type>::Compute(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LDLT<t_row, t_col, t_type>::Compute(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> LDLT<t_row, t_col, t_type>::GetMatrix() const
{
    return Matrix<t_row, t_col, t_type>(m_elem);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> LDLT<t_row, t_col, t_type>::GetMatrixL() const
{
    int i, j;
    t_type L[t_row * t_col] = {
        0,
    };

    for (i = 0; i < t_row; i++)
        L[i * t_col + i] = 1;

    for (i = 1; i < t_row; i++)
        for (j = 0; j < i; j++)
            L[i * t_col + j] = m_elem[i * t_col + j];

    return Matrix<t_row, t_col, t_type>(L);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> LDLT<t_row, t_col, t_type>::GetMatrixD() const
{
    int i;
    t_type D[t_row * t_col] = {
        0,
    };

    for (i = 0; i < t_row; i++)
        D[i * t_col + i] = m_elem[i * t_col + i];

    return Matrix<t_row, t_col, t_type>(D);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> LDLT<t_row, t_col, t_type>::GetMatrixU() const
{
    int i, j;
    t_type U[t_row * t_col] = {
        0,
    };

    for (i = 0; i < t_row; i++)
        U[i * t_col + i] = 1;

    for (i = 1; i < t_row; i++)
        for (j = 0; j < i; j++)
            U[i + t_col * j] = m_elem[i * t_col + j];

    return Matrix<t_row, t_col, t_type>(U);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t col>
inline int8_t LDLT<t_row, t_col, t_type>::Solve(const Matrix<t_row, col, t_type> &b, Matrix<t_col, col, t_type> &x)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
    // Ux = y is solved by backward substitution for x

    int i, k, j;
    t_type *pMi;

    if (!m_isOk)
        return -1;

    for (j = 0; j < col; j++)
    {
        /* Solve Lz = b */
        // Solve the unit lower triangular (forward substitution), here x is z
        for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
        {
            x.m_elem[i * col + j] = b.m_elem[i * col + j];

            for (k = 0; k < i; k++)
                x.m_elem[i * col + j] -= *(pMi + k) * x.m_elem[k * col + j];
        }

        /* Solve Dy = z */
        // Solve the diagonal matrix, here x is y
        for (i = 0, pMi = m_elem; i < t_row; i++, pMi += t_col)
        {
            x.m_elem[i * col + j] /= *(pMi + i);
        }

        /* Solve Ux = y */
        // Solve the unit upper triangular (backward substitution)
        for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
        {
            for (k = i + 1; k < t_col; k++)
                x.m_elem[i * col + j] -= *(pMi + k) * x.m_elem[k * col + j];
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LDLT<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, Vector<t_col, t_type> &x)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    t_type *pMi;

    if (!m_isOk)
        return -1;

    /* Solve Lz = b */
    // Solve the unit lower triangular (forward substitution), here x is z
    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Dy = z */
    // Solve the diagonal matrix, here x is y
    for (i = 0, pMi = m_elem; i < t_row; i++, pMi += t_col)
    {
        x.m_elem[i] /= *(pMi + i);
    }

    /* Solve Ux = y */
    // Solve the unit upper triangular (backward substitution)
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        for (k = i + 1; k < t_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t col>
inline Matrix<t_col, col, t_type> LDLT<t_row, t_col, t_type>::Solve(const Matrix<t_row, col, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
    // Ux = y is solved by backward substitution for x

    int i, k, j;
    t_type *pMi;
    t_type x[t_col * col] = {
        0,
    };

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return Matrix<t_col, col, t_type>(x);
    }

    for (j = 0; j < col; j++)
    {
        /* Solve Lz = b */
        // Solve the unit lower triangular (forward substitution), here x is z
        for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
        {
            x[i * col + j] = b.m_elem[i * col + j];

            for (k = 0; k < i; k++)
                x[i * col + j] -= *(pMi + k) * x[k * col + j];
        }

        /* Solve Dy = z */
        // Solve the diagonal matrix (divided by diagonal elements), here x is y
        for (i = 0, pMi = m_elem; i < t_row; i++, pMi += t_col)
        {
            x[i * col + j] /= *(pMi + i);
        }

        /* Solve Ux = y */
        // Solve the upper triangular (backward substitution)
        for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
        {
            for (k = i + 1; k < t_col; k++)
                x[i * col + j] -= *(pMi + k) * x[k * col + j];
        }
    }

    return Matrix<t_col, col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> LDLT<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LDUx = b
    // where L is a unit lower triangular matrix with an all diagonal element is 1
    //       D is a diagonal matrix
    //       U is upper triangular matrix or L^(T)
    // define DUx = z and Ux = y
    // Lz = b is solved by forward substitution for z
    // Dy = z is solved by 1 / dii for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    t_type *pMi;
    t_type x[t_col] = {
        0,
    };

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return Vector<t_col, t_type>(x);
    }

    /* Solve Lz = b */
    // Solve the unit lower triangular (forward substitution), here x is z
    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
    {
        x[i] = b.m_elem[i];

        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];
    }

    /* Solve Dy = z */
    // Solve the diagonal matrix (divided by diagonal elements), here x is y
    for (i = 0, pMi = m_elem; i < t_row; i++, pMi += t_col)
    {
        x[i] /= *(pMi + i);
    }

    /* Solve Ux = y */
    // Solve the upper triangular (backward substitution)
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        for (k = i + 1; k < t_col; k++)
            x[i] -= *(pMi + k) * x[k];
    }

    return Vector<t_col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LDLT<t_row, t_col, t_type>::Inverse(Matrix<t_row, t_col, t_type> &inv)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;

    if (!m_isOk)
        return -1;

    memcpy(inv.m_elem, m_elem, sizeof(t_type) * t_row * t_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv.m_elem + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv.m_elem; j < i; pMj += t_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + t_col; k < i; k++, pMk += t_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv.m_elem; j < t_col; j++, pMj += t_row)
    {
        for (i = j, pMi = pMj; i < t_row; pMi += t_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + t_col; k < t_row; k++, pMk += t_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LDLT<t_row, t_col, t_type>::Inverse(Matrix3<t_type, t_row, t_col> &inv)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;

    if (!m_isOk)
        return -1;

    memcpy(inv.m_elem, m_elem, sizeof(t_type) * t_row * t_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv.m_elem + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv.m_elem; j < i; pMj += t_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + t_col; k < i; k++, pMk += t_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv.m_elem; j < t_col; j++, pMj += t_row)
    {
        for (i = j, pMi = pMj; i < t_row; pMi += t_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + t_col; k < t_row; k++, pMk += t_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> LDLT<t_row, t_col, t_type>::Inverse(int8_t *isOk)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type inv[t_row * t_col] = {
        0,
    };

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return Matrix<t_row, t_col, t_type>();
    }

    memcpy(inv, m_elem, sizeof(t_type) * t_row * t_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += t_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + t_col; k < i; k++, pMk += t_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv; j < t_col; j++, pMj += t_row)
    {
        for (i = j, pMi = pMj; i < t_row; pMi += t_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + t_col; k < t_row; k++, pMk += t_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

    return Matrix<t_row, t_col, t_type>(inv);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t LDLT<t_row, t_col, t_type>::InverseArray(t_type *inv)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;

    if (!m_isOk)
        return -1;

    memcpy(inv, m_elem, sizeof(t_type) * t_row * t_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += t_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + t_col; k < i; k++, pMk += t_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv; j < t_col; j++, pMj += t_row)
    {
        for (i = j, pMi = pMj; i < t_row; pMi += t_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + t_col; k < t_row; k++, pMk += t_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type *LDLT<t_row, t_col, t_type>::InverseArray(int8_t *isOk)
{
    int i, j, k;
    t_type *pMi, *pMj, *pMk;
    t_type inv[t_row * t_col] = {
        0,
    };

    if (isOk)
        *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return inv;
    }

    memcpy(inv, m_elem, sizeof(t_type) * t_row * t_col);

    /* Calculate the inverse of a unit lower triangular matrix */
    // Invert the subdiagonal part of the matrix L, for row i
    // where the diagonal elements are assumed to be 1
    for (i = 1, pMi = inv + t_col; i < t_row; i++, pMi += t_col)
    {
        for (j = 0, pMj = inv; j < i; pMj += t_col, j++)
        {
            *(pMi + j) = -*(pMi + j);

            for (k = j + 1, pMk = pMj + t_col; k < i; k++, pMk += t_col)
                *(pMi + j) -= *(pMi + k) * *(pMk + j);
        }
    }

    /* Calculate the inverse of LDLT, inv(LT) * inv(D) * inv(L) */
    // inv(LDLT) is also positive definite symmetric matrix
    for (j = 0, pMj = inv; j < t_col; j++, pMj += t_row)
    {
        for (i = j, pMi = pMj; i < t_row; pMi += t_col, i++)
        {
            if (j == i)
                *(pMi + j) = 1 / *(pMi + i);
            else
                *(pMi + j) /= *(pMi + i);

            for (k = i + 1, pMk = pMi + t_col; k < t_row; k++, pMk += t_col)
                *(pMi + j) += *(pMk + i) * *(pMk + j) / *(pMk + k);

            *(pMj + i) = *(pMi + j);
        }
    }

    return inv;
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTLDLT_TPP_
