/*!
\file       dtNoPivLU.h
\brief      dtMath, LU Decomposition without pivoting(Doolittle form) class, A = LU
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTNO_PIV_LU_TPP_
#define DTMATH_DTNO_PIV_LU_TPP_

#include "dtNoPivLU.h"

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline NoPivLU<t_row, t_col, t_type>::NoPivLU() : m_isOk(0)
{
    memset(m_elem, 0, sizeof(t_type) * t_row * t_col);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline NoPivLU<t_row, t_col, t_type>::NoPivLU(const t_type *element, const size_t n_byte)
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
inline NoPivLU<t_row, t_col, t_type>::NoPivLU(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline NoPivLU<t_row, t_col, t_type>::NoPivLU(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t NoPivLU<t_row, t_col, t_type>::Compute()
{
    if (t_row != t_col)
    {
        m_isOk = 0;
        return -1;
    }

    int x, i, j, k;
    t_type *pMx, *pMi, *pMk;

    for (x = 0, pMx = m_elem; x < t_row; pMx += t_col, x++)
    {
        /* find the uppper triangular matrix element for row x */
        for (j = x; j < t_col; j++)
        {
            // U(x,j) = A(x,j) - sum of (L(x,k) * U(k,j))
            for (k = 0, pMk = m_elem; k < x; pMk += t_col, k++)
                *(pMx + j) -= *(pMx + k) * *(pMk + j); // in-place
        }

        if (std::abs(*(pMx + x)) <= std::numeric_limits<t_type>::epsilon())
        {
            m_isOk = 0;
            return -1; // matrix is singular, if Diagonal is 0
        }

        /* find the lower triangular matrix element for col x */
        for (i = x + 1, pMi = pMx + t_row; i < t_row; pMi += t_row, i++)
        {
            // L(i,x) = [A(i,x) - sum of (L(i,k) * U(k,x))] / Uxx
            for (k = 0, pMk = m_elem; k < x; pMk += t_row, k++)
                *(pMi + x) -= *(pMi + k) * *(pMk + x);
            *(pMi + x) /= *(pMx + x);
        }
    }

    m_isOk = 1;
    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t NoPivLU<t_row, t_col, t_type>::Compute(const t_type *element, const size_t n_byte)
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
inline int8_t NoPivLU<t_row, t_col, t_type>::Compute(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t NoPivLU<t_row, t_col, t_type>::Compute(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> NoPivLU<t_row, t_col, t_type>::GetMatrix() const
{
    return Matrix<t_row, t_col, t_type>(m_elem);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> NoPivLU<t_row, t_col, t_type>::GetMatrixL() const
{
    uint16_t i, j;
    t_type L[t_row * t_col]{0};

    /* Set diagonal elements as 1 */
    for (i = 0; i < t_row; i++)
    {
        L[i * (t_col + 1)] = 1;
    }

    /* Update remaining matrix from m_elem to L*/
    for (i = 1; i < t_row; i++)
        for (j = 0; j < i; j++)
            L[i * t_col + j] = m_elem[i * t_col + j];

    return Matrix<t_row, t_col, t_type>(L);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> NoPivLU<t_row, t_col, t_type>::GetMatrixU() const
{
    uint16_t i, j;
    t_type U[t_row * t_col]{0};

    for (i = 0; i < t_row; i++)
        for (j = i; j < t_col; j++)
            U[i * t_col + j] = m_elem[i * t_col + j];

    return Matrix<t_row, t_col, t_type>(U);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t NoPivLU<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, Vector<t_col, t_type> &x)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    if (!m_isOk)
    {
        assert(false && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, k;
    t_type *pMi;

    /* Solve Ly = b */
    // Solve the unit lower triangular (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
    {
        x.m_elem[i] = b.m_elem[i];
        for (k = 0; k < i; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];
    }

    /* Solve Ux = y */
    // Solve the upper triangular (backward substitution)
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        for (k = i + 1; k < t_col; k++)
            x.m_elem[i] -= *(pMi + k) * x.m_elem[k];

        // if (std::abs(*(pMi + k)) <= std::numeric_limits<t_type>::epsilon()) return -1; // The matrix U is singular

        x.m_elem[i] /= *(pMi + i);
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> NoPivLU<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LUx = b
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return Vector<t_col, t_type>();
    }

    int i, k;
    t_type *pMi;
    t_type x[t_col]{0};

    /* Solve Ly = b */
    // Unit Lower Triangular Solve (forward substitution)
    for (i = 0, pMi = m_elem; i < t_row; i++, pMi += t_col)
    {
        x[i] = b.m_elem[i];
        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];
    }

    /* Solve Ux = y */
    // Upper Triangular Solve (backward substitution)
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        for (k = i + 1; k < t_col; k++)
            x[i] -= *(pMi + k) * x[k];

        // if (std::abs(*(pU + k)) <= std::numeric_limits<t_type>::epsilon()) m_isOk = 0; // The matrix U is singular

        x[i] /= *(pMi + i);
    }

    if (isOk) *isOk = 1;
    return Vector<t_col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t NoPivLU<t_row, t_col, t_type>::Inverse(Matrix<t_row, t_col, t_type> &inv)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    t_type *pMi;
    t_type *pInvMi, *pInvMj, *pInvMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    for (i = 0; i < t_row; i++)
        invL[i * (t_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + t_col;
    pInvMi = invL + t_col;
    for (i = 1; i < t_row; i++, pMi += t_col, pInvMi += t_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += t_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + t_col;
            for (k = j + 1; k < i; k++, pInvMk += t_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < t_row; k++, pMi += (t_col + 1), pInvMk += (t_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + t_col * (t_row - 2);
    pInvMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, pMi -= t_col, pInvMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + t_col;
            for (k = i + 1; k <= j; k++, pInvMk += t_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < t_row; i++)
        for (j = 0; j < t_col; j++)
            for (k = 0; k < t_col; k++)
                inv.m_elem[i * t_col + j] += invU[i * t_col + k] * invL[k * t_col + j];

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t NoPivLU<t_row, t_col, t_type>::Inverse(Matrix3<t_type, t_row, t_col> &inv)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    t_type *pMi;
    t_type *pInvMi, *pInvMj, *pInvMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    for (i = 0; i < t_row; i++)
        invL[i * (t_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + t_col;
    pInvMi = invL + t_col;
    for (i = 1; i < t_row; i++, pMi += t_col, pInvMi += t_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += t_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + t_col;
            for (k = j + 1; k < i; k++, pInvMk += t_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < t_row; k++, pMi += (t_col + 1), pInvMk += (t_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + t_col * (t_row - 2);
    pInvMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, pMi -= t_col, pInvMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + t_col;
            for (k = i + 1; k <= j; k++, pInvMk += t_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < t_row; i++)
        for (j = 0; j < t_col; j++)
            for (k = 0; k < t_col; k++)
                inv.m_elem[i * t_col + j] += invU[i * t_col + k] * invL[k * t_col + j];

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> NoPivLU<t_row, t_col, t_type>::Inverse(int8_t *isOk)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return Matrix<t_row, t_col, t_type>();
    }

    int i, j, k;
    t_type *pMi;
    t_type *pInvMi, *pInvMj, *pInvMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    for (i = 0; i < t_row; i++)
        invL[i * (t_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + t_col;
    pInvMi = invL + t_col;
    for (i = 1; i < t_row; i++, pMi += t_col, pInvMi += t_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += t_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + t_col;
            for (k = j + 1; k < i; k++, pInvMk += t_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < t_row; k++, pMi += (t_col + 1), pInvMk += (t_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + t_col * (t_row - 2);
    pInvMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, pMi -= t_col, pInvMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + t_col;
            for (k = i + 1; k <= j; k++, pInvMk += t_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    memset(m_inv, 0, sizeof(t_type) * t_row * t_col);
    for (i = 0; i < t_row; i++)
        for (j = 0; j < t_col; j++)
            for (k = 0; k < t_col; k++)
                m_inv[i * t_col + j] += invU[i * t_col + k] * invL[k * t_col + j];

    if (isOk) *isOk = 1;
    return Matrix<t_row, t_col, t_type>(m_inv);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t NoPivLU<t_row, t_col, t_type>::InverseArray(t_type *inv)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    t_type *pMi;
    t_type *pInvMi, *pInvMj, *pInvMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    for (i = 0; i < t_row; i++)
        invL[i * (t_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + t_col;
    pInvMi = invL + t_col;
    for (i = 1; i < t_row; i++, pMi += t_col, pInvMi += t_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += t_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + t_col;
            for (k = j + 1; k < i; k++, pInvMk += t_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < t_row; k++, pMi += (t_col + 1), pInvMk += (t_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + t_col * (t_row - 2);
    pInvMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, pMi -= t_col, pInvMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + t_col;
            for (k = i + 1; k <= j; k++, pInvMk += t_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < t_row; i++)
        for (j = 0; j < t_col; j++)
            for (k = 0; k < t_col; k++)
                inv[i * t_col + j] += invU[i * t_col + k] * invL[k * t_col + j];

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type *NoPivLU<t_row, t_col, t_type>::InverseArray(int8_t *isOk)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return nullptr;
    }

    int i, j, k;
    t_type *pMi;
    t_type *pInvMi, *pInvMj, *pInvMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};
    memset(m_inv, 0, sizeof(t_type) * t_row * t_col);

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    for (i = 0; i < t_row; i++)
        invL[i * (t_col + 1)] = 1;

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    pMi = m_elem + t_col;
    pInvMi = invL + t_col;
    for (i = 1; i < t_row; i++, pMi += t_col, pInvMi += t_col)
    {
        pInvMj = invL;
        for (j = 0; j < i; j++, pInvMj += t_col)
        {
            *(pInvMi + j) = -*(pMi + j);

            pInvMk = pInvMj + t_col;
            for (k = j + 1; k < i; k++, pInvMk += t_col)
            {
                *(pInvMi + j) -= *(pMi + k) * *(pInvMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    pMi = m_elem;
    pInvMk = invU;
    for (k = 0; k < t_row; k++, pMi += (t_col + 1), pInvMk += (t_col + 1))
    {
        // if (std::abs(*pMi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *pInvMk = 1 / *pMi;
        *pInvMk = 1 / *pMi;
    }

    // Invert the remaining upper triangular matrix U.
    pMi = m_elem + t_col * (t_row - 2);
    pInvMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, pMi -= t_col, pInvMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            pInvMk = pInvMi + t_col;
            for (k = i + 1; k <= j; k++, pInvMk += t_col)
            {
                sum += *(pMi + k) * *(pInvMk + j);
            }
            *(pInvMi + j) = -*(pInvMi + i) * sum;
        }
    }

    /* Inv(A) = inv(U) * inv(L) */
    for (i = 0; i < t_row; i++)
        for (j = 0; j < t_col; j++)
            for (k = 0; k < t_col; k++)
                m_inv[i * t_col + j] += invU[i * t_col + k] * invL[k * t_col + j];

    if (isOk) *isOk = 1;
    return m_inv;
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTNO_PIV_LU_TPP_