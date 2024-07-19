/*!
\file       dtSVD.h
\brief      dtMath, Singular Value Decomposition solver class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTSVD_TPP_
#define DTMATH_DTSVD_TPP_

#include "dtSVD.h"

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline SVD<t_row, t_col, t_type>::SVD()
    : m_A(), m_U(), m_S(), m_V(), m_superDiagonal(), m_inv()
{
    m_isOk = 0;
    m_isInv = 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline SVD<t_row, t_col, t_type>::SVD(const t_type *element, const size_t n_byte)
    : m_U(), m_S(), m_V(), m_superDiagonal(), m_inv()
{
    if ((sizeof(t_type) * t_row * t_col) != n_byte)
    {
        m_isOk = 0;
        return;
    }

    if (t_row >= t_col)
    {
        memcpy(m_A, element, n_byte);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn();
        if (GivensReductionToDiagonalForm_Mxn())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_Mxn();
    }
    else
    { // Transposed A
        // for (uint16_t irow = 0; irow < t_row; ++irow)
        //     for (uint16_t icol = 0; icol < t_col; ++icol)
        //         m_A[icol * t_row + irow] = element[irow * t_col + icol];
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < t_row; ++irow)
        {
            for (cnt = t_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + t_row * icol] = element[irow * t_col + icol];
                m_A[irow + t_row * (icol + 1)] = element[irow * t_col + icol + 1];
                m_A[irow + t_row * (icol + 2)] = element[irow * t_col + icol + 2];
                m_A[irow + t_row * (icol + 3)] = element[irow * t_col + icol + 3];
            }

            for (cnt = t_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + t_row * icol] = element[irow * t_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN();
        if (GivensReductionToDiagonalForm_mxN())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_mxN();
    }

    m_isOk = 1;
    m_isInv = 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline SVD<t_row, t_col, t_type>::SVD(const Matrix<t_row, t_col, t_type> &m)
    : m_U(), m_S(), m_V(), m_superDiagonal(), m_inv()
{
    if (t_row >= t_col)
    {
        memcpy(m_A, m.m_elem, sizeof(t_type) * t_row * t_col);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn();
        if (GivensReductionToDiagonalForm_Mxn())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_Mxn();
    }
    else
    { // Transposed A
        // for (uint16_t irow = 0; irow < t_row; ++irow)
        //     for (uint16_t icol = 0; icol < t_col; ++icol)
        //         m_A[icol * t_row + irow] = m.m_elem[irow * t_col + icol];
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < t_row; ++irow)
        {
            for (cnt = t_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + t_row * icol] = m.m_elem[irow * t_col + icol];
                m_A[irow + t_row * (icol + 1)] = m.m_elem[irow * t_col + icol + 1];
                m_A[irow + t_row * (icol + 2)] = m.m_elem[irow * t_col + icol + 2];
                m_A[irow + t_row * (icol + 3)] = m.m_elem[irow * t_col + icol + 3];
            }

            for (cnt = t_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + t_row * icol] = m.m_elem[irow * t_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN();
        if (GivensReductionToDiagonalForm_mxN())
        {
            m_isOk = 0;
            return;
        }
        SortByDecreasingSingularValues_mxN();
    }

    m_isOk = 1;
    m_isInv = 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline SVD<t_row, t_col, t_type>::SVD(const Matrix3<t_type, t_row, t_col> &m)
    : m_U(), m_S(), m_V(), m_superDiagonal(), m_inv()
{
    memcpy(m_A, m.m_elem, sizeof(t_type) * t_row * t_col);

    // Compute SVD
    HouseholdersReductionToBidiagonalForm_Mxn();
    if (GivensReductionToDiagonalForm_Mxn())
    {
        m_isOk = 0;
        return;
    }
    SortByDecreasingSingularValues_Mxn();

    m_isOk = 1;
    m_isInv = 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t SVD<t_row, t_col, t_type>::Compute(const t_type *element, const size_t n_byte)
{
    if ((sizeof(t_type) * t_row * t_col) != n_byte)
    {
        m_isOk = 0;
        return -1;
    }

    if (t_row >= t_col)
    {
        memcpy(m_A, element, n_byte);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn();
        if (GivensReductionToDiagonalForm_Mxn())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_Mxn();
    }
    else
    { // Transposed A
        // for (uint16_t irow = 0; irow < t_row; ++irow)
        //     for (uint16_t icol = 0; icol < t_col; ++icol)
        //         m_A[icol * t_row + irow] = element[irow * t_col + icol];
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < t_row; ++irow)
        {
            for (cnt = t_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + t_row * icol] = element[irow * t_col + icol];
                m_A[irow + t_row * (icol + 1)] = element[irow * t_col + icol + 1];
                m_A[irow + t_row * (icol + 2)] = element[irow * t_col + icol + 2];
                m_A[irow + t_row * (icol + 3)] = element[irow * t_col + icol + 3];
            }

            for (cnt = t_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + t_row * icol] = element[irow * t_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN();
        if (GivensReductionToDiagonalForm_mxN())
        {
            m_isInv = 0;
            return -1;
        }
        SortByDecreasingSingularValues_mxN();
    }

    m_isOk = 1;
    m_isInv = 0;

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t SVD<t_row, t_col, t_type>::Compute(const Matrix<t_row, t_col, t_type> &m)
{
    if (t_row >= t_col)
    {
        memcpy(m_A, m.m_elem, sizeof(t_type) * t_row * t_col);

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_Mxn();
        if (GivensReductionToDiagonalForm_Mxn())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_Mxn();
    }
    else
    { // Transposed A
        // for (uint16_t irow = 0; irow < t_row; ++irow)
        //     for (uint16_t icol = 0; icol < t_col; ++icol)
        //         m_A[icol * t_row + irow] = m.m_elem[irow * t_col + icol];
        uint16_t cnt;
        uint16_t irow, icol;

        for (irow = 0; irow < t_row; ++irow)
        {
            for (cnt = t_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
            {
                m_A[irow + t_row * icol] = m.m_elem[irow * t_col + icol];
                m_A[irow + t_row * (icol + 1)] = m.m_elem[irow * t_col + icol + 1];
                m_A[irow + t_row * (icol + 2)] = m.m_elem[irow * t_col + icol + 2];
                m_A[irow + t_row * (icol + 3)] = m.m_elem[irow * t_col + icol + 3];
            }

            for (cnt = t_col % 4u; cnt > 0u; cnt--, icol++)
            {
                m_A[irow + t_row * icol] = m.m_elem[irow * t_col + icol];
            }
        }

        // Compute SVD
        HouseholdersReductionToBidiagonalForm_mxN();
        if (GivensReductionToDiagonalForm_mxN())
        {
            m_isOk = 0;
            return -1;
        }
        SortByDecreasingSingularValues_mxN();
    }

    m_isOk = 1;
    m_isInv = 0;

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t SVD<t_row, t_col, t_type>::Compute(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_A, m.m_elem, sizeof(t_type) * t_row * t_col);

    // Compute SVD
    HouseholdersReductionToBidiagonalForm_Mxn();
    if (GivensReductionToDiagonalForm_Mxn())
    {
        m_isOk = 0;
        return -1;
    }
    SortByDecreasingSingularValues_Mxn();

    m_isOk = 1;
    m_isInv = 0;

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_row, t_type> SVD<t_row, t_col, t_type>::GetMatrixU() const
{
    if (t_row > t_col)
    {
        t_type U[t_row * t_row]{0};

        for (int i = 0; i < t_row; i++)
            memcpy(&U[i * t_row], &m_U[i * t_col], sizeof(t_type) * t_col);

        return Matrix<t_row, t_row, t_type>(U);
    }
    else
    {
        return Matrix<t_row, t_row, t_type>(m_U);
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> SVD<t_row, t_col, t_type>::GetMatrixS() const
{
    return Matrix<t_row, t_col, t_type>('d', m_S, sizeof(m_S));
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_col, t_col, t_type> SVD<t_row, t_col, t_type>::GetMatrixV() const
{
    if (t_row < t_col)
    {
        t_type V[t_col * t_col]{0};

        for (int i = 0; i < t_col; i++)
            memcpy(&V[i * t_col], &m_V[i * t_row], sizeof(t_type) * t_row);

        return Matrix<t_col, t_col, t_type>(V);
    }
    else
    {
        return Matrix<t_col, t_col, t_type>(m_V);
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t col>
inline int8_t SVD<t_row, t_col, t_type>::Solve(const Matrix<t_row, col, t_type> &b, Matrix<t_col, col, t_type> &x, t_type tolerance)
{
    int i, j, k, c;
    t_type *pU, *pV;
    t_type sum;

    if (!m_isOk)
        return -1;

    if (t_row >= t_col)
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_col;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < t_col; i++, pV += t_col)
            {
                x.m_elem[i * col + c] = 0;
                for (j = 0; j < t_col; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < t_row; k++, pU += t_col)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x.m_elem[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_row;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < t_col; i++, pV += t_row)
            {
                x.m_elem[i * col + c] = 0;
                for (j = 0; j < t_row; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < t_row; k++, pU += t_row)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x.m_elem[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t SVD<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, Vector<t_col, t_type> &x, t_type tolerance)
{
    int i, j, k;
    t_type *pU, *pV;
    t_type sum;

    if (!m_isOk)
        return -1;

    if (t_row >= t_col)
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_col;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < t_col; i++, pV += t_col)
        {
            x.m_elem[i] = 0;
            for (j = 0; j < t_col; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < t_row; k++, pU += t_col)
                        sum += *(pU + j) * b.m_elem[k];

                    x.m_elem[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_row;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < t_col; i++, pV += t_row)
        {
            x.m_elem[i] = 0;
            for (j = 0; j < t_row; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < t_row; k++, pU += t_row)
                        sum += *(pU + j) * b.m_elem[k];

                    x.m_elem[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t col>
inline Matrix<t_col, col, t_type> SVD<t_row, t_col, t_type>::Solve(const Matrix<t_row, col, t_type> &b, int8_t *isOk, t_type tolerance)
{
    int i, j, k, c;
    t_type *pU, *pV;
    t_type sum;
    t_type x[t_col * col]{0};

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return Matrix<t_col, col, t_type>();
    }

    if (t_row >= t_col)
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_col;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < t_col; i++, pV += t_col)
            {
                x[i * col + c] = 0;
                for (j = 0; j < t_col; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < t_row; k++, pU += t_col)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_row;
        if (tolerance < sum)
            tolerance = sum;

        for (c = 0; c < col; c++)
        {
            for (i = 0, pV = m_V; i < t_col; i++, pV += t_row)
            {
                x[i * col + c] = 0;
                for (j = 0; j < t_row; j++)
                {
                    if (m_S[j] > tolerance)
                    {
                        for (k = 0, sum = 0, pU = m_U; k < t_row; k++, pU += t_row)
                            sum += *(pU + j) * b.m_elem[k * col + c];

                        x[i * col + c] += sum * *(pV + j) / m_S[j];
                    }
                }
            }
        }
    }

    if (isOk) *isOk = 1;

    return Matrix<t_col, col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> SVD<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, int8_t *isOk, t_type tolerance)
{
    int i, j, k;
    t_type *pU, *pV;
    t_type sum;
    t_type x[t_col]{0};

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return Vector<t_col, t_type>();
    }

    if (t_row >= t_col)
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_col;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < t_col; i++, pV += t_col)
        {
            x[i] = 0;
            for (j = 0; j < t_col; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < t_row; k++, pU += t_col)
                        sum += *(pU + j) * b.m_elem[k];

                    x[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }
    else
    {
        sum = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_row;
        if (tolerance < sum)
            tolerance = sum;

        for (i = 0, pV = m_V; i < t_col; i++, pV += t_row)
        {
            x[i] = 0;
            for (j = 0; j < t_row; j++)
            {
                if (m_S[j] > tolerance)
                {
                    for (k = 0, sum = 0, pU = m_U; k < t_row; k++, pU += t_row)
                        sum += *(pU + j) * b.m_elem[k];

                    x[i] += sum * *(pV + j) / m_S[j];
                }
            }
        }
    }

    if (isOk) *isOk = 1;

    return Vector<t_col, t_type>(x);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t SVD<t_row, t_col, t_type>::Inverse(Matrix<t_col, t_row, t_type> &inv, t_type tolerance)
{
    int i, j, k;
    t_type *pU, *pV, *pInvA;
    t_type tol;

    if (!m_isOk) return -1;

    if (m_isInv)
    {
        memcpy(inv.m_elem, m_inv, sizeof(m_inv));
        return 0;
    }

    if (t_row >= t_col)
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_col;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv.m_elem; i < t_col; i++, pV += t_col)
        {
            for (j = 0, pU = m_U; j < t_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < t_col; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }
    else
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_row;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv.m_elem; i < t_col; i++, pV += t_row)
        {
            for (j = 0, pU = m_U; j < t_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < t_row; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }

    memcpy(m_inv, inv.m_elem, sizeof(m_inv));
    m_isInv = 1;

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t SVD<t_row, t_col, t_type>::Inverse(Matrix3<t_type, t_col, t_row> &inv, t_type tolerance)
{
    int i, j, k;
    t_type *pU, *pV, *pInvA;
    t_type tol;

    if (!m_isOk) return -1;

    if (m_isInv)
    {
        memcpy(inv.m_elem, m_inv, sizeof(m_inv));
        return 0;
    }

    tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_col;
    if (tolerance < tol)
        tolerance = tol;

    for (i = 0, pV = m_V, pInvA = inv.m_elem; i < t_col; i++, pV += t_col)
    {
        for (j = 0, pU = m_U; j < t_row; j++, pInvA++)
        {
            for (k = 0, *pInvA = 0; k < t_col; k++, pU++)
                if (m_S[k] > tolerance)
                    *pInvA += *(pV + k) * *pU / m_S[k];
        }
    }

    memcpy(m_inv, inv.m_elem, sizeof(m_inv));
    m_isInv = 1;

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_col, t_row, t_type> SVD<t_row, t_col, t_type>::Inverse(int8_t *isOk, t_type tolerance)
{
    int i, j, k;
    t_type *pU, *pV, *pInvA;
    t_type tol;

    if (!m_isOk)
    {
        if (isOk) *isOk = 0;
        return Matrix<t_col, t_row, t_type>();
    }

    if (m_isInv)
    {
        return Matrix<t_col, t_row, t_type>(m_inv);
    }

    memset(m_inv, 0, sizeof(t_type) * t_row * t_col); // why must?

    if (t_row >= t_col)
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_col;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = m_inv; i < t_col; i++, pV += t_col)
        {
            for (j = 0, pU = m_U; j < t_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < t_col; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }
    else
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_row;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = m_inv; i < t_col; i++, pV += t_row)
        {
            for (j = 0, pU = m_U; j < t_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < t_row; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }

    if (isOk) *isOk = 1;
    m_isInv = 1;

    return Matrix<t_col, t_row, t_type>(m_inv);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t SVD<t_row, t_col, t_type>::InverseArray(t_type *inv, t_type tolerance)
{
    int i, j, k;
    t_type *pU, *pV, *pInvA;
    t_type tol;

    if (!m_isOk) return -1;

    if (m_isInv)
    {
        memcpy(inv, m_inv, sizeof(m_inv));
        return 0;
    }

    if (t_row >= t_col)
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_col;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv; i < t_col; i++, pV += t_col)
        {
            for (j = 0, pU = m_U; j < t_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < t_col; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }
    else
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_row;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = inv; i < t_col; i++, pV += t_row)
        {
            for (j = 0, pU = m_U; j < t_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < t_row; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }

    memcpy(m_inv, inv, sizeof(m_inv));
    m_isInv = 1;

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type *SVD<t_row, t_col, t_type>::InverseArray(int8_t *isOk, t_type tolerance)
{
    int i, j, k;
    t_type *pU, *pV, *pInvA;
    t_type tol;

    if (!m_isOk)
    {
        if (isOk) *isOk = 0;
        return m_inv;
    }

    if (m_isInv)
    {
        return m_inv;
    }

    if (t_row >= t_col)
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_col;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = m_inv; i < t_col; i++, pV += t_col)
        {
            for (j = 0, pU = m_U; j < t_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < t_col; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }
    else
    {
        tol = std::numeric_limits<t_type>::epsilon() * m_S[0] * (t_type)t_row;
        if (tolerance < tol)
            tolerance = tol;

        for (i = 0, pV = m_V, pInvA = m_inv; i < t_col; i++, pV += t_row)
        {
            for (j = 0, pU = m_U; j < t_row; j++, pInvA++)
            {
                for (k = 0, *pInvA = 0; k < t_row; k++, pU++)
                    if (m_S[k] > tolerance)
                        *pInvA += *(pV + k) * *pU / m_S[k];
            }
        }
    }

    if (isOk) *isOk = 1;
    m_isInv = 1;

    return m_inv;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void SVD<t_row, t_col, t_type>::HouseholdersReductionToBidiagonalForm_Mxn()
{
    int i, j, k, ip1;
    t_type s, s2, si, scale;
    t_type *pU, *pUi, *pV, *pVi;
    t_type halfNormSquared;

    // Copy A to U
    memcpy(m_U, m_A, sizeof(t_type) * t_row * t_col);

    m_S[0] = 0;
    s = 0;
    scale = 0;

    for (i = 0, pUi = m_U, ip1 = 1; i < t_col; pUi += t_col, i++, ip1++)
    {
        m_superDiagonal[i] = scale * s;

        // Perform Householder transform on columns.
        // Calculate the normed squared of the i-th column vector starting at row i.
        for (j = i, pU = pUi, scale = 0; j < t_row; j++, pU += t_col)
            scale += std::abs(*(pU + i));

        if (scale > std::numeric_limits<t_type>::epsilon())
        {
            for (j = i, pU = pUi, s2 = 0; j < t_row; j++, pU += t_col)
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

            for (j = ip1; j < t_col; j++)
            {
                for (k = i, si = 0, pU = pUi; k < t_row; k++, pU += t_col)
                    si += *(pU + i) * *(pU + j);

                si /= halfNormSquared;

                for (k = i, pU = pUi; k < t_row; k++, pU += t_col)
                    *(pU + j) += si * *(pU + i);
            }
        }

        for (j = i, pU = pUi; j < t_row; j++, pU += t_col)
            *(pU + i) *= scale;

        m_S[i] = s * scale;

        // Perform Householder transform on rows.
        // Calculate the normed squared of the i-th row vector starting at column i.
        s = 0;
        scale = 0;
        if (i >= t_row || i == (t_col - 1))
            continue;

        for (j = ip1; j < t_col; j++)
            scale += std::abs(*(pUi + j));

        if (scale > std::numeric_limits<t_type>::epsilon())
        {
            for (j = ip1, s2 = 0; j < t_col; j++)
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

            for (k = ip1; k < t_col; k++)
                m_superDiagonal[k] = *(pUi + k) / halfNormSquared;

            if (i < (t_row - 1))
            {
                for (j = ip1, pU = pUi + t_col; j < t_row; j++, pU += t_col)
                {
                    for (k = ip1, si = 0; k < t_col; k++)
                        si += *(pUi + k) * *(pU + k);

                    for (k = ip1; k < t_col; k++)
                        *(pU + k) += si * m_superDiagonal[k];
                }
            }
            for (k = ip1; k < t_col; k++)
                *(pUi + k) *= scale;
        }
    }

    /* Update V */
    pUi = m_U + t_col * (t_col - 2);
    pVi = m_V + t_col * (t_col - 1);
    *(pVi + t_col - 1) = 1;
    s = m_superDiagonal[t_col - 1];
    pVi -= t_col;

    for (i = t_col - 2, ip1 = t_col - 1; i >= 0; i--, pUi -= t_col, pVi -= t_col, ip1--)
    {
        if (std::abs(s) > std::numeric_limits<t_type>::epsilon())
        {
            pV = pVi + t_col;

            for (j = ip1; j < t_col; j++, pV += t_col)
                *(pV + i) = (*(pUi + j) / *(pUi + ip1)) / s;

            for (j = ip1; j < t_col; j++)
            {
                si = 0;

                for (k = ip1, pV = pVi + t_col; k < t_col; k++, pV += t_col)
                    si += *(pUi + k) * *(pV + j);

                for (k = ip1, pV = pVi + t_col; k < t_col; k++, pV += t_col)
                    *(pV + j) += si * *(pV + i);
            }
        }

        pV = pVi + t_col;
        for (j = ip1; j < t_col; j++, pV += t_col)
        {
            *(pVi + j) = 0;
            *(pV + i) = 0;
        }

        *(pVi + i) = 1;
        s = m_superDiagonal[i];
    }

    /* Update U */
    pUi = m_U + t_col * (t_col - 1);
    for (i = t_col - 1, ip1 = t_col; i >= 0; ip1 = i, i--, pUi -= t_col)
    {
        s = m_S[i];

        for (j = ip1; j < t_col; j++)
            *(pUi + j) = 0;

        if (std::abs(s) > std::numeric_limits<t_type>::epsilon())
        {
            for (j = ip1; j < t_col; j++)
            {
                si = 0;
                pU = pUi + t_col;

                for (k = ip1; k < t_row; k++, pU += t_col)
                    si += *(pU + i) * *(pU + j);

                si = (si / *(pUi + i)) / s;
                for (k = i, pU = pUi; k < t_row; k++, pU += t_col)
                    *(pU + j) += si * *(pU + i);
            }

            for (j = i, pU = pUi; j < t_row; j++, pU += t_col)
                *(pU + i) /= s;
        }
        else
        {
            for (j = i, pU = pUi; j < t_row; j++, pU += t_col)
                *(pU + i) = 0;
        }

        *(pUi + i) += 1;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t SVD<t_row, t_col, t_type>::GivensReductionToDiagonalForm_Mxn()
{
    t_type epsilon;
    t_type c, s;
    t_type f, g, h;
    t_type x, y, z;
    t_type *pU, *pV;
    int i, j, k, m;
    int rotation_test;
    int iteration_count;

    for (i = 0, x = 0; i < t_col; i++)
    {
        y = std::abs(m_S[i]) + std::abs(m_superDiagonal[i]);
        if (x < y)
            x = y;
    }

    epsilon = x * std::numeric_limits<t_type>::epsilon();

    for (k = t_col - 1; k >= 0; k--)
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

                    for (j = 0, pU = m_U; j < t_row; j++, pU += t_col)
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
                    for (j = 0, pV = m_V; j < t_col; j++, pV += t_col)
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

                    for (j = 0, pV = m_V; j < t_col; j++, pV += t_col)
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

                    for (j = 0, pU = m_U; j < t_row; j++, pU += t_col)
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

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void SVD<t_row, t_col, t_type>::SortByDecreasingSingularValues_Mxn()
{
    int i, j, maxIdx;
    t_type temp;
    t_type *pM1, *pM2;

    for (i = 0; i < t_col - 1; i++)
    {
        maxIdx = i;

        for (j = i + 1; j < t_col; j++)
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

        for (j = 0; j < t_row; j++, pM1 += t_col, pM2 += t_col)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }

        pM1 = m_V + maxIdx;
        pM2 = m_V + i;

        for (j = 0; j < t_col; j++, pM1 += t_col, pM2 += t_col)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void SVD<t_row, t_col, t_type>::HouseholdersReductionToBidiagonalForm_mxN()
{
    int i, j, k, ip1;
    t_type s, s2, si, scale;
    t_type *pV, *pVi, *pU, *pUi;
    t_type halfNormSquared;

    // Copy A to V
    memcpy(m_V, m_A, sizeof(t_type) * t_col * t_row);

    m_S[0] = 0;
    s = 0;
    scale = 0;

    for (i = 0, pVi = m_V, ip1 = 1; i < t_row; pVi += t_row, i++, ip1++)
    {
        m_superDiagonal[i] = scale * s;

        // Perform Householder transform on columns.
        // Calculate the normed squared of the i-th column vector starting at row i.
        for (j = i, pV = pVi, scale = 0; j < t_col; j++, pV += t_row)
            scale += std::abs(*(pV + i));

        if (scale > std::numeric_limits<t_type>::epsilon())
        {
            for (j = i, pV = pVi, s2 = 0; j < t_col; j++, pV += t_row)
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

            for (j = ip1; j < t_row; j++)
            {
                for (k = i, si = 0, pV = pVi; k < t_col; k++, pV += t_row)
                    si += *(pV + i) * *(pV + j);

                si /= halfNormSquared;

                for (k = i, pV = pVi; k < t_col; k++, pV += t_row)
                    *(pV + j) += si * *(pV + i);
            }
        }
        for (j = i, pV = pVi; j < t_col; j++, pV += t_row)
            *(pV + i) *= scale;
        m_S[i] = s * scale;

        // Perform Householder transform on rows.
        // Calculate the normed squared of the i-th row vector starting at column i.
        s = 0;
        scale = 0;
        if (i >= t_col || i == (t_row - 1))
            continue;

        for (j = ip1; j < t_row; j++)
            scale += std::abs(*(pVi + j));

        if (scale > std::numeric_limits<t_type>::epsilon())
        {
            for (j = ip1, s2 = 0; j < t_row; j++)
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

            for (k = ip1; k < t_row; k++)
                m_superDiagonal[k] = *(pVi + k) / halfNormSquared;

            if (i < (t_col - 1))
            {
                for (j = ip1, pV = pVi + t_row; j < t_col; j++, pV += t_row)
                {
                    for (k = ip1, si = 0; k < t_row; k++)
                        si += *(pVi + k) * *(pV + k);

                    for (k = ip1; k < t_row; k++)
                        *(pV + k) += si * m_superDiagonal[k];
                }
            }
            for (k = ip1; k < t_row; k++)
                *(pVi + k) *= scale;
        }
    }

    /* Update U */
    pVi = m_V + t_row * (t_row - 2);
    pUi = m_U + t_row * (t_row - 1);
    *(pUi + t_row - 1) = 1;
    s = m_superDiagonal[t_row - 1];
    pUi -= t_row;

    for (i = t_row - 2, ip1 = t_row - 1; i >= 0; i--, pVi -= t_row, pUi -= t_row, ip1--)
    {
        if (std::abs(s) > std::numeric_limits<t_type>::epsilon())
        {
            pU = pUi + t_row;

            for (j = ip1; j < t_row; j++, pU += t_row)
                *(pU + i) = (*(pVi + j) / *(pVi + ip1)) / s;

            for (j = ip1; j < t_row; j++)
            {
                si = 0;

                for (k = ip1, pU = pUi + t_row; k < t_row; k++, pU += t_row)
                    si += *(pVi + k) * *(pU + j);

                for (k = ip1, pU = pUi + t_row; k < t_row; k++, pU += t_row)
                    *(pU + j) += si * *(pU + i);
            }
        }

        pU = pUi + t_row;
        for (j = ip1; j < t_row; j++, pU += t_row)
        {
            *(pUi + j) = 0;
            *(pU + i) = 0;
        }

        *(pUi + i) = 1;
        s = m_superDiagonal[i];
    }

    /* Update V */
    pVi = m_V + t_row * (t_row - 1);
    for (i = t_row - 1, ip1 = t_row; i >= 0; ip1 = i, i--, pVi -= t_row)
    {
        s = m_S[i];

        for (j = ip1; j < t_row; j++)
            *(pVi + j) = 0;

        if (std::abs(s) > std::numeric_limits<t_type>::epsilon())
        {
            for (j = ip1; j < t_row; j++)
            {
                si = 0;
                pV = pVi + t_row;

                for (k = ip1; k < t_col; k++, pV += t_row)
                    si += *(pV + i) * *(pV + j);

                si = (si / *(pVi + i)) / s;
                for (k = i, pV = pVi; k < t_col; k++, pV += t_row)
                    *(pV + j) += si * *(pV + i);
            }

            for (j = i, pV = pVi; j < t_col; j++, pV += t_row)
                *(pV + i) /= s;
        }
        else
        {
            for (j = i, pV = pVi; j < t_col; j++, pV += t_row)
                *(pV + i) = 0;
        }

        *(pVi + i) += 1;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t SVD<t_row, t_col, t_type>::GivensReductionToDiagonalForm_mxN()
{
    t_type epsilon;
    t_type c, s;
    t_type f, g, h;
    t_type x, y, z;
    t_type *pV, *pU;
    int i, j, k, m;
    int rotation_test;
    int iteration_count;

    for (i = 0, x = 0; i < t_row; i++)
    {
        y = std::abs(m_S[i]) + std::abs(m_superDiagonal[i]);
        if (x < y)
            x = y;
    }

    epsilon = x * std::numeric_limits<t_type>::epsilon();

    for (k = t_row - 1; k >= 0; k--)
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

                    for (j = 0, pV = m_V; j < t_col; j++, pV += t_row)
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
                    for (j = 0, pU = m_U; j < t_row; j++, pU += t_row)
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

                    for (j = 0, pU = m_U; j < t_row; j++, pU += t_row)
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

                    for (j = 0, pV = m_V; j < t_col; j++, pV += t_row)
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

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void SVD<t_row, t_col, t_type>::SortByDecreasingSingularValues_mxN()
{
    int i, j, maxIdx;
    t_type temp;
    t_type *pM1, *pM2;

    for (i = 0; i < t_row - 1; i++)
    {
        maxIdx = i;

        for (j = i + 1; j < t_row; j++)
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

        for (j = 0; j < t_col; j++, pM1 += t_row, pM2 += t_row)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }

        pM1 = m_U + maxIdx;
        pM2 = m_U + i;

        for (j = 0; j < t_row; j++, pM1 += t_row, pM2 += t_row)
        {
            temp = *pM1;
            *pM1 = *pM2;
            *pM2 = temp;
        }
    }
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTSVD_TPP_
