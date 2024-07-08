/*!
\file       dtHouseholderQR.h
\brief      dtMath, QR Decomposition without pivoting(Householder method) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTQR_TPP_
#define DTMATH_DTQR_TPP_

#include "dtQR.h"

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline QR<t_row, t_col, t_type>::QR() : m_elem(), m_R(), m_Q(), m_isOk(0)
{
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline QR<t_row, t_col, t_type>::QR(const t_type *element, const size_t n_byte)
    : m_R(), m_Q()
{
    if ((sizeof(t_type) * t_row * t_col) != n_byte)
        m_isOk = 0;
    else
    {
        memcpy(m_elem, element, n_byte);
        Compute();
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline QR<t_row, t_col, t_type>::QR(const Matrix<t_row, t_col, t_type> &m)
    : m_R(), m_Q()
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline QR<t_row, t_col, t_type>::QR(const Matrix3<t_type, t_row, t_col> &m)
    : m_R(), m_Q()
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t QR<t_row, t_col, t_type>::Compute()
{
    // A = Q*R
    int col = (t_row > t_col) ? t_col : t_row;
    int i, j, k, jj;

    t_type H1n[t_row * t_row]; // size: mxm
    t_type Hn1[t_row * t_col] = {
        0,
    }; // Hn*...*H1*A, mxn
    t_type H[t_row * t_row];
    t_type u[t_row];
    t_type *pQi, *pQij;
    t_type *pRi, *pRij;
    t_type *pHi, *pHkj;
    t_type *pH1ni;
    t_type *pHn1kj;
    t_type *pUi, *pUij;
    t_type sqNorm;

    /* QR Decomposition */
    // R = Hn * ... * H1 * A, size:mxn
    // Q = H1 * ... * Hn, size:mxm
    // H = I - 2 * u * uT / uT * u, size:mxm
    // u = a + sign(a1) * norm(a) * e1
    // a is column vector of R or A
    // e1 = [1 0 .. 0]T

    // copy mat A to mat Hn1A
    memcpy(Hn1, m_elem, sizeof(t_type) * t_row * t_col);

    for (j = 0, pUi = Hn1; j < col; j++, pUi += t_col)
    {
        // calculate the vector u
        sqNorm = 0;
        for (i = j, pUij = pUi + j; i < t_row; i++, pUij += t_col)
        {
            u[i] = *(pUij);              // vector a (column of R or Hn1*A)
            sqNorm += *(pUij) * *(pUij); // squared norm of vector a
        }

        u[j] += (u[j] > 0) ? std::sqrt(sqNorm) : -std::sqrt(sqNorm);

        // uT*u, squared norm of vector u
        sqNorm = 0;
        for (i = j; i < t_row; i++)
            sqNorm += *(u + i) * *(u + i);

        if (sqNorm > std::numeric_limits<t_type>::epsilon())
        {
            // calculate the householder matrix H = I - 2*u*uT / uT*u
            for (i = j, pHi = H + j * t_row; i < t_row; i++, pHi += t_row)
            {
                for (k = j; k < t_row; k++)
                {
                    if (i == k)
                        *(pHi + k) = 1 - 2 * *(u + i) * *(u + k) / sqNorm;
                    else
                        *(pHi + k) = -2 * *(u + i) * *(u + k) / sqNorm;
                }
            }

            if (j == 0)
            {
                // Q = I*H1;
                memcpy(H1n, H, sizeof(t_type) * t_row * t_row);
                memcpy(m_Q, H, sizeof(t_type) * t_row * t_row);

                // R = H1*A(=Hn1)
                pRi = m_R;
                pHi = H;
                for (i = 0; i < t_row; i++, pRi += t_col, pHi += t_row)
                {
                    for (jj = 0; jj < t_col; jj++)
                    {
                        pHn1kj = Hn1 + jj;
                        pRij = pRi + jj;
                        for (k = 0; k < t_row; k++, pHn1kj += t_col)
                            *(pRij) += *(pHi + k) * *pHn1kj;
                    }
                }
                memcpy(Hn1, m_R, sizeof(t_type) * t_row * t_col);
            }
            else
            {
                // Q = (H1*...*Hn-1) * H
                pQi = m_Q;
                pH1ni = H1n;
                for (i = 0; i < t_row; i++, pQi += t_row, pH1ni += t_row)
                {
                    for (jj = j; jj < t_row; jj++)
                    {
                        pHkj = H + t_row * j + jj;
                        pQij = pQi + jj;
                        *pQij = 0;
                        for (k = j; k < t_row; k++, pHkj += t_row)
                            *(pQij) += *(pH1ni + k) * *pHkj;
                    }
                }
                memcpy(H1n, m_Q, sizeof(t_type) * t_row * t_row);

                // R = H * (Hn-1*...*H1*A)
                pRi = m_R + t_col * j;
                pHi = H + t_row * j;
                for (i = j; i < t_row; i++, pRi += t_col, pHi += t_row)
                {
                    for (jj = j; jj < t_col; jj++)
                    {
                        pHn1kj = Hn1 + t_col * j + jj;
                        pRij = pRi + jj;
                        *pRij = 0;
                        for (k = j; k < t_row; k++, pHn1kj += t_col)
                            *pRij += *(pHi + k) * *pHn1kj;
                    }
                }
                memcpy(Hn1, m_R, sizeof(t_type) * t_row * t_col);
            }
        }
    }

    m_isOk = 1;
    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t QR<t_row, t_col, t_type>::Compute(const t_type *element, const size_t n_byte)
{
    if ((sizeof(t_type) * t_row * t_col) != n_byte)
    {
        m_isOk = 0;
        return -1;
    }

    memcpy(m_elem, element, n_byte);
    memset(m_Q, 0, sizeof(t_type) * t_row * t_row);
    memset(m_R, 0, sizeof(t_type) * t_row * t_col);

    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t QR<t_row, t_col, t_type>::Compute(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    memset(m_Q, 0, sizeof(t_type) * t_row * t_row);
    memset(m_R, 0, sizeof(t_type) * t_row * t_col);

    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t QR<t_row, t_col, t_type>::Compute(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    memset(m_Q, 0, sizeof(t_type) * t_row * t_row);
    memset(m_R, 0, sizeof(t_type) * t_row * t_col);

    return Compute();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_row, t_type> QR<t_row, t_col, t_type>::GetMatrixQ() const
{
    return Matrix<t_row, t_row, t_type>(m_Q);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> QR<t_row, t_col, t_type>::GetMatrixR() const
{
    return Matrix<t_row, t_col, t_type>(m_R);
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTQR_TPP_