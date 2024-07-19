#pragma once

#ifndef DTMATH_DTFULL_PIV_LU_TPP_
#define DTMATH_DTFULL_PIV_LU_TPP_

#include "dtFullPivLU.h"

namespace dt
{
namespace Math
{

/** Initializes the FullPivLU object
 */

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline FullPivLU<t_row, t_col, t_type>::FullPivLU()
{
    memset(m_elem, 0, sizeof(t_type) * t_row * t_col);
    m_isOk = 0;
}

/** Initializes the FullPivLU object

* Function Arguments:
*
* element   -> vector (array) containing elements of the matrix
* n_byte    -> size of the matrix
*
* \return    -> 0 on success, -1 on failure
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline FullPivLU<t_row, t_col, t_type>::FullPivLU(const t_type *element, const size_t n_byte)
{
    if ((sizeof(t_type) * t_row * t_col) != n_byte) m_isOk = 0;
    else
    {
        memcpy(m_elem, element, n_byte);
        Compute();
    }
}

/** Initializes the CdtFullPivLU object

* Function Arguments:
*
* m   ->  dt::Matrix object
*
* \return    -> 0 on success, -1 on failure
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline FullPivLU<t_row, t_col, t_type>::FullPivLU(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

/** Initializes the CdtFullPivLU object

* Function Arguments:
*
* m   ->  dt::Matrix3 object
*
* \return    -> 0 on success, -1 on failure
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline FullPivLU<t_row, t_col, t_type>::FullPivLU(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    Compute();
}

/** Compute the PAQ = LU factorization with full pivoting

* \return    -> 0 on success, -1 on failure
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t FullPivLU<t_row, t_col, t_type>::Compute()
{
    if (t_row != t_col)
    {
        m_isOk = 0;
        return -1;
    }

    int i, j, k;
    t_type *pMi, *pMk, *p_pivotRow = nullptr;
    t_type max_row, max_col, absElem;
    t_type pivotRow[t_row];
    t_type temp_clm;
    int tempPerm;

    for (i = 0; i < t_row; i++) // filling the Q permutation vector as identity
        m_permCol[i] = i;

    for (i = 0, pMi = m_elem; i < t_row; pMi += t_row, i++)
    {
        // Pivoting /
        // find the pivot row
        m_pivot[i] = i;
        max_row = std::abs(*(pMi + i));

        for (j = i, pMk = pMi; j < t_row; j++, pMk += t_row)
        {
            for (k = i; k < t_row; k++) // Loop for checking maximum value in each element of row j
            {
                if (max_row < (absElem = std::abs(*(pMk + k))))
                {
                    max_row = absElem;
                    m_pivot[i] = j;
                    p_pivotRow = pMk; // pMk is pivot row
                }
            }
        }

        if (m_pivot[i] != i)
        {
            memcpy(pivotRow, p_pivotRow, sizeof(t_type) * t_col);
            memcpy(p_pivotRow, pMi, sizeof(t_type) * t_col);
            memcpy(pMi, pivotRow, sizeof(t_type) * t_col);
        }

        m_pivotCol[i] = i;                                    // Stores the column to be swapped temporarily as the ith column
        max_col = std::abs(*(pMi + i));                       // Temporarily stores the value of diagonal "i x i" as maximum value
        for (k = i, pMk = pMi + i; k < t_col - 1; k++, pMk++) // Loop for finding for highest value in the "i"th row starting from i x i value
        {
            if (max_col < (absElem = std::abs(*(pMk + 1)))) // Check for maximum value for each element in row "i" (Rows already arranged in previous step)
            {
                max_col = absElem; // Replace maximum value in "max_col"
                m_pivotCol[i] = k + 1;
            }
        }

        if (m_pivotCol[i] != i)
        {

            for (k = 0, pMk = m_elem; k < t_col; k++, pMk += t_col)
            {

                temp_clm = *(pMk + i);
                *(pMk + i) = *(pMk + m_pivotCol[i]); // These three lines swap the "k"th column with the column with the next maximum value
                *(pMk + m_pivotCol[i]) = temp_clm;
            }
            tempPerm = m_permCol[i];
            m_permCol[i] = m_permCol[m_pivotCol[i]]; // These three lines swap the "i"th row of the permutation vector
            m_permCol[m_pivotCol[i]] = tempPerm;
        }

        // matrix is singular, return error
        if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon())
        {
            m_isOk = 0;
            return -1;
        }

        // LU Decompostion using Gaussian Elimination //
        for (k = i + 1, pMk = pMi + t_col; k < t_row; pMk += t_col, k++)
        {
            // find the lower triangular matrix elements for column i.
            *(pMk + i) /= *(pMi + i);

            // update the upper triangular matrix for remaining matrix
            for (j = i + 1; j < t_col; j++)
                *(pMk + j) -= *(pMk + i) * *(pMi + j);
        }
    }

    m_isOk = 1;
    return 0;
}

/** Compute the LU factorization with full pivoting

* Function Arguments:
*
* element   -> vector (array) containing elements of the matrix
* n_byte    -> size of the matrix
*
* \return    -> 0 on success, -1 on failure
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t FullPivLU<t_row, t_col, t_type>::Compute(const t_type *element, const size_t n_byte)
{
    if ((sizeof(t_type) * t_row * t_col) != n_byte)
    {
        m_isOk = 0;
        return -1;
    }

    memcpy(m_elem, element, n_byte);
    return Compute();
}

/** Compute the LU factorization with full pivoting

* Function Arguments:
*
* m   ->  dt::Matrix object
*
* \return    -> 0 on success, -1 on failure
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t FullPivLU<t_row, t_col, t_type>::Compute(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    return Compute();
}

/** Compute the LU factorization with full pivoting

* Function Arguments:
*
* m   ->  dt::Matrix3 object
*
* \return    -> 0 on success, -1 on failure
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t FullPivLU<t_row, t_col, t_type>::Compute(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    return Compute();
}

/** Get the complete matrix after factorization with with full pivoting

* Get the complete solved matrix after the PAQ = LU factorization
*
* \return    -> LU matrix elements in a full matrix
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> FullPivLU<t_row, t_col, t_type>::GetMatrix() const
{
    return Matrix<t_row, t_col, t_type>(m_elem);
}

/** Get the L matrix after factorization with full pivoting

* Get the lower triangular matrix after the PAQ = LU factorization with full pivoting
*
* \return    -> L matrix in the PAQ = LU equation
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> FullPivLU<t_row, t_col, t_type>::GetMatrixL() const
{
    int i, j;
    t_type *pMi;
    t_type L[t_row * t_col]{0};

    for (i = 0, pMi = L; i < t_row; pMi += t_row, i++)
        *(pMi + i) = 1;

    // Update remaining matrix from m_elem to L//

    for (i = 1, pMi = L + t_row; i < t_row; pMi += t_row, i++)
    {
        for (j = 0; j < i; j++)
            *(pMi + j) = m_elem[i * t_col + j];
    }

    return Matrix<t_row, t_col, t_type>(L);
}

/** Get the U matrix after factorization with full pivoting

* Get the upper triangular matrix after the PAQ = LU factorization with full pivoting
*
* \return    -> U matrix in the PAQ = LU equation
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> FullPivLU<t_row, t_col, t_type>::GetMatrixU() const
{
    int i, j;
    t_type U[t_row * t_col]{0};
    t_type *pMi;

    for (i = 0, pMi = U; i < t_row; pMi += t_row, i++)
    {
        for (j = i; j < t_col; j++)
            *(pMi + j) = m_elem[i * t_col + j];
    }

    return Matrix<t_row, t_col, t_type>(U);
}

/** Get the P matrix

* Get the row permutation matrix 'P' in the PAQ = LU equation
*
* \return    -> P matrix in the PAQ = LU equation
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> FullPivLU<t_row, t_col, t_type>::GetMatrixP() const
{
    int i;
    t_type P[t_row * t_col]{0};
    t_type row[t_col]{0};
    t_type *pMi;

    for (i = 0, pMi = P; i < t_row; pMi += t_row, i++)
        *(pMi + i) = 1;

    for (i = 0; i < t_row; i++)
    {
        if (m_pivot[i] != i)
        {
            memcpy(row, &P[i * t_col], sizeof(t_type) * t_col);
            memcpy(&P[i * t_col], &P[m_pivot[i] * t_col], sizeof(t_type) * t_col);
            memcpy(&P[m_pivot[i] * t_col], row, sizeof(t_type) * t_col);
        }
    }

    return Matrix<t_row, t_col, t_type>(P);
}

/** Get the Q matrix

* Get the column permutation matrix 'Q' in the PAQ = LU equation
*
* \return    -> Q matrix in the PAQ = LU equation
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> FullPivLU<t_row, t_col, t_type>::GetMatrixQ() const
{
    int i, j;
    t_type Q[t_row * t_col]{0};
    t_type temp_clm;
    t_type *pMk;

    for (i = 0, pMk = Q; i < t_row; pMk += t_row, i++)
        *(pMk + i) = 1;

    for (i = 0; i < t_row; i++)
    {
        if (m_pivotCol[i] != i)
        {
            for (j = 0, pMk = Q; j < t_col; j++, pMk += t_col)
            {

                temp_clm = *(pMk + i);
                *(pMk + i) = *(pMk + m_pivotCol[i]); // These three lines swap the "j"th column with the m_pivotCol[i]th column
                *(pMk + m_pivotCol[i]) = temp_clm;
            }
        }
    }
    return Matrix<t_row, t_col, t_type>(Q);
}

/** Perform the P*A*Q*x = L*U*x = b operation

* Function Arguments:
*
* b ->  The 'b' vector in the L*U*x = b equation
* x ->  The 'x' vector in the L*U*x = b equation
*
* \return    -> 0 on success, -1 on failure
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t FullPivLU<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, Vector<t_col, t_type> &x)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    t_type *pMi;
    t_type tmp;
    t_type vb[t_row];
    t_type x_out[t_col]{0};

    memcpy(vb, b.m_elem, sizeof(vb));

    if (!m_isOk) return -1;

    /* Solve Ly = b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
    {
        if (m_pivot[i] != i)
        {
            tmp = vb[i];
            vb[i] = vb[m_pivot[i]];
            vb[m_pivot[i]] = tmp;
        }

        x_out[i] = vb[i];
        for (k = 0; k < i; k++)
            x_out[i] -= *(pMi + k) * x_out[k];
    }

    /* Solve Ux = y */
    // interchange the row of vector b along the original position
    // Solve the upper triangular (backward substitution)
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < t_col; k++)
            x_out[i] -= *(pMi + k) * x_out[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x_out[i] /= *(pMi + i);
        x.m_elem[m_permCol[i]] = x_out[i];
    }

    return 0;
}

/** Perform the P*A*Q*x = L*U*x = b operation

* Function Arguments:
*
* b     ->  The 'b' vector in the L*U*x = b equation
* isOk  ->  the Ok flag to indicate that the process completed properly
*
* \return    -> The 'x' vector in the L*U*x = b
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> FullPivLU<t_row, t_col, t_type>::Solve(const Vector<t_row, t_type> &b, int8_t *isOk)
{
    // Solve, Ax = LUx = b
    // where L is a lower triangular matrix with an all diagonal element is 1
    //       U is upper triangular matrix
    // define Ux = y
    // Ly = b is solved by forward substitution for y
    // Ux = y is solved by backward substitution for x

    int i, k;
    t_type *pMi;
    t_type x[t_col] = {
        0,
    },
           x_out[t_col];
    t_type tmp;
    t_type vb[t_row];

    memcpy(vb, b.m_elem, sizeof(vb));

    if (isOk) *isOk = 1;

    if (!m_isOk && isOk)
    {
        *isOk = 0;
        return Vector<t_col, t_type>();
    }

    /* Solve Ly =b */
    // interchange the row of vector b with the pivot order
    // Solve the unit lower triangular matrix for y (forward substitution), here x is y
    for (i = 0, pMi = m_elem; i < t_row; pMi += t_col, i++)
    {
        if (m_pivot[i] != i)
        {
            tmp = vb[i];
            vb[i] = vb[m_pivot[i]]; //  x = P * vb
            vb[m_pivot[i]] = tmp;
        }

        x[i] = vb[i];

        for (k = 0; k < i; k++)
            x[i] -= *(pMi + k) * x[k];
    }

    /* Solve Ux = y */
    // interchange the row of vector b along the original position
    // Solve the upper triangular (backward substitution)
    for (i = t_row - 1, pMi = m_elem + (t_row - 1) * t_col; i >= 0; i--, pMi -= t_col)
    {
        // if (m_pivot[i] != i)
        //{
        //     tmp = b.m_elem[i];
        //     b.m_elem[i] = b.m_elem[m_pivot[i]];
        //     b.m_elem[m_pivot[i]] = tmp;
        // }

        for (k = i + 1; k < t_col; k++)
            x[i] -= *(pMi + k) * x[k];

        // if (std::abs(*(pMi + i)) <= std::numeric_limits<t_type>::epsilon()) return -1;

        x[i] /= *(pMi + i);
        x_out[m_permCol[i]] = x[i]; //  x_out = Q * x
    }

    return Vector<t_col, t_type>(x_out);
}

/** Inverse of the A matrix
 *
 * Function Arguments:
 *
 * inv    -> the output inverse matrix in Matrix class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t FullPivLU<t_row, t_col, t_type>::Inverse(Matrix<t_row, t_col, t_type> &inv)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    uint16_t colIdx[t_col]; // column index for interchange the current col with the pivot col
    uint16_t tmpColIdx;
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};
    t_type QinvPAQ[t_row * t_col]{0};
    t_type *p_pivotRow;

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < t_row; i++)
    {
        invL[i * (t_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + t_col;
    p_invMi = invL + t_col;
    for (i = 1; i < t_row; i++, p_Mi += t_col, p_invMi += t_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += t_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + t_col;
            for (k = j + 1; k < i; k++, p_invMk += t_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < t_row; k++, p_Mi += (t_col + 1), p_invMk += (t_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + t_col * (t_row - 2);
    p_invMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, p_Mi -= t_col, p_invMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + t_col;
            for (k = i + 1; k <= j; k++, p_invMk += t_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }
    // int perm[t_col];
    // for (i = 0; i < t_row; i++) colIdx[i] = i;

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < t_row; i++)
    {
        for (j = 0; j < t_col; j++)
        {
            // if (m_pivot[j] != j)
            //{
            //     colIdx[j] = colIdx[m_pivot[j]];
            //     colIdx[m_pivot[j]] = j;
            // }
            if (m_pivot[j] != j)
            {

                tmpColIdx = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tmpColIdx;
            }

            for (k = 0; k < t_col; k++)
                inv.m_elem[i * t_col + colIdx[j]] += invU[i * t_col + k] * invL[k * t_col + j];

            colIdx[j] = j;
        }
    }

    for (i = 0, p_Mi = inv.m_elem; i < t_row; i++, p_Mi += t_col)
    {
        p_pivotRow = QinvPAQ + (m_permCol[i] * t_col);
        memcpy(p_pivotRow, p_Mi, sizeof(t_type) * t_col);
        // printf("m_pivot[%d] = %d, colIdx[%d] = %d \n", i, m_pivot[i], i, colIdx[i]);
    }

    memcpy(inv.m_elem, QinvPAQ, sizeof(t_type) * t_col * t_row);

    return 0;
}

/** Inverse of the A matrix
 *
 * Function Arguments:
 *
 * inv    -> the output inverse matrix in Matrix3 class
 *
 * \return    -> 0 on success, -1 on failure
 */
template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t FullPivLU<t_row, t_col, t_type>::Inverse(Matrix3<t_type, t_row, t_col> &inv)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        return -1;
    }

    int i, j, k;
    uint16_t colIdx[t_col]; // column index for interchange the current col with the pivot col
    uint16_t tmpColIdx;
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};
    t_type QinvPAQ[t_row * t_col]{0};
    t_type *p_pivotRow;

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < t_row; i++)
    {
        invL[i * (t_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + t_col;
    p_invMi = invL + t_col;
    for (i = 1; i < t_row; i++, p_Mi += t_col, p_invMi += t_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += t_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + t_col;
            for (k = j + 1; k < i; k++, p_invMk += t_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < t_row; k++, p_Mi += (t_col + 1), p_invMk += (t_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + t_col * (t_row - 2);
    p_invMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, p_Mi -= t_col, p_invMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + t_col;
            for (k = i + 1; k <= j; k++, p_invMk += t_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }
    // int perm[t_col];
    // for (i = 0; i < t_row; i++) colIdx[i] = i;

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < t_row; i++)
    {
        for (j = 0; j < t_col; j++)
        {
            // if (m_pivot[j] != j)
            //{
            //     colIdx[j] = colIdx[m_pivot[j]];
            //     colIdx[m_pivot[j]] = j;
            // }
            if (m_pivot[j] != j)
            {

                tmpColIdx = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tmpColIdx;
            }

            for (k = 0; k < t_col; k++)
                inv.m_elem[i * t_col + colIdx[j]] += invU[i * t_col + k] * invL[k * t_col + j];

            colIdx[j] = j;
        }
    }

    for (i = 0, p_Mi = inv.m_elem; i < t_row; i++, p_Mi += t_col)
    {
        p_pivotRow = QinvPAQ + (m_permCol[i] * t_col);
        memcpy(p_pivotRow, p_Mi, sizeof(t_type) * t_col);

        // printf("m_pivot[%d] = %d, colIdx[%d] = %d \n", i, m_pivot[i], i, colIdx[i]);
    }

    memcpy(inv.m_elem, QinvPAQ, sizeof(t_type) * t_col * t_row);

    return 0;
}

/** Inverse of the A matrix
 *
 * Function Arguments:
 *
 * isOk  ->  the Ok flag to indicate that the process completed properly
 *
 * \return    -> the inverse matrix in Matrix class
 */
template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> FullPivLU<t_row, t_col, t_type>::Inverse(int8_t *isOk)
{
    if (!m_isOk)
    {
        assert(m_isOk && "The matrix is not decomposed into LU");
        if (isOk) *isOk = 0;
        return Matrix<t_row, t_col, t_type>();
    }

    int i, j, k;
    uint16_t colIdx[t_col]; // column index for interchange the current col with the pivot col
    uint16_t tmpColIdx;
    t_type *p_Mi;
    t_type *p_invMi, *p_invMj, *p_invMk;
    t_type sum;
    t_type invL[t_row * t_col]{0};
    t_type invU[t_row * t_col]{0};
    t_type QinvPAQ[t_row * t_col]{0};
    t_type m_inv[t_row * t_col]{0};
    t_type *p_pivotRow;

    /* Initialization */
    // Set the diagonal elements of the lower triangular matrix as "1"
    // Initialize permutation vector.
    for (i = 0; i < t_row; i++)
    {
        invL[i * (t_col + 1)] = 1;
        colIdx[i] = i;
    }

    /* Inverse of Lower triangular matrix */
    // Invert the subdiagonal part of the matrix L row by row where
    // the diagonal elements are assumed to be 1.
    p_Mi = m_elem + t_col;
    p_invMi = invL + t_col;
    for (i = 1; i < t_row; i++, p_Mi += t_col, p_invMi += t_col)
    {
        p_invMj = invL;
        for (j = 0; j < i; j++, p_invMj += t_col)
        {
            *(p_invMi + j) = -*(p_Mi + j);

            p_invMk = p_invMj + t_col;
            for (k = j + 1; k < i; k++, p_invMk += t_col)
            {
                *(p_invMi + j) -= *(p_Mi + k) * *(p_invMk + j);
            }
        }
    }

    /* Inverse of Upper triangular matrix */
    // Invert the diagonal elements of the upper triangular matrix U.
    p_Mi = m_elem;
    p_invMk = invU;
    for (k = 0; k < t_row; k++, p_Mi += (t_col + 1), p_invMk += (t_col + 1))
    {
        // if (std::abs(*p_Mi) <= std::numeric_limits<t_type>::epsilon()) return -1;
        // else *p_invMk = 1 / *p_Mi;
        *p_invMk = 1 / *p_Mi;
    }

    // Invert the remaining upper triangular matrix U.
    p_Mi = m_elem + t_col * (t_row - 2);
    p_invMi = invU + t_col * (t_row - 2);
    for (i = t_row - 2; i >= 0; i--, p_Mi -= t_col, p_invMi -= t_col)
    {
        for (j = t_col - 1; j > i; j--)
        {
            sum = 0;
            p_invMk = p_invMi + t_col;
            for (k = i + 1; k <= j; k++, p_invMk += t_col)
            {
                sum += *(p_Mi + k) * *(p_invMk + j);
            }
            *(p_invMi + j) = -*(p_invMi + i) * sum;
        }
    }
    // int perm[t_col];
    // for (i = 0; i < t_row; i++) colIdx[i] = i;

    /* Inv(A) = inv(U) * inv(L) * P */
    for (i = 0; i < t_row; i++)
    {
        for (j = 0; j < t_col; j++)
        {
            // if (m_pivot[j] != j)
            //{
            //     colIdx[j] = colIdx[m_pivot[j]];
            //     colIdx[m_pivot[j]] = j;
            // }
            if (m_pivot[j] != j)
            {

                tmpColIdx = colIdx[j];
                colIdx[j] = colIdx[m_pivot[j]]; // These three lines swap the "i"th row of the permutation vector
                colIdx[m_pivot[j]] = tmpColIdx;
            }

            for (k = 0; k < t_col; k++)
                m_inv[i * t_col + colIdx[j]] += invU[i * t_col + k] * invL[k * t_col + j];

            colIdx[j] = j;
        }
    }

    for (i = 0, p_Mi = m_inv; i < t_row; i++, p_Mi += t_col)
    {
        p_pivotRow = QinvPAQ + (m_permCol[i] * t_col);
        memcpy(p_pivotRow, p_Mi, sizeof(t_type) * t_col);

        // printf("m_pivot[%d] = %d, colIdx[%d] = %d \n", i, m_pivot[i], i, colIdx[i]);
    }

    memcpy(m_inv, QinvPAQ, sizeof(t_type) * t_col * t_row);

    if (isOk) *isOk = 1;
    return Matrix<t_row, t_col, t_type>(m_inv);
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTNO_PIV_LU_TPP_
