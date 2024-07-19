/*!
\file       dtCscLLT.tpp
\brief      dtMath, Cholesky decomposition(L*L^T form) for Sparse Matrix Class
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       2024. 06. 11
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCSC_LLT_TPP_
#define DTMATH_DTCSC_LLT_TPP_

#include "dtCscLLT.h"

namespace dt
{
namespace Math
{

/** Initializes the  CscLLT object from a  CscMatrix object

* Function Arguments:
*
* order   ->  value '1' for AMD ordering; natural ordering for otherwise (default = 0)
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline CscLLT<m_row, m_col, m_type>::CscLLT(int order)
{
    m_elemNum = 0;
    memset(m_elem, 0, sizeof(m_type) * m_row * m_col);
    memset(m_rowIdx, 0, sizeof(int) * m_row * m_col);
    memset(m_colPtr, 0, sizeof(int) * (m_col + 1));

    if (order)
    {
        AMDFlag = 1;
    }
}

/** Initializes the  CscLLT object from a  CscMatrix object

* Function Arguments:
*
* m   ->   CscMatrix object
* order   ->  value '1' for AMD ordering; natural ordering for otherwise (default = 0)
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline CscLLT<m_row, m_col, m_type>::CscLLT(const CscMatrix<m_row, m_col, m_type> &m, int order)
{
    m_elemNum = m.m_elemNum;

    if (order)
    {
        int cnt = 0;
        Order.AMDOrder(m);

        AMDFlag = 1;
        P = Order.P;
        Pinv = Order.Pinv;

        for (int j = 0; j < m_col; j++)
        {
            m_colPtr[j] = cnt;
            for (int i = m.m_colPtr[P[j]]; i < m.m_colPtr[P[j] + 1]; i++)
            {
                m_rowIdx[cnt] = Pinv[m.m_rowIdx[i]];
                m_elem[cnt++] = m.m_elem[i];
            }
        }
        m_colPtr[m_col] = cnt;
    }
    else
    {
        memcpy(m_elem, m.m_elem, sizeof(m_type) * m_elemNum);
        memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
        memcpy(m_colPtr, m.m_colPtr, sizeof(int) * (m_col + 1));
    }

    etree();
    Compute();
}

/** Initializes the  CscLLT object from a  CscMatrix object

* Function Arguments:
*
* element   -> array containing the non-zero elements of the sparse matrix
* elemNum   -> number of non-zero elements of the sparse matrix
* rowIdx    -> array containing the row indices of the sparse matrix
* colPtr    -> array containing the column pointers of the sparse matrix
*
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline CscLLT<m_row, m_col, m_type>::CscLLT(const m_type *element, const int elemNum, const int *rowIdx, const int *colPtr)
{
    m_elemNum = elemNum;

    memcpy(m_elem, element, sizeof(m_type) * m_elemNum);
    memcpy(m_rowIdx, rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, colPtr, sizeof(int) * (m_col + 1));

    etree();
    Compute();
}

/** Computes the elimination tree for sparse Cholesky LLT factorization

* This function is called in the Compute function to compute the
* elimination tree (etree) of the resulting LT matrix. Though meant to be
* called within the Compute() function, it can be called separately as
* well.

*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline void CscLLT<m_row, m_col, m_type>::etree()
{

    int ancestor[m_col] = {
        0,
    };

    int i = 0;
    int inext = 0;

    for (int k = 0; k < m_col; k++)
    {
        m_parent[k] = -1;
        ancestor[k] = -1;

        for (int p = m_colPtr[k]; p < m_colPtr[k + 1]; p++)
        {
            i = m_rowIdx[p];
            for (; i != -1 && i < k; i = inext) /* traverse from i to k */
            {
                inext = ancestor[i];              /* inext = ancestor of i */
                ancestor[i] = k;                  /* path compression */
                if (inext == -1) m_parent[i] = k; /* no anc., parent is k */
            }
        }
    }
}

/** Computes the reach of the elimination tree for sparse Cholesky factorization

* This function is called in the Compute function to compute the reach of
* the elimination tree (etree) of the resulting LT matrix. The etree() function
* needs to be called first before calling this function. Though meant to be
* called within the Compute() function, it can be called separately as
* well.

* Function Arguments:

* k         -> kth column of the matrix
* parent    -> integer array containing the etree
* s         -> stack containing the pattern of the resulting Sparse LT Matrix
* w         -> matrix of flags for the computation of the etree
* colPtr    -> matrix of column pointers to be updated
* x         -> temporary array to store elements of the kth column of A Matrix

* \return    -> the top most index of the stack
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int CscLLT<m_row, m_col, m_type>::ereach(int k, int *s, int *w, int *colPtr, m_type *x) const
{
    int top = m_col;
    int cnt = 0;

    int i, len;
    w[k] = k;
    for (int p = m_colPtr[k]; p < m_colPtr[k + 1]; p++)
    {
        i = m_rowIdx[p];
        if (i > k) continue;
        x[i] = m_elem[p];
        for (len = 0; w[i] != k; i = m_parent[i])
        {
            s[len++] = i;
            w[i] = k;
            cnt++;
        }

        while (len > 0)
        {
            s[--top] = s[--len];
        }
    }
    colPtr[k + 1] = colPtr[k] + cnt + 1;
    return (top);
}

/** Computes the sparse numerical Cholesky factorization

* This function performs the numerical factorization to give the LT sparse
* matrix.

* \return    -> Zero
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t CscLLT<m_row, m_col, m_type>::Compute()
{
    int top = 0, top_1 = 0;
    int w[m_col] = {
        0,
    };
    int s[m_col] = {
        0,
    };
    int i = 0, j = 0, k = 0, p = 0;
    m_type x[m_row] = {
        0,
    };

    m_type d = 0;

    for (k = 0; k < m_col; k++)
    {
        top = ereach(k, s, w, colPtr_s, x);
        d = x[k];
        top_1 = top;
        for (; top < m_col; top++)
        {
            i = s[top];
            // U^T Solve for Cholesky
            for (j = colPtr_s[i]; j < colPtr_s[i + 1] - 1; j++)
            {
                x[i] -= elem_s[j] * x[rowIdx_s[j]];
            }
            x[i] = x[i] / elem_s[colPtr_s[i + 1] - 1];
            d -= x[i] * x[i];
            rowIdx_s[p] = i;
            elem_s[p++] = x[i];
        }
        // Reinitialization of x vector
        for (j = top_1; j < m_col; j++)
        {
            x[s[j]] = 0;
        }

        rowIdx_s[p] = k;
        if (d <= std::numeric_limits<m_type>::epsilon())
        {
            assert(false && "The matrix is not positive definite");
            return -1;
        }
        elem_s[p++] = sqrt(d);
        x[k] = 0;
    }
    elemNum_s = colPtr_s[m_col];

    return 0;
}

/** Computes the sparse numerical Cholesky factorization

* This function performs the numerical factorization to give the LT sparse
* matrix.
*
* Function Arguments:
*
* element   -> array containing non zero elements of A
* elemNum   -> num. of non zero elements of A
* rowIndex  -> array containing row indices of A
* colPointer-> array containing column pointers of A matrix
* flag      -> flag option for computing etree
*
* \return    -> Compute() function with Zero return
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t CscLLT<m_row, m_col, m_type>::Compute(const m_type *element, const int elemNumber, const int *rowIndex, const int *colPointer, int flag)
{
    m_elemNum = elemNumber;

    memcpy(m_elem, element, sizeof(m_type) * m_elemNum);
    memcpy(m_rowIdx, rowIndex, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, colPointer, sizeof(int) * (m_col + 1));

    if (flag == 0)
    {
        etree();
    }
    return Compute();
}

/** Computes the sparse numerical Cholesky factorization

* This function performs the numerical factorization to give the LT sparse
* matrix.
*
* Function Arguments:
*
* m     -> The input sparse matrix 'A'
* flag      -> flag option for computing etree
*
* \return    -> Compute() function with Zero return
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t CscLLT<m_row, m_col, m_type>::Compute(const CscMatrix<m_row, m_col, m_type> &m, int flag)
{
    m_elemNum = m.m_elemNum;
    int cnt = 0;
    if (AMDFlag)
    {

        if (flag == 0)
        {
            Order.AMDOrder(m);
            P = Order.P;
            Pinv = Order.Pinv;
        }

        for (int j = 0; j < m_col; j++)
        {
            m_colPtr[j] = cnt;
            for (int i = m.m_colPtr[P[j]]; i < m.m_colPtr[P[j] + 1]; i++)
            {
                m_rowIdx[cnt] = Pinv[m.m_rowIdx[i]];
                m_elem[cnt++] = m.m_elem[i];
            }
        }
        m_colPtr[m_col] = cnt;
    }
    else
    {
        memcpy(m_elem, m.m_elem, sizeof(m_type) * m_elemNum);
        memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
        memcpy(m_colPtr, m.m_colPtr, sizeof(int) * (m_col + 1));
    }

    if (flag == 0)
    {
        etree();
    }
    return Compute();
}

/** Gets the lower triangular matrix 'L'

* This function returns the lower triangular matrix 'L' found from the
* Cholesky factorization. This function performs a transpose of the LT
* found from the Compute() function.
*
* \return    -> The lower sparse triangular matrix 'L'
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline CscMatrix<m_row, m_col, m_type> CscLLT<m_row, m_col, m_type>::GetMatrixL() const
{
    m_type elem_tr[elemNum_s];
    memset(elem_tr, 0, sizeof(elem_tr));
    int colIdx[elemNum_s];
    memset(colIdx, 0, sizeof(colIdx));
    int rowPtr[m_row + 1]; // m_row + 1
    memset(rowPtr, 0, sizeof(rowPtr));
    int tmpPtr[m_row + 1];
    memset(tmpPtr, 0, sizeof(tmpPtr));

    if (!elemNum_s) return Matrix<m_col, m_row, m_type>();

    // compute number of non-zero entries per row
    for (uint16_t n = 0; n < elemNum_s; n++)
    {
        rowPtr[rowIdx_s[n]]++;
    }

    // cumsum the elemNum per row to get rowPtr[]
    for (uint16_t i = 0, cumsum = 0, temp; i < m_row; i++)
    {
        temp = rowPtr[i]; // number of non-zero entries per row
        rowPtr[i] = cumsum;
        cumsum += temp;
    }
    rowPtr[m_row] = elemNum_s;
    memcpy(tmpPtr, rowPtr, sizeof(rowPtr));

    // compute column index and data element
    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = colPtr_s[j], ii, jj; i < colPtr_s[j + 1]; i++)
        {
            ii = rowIdx_s[i];
            jj = tmpPtr[ii];

            colIdx[jj] = j;
            elem_tr[jj] = elem_s[i];

            tmpPtr[ii]++;
        }
    }

    return CscMatrix<m_col, m_row, m_type>(elem_tr, elemNum_s, colIdx, rowPtr);
}

/** Gets the upper triangular matrix 'U'

* This function returns the upper triangular matrix 'U' found from the
* Cholesky factorization. This function performs returns the LT found
* from the Compute() function.
*
* \return    -> The upper sparse triangular matrix 'U'
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline CscMatrix<m_col, m_row, m_type> CscLLT<m_row, m_col, m_type>::GetMatrixU() const
{
    return CscMatrix<m_row, m_col, m_type>(elem_s, elemNum_s, rowIdx_s, colPtr_s);
}

/** Computes the x = (L*LT)\b solve with dense vectors

* This function returns x in the Ax = b equation using U matrix
* from the compute() function. The solve function actually performs a
* x = (UT*U)\b solve and returns the dense x matrix.
*
* Function Arguments:
*
* b     -> The input 'b' dense vector from the Ax = b equation
* x     -> The output 'x' dense vector
*
* \return    -> Zero on completion
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline int8_t CscLLT<m_row, m_col, m_type>::Solve(const Vector<m_row, m_type> &b, Vector<m_col, m_type> &x)
{
    if (!m_elemNum) return -1;
    m_type x_perm[m_row];

    // if AMD option is selected
    if (AMDFlag)
    {
        // Permute vector
        for (int k = 0; k < m_col; k++)
        {
            x_perm[Pinv[k]] = b.m_elem[k];
        }
    }
    else
    {
        memcpy(x_perm, b.m_elem, sizeof(b.m_elem));
    }
    for (uint16_t j = 0; j < m_col; j++)
    {
        x_perm[j] = b.m_elem[j];
        for (uint16_t i = colPtr_s[j]; i < colPtr_s[j + 1] - 1; i++)
        {
            x_perm[j] -= elem_s[i] * x_perm[rowIdx_s[i]];
        }
        x_perm[j] /= elem_s[colPtr_s[j + 1] - 1];
    }

    for (int j = m_col - 1; j >= 0; j--)
    {
        x_perm[j] /= elem_s[colPtr_s[j + 1] - 1];
        for (int i = colPtr_s[j]; i < colPtr_s[j + 1] - 1; i++)
        {
            x_perm[rowIdx_s[i]] -= elem_s[i] * x_perm[j];
        }
    }
    // if AMD option is selected
    if (AMDFlag)
    {
        // Inverse permute vector
        for (int k = 0; k < m_col; k++)
        {
            x.m_elem[k] = x_perm[Pinv[k]];
        }
        return 0;
    }

    memcpy(x.m_elem, x_perm, sizeof(x_perm));
    return 0;
}

/** Computes the x = (L*LT)\b solve with dense vectors

* This function returns x in the Ax = b equation using U matrix
* from the compute() function. The solve function actually performs a
* x = (UT*U)\b solve and returns the dense x matrix.
*
* Function Arguments:
*
* b     -> The input 'b' dense vector from the Ax = b equation
*
* \return    -> The output 'x' dense vector
*/

template <uint16_t m_row, uint16_t m_col, typename m_type>
inline Vector<m_row, m_type> CscLLT<m_row, m_col, m_type>::Solve(const Vector<m_row, m_type> &b) const
{
    m_type vec[m_row] = {
        0,
    };
    m_type b_perm[m_row];
    m_type x_perm[m_row];

    if (!m_elemNum) return Vector<m_col, m_type>();

    if (AMDFlag)
    {
        for (int k = 0; k < m_col; k++)
        {
            b_perm[Pinv[k]] = b.m_elem[k];
        }
    }
    else
    {
        for (int k = 0; k < m_col; k++)
        {
            b_perm[k] = b.m_elem[k];
        }
    }

    for (uint16_t j = 0; j < m_col; j++)
    {
        vec[j] = b_perm[j];
        for (uint16_t i = colPtr_s[j]; i < colPtr_s[j + 1] - 1; i++)
        {
            vec[j] -= elem_s[i] * vec[rowIdx_s[i]];
        }
        vec[j] /= elem_s[colPtr_s[j + 1] - 1];
    }

    for (int j = m_col - 1; j >= 0; j--)
    {
        vec[j] /= elem_s[colPtr_s[j + 1] - 1];
        for (int i = colPtr_s[j]; i < colPtr_s[j + 1] - 1; i++)
        {
            vec[rowIdx_s[i]] -= elem_s[i] * vec[j];
        }
    }
    if (AMDFlag)
    {
        for (int k = 0; k < m_col; k++)
        {
            x_perm[k] = vec[Pinv[k]];
        }
        return Vector<m_row, m_type>(x_perm);
    }

    return Vector<m_row, m_type>(vec);
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTCSC_LLT_TPP_
