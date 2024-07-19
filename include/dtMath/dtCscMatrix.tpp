/*!
\file       dtCscMatrix.h
\brief      dtMath, Compressed Sparse Column Matrix (m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCSC_MATRIX_TPP_
#define DTMATH_DTCSC_MATRIX_TPP_

#include "dtCscMatrix.h"

#include <cassert>

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type>::CscMatrix(const t_type *element, const int elemNum, const int *rowIdx, const int *colPtr)
{
    m_elemNum = elemNum;
    memcpy(m_elem, element, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, colPtr, sizeof(m_colPtr));
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type>::CscMatrix()
    : m_elemNum(0), m_elem(), m_rowIdx(), m_colPtr()
{
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type>::CscMatrix(const t_type *element, const size_t n_byte)
    : m_elem(), m_rowIdx()
{
    assert(n_byte == sizeof(t_type) * t_row * t_col && "Check the element size");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    for (icol = 0; icol < t_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < t_row; irow++)
        {
            if (std::abs(element[irow * t_col + icol]) > tolerance)
            {
                m_elem[i] = element[irow * t_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type>::CscMatrix(const CscMatrix &m)
    : m_elem(), m_rowIdx()
{
    m_elemNum = m.m_elemNum;
    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(m_colPtr));
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type>::CscMatrix(const CscMatrix<0, 0, t_type> &m)
    : m_elem(), m_rowIdx()
{
    assert(t_row == m.m_row && "Row dimensions do not matched");
    assert(t_col == m.m_col && "Col dimensions do not matched");
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    m_elemNum = m.m_elemNum;
    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(m_colPtr));
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type>::CscMatrix(const Matrix<t_row, t_col, t_type> &m)
    : m_elem(), m_rowIdx()
{
    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    for (icol = 0; icol < t_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < t_row; irow++)
        {
            if (std::abs(m.m_elem[irow * t_col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * t_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type>::CscMatrix(const Matrix<0, 0, t_type> &m)
    : m_elem(), m_rowIdx()
{
    assert(t_row == m.m_row && "Row dimensions do not matched");
    assert(t_col == m.m_col && "Col dimensions do not matched");
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    for (icol = 0; icol < t_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < t_row; irow++)
        {
            if (std::abs(m.m_elem[irow * t_col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * t_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void CscMatrix<t_row, t_col, t_type>::SetElement(const t_type *element, const size_t n_byte)
{
    assert(n_byte == sizeof(t_type) * t_row * t_col && "Check the element size");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    for (icol = 0; icol < t_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < t_row; irow++)
        {
            if (std::abs(element[irow * t_col + icol]) > tolerance)
            {
                m_elem[i] = element[irow * t_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void CscMatrix<t_row, t_col, t_type>::SetElement(const Matrix<t_row, t_col, t_type> &m)
{
    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    for (icol = 0; icol < t_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < t_row; irow++)
        {
            if (std::abs(m.m_elem[irow * t_col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * t_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void CscMatrix<t_row, t_col, t_type>::SetElement(const Matrix3<t_type, t_row, t_col> &m)
{
    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    for (icol = 0; icol < t_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < t_row; irow++)
        {
            if (std::abs(m.m_elem[irow * t_col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * t_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void CscMatrix<t_row, t_col, t_type>::SetElement(const Matrix<0, 0, t_type> &m)
{
    assert(t_row == m.m_row && "Row dimensions do not matched");
    assert(t_col == m.m_col && "Col dimensions do not matched");
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    for (icol = 0; icol < t_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < t_row; irow++)
        {
            if (std::abs(m.m_elem[irow * t_col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * t_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline const t_type *const CscMatrix<t_row, t_col, t_type>::GetDataAddr() const
{
    return m_elem;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline const int *const CscMatrix<t_row, t_col, t_type>::GetRowIdx() const
{
    return m_rowIdx;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline const int *const CscMatrix<t_row, t_col, t_type>::GetColPtr() const
{
    return m_colPtr;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> CscMatrix<t_row, t_col, t_type>::GetRowVec(const uint16_t idxRow) const
{
    assert(t_row > idxRow && "Index out of range");

    t_type vec[t_col]{0};

    if (m_elemNum == 0) return Vector<t_col, t_type>();

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            if (m_rowIdx[i] == idxRow)
            {
                vec[j] = m_elem[i];
                break;
            }
        }
    }

    return Vector<t_col, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> CscMatrix<t_row, t_col, t_type>::GetColVec(const uint16_t idxCol) const
{
    assert(t_col > idxCol && "Index out of range");

    t_type vec[t_row]{0};

    if (m_elemNum == 0) return Vector<t_row, t_type>();

    for (uint16_t i = m_colPtr[idxCol]; i < m_colPtr[idxCol + 1]; i++)
    {
        vec[m_rowIdx[i]] = m_elem[i];
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<0, t_type> CscMatrix<t_row, t_col, t_type>::GetRowVec(const uint16_t idxRow, const uint16_t col) const
{
    assert(t_row > idxRow && "Index out of range");

    Vector<0, t_type> vec(t_col);

    if (m_elemNum == 0) return vec;

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            if (m_rowIdx[i] == idxRow)
            {
                vec.m_elem[j] = m_elem[i];
                break;
            }
        }
    }

    return vec;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<0, t_type> CscMatrix<t_row, t_col, t_type>::GetColVec(const uint16_t idxCol, const uint16_t row) const
{
    assert(t_col > idxCol && "Index out of range");

    Vector<0, t_type> vec(t_row);

    if (m_elemNum == 0) return vec;

    for (uint16_t i = m_colPtr[idxCol]; i < m_colPtr[idxCol + 1]; i++)
    {
        vec.m_elem[m_rowIdx[i]] = m_elem[i];
    }

    return vec;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t CscMatrix<t_row, t_col, t_type>::GetRowVec(const uint16_t idxRow, Vector<t_col, t_type> &v) const
{
    assert(t_row > idxRow && "Index out of range");

    memset(v.m_elem, 0, sizeof(t_type) * t_col);

    if (m_elemNum == 0) return 0;

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            if (m_rowIdx[i] == idxRow)
            {
                v.m_elem[j] = m_elem[i];
                break;
            }
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t CscMatrix<t_row, t_col, t_type>::GetColVec(const uint16_t idxCol, Vector<t_row, t_type> &v) const
{
    assert(t_col > idxCol && "Index out of range");

    memset(v.m_elem, 0, sizeof(t_type) * t_row);

    if (m_elemNum == 0) return 0;

    for (uint16_t i = m_colPtr[idxCol]; i < m_colPtr[idxCol + 1]; i++)
    {
        v.m_elem[m_rowIdx[i]] = m_elem[i];
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t CscMatrix<t_row, t_col, t_type>::GetRowVec(const uint16_t idxRow, Vector<0, t_type> &v) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col == v.m_row && "Check the dimensions");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    memset(v.m_elem, 0, sizeof(t_type) * t_col);

    if (m_elemNum == 0) return 0;

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            if (m_rowIdx[i] == idxRow)
            {
                v.m_elem[j] = m_elem[i];
                break;
            }
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t CscMatrix<t_row, t_col, t_type>::GetColVec(const uint16_t idxCol, Vector<0, t_type> &v) const
{
    assert(t_col > idxCol && "Index out of range");
    assert(t_row == v.m_row && "Check the dimensions");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    memset(v.m_elem, 0, sizeof(t_type) * t_row);

    if (m_elemNum == 0) return 0;

    for (uint16_t i = m_colPtr[idxCol]; i < m_colPtr[idxCol + 1]; i++)
    {
        v.m_elem[m_rowIdx[i]] = m_elem[i];
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> CscMatrix<t_row, t_col, t_type>::GetDenseMat() const
{
    t_type mat[t_row * t_col]{0};

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            mat[m_rowIdx[i] * t_col + j] = m_elem[i];
        }
    }

    return Matrix<t_row, t_col, t_type>(mat);
}

/** Returns the Upper Triangular matrix from the Original Matrix in sparse format

* \return    -> The upper triangular CSC matrix from original matrix
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> CscMatrix<t_row, t_col, t_type>::GetUpperTriangular() const
{
    t_type elem[m_elemNum];
    int elemNum;
    int rowIdx[m_elemNum] = {
        0,
    };
    int colPtr[t_col + 1] = {
        0,
    };
    int i, k = 0;

    if (!m_elemNum)
    {
        assert(false && "The sparse matrix is empty");
        return CscMatrix<t_row, t_col, t_type>();
    }

    if (t_col != t_row)
    {
        assert(false && "Row and column sizes are not equal");
        return CscMatrix<t_row, t_col, t_type>();
    }

    for (int j = 0; j < t_col; j++)
    {
        for (int p = m_colPtr[j]; p < m_colPtr[j + 1]; p++)
        {
            i = m_rowIdx[p];
            if (i > j) continue;
            rowIdx[k] = i;
            elem[k++] = m_elem[p];
        }
        colPtr[j + 1] = k;
    }
    elemNum = colPtr[t_col];

    return CscMatrix<t_row, t_col, t_type>(elem, elemNum, rowIdx, colPtr);
}

/** Returns the Lower Triangular matrix from the Original Matrix in sparse format

* \return    -> The lower triangular CSC matrix from original matrix
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> CscMatrix<t_row, t_col, t_type>::GetLowerTriangular() const
{

    t_type elem[m_elemNum];
    int elemNum;
    int rowIdx[m_elemNum] = {
        0,
    };
    int colPtr[t_col + 1] = {
        0,
    }; // t_row + 1
    int i, k = 0;

    if (!m_elemNum)
    {
        assert(false && "The sparse matrix is empty");
        return CscMatrix<t_row, t_col, t_type>();
    }

    if (t_col != t_row)
    {
        assert(false && "Row and column sizes are not equal");
        return CscMatrix<t_row, t_col, t_type>();
    }

    for (int j = 0; j < t_col; j++)
    {
        for (int p = m_colPtr[j]; p < m_colPtr[j + 1]; p++)
        {
            i = m_rowIdx[p];
            if (i < j) continue;
            rowIdx[k] = i;
            elem[k++] = m_elem[p];
        }
        colPtr[j + 1] = k;
    }
    elemNum = colPtr[t_col];

    return CscMatrix<t_row, t_col, t_type>(elem, elemNum, rowIdx, colPtr);
}

/** Compares the patterns of two sparse matrices; two matrix MUST of same size
 *
 * Function Arguments:
 *
 * b         -> The matrix whose pattern is to be compared with
 *
 * \return    -> 0 on success, -1 on failure
 */

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int CscMatrix<t_row, t_col, t_type>::ComparePattern(const CscMatrix<t_row, t_col, t_type> &b) const
{

    int i, j;

    if (m_elemNum != b.m_elemNum)
    {
        // printf("The matrix pattern does not match 1 \n");
        return -1;
    }

    for (j = 1; j < t_col + 1; j++)
    {
        if (m_colPtr[j] != b.m_colPtr[j])
        {
            // printf("The matrix pattern does not match 2 \n");
            return -1;
        }
    }

    for (j = 0; j < t_col + 1; j++)
    {
        for (i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            if (m_rowIdx[i] != b.m_rowIdx[i])
            {
                // printf("The matrix pattern does not match 3 \n");
                return -1;
            }
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_col, t_row, t_type> CscMatrix<t_row, t_col, t_type>::Transpose() const
{
    // Transpose ==  csc -> csr
    t_type elem[t_row * t_col];
    int colIdx[t_row * t_col]{0};
    int rowPtr[t_row + 1]{0}; // t_row + 1
    int tmpPtr[t_row + 1]{0};

    if (m_elemNum == 0) return Matrix<t_col, t_row, t_type>();

    // compute number of non-zero entries per row
    for (uint16_t n = 0; n < m_elemNum; n++)
    {
        rowPtr[m_rowIdx[n]]++;
    }

    // cumsum the elemNum per row to get rowPtr[]
    for (uint16_t i = 0, cumsum = 0, temp; i < t_row; i++)
    {
        temp = rowPtr[i]; // number of non-zero entries per row
        rowPtr[i] = cumsum;
        cumsum += temp;
    }

    rowPtr[t_row] = m_elemNum;
    memcpy(tmpPtr, rowPtr, sizeof(rowPtr));

    // compute column index and data element
    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j], jj; i < m_colPtr[j + 1]; i++)
        {
            jj = tmpPtr[m_rowIdx[i]];
            colIdx[jj] = j;
            elem[jj] = m_elem[i];
            tmpPtr[m_rowIdx[i]]++;
        }
    }

    return CscMatrix<t_col, t_row, t_type>(elem, m_elemNum, colIdx, rowPtr);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type CscMatrix<t_row, t_col, t_type>::GetNorm() const
{
    t_type sqSum = 0;

    for (uint16_t i = 0; i < m_elemNum; i++)
        sqSum += m_elem[i] * m_elem[i];

    return std::sqrt(sqSum);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type CscMatrix<t_row, t_col, t_type>::GetSqNorm() const
{
    t_type sqSum = 0;

    for (uint16_t i = 0; i < m_elemNum; i++)
        sqSum += m_elem[i] * m_elem[i];

    return sqSum;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type CscMatrix<t_row, t_col, t_type>::GetLpNorm(const int p) const
{
    t_type powSum = 0;

    for (uint16_t i = 0; i < m_elemNum; i++)
        powSum += std::pow(std::abs(m_elem[i]), (t_type)p);

    return std::pow(powSum, (t_type)1 / p);
}

/* Assignment operators */
template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> &CscMatrix<t_row, t_col, t_type>::operator=(const CscMatrix &m)
{
    m_elemNum = m.m_elemNum;
    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(m_colPtr));

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> &CscMatrix<t_row, t_col, t_type>::operator=(const CscMatrix<0, 0, t_type> &m)
{
    assert(t_row == m.m_row && "Row dimensions do not matched");
    assert(t_col == m.m_col && "Col dimensions do not matched");
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    m_elemNum = m.m_elemNum;
    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(m_colPtr));

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> &CscMatrix<t_row, t_col, t_type>::operator*=(const t_type s)
{
    uint16_t cnt, i = 0;

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        m_elem[i] *= s;
        m_elem[i + 2] *= s;
        m_elem[i + 1] *= s;
        m_elem[i + 3] *= s;
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] *= s;
    }

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> &CscMatrix<t_row, t_col, t_type>::operator/=(const t_type s)
{
    t_type scalar = s;
    uint16_t cnt, i = 0;

    if (std::abs(scalar) < std::numeric_limits<t_type>::epsilon())
    {
        if (scalar < 0) scalar = -std::numeric_limits<t_type>::epsilon();
        else scalar = std::numeric_limits<t_type>::epsilon();
    }

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        m_elem[i] /= scalar;
        m_elem[i + 2] /= scalar;
        m_elem[i + 1] /= scalar;
        m_elem[i + 3] /= scalar;
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] /= scalar;
    }

    return (*this);
}

/* Arithmetic operators */
template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> CscMatrix<t_row, t_col, t_type>::operator-() const
{
    t_type elem[m_elemNum];
    uint16_t cnt, i = 0;

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        elem[i] = -m_elem[i];
        elem[i + 2] = -m_elem[i + 2];
        elem[i + 1] = -m_elem[i + 1];
        elem[i + 3] = -m_elem[i + 3];
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        elem[i] = -m_elem[i];
    }

    return CscMatrix(elem, m_elemNum, m_rowIdx, m_colPtr);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> CscMatrix<t_row, t_col, t_type>::operator*(const t_type s) const
{
    t_type elem[m_elemNum];
    uint16_t cnt, i = 0;

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        elem[i] = m_elem[i] * s;
        elem[i + 2] = m_elem[i + 2] * s;
        elem[i + 1] = m_elem[i + 1] * s;
        elem[i + 3] = m_elem[i + 3] * s;
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        elem[i] = m_elem[i] * s;
    }

    return CscMatrix(elem, m_elemNum, m_rowIdx, m_colPtr);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> CscMatrix<t_row, t_col, t_type>::operator/(const t_type s) const
{
    t_type elem[m_elemNum];
    uint16_t cnt, i = 0;
    t_type scalar = s;

    if (std::abs(scalar) < std::numeric_limits<t_type>::epsilon())
    {
        if (scalar < 0) scalar = -std::numeric_limits<t_type>::epsilon();
        else scalar = std::numeric_limits<t_type>::epsilon();
    }

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        elem[i] = m_elem[i] / scalar;
        elem[i + 2] = m_elem[i + 2] / scalar;
        elem[i + 1] = m_elem[i + 1] / scalar;
        elem[i + 3] = m_elem[i + 3] / scalar;
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        elem[i] = m_elem[i] / scalar;
    }

    return CscMatrix(elem, m_elemNum, m_rowIdx, m_colPtr);
}

/** Addition operator between two sparse matrices (A + B)

* \return    -> Sparse matrix after the operation is performed
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> CscMatrix<t_row, t_col, t_type>::operator+(const CscMatrix<t_row, t_col, t_type> &m) const
{
    t_type x[t_row] = {
        0,
    };
    t_type elem[m_elemNum + m.m_elemNum] = {
        0,
    };
    int rowIdx[m_elemNum + m.m_elemNum] = {
        0,
    };
    int colPtr[t_col + 1] = {
        0,
    };
    int elemNum = 0;
    int w[t_row] = {
        0,
    };
    /*
        if (t_row != m.t_row || t_col != m.t_col)
        {
            assert(false && "A and B matrices are not of same size");
            return CscMatrix<t_row, t_col, t_type>();
        }
    */
    for (int i = 0; i < t_col; i++)
    {
        // Scatter Function
        for (uint16_t j = m_colPtr[i]; j < m_colPtr[i + 1]; j++)
        {
            x[m_rowIdx[j]] = m_elem[j];

            w[m_rowIdx[j]] = i + 1;
            rowIdx[elemNum] = m_rowIdx[j];
            elemNum++;
        }
        for (uint16_t j = m.m_colPtr[i]; j < m.m_colPtr[i + 1]; j++)
        {
            if (w[m.m_rowIdx[j]] < i + 1)
            {
                w[m.m_rowIdx[j]] = i + 1;
                rowIdx[elemNum] = m.m_rowIdx[j];
                elemNum++;
            }
            x[m.m_rowIdx[j]] += m.m_elem[j];
        }
        // Gather Function
        for (uint16_t k = colPtr[i]; k < elemNum; k++)
        {
            elem[k] = x[rowIdx[k]];
            x[rowIdx[k]] = 0;
        }
        colPtr[i + 1] = elemNum;
    }
    return CscMatrix(elem, elemNum, rowIdx, colPtr);
}

/** Subtraction operator between two sparse matrices (A - B)

* \return    -> Sparse matrix after the operation is performed
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> CscMatrix<t_row, t_col, t_type>::operator-(const CscMatrix<t_row, t_col, t_type> &m) const
{
    t_type x[t_row] = {
        0,
    };
    t_type elem[m_elemNum + m.m_elemNum] = {
        0,
    };
    int rowIdx[m_elemNum + m.m_elemNum] = {
        0,
    };
    int colPtr[t_col + 1] = {
        0,
    };
    int elemNum = 0;
    int w[t_row] = {
        0,
    };
    /*
    if (t_row != m.t_row || t_col != m.t_col)
    {
        assert(false && "A and B matrices are not of same size");
        return CscMatrix<t_row, t_col, t_type>();
    }
    */
    for (int i = 0; i < t_col; i++)
    {
        for (uint16_t j = m_colPtr[i]; j < m_colPtr[i + 1]; j++)
        {
            // Scatter Function
            x[m_rowIdx[j]] = m_elem[j];

            w[m_rowIdx[j]] = i + 1;
            rowIdx[elemNum] = m_rowIdx[j];
            elemNum++;
        }
        for (uint16_t j = m.m_colPtr[i]; j < m.m_colPtr[i + 1]; j++)
        {
            if (w[m.m_rowIdx[j]] < i + 1)
            {
                w[m.m_rowIdx[j]] = i + 1;
                rowIdx[elemNum] = m.m_rowIdx[j];
                elemNum++;
            }
            x[m.m_rowIdx[j]] -= m.m_elem[j];
        }
        // Gather Function
        for (uint16_t k = colPtr[i]; k < elemNum; k++)
        {
            elem[k] = x[rowIdx[k]];
            x[rowIdx[k]] = 0;
        }
        colPtr[i + 1] = elemNum;
    }
    return CscMatrix(elem, elemNum, rowIdx, colPtr);
}

/** Multiplication operator between two sparse matrices (A * B)

* \return    -> Sparse matrix after the operation is performed
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t col>
inline CscMatrix<t_row, col, t_type> CscMatrix<t_row, t_col, t_type>::operator*(const CscMatrix<t_col, col, t_type> &m) const
{
    t_type x[t_row] = {
        0,
    };
    t_type elem[t_row * t_col] = {
        0,
    };
    int rowIdx[t_row * t_col] = {
        0,
    };
    int colPtr[col + 1] = {
        0,
    };
    int elemNum = 0;
    int w[t_row] = {
        0,
    };

    for (int i = 0; i < col; i++)
    {
        // Scatter Function
        for (uint16_t j = m.m_colPtr[i]; j < m.m_colPtr[i + 1]; j++)
        {
            for (uint16_t k = m_colPtr[m.m_rowIdx[j]]; k < m_colPtr[m.m_rowIdx[j] + 1]; k++)
            {
                if (w[m_rowIdx[k]] < i + 1)
                {
                    w[m_rowIdx[k]] = i + 1;
                    rowIdx[elemNum] = m_rowIdx[k];
                    elemNum++;
                }
                x[m_rowIdx[k]] += m_elem[k] * m.m_elem[j];
                // printf("%7.2f\t", x[m_rowIdx[k]]);
            }
        }
        // Gather Function
        for (uint16_t k = colPtr[i]; k < elemNum; k++)
        {
            elem[k] = x[rowIdx[k]];
            x[rowIdx[k]] = 0;
        }
        colPtr[i + 1] = elemNum;
    }
    return CscMatrix<t_row, col, t_type>(elem, elemNum, rowIdx, colPtr);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> CscMatrix<t_row, t_col, t_type>::operator*(const Vector<t_col, t_type> &v) const
{
    t_type vec[t_row]{0};

    if (!m_elemNum) return Vector<t_row, t_type>();

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[m_rowIdx[i]] += m_elem[i] * v.m_elem[j];
        }
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> CscMatrix<t_row, t_col, t_type>::operator*(const Vector<0, t_type> &v) const
{
    assert(t_col == v.m_row && "Check the dimensions");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    t_type vec[t_row]{0};

    if (!m_elemNum) return Vector<t_row, t_type>();

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[m_rowIdx[i]] += m_elem[i] * v.m_elem[j];
        }
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> CscMatrix<t_row, t_col, t_type>::operator*(const Vector3<t_type, t_col> &v) const
{
    t_type vec[t_row]{0};

    if (!m_elemNum) return Vector<t_row, t_type>();

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[m_rowIdx[i]] += m_elem[i] * v.m_elem[j];
        }
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> CscMatrix<t_row, t_col, t_type>::operator*(const Vector4<t_type, t_col> &v) const
{
    t_type vec[t_row]{0};

    if (!m_elemNum) return Vector<t_row, t_type>();

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[m_rowIdx[i]] += m_elem[i] * v.m_elem[j];
        }
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> CscMatrix<t_row, t_col, t_type>::operator*(const Vector6<t_type, t_col> &v) const
{
    t_type vec[t_row]{0};

    if (!m_elemNum) return Vector<t_row, t_type>();

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[m_rowIdx[i]] += m_elem[i] * v.m_elem[j];
        }
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> CscMatrix<t_row, t_col, t_type>::TposeVec(const Vector<t_row, t_type> &v) const
{
    t_type vec[t_col]{0};

    if (!m_elemNum) return Vector<t_col, t_type>();

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[j] += m_elem[i] * v.m_elem[m_rowIdx[i]];
        }
    }

    return Vector<t_col, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> CscMatrix<t_row, t_col, t_type>::TposeVec(const Vector<0, t_type> &v) const
{
    assert(t_row == v.m_row && "Check the dimensions");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    t_type vec[t_col]{0};

    if (!m_elemNum) return Vector<t_col, t_type>();

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[j] += m_elem[i] * v.m_elem[m_rowIdx[i]];
        }
    }

    return Vector<t_col, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> CscMatrix<t_row, t_col, t_type>::TposeVec(const Vector3<t_type, t_row> &v) const
{
    t_type vec[t_col]{0};

    if (!m_elemNum) return Vector<t_col, t_type>();

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[j] += m_elem[i] * v.m_elem[m_rowIdx[i]];
        }
    }

    return Vector<t_col, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> CscMatrix<t_row, t_col, t_type>::TposeVec(const Vector4<t_type, t_row> &v) const
{
    t_type vec[t_col]{0};

    if (!m_elemNum) return Vector<t_col, t_type>();

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[j] += m_elem[i] * v.m_elem[m_rowIdx[i]];
        }
    }

    return Vector<t_col, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_col, t_type> CscMatrix<t_row, t_col, t_type>::TposeVec(const Vector6<t_type, t_row> &v) const
{
    t_type vec[t_col]{0};

    if (!m_elemNum) return Vector<t_col, t_type>();

    for (uint16_t j = 0; j < t_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec[j] += m_elem[i] * v.m_elem[m_rowIdx[i]];
        }
    }

    return Vector<t_col, t_type>(vec);
}

/** Lower triangular solve on a dense vector (x = L\b)

* This returns a dense column vector after triangular solve on a dense vector
*
* Function Arguments:
*
* v         -> the dense column vector 'b' before the solve
*
* \return    -> the dense column vector 'x' after the solve
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> CscMatrix<t_row, t_col, t_type>::LowerTriangularSolve(const Vector<t_row, t_type> &v) const
{
    t_type vec[t_row] = {
        0,
    };
    memcpy(vec, v.m_elem, sizeof(v.m_elem));

    if (!m_elemNum)
    {
        assert(false && "The matrix is empty");
        return Vector<t_col, t_type>();
    }

    for (uint16_t j = 0; j < t_row; j++)
    {
        vec[j] = vec[j] / m_elem[m_colPtr[j]];
        for (uint16_t i = m_colPtr[j] + 1; i < m_colPtr[j + 1]; i++)
        {
            vec[m_rowIdx[i]] -= m_elem[i] * vec[j];
        }
    }

    return Vector<t_col, t_type>(vec);
}

/** LT triangular solve on a dense vector (x = L^T\b)

* This returns a dense column vector after triangular solve on a dense vector
*
* Function Arguments:
*
* v         -> the dense column vector 'b' before the solve
*
* \return    -> the dense column vector 'x' after the solve
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> CscMatrix<t_row, t_col, t_type>::LTTriangularSolve(const Vector<t_row, t_type> &v) const
{
    t_type vec[t_row] = {
        0,
    };
    memcpy(vec, v.m_elem, sizeof(v.m_elem));

    if (!m_elemNum)
    {
        assert(false && "The matrix is empty");
        return Vector<t_col, t_type>();
    }

    for (int j = t_row - 1; j >= 0; j--)
    {
        for (int i = m_colPtr[j] + 1; i < m_colPtr[j + 1]; i++)
        {
            vec[j] -= m_elem[i] * vec[m_rowIdx[i]];
        }
        vec[j] = vec[j] / m_elem[m_colPtr[j]];
    }
    return Vector<t_col, t_type>(vec);
}

/** Upper triangular solve on a dense vector (x = U\b)

* This returns a dense column vector after triangular solve on a dense vector
*
* Function Arguments:
*
* v         -> the dense column vector 'b' before the solve
*
* \return    -> the dense column vector 'x' after the solve
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> CscMatrix<t_row, t_col, t_type>::UpperTriangularSolve(const Vector<t_row, t_type> &v) const
{
    t_type vec[t_row] = {
        0,
    };
    memcpy(vec, v.m_elem, sizeof(v.m_elem));

    if (!m_elemNum)
    {
        assert(false && "The matrix is empty");
        return Vector<t_col, t_type>();
    }

    for (int j = t_row - 1; j >= 0; j--)
    {
        vec[j] = vec[j] / m_elem[m_colPtr[j + 1] - 1];
        for (int i = m_colPtr[j]; i < m_colPtr[j + 1] - 1; i++)
        {
            vec[m_rowIdx[i]] -= m_elem[i] * vec[j];
        }
    }
    return Vector<t_col, t_type>(vec);
}

/** UT triangular solve on a dense vector (x = U^T\b)

* This returns a dense column vector after triangular solve on a dense vector
*
* Function Arguments:
*
* v         -> the dense column vector 'b' before the solve
*
* \return    -> the dense column vector 'x' after the solve
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> CscMatrix<t_row, t_col, t_type>::UTTriangularSolve(const Vector<t_row, t_type> &v) const
{
    t_type vec[t_row] = {
        0,
    };
    memcpy(vec, v.m_elem, sizeof(v.m_elem));

    if (!m_elemNum)
    {
        assert(false && "The matrix is empty");
        return Vector<t_col, t_type>();
    }

    for (int j = 0; j < t_col; j++)
    {
        for (int i = m_colPtr[j]; i < m_colPtr[j + 1] - 1; i++)
        {
            vec[j] -= m_elem[i] * vec[m_rowIdx[i]];
        }
        vec[j] = vec[j] / m_elem[m_colPtr[j + 1] - 1];
    }
    return Vector<t_col, t_type>(vec);
}

/** Permute rows of a sparse matrix using a permutation vector

* Function Arguments:
*
* p     -> row permutation vector
*
* \return    -> row permuted sparse matrix
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> CscMatrix<t_row, t_col, t_type>::PermuteRow(const Vector<t_row, int> &p) const
{
    t_type elem[t_row * t_col];
    int i;
    int rowIdx[t_row * t_col] = {
        0,
    };

    int p_inv[t_row] = {
        0,
    };

    if (!m_elemNum)
    {
        assert(false && "The matrix is empty");
        return CscMatrix<t_row, t_col, t_type>();
    }

    // evaluate P^(-1) or P^T

    for (i = 0; i < t_row; i++)
    {
        p_inv[p.m_elem[i]] = i;
        assert(p.m_elem[i] >= 0 && p.m_elem[i] < t_row && "Check entries of permutation vector");
    }

    for (i = 0; i < t_col; i++)
    {
        for (int j = m_colPtr[i]; j < m_colPtr[i + 1]; j++)
        {
            rowIdx[j] = p_inv[m_rowIdx[j]];
            elem[j] = m_elem[j];
        }
    }
    return CscMatrix<t_row, t_col, t_type>(elem, m_elemNum, rowIdx, m_colPtr);
}

/** Permute columns of a sparse matrix using a permutation vector

* Function Arguments:
*
* p     -> column permutation vector
*
* \return    -> column permuted sparse matrix
*/

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> CscMatrix<t_row, t_col, t_type>::PermuteCol(const Vector<t_col, int> &p) const
{
    t_type elem[t_row * t_col];
    int rowIdx[t_row * t_col] = {
        0,
    };
    int colIdx[t_row * t_col] = {
        0,
    };
    int colPtr[t_col + 1] = {
        0,
    };

    if (!m_elemNum)
    {
        assert(false && "The matrix is empty");
        return CscMatrix<t_row, t_col, t_type>();
    }

    int temp = 0, cnt = 0, cPtr = 0;

    for (int i = 0; i < t_col; i++)
    {
        for (int j = m_colPtr[i]; j < m_colPtr[i + 1]; j++)
        {
            colIdx[j] = i;
        }
    }
    for (int i = 0; i < t_col; i++)
    {
        temp = colIdx[m_colPtr[i]];
        cPtr = p.m_elem[temp];
        for (int j = m_colPtr[cPtr]; j < m_colPtr[cPtr + 1]; j++)
        {
            elem[cnt] = m_elem[j];
            rowIdx[cnt] = m_rowIdx[j];
            cnt++;
        }
        colPtr[i + 1] = cnt;
    }
    return CscMatrix<t_row, t_col, t_type>(elem, m_elemNum, rowIdx, colPtr);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void CscMatrix<t_row, t_col, t_type>::Print(const char endChar)
{
    printf("dat vec =");
    for (uint16_t i = 0; i < m_elemNum; i++)
    {
        printf("%7.2f ", (t_type)(m_elem[i]));
    }
    printf("\nvec num = %d\n", m_elemNum);

    printf("row idx = ");
    for (uint16_t i = 0; i < m_elemNum; i++)
    {
        printf("%d ", m_rowIdx[i]);
    }
    printf("\n");

    printf("col ptr = ");
    for (uint16_t i = 0; i < t_col + 1; i++)
    {
        printf("%d ", m_colPtr[i]);
    }
    printf("\n%c", endChar);
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTCSC_MATRIX_TPP_
