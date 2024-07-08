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
