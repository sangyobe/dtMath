/*!
\file       dtCscMatrix0.tpp
\brief      dtMath, Dynamic Memory Allocation Compressed Sparse Column Matrix (m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 04. 22
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCSC_MATRIX0_TPP_
#define DTMATH_DTCSC_MATRIX0_TPP_

#include "dtCscMatrix0.h"
#include "dtMathMem.h"

#include <cassert>

namespace dt
{
namespace Math
{

// template <typename t_type>
// inline CscMatrix<0, 0, t_type>::CscMatrix(const t_type *element, const int elemNum, const int *rowIdx, const int *colPtr)
// {
//     m_elemNum = elemNum;
//     memcpy(m_elem, element, sizeof(t_type) * m_elemNum);
//     memcpy(m_rowIdx, rowIdx, sizeof(m_rowIdx));
//     memcpy(m_colPtr, colPtr, sizeof(m_colPtr));
// }

template <typename t_type>
inline CscMatrix<0, 0, t_type>::CscMatrix()
    : m_row(0), m_col(0), m_size(0), m_elemNum(0), m_elem(nullptr), m_rowIdx(nullptr), m_colPtr(nullptr), m_compressed(false)
{
}

template <typename t_type>
inline CscMatrix<0, 0, t_type>::CscMatrix(const uint16_t row, const uint16_t col)
    : m_row(row), m_col(col), m_size(row * col), m_elemNum(0), m_compressed(false)
{
    m_elem = MemAlloc<t_type>(m_size);
    m_rowIdx = MemAlloc<int>(m_size);
    m_colPtr = MemAllocZeroInit<int>(m_col + 1);
}

template <typename t_type>
inline CscMatrix<0, 0, t_type>::CscMatrix(const uint16_t row, const uint16_t col, const t_type *element)
    : m_row(row), m_col(col), m_size(row * col), m_compressed(false)
{
    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    m_elem = MemAlloc<t_type>(m_size);
    m_rowIdx = MemAlloc<int>(m_size);
    m_colPtr = MemAlloc<int>(m_col + 1);

    for (icol = 0; icol < m_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < m_row; irow++)
        {
            if (std::abs(element[irow * m_col + icol]) > tolerance)
            {
                m_elem[i] = element[irow * m_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
inline CscMatrix<0, 0, t_type>::CscMatrix(const CscMatrix &m)
    : m_row(m.m_row), m_col(m.m_col), m_size(m.m_size), m_elemNum(m.m_elemNum), m_compressed(m.m_compressed)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    if (m_compressed)
    {
        m_elem = MemAlloc<t_type>(m_elemNum);
        m_rowIdx = MemAlloc<int>(m_elemNum);
        m_colPtr = MemAlloc<int>(m_col + 1);
    }
    else
    {
        m_elem = MemAlloc<t_type>(m_size);
        m_rowIdx = MemAlloc<int>(m_size);
        m_colPtr = MemAlloc<int>(m_col + 1);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(int) * (m_col + 1));
}

template <typename t_type>
inline CscMatrix<0, 0, t_type>::CscMatrix(const Matrix<0, 0, t_type> &m)
    : m_row(m.m_row), m_col(m.m_col), m_size(m.m_size), m_compressed(false)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    m_elem = MemAlloc<t_type>(m_size);
    m_rowIdx = MemAlloc<int>(m_size);
    m_colPtr = MemAlloc<int>(m_col + 1);

    for (icol = 0; icol < m_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < m_row; irow++)
        {
            if (std::abs(m.m_elem[irow * m_col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * m_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline CscMatrix<0, 0, t_type>::CscMatrix(const CscMatrix<row, col, t_type> &m)
    : m_row(row), m_col(col), m_size(row * col), m_elemNum(m.m_elemNum), m_compressed(false)
{
    m_elem = MemAlloc<t_type>(m_size);
    m_rowIdx = MemAlloc<int>(m_size);
    m_colPtr = MemAlloc<int>(col + 1);
    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(int) * (col + 1));
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline CscMatrix<0, 0, t_type>::CscMatrix(const Matrix<row, col, t_type> &m)
    : m_row(row), m_col(col), m_size(row * col), m_compressed(false)
{
    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    m_elem = MemAlloc<t_type>(m_size);
    m_rowIdx = MemAlloc<int>(m_size);
    m_colPtr = MemAlloc<int>(m_col + 1);

    for (icol = 0; icol < col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < row; irow++)
        {
            if (std::abs(m.m_elem[irow * col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
inline CscMatrix<0, 0, t_type>::~CscMatrix()
{
    if (m_elem)
    {
        MemFree<t_type>(m_elem);
        MemFree<int>(m_rowIdx);
        MemFree<int>(m_colPtr);
        m_elem = nullptr;
        m_rowIdx = nullptr;
        m_colPtr = nullptr;
    }
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::NewSize(const uint16_t row, const uint16_t col)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    if (m_elem) Release();

    m_row = row;
    m_col = col;
    m_size = row * col;
    m_elemNum = 0;
    m_compressed = false;

    m_elem = MemAlloc<t_type>(m_size);
    m_rowIdx = MemAlloc<int>(m_size);
    m_colPtr = MemAllocZeroInit<int>(m_col + 1);
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::NewSize(const uint16_t row, const uint16_t col, const t_type *element)
{
    assert(m_elem == nullptr && "Memory has been allocated");
    assert(element != nullptr && "element is nullptr");

    if (m_elem) Release();

    m_row = row;
    m_col = col;
    m_size = row * col;
    m_compressed = false;

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    m_elem = MemAlloc<t_type>(m_size);
    m_rowIdx = MemAlloc<int>(m_size);
    m_colPtr = MemAlloc<int>(m_col + 1);

    for (icol = 0; icol < m_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < m_row; irow++)
        {
            if (std::abs(element[irow * m_col + icol]) > tolerance)
            {
                m_elem[i] = element[irow * m_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::NewSize(const CscMatrix &m)
{
    assert(m_elem == nullptr && "Memory has been allocated");
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    if (m_elem) Release();

    m_row = m.m_row;
    m_col = m.m_col;
    m_size = m.m_size;
    m_elemNum = m.m_elemNum;
    m_compressed = m.m_compressed;

    if (m_compressed)
    {
        m_elem = MemAlloc<t_type>(m_elemNum);
        m_rowIdx = MemAlloc<int>(m_elemNum);
        m_colPtr = MemAlloc<int>(m_col + 1);
    }
    else
    {
        m_elem = MemAlloc<t_type>(m_size);
        m_rowIdx = MemAlloc<int>(m_size);
        m_colPtr = MemAlloc<int>(m_col + 1);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(int) * (m_col + 1));
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::NewSize(const Matrix<0, 0, t_type> &m)
{
    assert(m_elem == nullptr && "Memory has been allocated");
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    if (m_elem) Release();

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    m_row = m.m_row;
    m_col = m.m_col;
    m_size = m.m_size;
    m_compressed = false;

    m_elem = MemAlloc<t_type>(m_size);
    m_rowIdx = MemAlloc<int>(m_size);
    m_colPtr = MemAlloc<int>(m_col + 1);

    for (icol = 0; icol < m_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < m_row; irow++)
        {
            if (std::abs(m.m_elem[irow * m_col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * m_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline void CscMatrix<0, 0, t_type>::NewSize(const CscMatrix<row, col, t_type> &m)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    if (m_elem) Release();

    m_row = row;
    m_col = col;
    m_size = row * col;
    m_elemNum = m.m_elemNum;
    m_compressed = false;

    m_elem = MemAlloc<t_type>(m_size);
    m_rowIdx = MemAlloc<int>(m_size);
    m_colPtr = MemAlloc<int>(col + 1);
    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(int) * (col + 1));
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline void CscMatrix<0, 0, t_type>::NewSize(const Matrix<row, col, t_type> &m)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    if (m_elem) Release();

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    m_row = row;
    m_col = col;
    m_size = row * col;
    m_compressed = false;

    m_elem = MemAlloc<t_type>(m_size);
    m_rowIdx = MemAlloc<int>(m_size);
    m_colPtr = MemAlloc<int>(col + 1);

    for (icol = 0; icol < col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < row; irow++)
        {
            if (std::abs(m.m_elem[irow * col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::ReSize(const uint16_t row, const uint16_t col)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    if (m_row != row || m_col != col || m_compressed)
    {
        m_row = row;
        m_col = col;
        m_size = row * col;
        m_elemNum = 0;
        m_compressed = false;
        m_elem = MemReAlloc<t_type>(m_elem, m_size);
        m_rowIdx = MemReAlloc<int>(m_rowIdx, m_size);
        m_colPtr = MemReAlloc<int>(m_colPtr, m_col + 1);
    }
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::ReSize(const uint16_t row, const uint16_t col, const t_type *element)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(element != nullptr && "element is nullptr");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    if (m_row != row || m_col != col || m_compressed)
    {
        m_row = row;
        m_col = col;
        m_size = row * col;
        m_compressed = false;
        m_elem = MemReAlloc<t_type>(m_elem, m_size);
        m_rowIdx = MemReAlloc<int>(m_rowIdx, m_size);
        m_colPtr = MemReAlloc<int>(m_colPtr, m_col + 1);
    }

    for (icol = 0; icol < m_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < m_row; irow++)
        {
            if (std::abs(element[irow * m_col + icol]) > tolerance)
            {
                m_elem[i] = element[irow * m_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::ReSize(const CscMatrix &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    m_row = m.m_row;
    m_col = m.m_col;
    m_size = m.m_size;
    m_elemNum = m.m_elemNum;
    m_compressed = m.m_compressed;

    if (m_compressed)
    {
        m_elem = MemReAlloc<t_type>(m_elem, m_elemNum);
        m_rowIdx = MemReAlloc<int>(m_rowIdx, m_elemNum);
        m_colPtr = MemReAlloc<int>(m_colPtr, m_col + 1);
    }
    else
    {
        m_elem = MemReAlloc<t_type>(m_elem, m_size);
        m_rowIdx = MemReAlloc<int>(m_rowIdx, m_size);
        m_colPtr = MemReAlloc<int>(m_colPtr, m_col + 1);
    }

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(int) * (m_col + 1));
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::ReSize(const Matrix<0, 0, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    if (m_row != m.m_row || m_col != m.m_col || m_compressed)
    {
        m_row = m.m_row;
        m_col = m.m_col;
        m_size = m.m_size;
        m_compressed = false;
        m_elem = MemReAlloc<t_type>(m_elem, m_size);
        m_rowIdx = MemReAlloc<int>(m_rowIdx, m_size);
        m_colPtr = MemReAlloc<int>(m_colPtr, m_col + 1);
    }

    for (icol = 0; icol < m_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < m_row; irow++)
        {
            if (std::abs(m.m_elem[irow * m_col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * m_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline void CscMatrix<0, 0, t_type>::ReSize(const CscMatrix<row, col, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    if (m_row != row || m_col != col || m_compressed)
    {
        m_row = row;
        m_col = col;
        m_size = row * col;
        m_compressed = false;
        m_elem = MemReAlloc<t_type>(m_elem, m_size);
        m_rowIdx = MemReAlloc<int>(m_rowIdx, m_size);
        m_colPtr = MemReAlloc<int>(m_colPtr, col + 1);
    }

    m_elemNum = m.m_elemNum;

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(int) * (col + 1));
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline void CscMatrix<0, 0, t_type>::ReSize(const Matrix<row, col, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    if (m_row != row || m_col != col || m_compressed)
    {
        m_row = row;
        m_col = col;
        m_size = row * col;
        m_compressed = false;
        m_elem = MemReAlloc<t_type>(m_elem, m_size);
        m_rowIdx = MemReAlloc<int>(m_rowIdx, m_size);
        m_colPtr = MemReAlloc<int>(m_colPtr, col + 1);
    }

    for (icol = 0; icol < col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < row; irow++)
        {
            if (std::abs(m.m_elem[irow * col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::Compress()
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    m_elem = MemReAlloc<t_type>(m_elem, m_elemNum);
    m_rowIdx = MemReAlloc<int>(m_rowIdx, m_elemNum);
    m_compressed = true;
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::Release()
{
    if (m_elem)
    {
        MemFree<t_type>(m_elem);
        MemFree<int>(m_rowIdx);
        MemFree<int>(m_colPtr);
        m_elem = nullptr;
        m_rowIdx = nullptr;
        m_colPtr = nullptr;
        m_row = 0;
        m_col = 0;
        m_size = 0;
        m_elemNum = 0;
        m_compressed = false;
    }
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::SetElement(const t_type *element, const size_t n_byte)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(n_byte == sizeof(t_type) * m_row * m_col && "Check the element size");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    if (m_compressed) ReSize(m_row, m_col);

    for (icol = 0; icol < m_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < m_row; irow++)
        {
            if (std::abs(element[irow * m_col + icol]) > tolerance)
            {
                m_elem[i] = element[irow * m_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::SetElement(const t_type *element, const uint16_t row, const uint16_t col)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(element != nullptr && "element is nullptr");
    assert(m_row == row && "Row dimensions do not matched");
    assert(m_col == col && "Col dimensions do not matched");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    if (m_compressed) ReSize(m_row, m_col);

    for (icol = 0; icol < m_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < m_row; irow++)
        {
            if (std::abs(element[irow * m_col + icol]) > tolerance)
            {
                m_elem[i] = element[irow * m_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::SetElement(const Matrix<0, 0, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == m.m_row && "Check the row size");
    assert(m_col == m.m_col && "Check the col size");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    if (m_compressed) ReSize(m_row, m_col);

    for (icol = 0; icol < m_col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < m_row; irow++)
        {
            if (std::abs(m.m_elem[irow * m_col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * m_col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::SetElement(const Matrix3<t_type, 3, 3> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Check the row size");
    assert(m_col == 3 && "Check the col size");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    if (m_compressed) ReSize(3, 3);

    for (icol = 0; icol < 3; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < 3; irow++)
        {
            if (std::abs(m.m_elem[irow * 3 + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * 3 + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline void CscMatrix<0, 0, t_type>::SetElement(const Matrix<row, col, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Check the row size");
    assert(m_col == col && "Check the col size");

    uint16_t icol;
    uint16_t i = 0;
    t_type tolerance = std::numeric_limits<t_type>::epsilon();

    if (m_compressed) ReSize(row, col);

    for (icol = 0; icol < col; icol++)
    {
        m_colPtr[icol] = i;
        for (uint16_t irow = 0; irow < row; irow++)
        {
            if (std::abs(m.m_elem[irow * col + icol]) > tolerance)
            {
                m_elem[i] = m.m_elem[irow * col + icol];
                m_rowIdx[i++] = irow;
            }
        }
    }

    m_colPtr[icol] = i;
    m_elemNum = i;
}

template <typename t_type>
inline const t_type *const CscMatrix<0, 0, t_type>::GetDataAddr() const
{
    return m_elem;
}

template <typename t_type>
inline const int *const CscMatrix<0, 0, t_type>::GetRowIdx() const
{
    return m_rowIdx;
}

template <typename t_type>
inline const int *const CscMatrix<0, 0, t_type>::GetColPtr() const
{
    return m_colPtr;
}

template <typename t_type>
inline size_t CscMatrix<0, 0, t_type>::GetMemSize() const
{
    return MemSize<t_type>(m_elem) + MemSize<int>(m_rowIdx) + MemSize<int>(m_colPtr);
}

template <typename t_type>
template <uint16_t col>
inline Vector<col, t_type> CscMatrix<0, 0, t_type>::GetRowVec(const uint16_t idxRow) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");
    assert(m_col == col && "Col dimensions do not matched");

    t_type vec[col]{0};

    if (m_elemNum == 0) return Vector<col, t_type>();

    for (uint16_t j = 0; j < m_col; j++)
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

    return Vector<col, t_type>(vec);
}

template <typename t_type>
template <uint16_t row>
inline Vector<row, t_type> CscMatrix<0, 0, t_type>::GetColVec(const uint16_t idxCol) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col > idxCol && "Index out of range");
    assert(m_row == row && "Row dimensions do not matched");

    t_type vec[row]{0};

    if (m_elemNum == 0) return Vector<row, t_type>();

    for (uint16_t i = m_colPtr[idxCol]; i < m_colPtr[idxCol + 1]; i++)
    {
        vec[m_rowIdx[i]] = m_elem[i];
    }

    return Vector<row, t_type>(vec);
}

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::GetRowVec(const uint16_t idxRow, const uint16_t col) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");
    assert(m_col == col && "Col dimensions do not matched");

    Vector<0, t_type> vec(col);

    if (m_elemNum == 0) return vec;

    for (uint16_t j = 0; j < m_col; j++)
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

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::GetColVec(const uint16_t idxCol, const uint16_t row) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col > idxCol && "Index out of range");
    assert(m_row == row && "Row dimensions do not matched");

    Vector<0, t_type> vec(row);

    if (m_elemNum == 0) return vec;

    for (uint16_t i = m_colPtr[idxCol]; i < m_colPtr[idxCol + 1]; i++)
    {
        vec.m_elem[m_rowIdx[i]] = m_elem[i];
    }

    return vec;
}

template <typename t_type>
template <uint16_t col>
inline int8_t CscMatrix<0, 0, t_type>::GetRowVec(const uint16_t idxRow, Vector<col, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");
    assert(m_col == col && "Col dimensions do not matched");

    memset(v.m_elem, 0, sizeof(t_type) * col);

    if (m_elemNum == 0) return 0;

    for (uint16_t j = 0; j < m_col; j++)
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

template <typename t_type>
template <uint16_t row>
inline int8_t CscMatrix<0, 0, t_type>::GetColVec(const uint16_t idxCol, Vector<row, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col > idxCol && "Index out of range");
    assert(m_row == row && "Row dimensions do not matched");

    memset(v.m_elem, 0, sizeof(t_type) * row);

    if (m_elemNum == 0) return 0;

    for (uint16_t i = m_colPtr[idxCol]; i < m_colPtr[idxCol + 1]; i++)
    {
        v.m_elem[m_rowIdx[i]] = m_elem[i];
    }

    return 0;
}

template <typename t_type>
inline int8_t CscMatrix<0, 0, t_type>::GetRowVec(const uint16_t idxRow, Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");
    assert(m_col == v.m_row && "Col dimensions do not matched");

    memset(v.m_elem, 0, sizeof(t_type) * v.m_row);

    if (m_elemNum == 0) return 0;

    for (uint16_t j = 0; j < m_col; j++)
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

template <typename t_type>
inline int8_t CscMatrix<0, 0, t_type>::GetColVec(const uint16_t idxCol, Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(m_col > idxCol && "Index out of range");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    memset(v.m_elem, 0, sizeof(t_type) * v.m_row);

    if (m_elemNum == 0) return 0;

    for (uint16_t i = m_colPtr[idxCol]; i < m_colPtr[idxCol + 1]; i++)
    {
        v.m_elem[m_rowIdx[i]] = m_elem[i];
    }

    return 0;
}

template <typename t_type>
inline Matrix<0, 0, t_type> CscMatrix<0, 0, t_type>::GetDenseMat() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    Matrix<0, 0, t_type> mat(m_row, m_col);

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            mat.m_elem[m_rowIdx[i] * m_col + j] = m_elem[i];
        }
    }

    return mat;
}

template <typename t_type>
inline CscMatrix<0, 0, t_type> CscMatrix<0, 0, t_type>::Transpose() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    // Transpose ==  csc -> csr
    CscMatrix<0, 0, t_type> mat(m_col, m_row);
    int *tmpPtr = MemAlloc<int>(m_row + 1);

    if (m_elemNum == 0) return mat;

    // compute number of non-zero entries per row
    for (uint16_t n = 0; n < m_elemNum; n++)
    {
        mat.m_colPtr[m_rowIdx[n]]++;
    }

    // cumsum the elemNum per row to get rowPtr[]
    for (uint16_t i = 0, cumsum = 0, temp; i < m_row; i++)
    {
        temp = mat.m_colPtr[i]; // number of non-zero entries per row
        mat.m_colPtr[i] = cumsum;
        cumsum += temp;
    }

    mat.m_colPtr[m_row] = m_elemNum;
    mat.m_elemNum = m_elemNum;
    memcpy(tmpPtr, mat.m_colPtr, sizeof(int) * (m_row + 1));

    // compute column index and data element
    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j], jj; i < m_colPtr[j + 1]; i++)
        {
            jj = tmpPtr[m_rowIdx[i]];
            mat.m_rowIdx[jj] = j;
            mat.m_elem[jj] = m_elem[i];
            tmpPtr[m_rowIdx[i]]++;
        }
    }

    MemFree<int>(tmpPtr);

    return mat;
}

template <typename t_type>
inline t_type CscMatrix<0, 0, t_type>::GetNorm() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    t_type sqSum = 0;

    for (uint16_t i = 0; i < m_elemNum; i++)
        sqSum += m_elem[i] * m_elem[i];

    return std::sqrt(sqSum);
}

template <typename t_type>
inline t_type CscMatrix<0, 0, t_type>::GetSqNorm() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    t_type sqSum = 0;

    for (uint16_t i = 0; i < m_elemNum; i++)
        sqSum += m_elem[i] * m_elem[i];

    return sqSum;
}

template <typename t_type>
inline t_type CscMatrix<0, 0, t_type>::GetLpNorm(const int p) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    t_type powSum = 0;

    for (uint16_t i = 0; i < m_elemNum; i++)
        powSum += std::pow(std::abs(m_elem[i]), (t_type)p);

    return std::pow(powSum, (t_type)1 / p);
}

/* Assignment operators */
template <typename t_type>
template <uint16_t row, uint16_t col>
inline CscMatrix<0, 0, t_type> &CscMatrix<0, 0, t_type>::operator=(const CscMatrix<row, col, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");
    assert(m_col == col && "Col dimensions do not matched");

    if (m_compressed) ReSize(row, col);

    m_elemNum = m.m_elemNum;
    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(int) * (col + 1));

    return (*this);
}

template <typename t_type>
inline CscMatrix<0, 0, t_type> &CscMatrix<0, 0, t_type>::operator=(const CscMatrix<0, 0, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == m.m_row && "Row dimensions do not matched");
    assert(m_col == m.m_col && "Col dimensions do not matched");

    m_elemNum = m.m_elemNum;

    if (m.m_compressed == true) Compress();
    else ReSize(m_row, m_col);

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, m.m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(m_colPtr, m.m_colPtr, sizeof(int) * (m_col + 1));

    return (*this);
}

template <typename t_type>
inline CscMatrix<0, 0, t_type> &CscMatrix<0, 0, t_type>::operator*=(const t_type s)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

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

template <typename t_type>
inline CscMatrix<0, 0, t_type> &CscMatrix<0, 0, t_type>::operator/=(const t_type s)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

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
template <typename t_type>
inline CscMatrix<0, 0, t_type> CscMatrix<0, 0, t_type>::operator-() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    CscMatrix<0, 0, t_type> mat(m_row, m_col);

    mat.m_elemNum = m_elemNum;
    memcpy(mat.m_rowIdx, m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(mat.m_colPtr, m_colPtr, sizeof(int) * (m_col + 1));

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        mat.m_elem[i] = -m_elem[i];
        mat.m_elem[i + 2] = -m_elem[i + 2];
        mat.m_elem[i + 1] = -m_elem[i + 1];
        mat.m_elem[i + 3] = -m_elem[i + 3];
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = -m_elem[i];
    }

    return mat;
}

template <typename t_type>
inline CscMatrix<0, 0, t_type> CscMatrix<0, 0, t_type>::operator*(const t_type s) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    CscMatrix<0, 0, t_type> mat(m_row, m_col);

    mat.m_elemNum = m_elemNum;
    memcpy(mat.m_rowIdx, m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(mat.m_colPtr, m_colPtr, sizeof(int) * (m_col + 1));

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        mat.m_elem[i] = m_elem[i] * s;
        mat.m_elem[i + 2] = m_elem[i + 2] * s;
        mat.m_elem[i + 1] = m_elem[i + 1] * s;
        mat.m_elem[i + 3] = m_elem[i + 3] * s;
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = m_elem[i] * s;
    }

    return mat;
}

template <typename t_type>
inline CscMatrix<0, 0, t_type> CscMatrix<0, 0, t_type>::operator/(const t_type s) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    t_type scalar = s;
    CscMatrix<0, 0, t_type> mat(m_row, m_col);

    mat.m_elemNum = m_elemNum;
    memcpy(mat.m_rowIdx, m_rowIdx, sizeof(int) * m_elemNum);
    memcpy(mat.m_colPtr, m_colPtr, sizeof(int) * (m_col + 1));

    if (std::abs(scalar) < std::numeric_limits<t_type>::epsilon())
    {
        if (scalar < 0) scalar = -std::numeric_limits<t_type>::epsilon();
        else scalar = std::numeric_limits<t_type>::epsilon();
    }

    for (cnt = m_elemNum >> 2u; cnt > 0u; cnt--, i += 4)
    {
        mat.m_elem[i] = m_elem[i] / scalar;
        mat.m_elem[i + 2] = m_elem[i + 2] / scalar;
        mat.m_elem[i + 1] = m_elem[i + 1] / scalar;
        mat.m_elem[i + 3] = m_elem[i + 3] / scalar;
    }

    for (cnt = m_elemNum % 4u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = m_elem[i] / scalar;
    }

    return mat;
}

template <typename t_type>
template <uint16_t col>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::operator*(const Vector<col, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == col && "Col dimensions do not matched");

    Vector<0, t_type> vec(m_row);

    if (m_elemNum == 0) return vec;

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec.m_elem[m_rowIdx[i]] += m_elem[i] * v.m_elem[j];
        }
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::operator*(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == v.m_row && "Col dimensions do not matched");

    Vector<0, t_type> vec(m_row);

    if (m_elemNum == 0) return vec;

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec.m_elem[m_rowIdx[i]] += m_elem[i] * v.m_elem[j];
        }
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::operator*(const Vector3<t_type, 3> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == 3 && "Col dimensions do not matched");

    Vector<0, t_type> vec(m_row);

    if (m_elemNum == 0) return vec;

    for (uint16_t j = 0; j < 3; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec.m_elem[m_rowIdx[i]] += m_elem[i] * v.m_elem[j];
        }
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::operator*(const Vector4<t_type, 4> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == 4 && "Col dimensions do not matched");

    Vector<0, t_type> vec(m_row);

    if (m_elemNum == 0) return vec;

    for (uint16_t j = 0; j < 4; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec.m_elem[m_rowIdx[i]] += m_elem[i] * v.m_elem[j];
        }
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::operator*(const Vector6<t_type, 6> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == 6 && "Col dimensions do not matched");

    Vector<0, t_type> vec(m_row);

    if (m_elemNum == 0) return vec;

    for (uint16_t j = 0; j < 6; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec.m_elem[m_rowIdx[i]] += m_elem[i] * v.m_elem[j];
        }
    }

    return vec;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::TposeVec(const Vector<row, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    Vector<0, t_type> vec(m_col);

    if (!m_elemNum) return vec;

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec.m_elem[j] += m_elem[i] * v.m_elem[m_rowIdx[i]];
        }
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::TposeVec(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    Vector<0, t_type> vec(m_col);

    if (!m_elemNum) return vec;

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec.m_elem[j] += m_elem[i] * v.m_elem[m_rowIdx[i]];
        }
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::TposeVec(const Vector3<t_type, 3> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");

    Vector<0, t_type> vec(m_col);

    if (!m_elemNum) return vec;

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec.m_elem[j] += m_elem[i] * v.m_elem[m_rowIdx[i]];
        }
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::TposeVec(const Vector4<t_type, 4> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 4 && "Row dimensions do not matched");

    Vector<0, t_type> vec(m_col);

    if (!m_elemNum) return vec;

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec.m_elem[j] += m_elem[i] * v.m_elem[m_rowIdx[i]];
        }
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::TposeVec(const Vector6<t_type, 6> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 6 && "Row dimensions do not matched");

    Vector<0, t_type> vec(m_col);

    if (!m_elemNum) return vec;

    for (uint16_t j = 0; j < m_col; j++)
    {
        for (uint16_t i = m_colPtr[j]; i < m_colPtr[j + 1]; i++)
        {
            vec.m_elem[j] += m_elem[i] * v.m_elem[m_rowIdx[i]];
        }
    }

    return vec;
}

template <typename t_type>
inline void CscMatrix<0, 0, t_type>::Print(const char endChar)
{
    printf("dat vec =");
    for (uint16_t i = 0; i < m_elemNum; i++)
    {
        printf("%7.2f ", (t_type)(m_elem[i]));
    }
    printf("\ndat num = %d\n", m_elemNum);

    printf("row idx = ");
    for (uint16_t i = 0; i < m_elemNum; i++)
    {
        printf("%d ", m_rowIdx[i]);
    }
    printf("\n");

    printf("col ptr = ");
    for (uint16_t i = 0; i < m_col + 1; i++)
    {
        printf("%d ", m_colPtr[i]);
    }
    printf("\n%c", endChar);
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTCSC_MATRIX0_TPP_
