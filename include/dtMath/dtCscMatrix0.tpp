/*!
\file       dtCscMatrix0.tpp
\brief      dtMath, Dynamic Memory Allocation Compressed Sparse Column Matrix (m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 06. 10
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

template <typename t_type>
inline CscMatrix<0, 0, t_type>::CscMatrix(const uint16_t row, const uint16_t col, const t_type *element, const int elemNum, const int *rowIdx, const int *colPtr)
    : m_row(row),
      m_col(col), m_size(row * col)
{

    m_elem = MemAlloc<t_type>(elemNum);
    m_rowIdx = MemAlloc<int>(elemNum);
    m_colPtr = MemAlloc<int>(col + 1);

    m_elemNum = elemNum;
    memcpy(m_elem, element, sizeof(t_type) * m_elemNum);
    memcpy(m_rowIdx, rowIdx, sizeof(int) * elemNum);
    memcpy(m_colPtr, colPtr, sizeof(int) * (col + 1));
}

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

/** Returns the Upper Triangular matrix from the Original Matrix in sparse format

* \return    -> The upper triangular CSC matrix from original matrix
*/

template <typename t_type>
inline CscMatrix<0, 0, t_type> CscMatrix<0, 0, t_type>::GetUpperTriangular() const
{
    t_type elem[m_elemNum];
    int elemNum;
    int rowIdx[m_elemNum];
    memset(rowIdx, 0, sizeof(rowIdx));

    int colPtr[m_col + 1];
    memset(colPtr, 0, sizeof(colPtr));

    int i, k = 0;

    if (!m_elemNum)
    {
        assert(false && "The sparse matrix is empty");
        return CscMatrix<0, 0, t_type>(m_row, m_col);
    }

    if (m_col != m_row)
    {
        assert(false && "Row and column sizes are not equal");
        return CscMatrix<0, 0, t_type>(m_row, m_col);
    }

    for (int j = 0; j < m_col; j++)
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
    elemNum = colPtr[m_col];

    return CscMatrix<0, 0, t_type>(m_row, m_col, elem, elemNum, rowIdx, colPtr);
}

/** Returns the Lower Triangular matrix from the Original Matrix in sparse format

* \return    -> The lower triangular CSC matrix from original matrix
*/

template <typename t_type>
inline CscMatrix<0, 0, t_type> CscMatrix<0, 0, t_type>::GetLowerTriangular() const
{
    t_type elem[m_elemNum];
    int elemNum;
    int rowIdx[m_elemNum];
    memset(rowIdx, 0, sizeof(rowIdx));

    int colPtr[m_col + 1];
    memset(colPtr, 0, sizeof(colPtr));
    int i, k = 0;

    if (!m_elemNum)
    {
        assert(false && "The sparse matrix is empty");
        return CscMatrix<0, 0, t_type>(m_row, m_col);
    }

    if (m_col != m_row)
    {
        assert(false && "Row and column sizes are not equal");
        return CscMatrix<0, 0, t_type>(m_row, m_col);
    }

    for (int j = 0; j < m_col; j++)
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
    elemNum = colPtr[m_col];

    return CscMatrix<0, 0, t_type>(m_row, m_col, elem, elemNum, rowIdx, colPtr);
}

/** Compares the patterns of two sparse matrices; two matrix MUST of same size
 *
 * Function Arguments:
 *
 * b         -> The matrix whose pattern is to be compared with
 *
 * \return    -> 0 on success, -1 on failure
 */

template <typename t_type>
inline int CscMatrix<0, 0, t_type>::ComparePattern(const CscMatrix<0, 0, t_type> &b) const
{

    int i, j;

    if (m_elemNum != b.m_elemNum)
    {
        // printf("No. of elements is not equal \n");
        return -1;
    }

    for (j = 0; j < m_col + 1; j++)
    {
        if (m_colPtr[j] != b.m_colPtr[j])
        {
            printf("The column pointers array does not match \n");
            return -1;
        }
    }

    for (j = 0; j < m_elemNum; j++)
    {
        if (m_rowIdx[j] != b.m_rowIdx[j])
        {
            printf("The row indices do not match \n");
            return -1;
        }
    }

    return 0;
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

/** Addition operator between two sparse matrices (A + B)

* \return    -> Sparse matrix after the operation is performed
*/

template <typename t_type>
inline CscMatrix<0, 0, t_type> CscMatrix<0, 0, t_type>::operator+(const CscMatrix<0, 0, t_type> &m) const
{
    t_type x[m_row];
    memset(x, 0, sizeof(x));
    t_type elem[m_row * m.m_col];
    memset(elem, 0, sizeof(elem));
    int rowIdx[m_row * m.m_col];
    memset(rowIdx, 0, sizeof(rowIdx));
    int colPtr[m.m_col + 1];
    memset(colPtr, 0, sizeof(colPtr));
    int elemNum = 0;
    int w[m_row];
    memset(w, 0, sizeof(w));

    if (m_row != m.m_row || m_col != m.m_col)
    {
        assert(false && "A and B matrices are not of same size");
        return CscMatrix<0, 0, t_type>(m_row, m_col);
    }

    for (int i = 0; i < m_col; i++)
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
    return CscMatrix(m_row, m_col, elem, elemNum, rowIdx, colPtr);
}

/** Subtraction operator between two sparse matrices (A - B)

* \return    -> Sparse matrix after the operation is performed
*/

template <typename t_type>
inline CscMatrix<0, 0, t_type> CscMatrix<0, 0, t_type>::operator-(const CscMatrix<0, 0, t_type> &m) const
{
    t_type x[m_row];
    memset(x, 0, sizeof(x));
    t_type elem[m_row * m.m_col];
    memset(elem, 0, sizeof(elem));
    int rowIdx[m_row * m.m_col];
    memset(rowIdx, 0, sizeof(rowIdx));
    int colPtr[m.m_col + 1];
    memset(colPtr, 0, sizeof(colPtr));
    int elemNum = 0;
    int w[m_row];
    memset(w, 0, sizeof(w));

    if (m_row != m.m_row || m_col != m.m_col)
    {
        assert(false && "A and B matrices are not of same size");
        return CscMatrix<0, 0, t_type>(m_row, m_col);
    }

    for (int i = 0; i < m_col; i++)
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
    return CscMatrix(m_row, m_col, elem, elemNum, rowIdx, colPtr);
}

/** Multiplication operator between two sparse matrices (A * B)

* \return    -> Sparse matrix after the operation is performed
*/

template <typename t_type>
inline CscMatrix<0, 0, t_type> CscMatrix<0, 0, t_type>::operator*(const CscMatrix<0, 0, t_type> &m) const
{
    t_type x[m_row];
    memset(x, 0, sizeof(x));
    t_type elem[m_row * m.m_col];
    memset(elem, 0, sizeof(elem));
    int rowIdx[m_row * m.m_col];
    memset(rowIdx, 0, sizeof(rowIdx));
    int colPtr[m.m_col + 1];
    memset(colPtr, 0, sizeof(colPtr));
    int elemNum = 0;
    int w[m_row];
    memset(w, 0, sizeof(w));

    if (m_col != m.m_row)
    {
        assert(false && "Columns of Matrix A and Rows of Matrix B are unequal");
        return CscMatrix<0, 0, t_type>(m_row, m.m_col);
    }

    for (int i = 0; i < m.m_col; i++)
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
    return CscMatrix<0, 0, t_type>(m_row, m.m_col, elem, elemNum, rowIdx, colPtr);
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

/** Lower triangular solve on a dense vector (x = L\b)

* This returns a dense column vector after triangular solve on a dense vector
*
* Function Arguments:
*
* v         -> the dense column vector 'b' before the solve
*
* \return    -> the dense column vector 'x' after the solve
*/

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::LowerTriangularSolve(const Vector<0, t_type> &v)
{
    if (m_row != m_col || v.m_row != m_row)
    {
        assert(false && "Check matrix and vector dimensions (L -> square)");
        return Vector<0, t_type>(m_col);
    }

    if (!m_elemNum)
    {
        assert(false && "The CSC Matrix is empty");
        return Vector<0, t_type>(m_col);
    }

    t_type vec[m_row];

    memcpy(vec, v.m_elem, m_row * sizeof(t_type));

    for (uint16_t j = 0; j < m_row; j++)
    {
        vec[j] = vec[j] / m_elem[m_colPtr[j]];
        for (uint16_t i = m_colPtr[j] + 1; i < m_colPtr[j + 1]; i++)
        {
            vec[m_rowIdx[i]] -= m_elem[i] * vec[j];
        }
    }

    return Vector<0, t_type>(m_col, vec);
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

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::LTTriangularSolve(const Vector<0, t_type> &v)
{
    if (m_row != m_col || v.m_row != m_row)
    {
        assert(false && "Check matrix and vector dimensions (L -> square)");
        return Vector<0, t_type>(m_col);
    }

    if (!m_elemNum)
    {
        assert(false && "The matrix is empty");
        return Vector<0, t_type>(m_col);
    }

    t_type vec[m_row];
    memcpy(vec, v.m_elem, m_row * sizeof(t_type));

    for (int j = m_row - 1; j >= 0; j--)
    {
        for (int i = m_colPtr[j] + 1; i < m_colPtr[j + 1]; i++)
        {
            vec[j] -= m_elem[i] * vec[m_rowIdx[i]];
        }
        vec[j] = vec[j] / m_elem[m_colPtr[j]];
    }
    return Vector<0, t_type>(m_col, vec);
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

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::UpperTriangularSolve(const Vector<0, t_type> &v)
{
    if (m_row != m_col || v.m_row != m_row)
    {
        assert(false && "Check matrix and vector dimensions (L -> square)");
        return Vector<0, t_type>(m_col);
    }

    if (!m_elemNum)
    {
        assert(false && "The matrix is empty");
        return Vector<0, t_type>(m_col);
    }
    t_type vec[m_row];
    memcpy(vec, v.m_elem, m_row * sizeof(t_type));

    for (int j = m_row - 1; j >= 0; j--)
    {
        vec[j] = vec[j] / m_elem[m_colPtr[j + 1] - 1];
        for (int i = m_colPtr[j]; i < m_colPtr[j + 1] - 1; i++)
        {
            vec[m_rowIdx[i]] -= m_elem[i] * vec[j];
        }
    }
    return Vector<0, t_type>(m_row, vec);
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

template <typename t_type>
inline Vector<0, t_type> CscMatrix<0, 0, t_type>::UTTriangularSolve(const Vector<0, t_type> &v)
{
    if (m_row != m_col || v.m_row != m_row)
    {
        assert(false && "Check matrix and vector dimensions (L -> square)");
        return Vector<0, t_type>(m_col);
    }
    if (!m_elemNum)
    {
        assert(false && "The matrix is empty");
        return Vector<0, t_type>(m_col);
    }
    t_type vec[m_row];
    memcpy(vec, v.m_elem, m_row * sizeof(t_type));

    for (int j = 0; j < m_col; j++)
    {
        for (int i = m_colPtr[j]; i < m_colPtr[j + 1] - 1; i++)
        {
            vec[j] -= m_elem[i] * vec[m_rowIdx[i]];
        }
        vec[j] = vec[j] / m_elem[m_colPtr[j + 1] - 1];
    }
    return Vector<0, t_type>(m_col, vec);
}

/** Permute rows of a sparse matrix using a permutation vector

* Function Arguments:
*
* p     -> row permutation vector
*
* \return    -> row permuted sparse matrix
*/

template <typename t_type>
template <uint16_t row>
inline CscMatrix<0, 0, t_type> CscMatrix<0, 0, t_type>::PermuteRow(const Vector<row, int> &p) const
{
    t_type elem[m_elemNum];
    int i;
    int rowIdx[m_elemNum];

    int p_inv[m_row];

    if (!m_elemNum)
    {
        assert(false && "The matrix is empty");
        return CscMatrix<0, 0, t_type>(m_row, m_col);
    }
    if (row != m_row)
    {
        assert(false && "The rows of matrix and rows of vector are unequal");
        return CscMatrix<0, 0, t_type>(m_row, m_col);
    }

    // evaluate P^(-1) or P^T

    for (i = 0; i < m_row; i++)
    {
        p_inv[p.m_elem[i]] = i;
        assert(p.m_elem[i] >= 0 && p.m_elem[i] < m_row && "Check entries of permutation vector");
    }

    for (i = 0; i < m_col; i++)
    {
        for (int j = m_colPtr[i]; j < m_colPtr[i + 1]; j++)
        {
            rowIdx[j] = p_inv[m_rowIdx[j]];
            elem[j] = m_elem[j];
        }
    }
    return CscMatrix<0, 0, t_type>(m_row, m_col, elem, m_elemNum, rowIdx, m_colPtr);
}

/** Permute columns of a sparse matrix using a permutation vector

* Function Arguments:
*
* p     -> column permutation vector
*
* \return    -> column permuted sparse matrix
*/

template <typename t_type>
template <uint16_t col>
inline CscMatrix<0, 0, t_type> CscMatrix<0, 0, t_type>::PermuteCol(const Vector<col, int> &p) const
{

    t_type elem[m_elemNum];
    int rowIdx[m_elemNum];
    int colIdx[m_elemNum];
    int colPtr[m_col + 1];
    colPtr[0] = 0;

    if (!m_elemNum)
    {
        assert(false && "The matrix is empty");
        return CscMatrix<0, 0, t_type>(m_row, m_col);
    }
    if (col != m_col)
    {
        assert(false && "The columns of matrix and rows of vector are unequal");
        return CscMatrix<0, 0, t_type>(m_row, m_col);
    }

    int temp = 0, cnt = 0, cPtr = 0;

    for (int i = 0; i < m_col; i++)
    {
        for (int j = m_colPtr[i]; j < m_colPtr[i + 1]; j++)
        {
            colIdx[j] = i;
        }
    }
    for (int i = 0; i < m_col; i++)
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
    return CscMatrix<0, 0, t_type>(m_row, m_col, elem, m_elemNum, rowIdx, colPtr);
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
