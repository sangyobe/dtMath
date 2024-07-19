/*!
\file       dtMatrix0.tpp
\brief      dtMath, Dynamic Memory Alloction General Matrix(m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 04. 01
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATRIX0_TPP_
#define DTMATH_DTMATRIX0_TPP_

#include "dtMathMem.h"
#include "dtMatrix0.h"

#include <cassert>

namespace dt
{
namespace Math
{

template <typename t_type>
inline Matrix<0, 0, t_type>::Matrix() : m_row(0), m_col(0), m_size(0), m_elem(nullptr)
{
}

template <typename t_type>
inline Matrix<0, 0, t_type>::Matrix(const uint16_t row, const uint16_t col)
    : m_row(row), m_col(col), m_size(row * col)
{
    m_elem = MemAllocZeroInit<t_type>(m_size);
}

template <typename t_type>
inline Matrix<0, 0, t_type>::Matrix(const uint16_t row, const uint16_t col, const t_type *element)
    : m_row(row), m_col(col), m_size(row * col)
{
    m_elem = MemAlloc<t_type>(m_size);
    memcpy(m_elem, element, sizeof(t_type) * m_size);
}

template <typename t_type>
inline Matrix<0, 0, t_type>::Matrix(const Matrix<0, 0, t_type> &m)
{
    m_row = m.m_row;
    m_col = m.m_col;
    m_size = m_row * m_col;
    m_elem = MemAlloc<t_type>(m_size);
    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Matrix<0, 0, t_type>::Matrix(const Matrix<row, col, t_type> &m)
{
    m_row = row;
    m_col = col;
    m_size = row * col;
    m_elem = MemAlloc<t_type>(row * col);
    memcpy(m_elem, m.m_elem, sizeof(t_type) * row * col);
}

template <typename t_type>
inline Matrix<0, 0, t_type>::~Matrix()
{
    if (m_elem)
    {
        MemFree<t_type>(m_elem);
        m_elem = nullptr;
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::NewSize(const uint16_t row, const uint16_t col)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    if (m_elem) Release();

    m_row = row;
    m_col = col;
    m_size = row * col;
    m_elem = MemAllocZeroInit<t_type>(m_size);
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::NewSize(const uint16_t row, const uint16_t col, const t_type *element)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    if (m_elem) Release();

    m_row = row;
    m_col = col;
    m_size = row * col;
    m_elem = MemAlloc<t_type>(m_size);
    memcpy(m_elem, element, sizeof(t_type) * m_size);
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::NewSize(const Matrix<0, 0, t_type> &m)
{
    assert(m_elem == nullptr && "Memory has been allocated");
    assert(m.m_elem != nullptr && "Argument Matrix has not been allocated");

    if (m_elem) Release();

    m_row = m.m_row;
    m_col = m.m_col;
    m_size = m_row * m_col;
    m_elem = MemAlloc<t_type>(m_size);
    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline void Matrix<0, 0, t_type>::NewSize(const Matrix<row, col, t_type> &m)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    if (m_elem) Release();

    m_row = row;
    m_col = col;
    m_size = row * col;
    m_elem = MemAlloc<t_type>(row * col);
    memcpy(m_elem, m.m_elem, sizeof(t_type) * row * col);
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::ReSize(const uint16_t row, const uint16_t col)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    m_row = row;
    m_col = col;
    m_size = row * col;
    m_elem = MemReAlloc<t_type>(m_elem, m_size);
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::ReSize(const uint16_t row, const uint16_t col, const t_type *element)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(element != nullptr && "Argument is empty");

    m_row = row;
    m_col = col;
    m_size = row * col;
    m_elem = MemReAlloc<t_type>(m_elem, m_size);
    memcpy(m_elem, element, sizeof(t_type) * m_size);
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::ReSize(const Matrix<0, 0, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_elem != nullptr && "Argument Matrix has not been allocated");

    m_row = m.m_row;
    m_col = m.m_col;
    m_size = m_row * m_col;
    m_elem = MemReAlloc<t_type>(m_elem, m_size);
    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline void Matrix<0, 0, t_type>::ReSize(const Matrix<row, col, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    m_row = row;
    m_col = col;
    m_size = row * col;
    m_elem = MemReAlloc<t_type>(m_elem, row * col);
    memcpy(m_elem, m.m_elem, sizeof(t_type) * row * col);
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::Release()
{
    if (m_elem)
    {
        MemFree<t_type>(m_elem);
        m_elem = nullptr;
        m_row = 0;
        m_col = 0;
        m_size = 0;
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetZero()
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    memset(m_elem, 0, sizeof(t_type) * m_size);
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetIdentity()
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    memset(m_elem, 0, sizeof(t_type) * m_size);
    uint16_t num = (m_row > m_col) ? m_col : m_row;
    uint16_t offset = m_col + 1;

    for (uint16_t i = 0; i < num; i++)
        m_elem[i * offset] = 1;
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetDiagonal(const t_type *element, const size_t n_byte)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t num = (m_row > m_col) ? m_col : m_row;
    uint16_t offset = m_col + 1;
    uint16_t elemNum = (uint16_t)(n_byte / sizeof(t_type));

    num = (num > elemNum) ? elemNum : num;

    for (uint16_t i = 0; i < num; i++)
        m_elem[i * offset] = element[i];
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetFill(const t_type value)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 2u; cnt > 0u; cnt--, i += 4)
    {
        m_elem[i] = value;
        m_elem[i + 1] = value;
        m_elem[i + 2] = value;
        m_elem[i + 3] = value;
    }

    for (cnt = m_size % 4u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] = value;
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetFillRow(const uint16_t idxRow, const t_type value)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");

    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = m_col >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        m_elem[idxRow * m_col + icol] = value;
        m_elem[idxRow * m_col + icol + 1] = value;
        m_elem[idxRow * m_col + icol + 2] = value;
        m_elem[idxRow * m_col + icol + 3] = value;
    }

    for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
    {
        m_elem[idxRow * m_col + icol] = value;
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetFillCol(const uint16_t idxCol, const t_type value)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxCol < m_col && "Index out of range");

    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * m_col + idxCol] = value;
        m_elem[(irow + 1) * m_col + idxCol] = value;
        m_elem[(irow + 2) * m_col + idxCol] = value;
        m_elem[(irow + 3) * m_col + idxCol] = value;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * m_col + idxCol] = value;
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetElement(const t_type *element, const size_t n_byte)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    size_t matSz = sizeof(t_type) * m_size;

    if (matSz > n_byte) memcpy(m_elem, element, n_byte);
    else memcpy(m_elem, element, matSz);
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetElement(const t_type *element, const uint16_t row, const uint16_t col)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");
    assert(m_col == col && "Col dimensions do not matched");

    if (m_size > row * col) memcpy(m_elem, element, sizeof(t_type) * row * col);
    else memcpy(m_elem, element, sizeof(t_type) * m_size);
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline void Matrix<0, 0, t_type>::SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix<row, col, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");
    assert(idxCol < m_col && "Index out of range");
    assert((m_row - idxRow) >= row && "Check the matrix size or row index");
    assert((m_col - idxCol) >= col && "Check the matrix size or col index");

    for (uint16_t irow = 0; irow < row; ++irow)
        for (uint16_t icol = 0; icol < col; ++icol)
            m_elem[(irow + idxRow) * m_col + idxCol + icol] = m.m_elem[irow * col + icol];
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix<0, 0, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");
    assert(idxCol < m_col && "Index out of range");
    assert((m_row - idxRow) >= m.m_row && "Check the matrix size or row index");
    assert((m_col - idxCol) >= m.m_col && "Check the matrix size or col index");

    for (uint16_t irow = 0; irow < m.m_row; ++irow)
        for (uint16_t icol = 0; icol < m.m_col; ++icol)
            m_elem[(irow + idxRow) * m_col + idxCol + icol] = m.m_elem[irow * m.m_col + icol];
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix3<t_type, 3, 3> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");
    assert(idxCol < m_col && "Index out of range");
    assert((m_row - idxRow) >= m.m_row && "Check the matrix size or row index");
    assert((m_col - idxCol) >= m.m_col && "Check the matrix size or col index");

    m_elem[idxRow * m_col + idxCol + 0] = m.m_elem[0];
    m_elem[idxRow * m_col + idxCol + 1] = m.m_elem[1];
    m_elem[idxRow * m_col + idxCol + 2] = m.m_elem[2];

    m_elem[(1 + idxRow) * m_col + idxCol + 0] = m.m_elem[3];
    m_elem[(1 + idxRow) * m_col + idxCol + 1] = m.m_elem[4];
    m_elem[(1 + idxRow) * m_col + idxCol + 2] = m.m_elem[5];

    m_elem[(2 + idxRow) * m_col + idxCol + 0] = m.m_elem[6];
    m_elem[(2 + idxRow) * m_col + idxCol + 1] = m.m_elem[7];
    m_elem[(2 + idxRow) * m_col + idxCol + 2] = m.m_elem[8];
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Rotation<t_type, 3, 3> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");
    assert(idxCol < m_col && "Index out of range");
    assert((m_row - idxRow) >= m.m_row && "Check the matrix size or row index");
    assert((m_col - idxCol) >= m.m_col && "Check the matrix size or col index");

    m_elem[idxRow * m_col + idxCol + 0] = m.m_elem[0];
    m_elem[idxRow * m_col + idxCol + 1] = m.m_elem[1];
    m_elem[idxRow * m_col + idxCol + 2] = m.m_elem[2];

    m_elem[(1 + idxRow) * m_col + idxCol + 0] = m.m_elem[3];
    m_elem[(1 + idxRow) * m_col + idxCol + 1] = m.m_elem[4];
    m_elem[(1 + idxRow) * m_col + idxCol + 2] = m.m_elem[5];

    m_elem[(2 + idxRow) * m_col + idxCol + 0] = m.m_elem[6];
    m_elem[(2 + idxRow) * m_col + idxCol + 1] = m.m_elem[7];
    m_elem[(2 + idxRow) * m_col + idxCol + 2] = m.m_elem[8];
}

template <typename t_type>
template <uint16_t col>
inline void Matrix<0, 0, t_type>::SetRowVec(const uint16_t idxRow, const Vector<col, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");

    uint16_t maxCol = (m_col < col) ? m_col : col;
    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = maxCol >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        m_elem[idxRow * m_col + icol] = v.m_elem[icol];
        m_elem[idxRow * m_col + icol + 1] = v.m_elem[icol + 1];
        m_elem[idxRow * m_col + icol + 2] = v.m_elem[icol + 2];
        m_elem[idxRow * m_col + icol + 3] = v.m_elem[icol + 3];
    }

    for (cnt = maxCol % 4u; cnt > 0u; cnt--, icol++)
    {
        m_elem[idxRow * m_col + icol] = v.m_elem[icol];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetRowVec(const uint16_t idxRow, const Vector<0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t maxCol = (m_col < v.m_row) ? m_col : v.m_row;
    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = maxCol >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        m_elem[idxRow * m_col + icol] = v.m_elem[icol];
        m_elem[idxRow * m_col + icol + 1] = v.m_elem[icol + 1];
        m_elem[idxRow * m_col + icol + 2] = v.m_elem[icol + 2];
        m_elem[idxRow * m_col + icol + 3] = v.m_elem[icol + 3];
    }

    for (cnt = maxCol % 4u; cnt > 0u; cnt--, icol++)
    {
        m_elem[idxRow * m_col + icol] = v.m_elem[icol];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetRowVec(const uint16_t idxRow, const Vector3<t_type, 3> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");

    if (m_col >= 3)
    {
        m_elem[m_col * idxRow] = v.m_elem[0];
        m_elem[m_col * idxRow + 1] = v.m_elem[1];
        m_elem[m_col * idxRow + 2] = v.m_elem[2];
    }
    else
    {
        for (uint16_t icol = 0; icol < m_col; icol++)
            m_elem[idxRow * m_col + icol] = v.m_elem[icol];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetRowVec(const uint16_t idxRow, const Vector4<t_type, 4> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");

    if (m_col >= 4)
    {
        m_elem[m_col * idxRow] = v.m_elem[0];
        m_elem[m_col * idxRow + 1] = v.m_elem[1];
        m_elem[m_col * idxRow + 2] = v.m_elem[2];
        m_elem[m_col * idxRow + 3] = v.m_elem[3];
    }
    else
    {
        for (uint16_t icol = 0; icol < m_col; icol++)
            m_elem[idxRow * m_col + icol] = v.m_elem[icol];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetRowVec(const uint16_t idxRow, const Vector6<t_type, 6> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");

    if (m_col >= 6)
    {
        m_elem[m_col * idxRow] = v.m_elem[0];
        m_elem[m_col * idxRow + 1] = v.m_elem[1];
        m_elem[m_col * idxRow + 2] = v.m_elem[2];
        m_elem[m_col * idxRow + 3] = v.m_elem[3];
        m_elem[m_col * idxRow + 4] = v.m_elem[4];
        m_elem[m_col * idxRow + 5] = v.m_elem[5];
    }
    else
    {
        for (uint16_t icol = 0; icol < m_col; icol++)
            m_elem[idxRow * m_col + icol] = v.m_elem[icol];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetRowVec(const uint16_t idxRow, const t_type *v, const size_t n_byte)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");

    uint16_t col = (uint16_t)(n_byte / sizeof(t_type));
    uint16_t maxCol = (m_col < col) ? m_col : col;
    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = maxCol >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        m_elem[idxRow * m_col + icol] = v[icol];
        m_elem[idxRow * m_col + icol + 1] = v[icol + 1];
        m_elem[idxRow * m_col + icol + 2] = v[icol + 2];
        m_elem[idxRow * m_col + icol + 3] = v[icol + 3];
    }

    for (cnt = maxCol % 4u; cnt > 0u; cnt--, icol++)
    {
        m_elem[idxRow * m_col + icol] = v[icol];
    }
}

template <typename t_type>
template <uint16_t row>
inline void Matrix<0, 0, t_type>::SetColVec(const uint16_t idxCol, const Vector<row, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxCol < m_col && "Index out of range");

    uint16_t maxRow = (m_row < row) ? m_row : row;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = maxRow >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * m_col + idxCol] = v.m_elem[irow];
        m_elem[(irow + 1) * m_col + idxCol] = v.m_elem[irow + 1];
        m_elem[(irow + 2) * m_col + idxCol] = v.m_elem[irow + 2];
        m_elem[(irow + 3) * m_col + idxCol] = v.m_elem[irow + 3];
    }

    for (cnt = maxRow % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * m_col + idxCol] = v.m_elem[irow];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetColVec(const uint16_t idxCol, const Vector<0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxCol < m_col && "Index out of range");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t maxRow = (m_row < v.m_row) ? m_row : v.m_row;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = maxRow >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * m_col + idxCol] = v.m_elem[irow];
        m_elem[(irow + 1) * m_col + idxCol] = v.m_elem[irow + 1];
        m_elem[(irow + 2) * m_col + idxCol] = v.m_elem[irow + 2];
        m_elem[(irow + 3) * m_col + idxCol] = v.m_elem[irow + 3];
    }

    for (cnt = maxRow % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * m_col + idxCol] = v.m_elem[irow];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetColVec(const uint16_t idxCol, const Vector3<t_type, 3> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxCol < m_col && "Index out of range");

    if (m_row >= 3)
    {
        m_elem[idxCol] = v.m_elem[0];
        m_elem[m_col + idxCol] = v.m_elem[1];
        m_elem[m_col * 2 + idxCol] = v.m_elem[2];
    }
    else
    {
        for (uint16_t irow = 0; irow < m_row; irow++)
            m_elem[irow * m_col + idxCol] = v.m_elem[irow];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetColVec(const uint16_t idxCol, const Vector4<t_type, 4> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxCol < m_col && "Index out of range");

    if (m_row >= 4)
    {
        m_elem[idxCol] = v.m_elem[0];
        m_elem[m_col + idxCol] = v.m_elem[1];
        m_elem[m_col * 2 + idxCol] = v.m_elem[2];
        m_elem[m_col * 3 + idxCol] = v.m_elem[3];
    }
    else
    {
        for (uint16_t irow = 0; irow < m_row; irow++)
            m_elem[irow * m_col + idxCol] = v.m_elem[irow];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetColVec(const uint16_t idxCol, const Vector6<t_type, 6> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxCol < m_col && "Index out of range");

    if (m_row >= 6)
    {
        m_elem[idxCol] = v.m_elem[0];
        m_elem[m_col + idxCol] = v.m_elem[1];
        m_elem[m_col * 2 + idxCol] = v.m_elem[2];
        m_elem[m_col * 3 + idxCol] = v.m_elem[3];
        m_elem[m_col * 4 + idxCol] = v.m_elem[4];
        m_elem[m_col * 5 + idxCol] = v.m_elem[5];
    }
    else
    {
        for (uint16_t irow = 0; irow < m_row; irow++)
            m_elem[irow * m_col + idxCol] = v.m_elem[irow];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetColVec(const uint16_t idxCol, const t_type *v, const size_t n_byte)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxCol < m_col && "Index out of range");

    uint16_t row = n_byte / sizeof(t_type);
    uint16_t maxRow = (m_row < row) ? m_row : row;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = maxRow >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * m_col + idxCol] = v[irow];
        m_elem[(irow + 1) * m_col + idxCol] = v[irow + 1];
        m_elem[(irow + 2) * m_col + idxCol] = v[irow + 2];
        m_elem[(irow + 3) * m_col + idxCol] = v[irow + 3];
    }

    for (cnt = maxRow % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * m_col + idxCol] = v[irow];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow1 < m_row && "Index out of range");
    assert(idxRow2 < m_row && "Index out of range");

    t_type tmpVec[m_col];
    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = m_col >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        tmpVec[icol] = m_elem[idxRow1 * m_col + icol];
        m_elem[idxRow1 * m_col + icol] = m_elem[idxRow2 * m_col + icol];
        m_elem[idxRow2 * m_col + icol] = tmpVec[icol];

        tmpVec[icol + 1] = m_elem[idxRow1 * m_col + icol + 1];
        m_elem[idxRow1 * m_col + icol + 1] = m_elem[idxRow2 * m_col + icol + 1];
        m_elem[idxRow2 * m_col + icol + 1] = tmpVec[icol + 1];

        tmpVec[icol + 2] = m_elem[idxRow1 * m_col + icol + 2];
        m_elem[idxRow1 * m_col + icol + 2] = m_elem[idxRow2 * m_col + icol + 2];
        m_elem[idxRow2 * m_col + icol + 2] = tmpVec[icol + 2];

        tmpVec[icol + 3] = m_elem[idxRow1 * m_col + icol + 3];
        m_elem[idxRow1 * m_col + icol + 3] = m_elem[idxRow2 * m_col + icol + 3];
        m_elem[idxRow2 * m_col + icol + 3] = tmpVec[icol + 3];
    }

    for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
    {
        tmpVec[icol] = m_elem[idxRow1 * m_col + icol];
        m_elem[idxRow1 * m_col + icol] = m_elem[idxRow2 * m_col + icol];
        m_elem[idxRow2 * m_col + icol] = tmpVec[icol];
    }
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxCol1 < m_col && "Index out of range");
    assert(idxCol2 < m_col && "Index out of range");

    t_type tmpVec[m_row];
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        tmpVec[irow] = m_elem[irow * m_col + idxCol1];
        m_elem[irow * m_col + idxCol1] = m_elem[irow * m_col + idxCol2];
        m_elem[irow * m_col + idxCol2] = tmpVec[irow];

        tmpVec[irow + 1] = m_elem[(irow + 1) * m_col + idxCol1];
        m_elem[(irow + 1) * m_col + idxCol1] = m_elem[(irow + 1) * m_col + idxCol2];
        m_elem[(irow + 1) * m_col + idxCol2] = tmpVec[irow + 1];

        tmpVec[irow + 2] = m_elem[(irow + 2) * m_col + idxCol1];
        m_elem[(irow + 2) * m_col + idxCol1] = m_elem[(irow + 2) * m_col + idxCol2];
        m_elem[(irow + 2) * m_col + idxCol2] = tmpVec[irow + 2];

        tmpVec[irow + 3] = m_elem[(irow + 3) * m_col + idxCol1];
        m_elem[(irow + 3) * m_col + idxCol1] = m_elem[(irow + 3) * m_col + idxCol2];
        m_elem[(irow + 3) * m_col + idxCol2] = tmpVec[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        tmpVec[irow] = m_elem[irow * m_col + idxCol1];
        m_elem[irow * m_col + idxCol1] = m_elem[irow * m_col + idxCol2];
        m_elem[irow * m_col + idxCol2] = tmpVec[irow];
    }
}

template <typename t_type>
inline const t_type *const Matrix<0, 0, t_type>::GetElementsAddr() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    return m_elem;
}

template <typename t_type>
inline size_t Matrix<0, 0, t_type>::GetMemSize() const
{
    return MemSize<t_type>(m_elem);
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Matrix<row, col, t_type> Matrix<0, 0, t_type>::GetBlock(const uint16_t idxRow, const uint16_t idxCol)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");
    assert(idxCol < m_col && "Index out of range");

    uint16_t irow, icol, cnt;
    t_type elem[row * col]{0};
    uint16_t rowSize = m_row - idxRow;
    uint16_t colSize = m_col - idxCol;

    if (idxRow >= m_row) return Matrix<row, col, t_type>(elem, sizeof(t_type) * row * col); // Index out of range, return zero mat
    if (idxCol >= m_col) return Matrix<row, col, t_type>(elem, sizeof(t_type) * row * col); // Index out of range, return zero mat
    if (rowSize > row) rowSize = row;
    if (colSize > col) colSize = col;

    for (irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            elem[irow * col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
            elem[irow * col + icol + 1] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 1];
            elem[irow * col + icol + 2] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 2];
            elem[irow * col + icol + 3] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            elem[irow * col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
        }
    }

    return Matrix<row, col, t_type>(elem);
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::GetBlock(const uint16_t idxRow, const uint16_t idxCol, const uint16_t row, const uint16_t col)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");
    assert(m_col > idxCol && "Index out of range");

    uint16_t irow, icol, cnt;
    uint16_t rowSize = m_row - idxRow;
    uint16_t colSize = m_col - idxCol;
    Matrix<0, 0, t_type> m(row, col);

    if (rowSize > row) rowSize = row;
    if (colSize > col) colSize = col;

    for (irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            m.m_elem[irow * col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
            m.m_elem[irow * col + icol + 1] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 1];
            m.m_elem[irow * col + icol + 2] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 2];
            m.m_elem[irow * col + icol + 3] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            m.m_elem[irow * col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
        }
    }

    return m;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t Matrix<0, 0, t_type>::GetBlock(const uint16_t idxRow, const uint16_t idxCol, Matrix<row, col, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(idxRow < m_row && "Index out of range");
    assert(idxCol < m_col && "Index out of range");

    uint16_t irow, icol, cnt;
    uint16_t rowSize = m_row - idxRow;
    uint16_t colSize = m_col - idxCol;

    if (idxRow >= m_row) return -1; // Index out of range
    if (idxCol >= m_col) return -1; // Index out of range
    if (rowSize > row) rowSize = row;
    if (colSize > col) colSize = col;

    for (irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            m.m_elem[irow * col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
            m.m_elem[irow * col + icol + 1] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 1];
            m.m_elem[irow * col + icol + 2] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 2];
            m.m_elem[irow * col + icol + 3] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            m.m_elem[irow * col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
        }
    }

    return 0;
}

template <typename t_type>
inline int8_t Matrix<0, 0, t_type>::GetBlock(const uint16_t idxRow, const uint16_t idxCol, Matrix<0, 0, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");
    assert(m_col > idxCol && "Index out of range");

    uint16_t irow, icol, cnt;
    uint16_t rowSize = m_row - idxRow;
    uint16_t colSize = m_col - idxCol;

    if (idxRow >= m_row) return -1; // Index out of range
    if (idxCol >= m_col) return -1; // Index out of range
    if (rowSize > m.m_row) rowSize = m.m_row;
    if (colSize > m.m_col) colSize = m.m_col;

    for (irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            m.m_elem[irow * m.m_col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
            m.m_elem[irow * m.m_col + icol + 1] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 1];
            m.m_elem[irow * m.m_col + icol + 2] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 2];
            m.m_elem[irow * m.m_col + icol + 3] = m_elem[(irow + idxRow) * m_col + idxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            m.m_elem[irow * m.m_col + icol] = m_elem[(irow + idxRow) * m_col + idxCol + icol];
        }
    }

    return 0;
}

template <typename t_type>
template <uint16_t col>
inline Vector<col, t_type> Matrix<0, 0, t_type>::GetRowVec(const uint16_t idxRow) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");
    // assert(m_col == col && "The sizes of return vector and row vector are different");

    t_type vec[col]{0};
    uint16_t cnt;
    uint16_t icol = 0;
    uint16_t colSize = (m_col < col) ? m_col : col;

    for (cnt = colSize >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        vec[icol] = m_elem[idxRow * m_col + icol];
        vec[icol + 1] = m_elem[idxRow * m_col + icol + 1];
        vec[icol + 2] = m_elem[idxRow * m_col + icol + 2];
        vec[icol + 3] = m_elem[idxRow * m_col + icol + 3];
    }

    for (cnt = colSize % 4u; cnt > 0u; cnt--, icol++)
    {
        vec[icol] = m_elem[idxRow * m_col + icol];
    }

    return Vector<col, t_type>(vec, sizeof(t_type) * col);
}

template <typename t_type>
template <uint16_t row>
inline Vector<row, t_type> Matrix<0, 0, t_type>::GetColVec(const uint16_t idxCol) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col > idxCol && "Index out of range");
    // assert(m_row == row && "The sizes of return vector and row vector are different");

    t_type vec[row]{0};
    uint16_t cnt;
    uint16_t irow = 0;
    uint16_t rowSize = (m_row < row) ? m_row : row;

    for (cnt = rowSize >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow * m_col + idxCol];
        vec[irow + 1] = m_elem[(irow + 1) * m_col + idxCol];
        vec[irow + 2] = m_elem[(irow + 2) * m_col + idxCol];
        vec[irow + 3] = m_elem[(irow + 3) * m_col + idxCol];
    }

    for (cnt = rowSize % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow * m_col + idxCol];
    }

    return Vector<row, t_type>(vec, sizeof(t_type) * row);
}

template <typename t_type>
inline Vector<0, t_type> Matrix<0, 0, t_type>::GetRowVec(const uint16_t idxRow, const uint16_t col) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");
    // assert(m_col == col && "The sizes of return vector and row vector are different");

    uint16_t cnt;
    uint16_t icol = 0;
    uint16_t colSize = (m_col < col) ? m_col : col;
    Vector<0, t_type> vec(col);

    for (cnt = colSize >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        vec.m_elem[icol] = m_elem[idxRow * m_col + icol];
        vec.m_elem[icol + 1] = m_elem[idxRow * m_col + icol + 1];
        vec.m_elem[icol + 2] = m_elem[idxRow * m_col + icol + 2];
        vec.m_elem[icol + 3] = m_elem[idxRow * m_col + icol + 3];
    }

    for (cnt = colSize % 4u; cnt > 0u; cnt--, icol++)
    {
        vec.m_elem[icol] = m_elem[idxRow * m_col + icol];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Matrix<0, 0, t_type>::GetColVec(const uint16_t idxCol, const uint16_t row) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col > idxCol && "Index out of range");
    // assert(m_row == row && "The sizes of return vector and row vector are different");

    uint16_t cnt;
    uint16_t irow = 0;
    uint16_t rowSize = (m_row < row) ? m_row : row;
    Vector<0, t_type> vec(row);

    for (cnt = rowSize >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow * m_col + idxCol];
        vec.m_elem[irow + 1] = m_elem[(irow + 1) * m_col + idxCol];
        vec.m_elem[irow + 2] = m_elem[(irow + 2) * m_col + idxCol];
        vec.m_elem[irow + 3] = m_elem[(irow + 3) * m_col + idxCol];
    }

    for (cnt = rowSize % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow * m_col + idxCol];
    }

    return vec;
}

template <typename t_type>
template <uint16_t col>
inline int8_t Matrix<0, 0, t_type>::GetRowVec(const uint16_t idxRow, Vector<col, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");
    // assert(m_col == col && "The sizes of return vector and row vector are different");

    uint16_t cnt;
    uint16_t icol = 0;
    uint16_t colSize = (m_col < col) ? m_col : col;

    for (cnt = colSize >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        v.m_elem[icol] = m_elem[idxRow * m_col + icol];
        v.m_elem[icol + 1] = m_elem[idxRow * m_col + icol + 1];
        v.m_elem[icol + 2] = m_elem[idxRow * m_col + icol + 2];
        v.m_elem[icol + 3] = m_elem[idxRow * m_col + icol + 3];
    }

    for (cnt = colSize % 4u; cnt > 0u; cnt--, icol++)
    {
        v.m_elem[icol] = m_elem[idxRow * m_col + icol];
    }

    return 0;
}

template <typename t_type>
template <uint16_t row>
inline int8_t Matrix<0, 0, t_type>::GetColVec(const uint16_t idxCol, Vector<row, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col > idxCol && "Index out of range");
    // assert(m_row >= row && "The sizes of return vector and row vector are different");

    uint16_t cnt;
    uint16_t irow = 0;
    uint16_t rowSize = (m_row < row) ? m_row : row;

    for (cnt = rowSize >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        v.m_elem[irow] = m_elem[irow * m_col + idxCol];
        v.m_elem[irow + 1] = m_elem[(irow + 1) * m_col + idxCol];
        v.m_elem[irow + 2] = m_elem[(irow + 2) * m_col + idxCol];
        v.m_elem[irow + 3] = m_elem[(irow + 3) * m_col + idxCol];
    }

    for (cnt = rowSize % 4u; cnt > 0u; cnt--, irow++)
    {
        v.m_elem[irow] = m_elem[irow * m_col + idxCol];
    }

    return 0;
}

template <typename t_type>
inline int8_t Matrix<0, 0, t_type>::GetRowVec(const uint16_t idxRow, Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    // assert(m_col == v.m_row && "The sizes of return vector and row vector are different");

    uint16_t cnt;
    uint16_t icol = 0;
    uint16_t colSize = (m_col < v.m_row) ? m_col : v.m_row;

    for (cnt = colSize >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        v.m_elem[icol] = m_elem[idxRow * m_col + icol];
        v.m_elem[icol + 1] = m_elem[idxRow * m_col + icol + 1];
        v.m_elem[icol + 2] = m_elem[idxRow * m_col + icol + 2];
        v.m_elem[icol + 3] = m_elem[idxRow * m_col + icol + 3];
    }

    for (cnt = colSize % 4u; cnt > 0u; cnt--, icol++)
    {
        v.m_elem[icol] = m_elem[idxRow * m_col + icol];
    }

    return 0;
}

template <typename t_type>
inline int8_t Matrix<0, 0, t_type>::GetColVec(const uint16_t idxCol, Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col > idxCol && "Index out of range");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    // assert(m_row >= v.m_row && "The sizes of return vector and row vector are different");

    uint16_t cnt;
    uint16_t irow = 0;
    uint16_t rowSize = (m_row < v.m_row) ? m_row : v.m_row;

    for (cnt = rowSize >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        v.m_elem[irow] = m_elem[irow * m_col + idxCol];
        v.m_elem[irow + 1] = m_elem[(irow + 1) * m_col + idxCol];
        v.m_elem[irow + 2] = m_elem[(irow + 2) * m_col + idxCol];
        v.m_elem[irow + 3] = m_elem[(irow + 3) * m_col + idxCol];
    }

    for (cnt = rowSize % 4u; cnt > 0u; cnt--, irow++)
    {
        v.m_elem[irow] = m_elem[irow * m_col + idxCol];
    }

    return 0;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline CscMatrix<row, col, t_type> Matrix<0, 0, t_type>::GetCscMat() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Check the row size");
    assert(m_col == col && "Check the col size");

    return CscMatrix<row, col, t_type>(*this);
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t Matrix<0, 0, t_type>::GetCscMat(CscMatrix<row, col, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Check the row size");
    assert(m_col == col && "Check the col size");

    m.SetElement(m_elem, sizeof(t_type) * row * col);
    return 0;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::Transpose() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    Matrix<0, 0, t_type> m(m_col, m_row);
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < m_row; ++irow)
    {
        for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
        {
            m.m_elem[irow + m_row * icol] = m_elem[irow * m_col + icol];
            m.m_elem[irow + m_row * (icol + 1)] = m_elem[irow * m_col + icol + 1];
            m.m_elem[irow + m_row * (icol + 2)] = m_elem[irow * m_col + icol + 2];
            m.m_elem[irow + m_row * (icol + 3)] = m_elem[irow * m_col + icol + 3];
        }

        for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
        {
            m.m_elem[irow + m_row * icol] = m_elem[irow * m_col + icol];
        }
    }

    return m;
}

template <typename t_type>
inline t_type Matrix<0, 0, t_type>::Trace() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    int cnt, i;
    uint16_t num = (m_row > m_col) ? m_col : m_row;
    uint16_t offset = m_col + 1;
    t_type sum = 0;

    // for (uint16_t i = 0; i < num; i++)
    //     sum += m_elem[i * offset];
    for (cnt = num >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        sum += m_elem[i * offset];
        sum += m_elem[(i + 1) * offset];
        sum += m_elem[(i + 2) * offset];
        sum += m_elem[(i + 3) * offset];
    }

    for (cnt = num % 4; cnt > 0; cnt--, i++)
    {
        sum += m_elem[i * offset];
    }

    return sum;
}

template <typename t_type>
inline t_type Matrix<0, 0, t_type>::GetNorm() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    t_type sqSum = 0;

    for (cnt = m_size >> 2u; cnt > 0; cnt--, i += 4)
    {
        sqSum += m_elem[i] * m_elem[i];
        sqSum += m_elem[i + 1] * m_elem[i + 1];
        sqSum += m_elem[i + 2] * m_elem[i + 2];
        sqSum += m_elem[i + 3] * m_elem[i + 3];
    }

    for (cnt = m_size % 4u; cnt > 0u; cnt--, i++)
    {
        sqSum += m_elem[i] * m_elem[i];
    }

    return std::sqrt(sqSum);
}

template <typename t_type>
inline t_type Matrix<0, 0, t_type>::GetSqNorm() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    t_type sqSum = 0;

    for (cnt = m_size >> 2u; cnt > 0u; cnt--, i += 4)
    {
        sqSum += m_elem[i] * m_elem[i];
        sqSum += m_elem[i + 1] * m_elem[i + 1];
        sqSum += m_elem[i + 2] * m_elem[i + 2];
        sqSum += m_elem[i + 3] * m_elem[i + 3];
    }

    for (cnt = m_size % 4u; cnt > 0u; cnt--, i++)
    {
        sqSum += m_elem[i] * m_elem[i];
    }

    return sqSum;
}

template <typename t_type>
inline t_type Matrix<0, 0, t_type>::GetLpNorm(const int p) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    t_type powSum = 0;

    for (cnt = m_size >> 2u; cnt > 0u; cnt--, i += 4)
    {
        powSum += std::pow(std::abs(m_elem[i]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 1]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 2]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 3]), (t_type)p);
    }

    for (cnt = m_size % 4u; cnt > 0u; cnt--, i++)
    {
        powSum += std::pow(std::abs(m_elem[i]), (t_type)p);
    }

    return std::pow(powSum, (t_type)1 / p);
}

template <typename t_type>
inline t_type Matrix<0, 0, t_type>::Determinant() const
{
    return PartialPivLU<0, 0, t_type>(m_row, m_col, m_elem).Determinant();
}

template <typename t_type>
inline NoPivLU<0, 0, t_type> Matrix<0, 0, t_type>::GetNoPivLU() const
{
    return NoPivLU<0, 0, t_type>(*this);
}

template <typename t_type>
inline PartialPivLU<0, 0, t_type> Matrix<0, 0, t_type>::GetPartialPivLU() const
{
    return PartialPivLU<0, 0, t_type>(m_row, m_col, m_elem);
}

template <typename t_type>
inline FullPivLU<0, 0, t_type> Matrix<0, 0, t_type>::GetFullPivLU() const
{
    return FullPivLU<0, 0, t_type>(*this);
}

template <typename t_type>
inline LLT<0, 0, t_type> Matrix<0, 0, t_type>::GetLLT() const
{
    return LLT<0, 0, t_type>(*this);
}

template <typename t_type>
inline LDLT<0, 0, t_type> Matrix<0, 0, t_type>::GetLDLT() const
{
    return LDLT<0, 0, t_type>(*this);
}

// template <typename t_type>
// inline QR<m_row, m_col, t_type> Matrix<0, 0, t_type>::GetQR() const
// {
//     return QR<m_row, m_col, t_type>(*this);
// }

template <typename t_type>
inline SVD<0, 0, t_type> Matrix<0, 0, t_type>::GetSVD() const
{
    return SVD<0, 0, t_type>(*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::Inv(int8_t *isOk) const
{
    return PartialPivLU<0, 0, t_type>(m_row, m_col, m_elem).Inverse(isOk);
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::FInv(int8_t *isOk) const
{
    return FullPivLU<0, 0, t_type>(*this).Inverse(isOk);
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::PInv(int8_t *isOk, t_type tolerance) const
{
    return SVD<0, 0, t_type>(*this).Inverse(isOk, tolerance);
}

/* Member access operators */
template <typename t_type>
inline t_type &Matrix<0, 0, t_type>::operator()(uint16_t irow, uint16_t icol)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > irow && "Index out of range");
    assert(m_col > icol && "Index out of range");

    return m_elem[irow * m_col + icol];
}

template <typename t_type>
inline const t_type &Matrix<0, 0, t_type>::operator()(uint16_t irow, uint16_t icol) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > irow && "Index out of range");
    assert(m_col > icol && "Index out of range");

    return m_elem[irow * m_col + icol];
}

/* Assignment operators */
template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator=(const Matrix<0, 0, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == m.m_row && "Row dimensions do not matched");
    assert(m_col == m.m_col && "Col dimensions do not matched");

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator+=(const Matrix<0, 0, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == m.m_row && "Row dimensions do not matched");
    assert(m_col == m.m_col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] += m.m_elem[i];
        m_elem[i + 2] += m.m_elem[i + 2];
        m_elem[i + 4] += m.m_elem[i + 4];
        m_elem[i + 6] += m.m_elem[i + 6];
        m_elem[i + 1] += m.m_elem[i + 1];
        m_elem[i + 3] += m.m_elem[i + 3];
        m_elem[i + 5] += m.m_elem[i + 5];
        m_elem[i + 7] += m.m_elem[i + 7];
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] += m.m_elem[i];
    }

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator-=(const Matrix<0, 0, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == m.m_row && "Row dimensions do not matched");
    assert(m_col == m.m_col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] -= m.m_elem[i];
        m_elem[i + 2] -= m.m_elem[i + 2];
        m_elem[i + 4] -= m.m_elem[i + 4];
        m_elem[i + 6] -= m.m_elem[i + 6];
        m_elem[i + 1] -= m.m_elem[i + 1];
        m_elem[i + 3] -= m.m_elem[i + 3];
        m_elem[i + 5] -= m.m_elem[i + 5];
        m_elem[i + 7] -= m.m_elem[i + 7];
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] -= m.m_elem[i];
    }

    return (*this);
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator=(const Matrix<row, col, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");
    assert(m_col == col && "Col dimensions do not matched");

    memcpy(m_elem, m.m_elem, sizeof(t_type) * m_size);

    return (*this);
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator+=(const Matrix<row, col, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");
    assert(m_col == col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] += m.m_elem[i];
        m_elem[i + 2] += m.m_elem[i + 2];
        m_elem[i + 4] += m.m_elem[i + 4];
        m_elem[i + 6] += m.m_elem[i + 6];
        m_elem[i + 1] += m.m_elem[i + 1];
        m_elem[i + 3] += m.m_elem[i + 3];
        m_elem[i + 5] += m.m_elem[i + 5];
        m_elem[i + 7] += m.m_elem[i + 7];
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] += m.m_elem[i];
    }

    return (*this);
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator-=(const Matrix<row, col, t_type> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");
    assert(m_col == col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] -= m.m_elem[i];
        m_elem[i + 2] -= m.m_elem[i + 2];
        m_elem[i + 4] -= m.m_elem[i + 4];
        m_elem[i + 6] -= m.m_elem[i + 6];
        m_elem[i + 1] -= m.m_elem[i + 1];
        m_elem[i + 3] -= m.m_elem[i + 3];
        m_elem[i + 5] -= m.m_elem[i + 5];
        m_elem[i + 7] -= m.m_elem[i + 7];
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] -= m.m_elem[i];
    }

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator=(const Matrix3<t_type, 3, 3> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(m_col == 3 && "Col dimensions do not matched");

    memcpy(m_elem, m.m_elem, sizeof(t_type) * 9);

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator+=(const Matrix3<t_type, 3, 3> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(m_col == 3 && "Col dimensions do not matched");

    m_elem[0] += m.m_elem[0];
    m_elem[1] += m.m_elem[1];
    m_elem[2] += m.m_elem[2];
    m_elem[3] += m.m_elem[3];
    m_elem[4] += m.m_elem[4];
    m_elem[5] += m.m_elem[5];
    m_elem[6] += m.m_elem[6];
    m_elem[7] += m.m_elem[7];
    m_elem[8] += m.m_elem[8];

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator-=(const Matrix3<t_type, 3, 3> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(m_col == 3 && "Col dimensions do not matched");

    m_elem[0] -= m.m_elem[0];
    m_elem[1] -= m.m_elem[1];
    m_elem[2] -= m.m_elem[2];
    m_elem[3] -= m.m_elem[3];
    m_elem[4] -= m.m_elem[4];
    m_elem[5] -= m.m_elem[5];
    m_elem[6] -= m.m_elem[6];
    m_elem[7] -= m.m_elem[7];
    m_elem[8] -= m.m_elem[8];

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator=(const Rotation<t_type, 3, 3> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(m_col == 3 && "Col dimensions do not matched");

    memcpy(m_elem, m.m_elem, sizeof(t_type) * 9);

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator+=(const Rotation<t_type, 3, 3> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(m_col == 3 && "Col dimensions do not matched");

    m_elem[0] += m.m_elem[0];
    m_elem[1] += m.m_elem[1];
    m_elem[2] += m.m_elem[2];
    m_elem[3] += m.m_elem[3];
    m_elem[4] += m.m_elem[4];
    m_elem[5] += m.m_elem[5];
    m_elem[6] += m.m_elem[6];
    m_elem[7] += m.m_elem[7];
    m_elem[8] += m.m_elem[8];

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator-=(const Rotation<t_type, 3, 3> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(m_col == 3 && "Col dimensions do not matched");

    m_elem[0] -= m.m_elem[0];
    m_elem[1] -= m.m_elem[1];
    m_elem[2] -= m.m_elem[2];
    m_elem[3] -= m.m_elem[3];
    m_elem[4] -= m.m_elem[4];
    m_elem[5] -= m.m_elem[5];
    m_elem[6] -= m.m_elem[6];
    m_elem[7] -= m.m_elem[7];
    m_elem[8] -= m.m_elem[8];

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator=(const Transform<t_type, 4, 4> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 4 && "Row dimensions do not matched");
    assert(m_col == 4 && "Col dimensions do not matched");

    m_elem[0] = m.m_R.m_elem[0];
    m_elem[1] = m.m_R.m_elem[1];
    m_elem[2] = m.m_R.m_elem[2];
    m_elem[3] = m.m_p.m_elem[0];
    m_elem[4] = m.m_R.m_elem[3];
    m_elem[5] = m.m_R.m_elem[4];
    m_elem[6] = m.m_R.m_elem[5];
    m_elem[7] = m.m_p.m_elem[1];
    m_elem[8] = m.m_R.m_elem[6];
    m_elem[9] = m.m_R.m_elem[7];
    m_elem[10] = m.m_R.m_elem[8];
    m_elem[11] = m.m_p.m_elem[2];
    m_elem[12] = 0;
    m_elem[13] = 0;
    m_elem[14] = 0;
    m_elem[15] = 1;

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator+=(const Transform<t_type, 4, 4> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 4 && "Row dimensions do not matched");
    assert(m_col == 4 && "Col dimensions do not matched");

    m_elem[0] += m.m_R.m_elem[0];
    m_elem[1] += m.m_R.m_elem[1];
    m_elem[2] += m.m_R.m_elem[2];
    m_elem[3] += m.m_p.m_elem[0];
    m_elem[4] += m.m_R.m_elem[3];
    m_elem[5] += m.m_R.m_elem[4];
    m_elem[6] += m.m_R.m_elem[5];
    m_elem[7] += m.m_p.m_elem[1];
    m_elem[8] += m.m_R.m_elem[6];
    m_elem[9] += m.m_R.m_elem[7];
    m_elem[10] += m.m_R.m_elem[8];
    m_elem[11] += m.m_p.m_elem[2];
    m_elem[15] += 1;

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator-=(const Transform<t_type, 4, 4> &m)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 4 && "Row dimensions do not matched");
    assert(m_col == 4 && "Col dimensions do not matched");

    m_elem[0] -= m.m_R.m_elem[0];
    m_elem[1] -= m.m_R.m_elem[1];
    m_elem[2] -= m.m_R.m_elem[2];
    m_elem[3] -= m.m_p.m_elem[0];
    m_elem[4] -= m.m_R.m_elem[3];
    m_elem[5] -= m.m_R.m_elem[4];
    m_elem[6] -= m.m_R.m_elem[5];
    m_elem[7] -= m.m_p.m_elem[1];
    m_elem[8] -= m.m_R.m_elem[6];
    m_elem[9] -= m.m_R.m_elem[7];
    m_elem[10] -= m.m_R.m_elem[8];
    m_elem[11] -= m.m_p.m_elem[2];
    m_elem[15] -= 1;

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator=(const t_type s)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] = s;
        m_elem[i + 2] = s;
        m_elem[i + 4] = s;
        m_elem[i + 6] = s;
        m_elem[i + 1] = s;
        m_elem[i + 3] = s;
        m_elem[i + 5] = s;
        m_elem[i + 7] = s;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] = s;
    }

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator+=(const t_type s)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] += s;
        m_elem[i + 2] += s;
        m_elem[i + 4] += s;
        m_elem[i + 6] += s;
        m_elem[i + 1] += s;
        m_elem[i + 3] += s;
        m_elem[i + 5] += s;
        m_elem[i + 7] += s;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] += s;
    }

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator-=(const t_type s)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] -= s;
        m_elem[i + 2] -= s;
        m_elem[i + 4] -= s;
        m_elem[i + 6] -= s;
        m_elem[i + 1] -= s;
        m_elem[i + 3] -= s;
        m_elem[i + 5] -= s;
        m_elem[i + 7] -= s;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] -= s;
    }

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator*=(const t_type s)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] *= s;
        m_elem[i + 2] *= s;
        m_elem[i + 4] *= s;
        m_elem[i + 6] *= s;
        m_elem[i + 1] *= s;
        m_elem[i + 3] *= s;
        m_elem[i + 5] *= s;
        m_elem[i + 7] *= s;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] *= s;
    }

    return (*this);
}

template <typename t_type>
inline Matrix<0, 0, t_type> &Matrix<0, 0, t_type>::operator/=(const t_type s)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    t_type scalar = s;

    if (std::abs(scalar) < std::numeric_limits<t_type>::epsilon())
    {
        if (scalar < 0)
            scalar = -std::numeric_limits<t_type>::epsilon();
        else
            scalar = std::numeric_limits<t_type>::epsilon();
    }

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] /= scalar;
        m_elem[i + 2] /= scalar;
        m_elem[i + 4] /= scalar;
        m_elem[i + 6] /= scalar;
        m_elem[i + 1] /= scalar;
        m_elem[i + 3] /= scalar;
        m_elem[i + 5] /= scalar;
        m_elem[i + 7] /= scalar;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] /= scalar;
    }

    return (*this);
}

template <typename t_type>
inline CommaInit<0, t_type> Matrix<0, 0, t_type>::operator<<(const t_type s)
{
    m_elem[0] = s;
    return CommaInit<0, t_type>(m_elem, m_size);
}

/* Arithmetic operators */
template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator-() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    Matrix<0, 0, t_type> mat(m_row, m_col);

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat.m_elem[i] = -m_elem[i];
        mat.m_elem[i + 2] = -m_elem[i + 2];
        mat.m_elem[i + 4] = -m_elem[i + 4];
        mat.m_elem[i + 6] = -m_elem[i + 6];
        mat.m_elem[i + 1] = -m_elem[i + 1];
        mat.m_elem[i + 3] = -m_elem[i + 3];
        mat.m_elem[i + 5] = -m_elem[i + 5];
        mat.m_elem[i + 7] = -m_elem[i + 7];
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = -m_elem[i];
    }

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator+(const Matrix<0, 0, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == m.m_row && "Row dimensions do not matched");
    assert(m_col == m.m_col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;
    Matrix<0, 0, t_type> mat(m_row, m_col);

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat.m_elem[i] = m_elem[i] + m.m_elem[i];
        mat.m_elem[i + 2] = m_elem[i + 2] + m.m_elem[i + 2];
        mat.m_elem[i + 4] = m_elem[i + 4] + m.m_elem[i + 4];
        mat.m_elem[i + 6] = m_elem[i + 6] + m.m_elem[i + 6];
        mat.m_elem[i + 1] = m_elem[i + 1] + m.m_elem[i + 1];
        mat.m_elem[i + 3] = m_elem[i + 3] + m.m_elem[i + 3];
        mat.m_elem[i + 5] = m_elem[i + 5] + m.m_elem[i + 5];
        mat.m_elem[i + 7] = m_elem[i + 7] + m.m_elem[i + 7];
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = m_elem[i] + m.m_elem[i];
    }

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator-(const Matrix<0, 0, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == m.m_row && "Row dimensions do not matched");
    assert(m_col == m.m_col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;
    Matrix<0, 0, t_type> mat(m_row, m_col);

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat.m_elem[i] = m_elem[i] - m.m_elem[i];
        mat.m_elem[i + 2] = m_elem[i + 2] - m.m_elem[i + 2];
        mat.m_elem[i + 4] = m_elem[i + 4] - m.m_elem[i + 4];
        mat.m_elem[i + 6] = m_elem[i + 6] - m.m_elem[i + 6];
        mat.m_elem[i + 1] = m_elem[i + 1] - m.m_elem[i + 1];
        mat.m_elem[i + 3] = m_elem[i + 3] - m.m_elem[i + 3];
        mat.m_elem[i + 5] = m_elem[i + 5] - m.m_elem[i + 5];
        mat.m_elem[i + 7] = m_elem[i + 7] - m.m_elem[i + 7];
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = m_elem[i] - m.m_elem[i];
    }

    return mat;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator+(const Matrix<row, col, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");
    assert(m_col == col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;
    Matrix<0, 0, t_type> mat(m_row, m_col);

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat.m_elem[i] = m_elem[i] + m.m_elem[i];
        mat.m_elem[i + 2] = m_elem[i + 2] + m.m_elem[i + 2];
        mat.m_elem[i + 4] = m_elem[i + 4] + m.m_elem[i + 4];
        mat.m_elem[i + 6] = m_elem[i + 6] + m.m_elem[i + 6];
        mat.m_elem[i + 1] = m_elem[i + 1] + m.m_elem[i + 1];
        mat.m_elem[i + 3] = m_elem[i + 3] + m.m_elem[i + 3];
        mat.m_elem[i + 5] = m_elem[i + 5] + m.m_elem[i + 5];
        mat.m_elem[i + 7] = m_elem[i + 7] + m.m_elem[i + 7];
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = m_elem[i] + m.m_elem[i];
    }

    return mat;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator-(const Matrix<row, col, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");
    assert(m_col == col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;
    Matrix<0, 0, t_type> mat(m_row, m_col);

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat.m_elem[i] = m_elem[i] - m.m_elem[i];
        mat.m_elem[i + 2] = m_elem[i + 2] - m.m_elem[i + 2];
        mat.m_elem[i + 4] = m_elem[i + 4] - m.m_elem[i + 4];
        mat.m_elem[i + 6] = m_elem[i + 6] - m.m_elem[i + 6];
        mat.m_elem[i + 1] = m_elem[i + 1] - m.m_elem[i + 1];
        mat.m_elem[i + 3] = m_elem[i + 3] - m.m_elem[i + 3];
        mat.m_elem[i + 5] = m_elem[i + 5] - m.m_elem[i + 5];
        mat.m_elem[i + 7] = m_elem[i + 7] - m.m_elem[i + 7];
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = m_elem[i] - m.m_elem[i];
    }

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator+(const Matrix3<t_type, 3, 3> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(m_col == 3 && "Col dimensions do not matched");

    Matrix<0, 0, t_type> mat(m_row, m_col);

    mat.m_elem[0] = m_elem[0] + m.m_elem[0];
    mat.m_elem[1] = m_elem[1] + m.m_elem[1];
    mat.m_elem[2] = m_elem[2] + m.m_elem[2];
    mat.m_elem[3] = m_elem[3] + m.m_elem[3];
    mat.m_elem[4] = m_elem[4] + m.m_elem[4];
    mat.m_elem[5] = m_elem[5] + m.m_elem[5];
    mat.m_elem[6] = m_elem[6] + m.m_elem[6];
    mat.m_elem[7] = m_elem[7] + m.m_elem[7];
    mat.m_elem[8] = m_elem[8] + m.m_elem[8];

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator-(const Matrix3<t_type, 3, 3> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(m_col == 3 && "Col dimensions do not matched");

    Matrix<0, 0, t_type> mat(m_row, m_col);

    mat.m_elem[0] = m_elem[0] - m.m_elem[0];
    mat.m_elem[1] = m_elem[1] - m.m_elem[1];
    mat.m_elem[2] = m_elem[2] - m.m_elem[2];
    mat.m_elem[3] = m_elem[3] - m.m_elem[3];
    mat.m_elem[4] = m_elem[4] - m.m_elem[4];
    mat.m_elem[5] = m_elem[5] - m.m_elem[5];
    mat.m_elem[6] = m_elem[6] - m.m_elem[6];
    mat.m_elem[7] = m_elem[7] - m.m_elem[7];
    mat.m_elem[8] = m_elem[8] - m.m_elem[8];

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator+(const Rotation<t_type, 3, 3> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(m_col == 3 && "Col dimensions do not matched");

    Matrix<0, 0, t_type> mat(m_row, m_col);

    mat.m_elem[0] = m_elem[0] + m.m_elem[0];
    mat.m_elem[1] = m_elem[1] + m.m_elem[1];
    mat.m_elem[2] = m_elem[2] + m.m_elem[2];
    mat.m_elem[3] = m_elem[3] + m.m_elem[3];
    mat.m_elem[4] = m_elem[4] + m.m_elem[4];
    mat.m_elem[5] = m_elem[5] + m.m_elem[5];
    mat.m_elem[6] = m_elem[6] + m.m_elem[6];
    mat.m_elem[7] = m_elem[7] + m.m_elem[7];
    mat.m_elem[8] = m_elem[8] + m.m_elem[8];

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator-(const Rotation<t_type, 3, 3> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(m_col == 3 && "Col dimensions do not matched");

    Matrix<0, 0, t_type> mat(m_row, m_col);

    mat.m_elem[0] = m_elem[0] - m.m_elem[0];
    mat.m_elem[1] = m_elem[1] - m.m_elem[1];
    mat.m_elem[2] = m_elem[2] - m.m_elem[2];
    mat.m_elem[3] = m_elem[3] - m.m_elem[3];
    mat.m_elem[4] = m_elem[4] - m.m_elem[4];
    mat.m_elem[5] = m_elem[5] - m.m_elem[5];
    mat.m_elem[6] = m_elem[6] - m.m_elem[6];
    mat.m_elem[7] = m_elem[7] - m.m_elem[7];
    mat.m_elem[8] = m_elem[8] - m.m_elem[8];

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator+(const Transform<t_type, 4, 4> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 4 && "Row dimensions do not matched");
    assert(m_col == 4 && "Col dimensions do not matched");

    Matrix<0, 0, t_type> mat(m_row, m_col);

    mat.m_elem[0] = m_elem[0] + m.m_R.m_elem[0];
    mat.m_elem[1] = m_elem[1] + m.m_R.m_elem[1];
    mat.m_elem[2] = m_elem[2] + m.m_R.m_elem[2];
    mat.m_elem[3] = m_elem[3] + m.m_p.m_elem[0];
    mat.m_elem[4] = m_elem[4] + m.m_R.m_elem[3];
    mat.m_elem[5] = m_elem[5] + m.m_R.m_elem[4];
    mat.m_elem[6] = m_elem[6] + m.m_R.m_elem[5];
    mat.m_elem[7] = m_elem[7] + m.m_p.m_elem[1];
    mat.m_elem[8] = m_elem[8] + m.m_R.m_elem[6];
    mat.m_elem[9] = m_elem[9] + m.m_R.m_elem[7];
    mat.m_elem[10] = m_elem[10] + m.m_R.m_elem[8];
    mat.m_elem[11] = m_elem[11] + m.m_p.m_elem[2];
    mat.m_elem[12] = m_elem[12];
    mat.m_elem[13] = m_elem[13];
    mat.m_elem[14] = m_elem[14];
    mat.m_elem[15] = m_elem[15] + 1;

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator-(const Transform<t_type, 4, 4> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 4 && "Row dimensions do not matched");
    assert(m_col == 4 && "Col dimensions do not matched");

    Matrix<0, 0, t_type> mat(m_row, m_col);

    mat.m_elem[0] = m_elem[0] - m.m_R.m_elem[0];
    mat.m_elem[1] = m_elem[1] - m.m_R.m_elem[1];
    mat.m_elem[2] = m_elem[2] - m.m_R.m_elem[2];
    mat.m_elem[3] = m_elem[3] - m.m_p.m_elem[0];
    mat.m_elem[4] = m_elem[4] - m.m_R.m_elem[3];
    mat.m_elem[5] = m_elem[5] - m.m_R.m_elem[4];
    mat.m_elem[6] = m_elem[6] - m.m_R.m_elem[5];
    mat.m_elem[7] = m_elem[7] - m.m_p.m_elem[1];
    mat.m_elem[8] = m_elem[8] - m.m_R.m_elem[6];
    mat.m_elem[9] = m_elem[9] - m.m_R.m_elem[7];
    mat.m_elem[10] = m_elem[10] - m.m_R.m_elem[8];
    mat.m_elem[11] = m_elem[11] - m.m_p.m_elem[2];
    mat.m_elem[12] = m_elem[12];
    mat.m_elem[13] = m_elem[13];
    mat.m_elem[14] = m_elem[14];
    mat.m_elem[15] = m_elem[15] - 1;

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator+(const t_type s) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    Matrix<0, 0, t_type> mat(m_row, m_col);

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat.m_elem[i] = m_elem[i] + s;
        mat.m_elem[i + 2] = m_elem[i + 2] + s;
        mat.m_elem[i + 4] = m_elem[i + 4] + s;
        mat.m_elem[i + 6] = m_elem[i + 6] + s;
        mat.m_elem[i + 1] = m_elem[i + 1] + s;
        mat.m_elem[i + 3] = m_elem[i + 3] + s;
        mat.m_elem[i + 5] = m_elem[i + 5] + s;
        mat.m_elem[i + 7] = m_elem[i + 7] + s;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = m_elem[i] + s;
    }

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator-(const t_type s) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    Matrix<0, 0, t_type> mat(m_row, m_col);

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat.m_elem[i] = m_elem[i] - s;
        mat.m_elem[i + 2] = m_elem[i + 2] - s;
        mat.m_elem[i + 4] = m_elem[i + 4] - s;
        mat.m_elem[i + 6] = m_elem[i + 6] - s;
        mat.m_elem[i + 1] = m_elem[i + 1] - s;
        mat.m_elem[i + 3] = m_elem[i + 3] - s;
        mat.m_elem[i + 5] = m_elem[i + 5] - s;
        mat.m_elem[i + 7] = m_elem[i + 7] - s;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = m_elem[i] - s;
    }

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator*(const t_type s) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    Matrix<0, 0, t_type> mat(m_row, m_col);

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat.m_elem[i] = m_elem[i] * s;
        mat.m_elem[i + 2] = m_elem[i + 2] * s;
        mat.m_elem[i + 4] = m_elem[i + 4] * s;
        mat.m_elem[i + 6] = m_elem[i + 6] * s;
        mat.m_elem[i + 1] = m_elem[i + 1] * s;
        mat.m_elem[i + 3] = m_elem[i + 3] * s;
        mat.m_elem[i + 5] = m_elem[i + 5] * s;
        mat.m_elem[i + 7] = m_elem[i + 7] * s;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = m_elem[i] * s;
    }

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator/(const t_type s) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    t_type scalar = s;
    Matrix<0, 0, t_type> mat(m_row, m_col);

    if (std::abs(scalar) < std::numeric_limits<t_type>::epsilon())
    {
        if (scalar < 0)
            scalar = -std::numeric_limits<t_type>::epsilon();
        else
            scalar = std::numeric_limits<t_type>::epsilon();
    }

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat.m_elem[i] = m_elem[i] / scalar;
        mat.m_elem[i + 2] = m_elem[i + 2] / scalar;
        mat.m_elem[i + 4] = m_elem[i + 4] / scalar;
        mat.m_elem[i + 6] = m_elem[i + 6] / scalar;
        mat.m_elem[i + 1] = m_elem[i + 1] / scalar;
        mat.m_elem[i + 3] = m_elem[i + 3] / scalar;
        mat.m_elem[i + 5] = m_elem[i + 5] / scalar;
        mat.m_elem[i + 7] = m_elem[i + 7] / scalar;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = m_elem[i] / scalar;
    }

    return mat;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator*(const Matrix<row, col, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == row && "Dimensions do not matched");

    uint16_t cnt;
    uint16_t irow, icol, i;
    Matrix<0, 0, t_type> mat(m_row, col);

    for (irow = 0; irow < m_row; ++irow)
    {
        for (icol = 0; icol < col; ++icol)
        {
            for (cnt = row >> 2u, i = 0; cnt > 0u; cnt--, i += 4u)
            {
                mat.m_elem[irow * col + icol] += m_elem[irow * row + i] * m.m_elem[i * col + icol];
                mat.m_elem[irow * col + icol] += m_elem[irow * row + i + 1] * m.m_elem[(i + 1) * col + icol];
                mat.m_elem[irow * col + icol] += m_elem[irow * row + i + 2] * m.m_elem[(i + 2) * col + icol];
                mat.m_elem[irow * col + icol] += m_elem[irow * row + i + 3] * m.m_elem[(i + 3) * col + icol];
            }

            for (cnt = row % 4u; cnt > 0u; cnt--, i++)
            {
                mat.m_elem[irow * col + icol] += m_elem[irow * row + i] * m.m_elem[i * col + icol];
            }
        }
    }

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator*(const Matrix<0, 0, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == m.m_row && "Dimensions do not matched");

    uint16_t cnt;
    uint16_t irow, icol, i;
    Matrix<0, 0, t_type> mat(m_row, m.m_col);

    for (irow = 0; irow < m_row; ++irow)
    {
        for (icol = 0; icol < m.m_col; ++icol)
        {
            for (cnt = m_col >> 2u, i = 0; cnt > 0u; cnt--, i += 4u)
            {
                mat.m_elem[irow * m.m_col + icol] += m_elem[irow * m_col + i] * m.m_elem[i * m.m_col + icol];
                mat.m_elem[irow * m.m_col + icol] += m_elem[irow * m_col + i + 1] * m.m_elem[(i + 1) * m.m_col + icol];
                mat.m_elem[irow * m.m_col + icol] += m_elem[irow * m_col + i + 2] * m.m_elem[(i + 2) * m.m_col + icol];
                mat.m_elem[irow * m.m_col + icol] += m_elem[irow * m_col + i + 3] * m.m_elem[(i + 3) * m.m_col + icol];
            }

            for (cnt = m_col % 4u; cnt > 0u; cnt--, i++)
            {
                mat.m_elem[irow * m.m_col + icol] += m_elem[irow * m_col + i] * m.m_elem[i * m.m_col + icol];
            }
        }
    }

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator*(const Matrix3<t_type, 3, 3> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == 3 && "Dimensions do not matched");

    Matrix<0, 0, t_type> mat(m_row, 3);

    for (uint16_t irow = 0; irow < m_row; ++irow)
    {
        for (uint16_t icol = 0; icol < 3; ++icol)
        {
            mat.m_elem[irow * 3 + icol] += m_elem[irow * 3] * m.m_elem[icol];
            mat.m_elem[irow * 3 + icol] += m_elem[irow * 3 + 1] * m.m_elem[3 + icol];
            mat.m_elem[irow * 3 + icol] += m_elem[irow * 3 + 2] * m.m_elem[6 + icol];
        }
    }

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator*(const Rotation<t_type, 3, 3> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == 3 && "Dimensions do not matched");

    Matrix<0, 0, t_type> mat(m_row, 3);

    for (uint16_t irow = 0; irow < m_row; ++irow)
    {
        for (uint16_t icol = 0; icol < 3; ++icol)
        {
            mat.m_elem[irow * 3 + icol] += m_elem[irow * 3] * m.m_elem[icol];
            mat.m_elem[irow * 3 + icol] += m_elem[irow * 3 + 1] * m.m_elem[3 + icol];
            mat.m_elem[irow * 3 + icol] += m_elem[irow * 3 + 2] * m.m_elem[6 + icol];
        }
    }

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator*(const Transform<t_type, 4, 4> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == 4 && "Dimensions do not matched");

    Matrix<0, 0, t_type> mat(m_row, 4);

    for (uint16_t irow = 0; irow < m_row * 4; irow += 4)
    {
        for (uint16_t icol = 0; icol < 3; ++icol)
        {
            mat.m_elem[irow + icol] =
                m_elem[irow] * m.m_R.m_elem[icol] +
                m_elem[irow + 1] * m.m_R.m_elem[icol + 3] +
                m_elem[irow + 2] * m.m_R.m_elem[icol + 6];
        }
        mat.m_elem[irow + 3] =
            m_elem[irow] * m.m_p.m_elem[0] +
            m_elem[irow + 1] * m.m_p.m_elem[1] +
            m_elem[irow + 2] * m.m_p.m_elem[2] +
            m_elem[irow + 3];
    }

    return mat;
}

template <typename t_type>
inline Vector<0, t_type> Matrix<0, 0, t_type>::operator*(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_col && "Dimensions do not matched");

    uint16_t cnt;
    uint16_t irow, icol;
    Vector<0, t_type> vec(m_row);

    for (irow = 0; irow < m_row; ++irow)
    {
        for (cnt = m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
        {
            vec.m_elem[irow] += m_elem[irow * m_col + icol] * v.m_elem[icol];
            vec.m_elem[irow] += m_elem[irow * m_col + icol + 1] * v.m_elem[icol + 1];
            vec.m_elem[irow] += m_elem[irow * m_col + icol + 2] * v.m_elem[icol + 2];
            vec.m_elem[irow] += m_elem[irow * m_col + icol + 3] * v.m_elem[icol + 3];
        }

        for (cnt = m_col % 4u; cnt > 0u; cnt--, icol++)
            vec.m_elem[irow] += m_elem[irow * m_col + icol] * v.m_elem[icol];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Matrix<0, 0, t_type>::operator*(const Vector3<t_type, 3> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == 3 && "Dimensions do not matched");

    Vector<0, t_type> vec(m_row);

    for (uint16_t irow = 0; irow < m_row; ++irow)
    {
        vec.m_elem[irow] = m_elem[irow * 3] * v.m_elem[0];
        vec.m_elem[irow] += m_elem[irow * 3 + 1] * v.m_elem[1];
        vec.m_elem[irow] += m_elem[irow * 3 + 2] * v.m_elem[2];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Matrix<0, 0, t_type>::operator*(const Vector4<t_type, 4> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == 4 && "Dimensions do not matched");

    Vector<0, t_type> vec(m_row);

    for (uint16_t irow = 0; irow < m_row; ++irow)
    {
        vec.m_elem[irow] = m_elem[irow * 4] * v.m_elem[0];
        vec.m_elem[irow] += m_elem[irow * 4 + 1] * v.m_elem[1];
        vec.m_elem[irow] += m_elem[irow * 4 + 2] * v.m_elem[2];
        vec.m_elem[irow] += m_elem[irow * 4 + 3] * v.m_elem[3];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Matrix<0, 0, t_type>::operator*(const Vector6<t_type, 6> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == 6 && "Dimensions do not matched");

    Vector<0, t_type> vec(m_row);

    for (uint16_t irow = 0; irow < m_row; ++irow)
    {
        vec.m_elem[irow] = m_elem[irow * 6] * v.m_elem[0];
        vec.m_elem[irow] += m_elem[irow * 6 + 1] * v.m_elem[1];
        vec.m_elem[irow] += m_elem[irow * 6 + 2] * v.m_elem[2];
        vec.m_elem[irow] += m_elem[irow * 6 + 3] * v.m_elem[3];
        vec.m_elem[irow] += m_elem[irow * 6 + 4] * v.m_elem[4];
        vec.m_elem[irow] += m_elem[irow * 6 + 5] * v.m_elem[5];
    }

    return vec;
}

template <typename t_type>
template <uint16_t col>
inline Vector<0, t_type> Matrix<0, 0, t_type>::operator*(const Vector<col, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_col == col && "Dimensions do not matched");

    uint16_t cnt;
    uint16_t irow, icol;
    Vector<0, t_type> vec(m_row);

    for (irow = 0; irow < m_row; ++irow)
    {
        for (cnt = col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
        {
            vec.m_elem[irow] += m_elem[irow * col + icol] * v.m_elem[icol];
            vec.m_elem[irow] += m_elem[irow * col + icol + 1] * v.m_elem[icol + 1];
            vec.m_elem[irow] += m_elem[irow * col + icol + 2] * v.m_elem[icol + 2];
            vec.m_elem[irow] += m_elem[irow * col + icol + 3] * v.m_elem[icol + 3];
        }

        for (cnt = col % 4u; cnt > 0u; cnt--, icol++)
            vec.m_elem[irow] += m_elem[irow * col + icol] * v.m_elem[icol];
    }

    return vec;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator&(const Vector<0, t_type> &v) const
{ // matrix3 * [v]x, []x is skew-symmetric matrix
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "This method is only for 3 x 3 matrix");
    assert(m_col == 3 && "This method is only for 3 x 3 matrix");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == 3 && "This method is only for 3 x 3 matrix");

    Matrix<0, 0, t_type> mat(3, 3);

    mat.m_elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    mat.m_elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    mat.m_elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    mat.m_elem[3] = m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1];
    mat.m_elem[4] = m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2];
    mat.m_elem[5] = m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0];

    mat.m_elem[6] = m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1];
    mat.m_elem[7] = m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2];
    mat.m_elem[8] = m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0];

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator&(const Vector<3, t_type> &v) const
{ // matrix3 * [v]x, []x is skew-symmetric matrix
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "This method is only for 3 x 3 matrix");
    assert(m_col == 3 && "This method is only for 3 x 3 matrix");

    Matrix<0, 0, t_type> mat(3, 3);

    mat.m_elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    mat.m_elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    mat.m_elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    mat.m_elem[3] = m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1];
    mat.m_elem[4] = m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2];
    mat.m_elem[5] = m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0];

    mat.m_elem[6] = m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1];
    mat.m_elem[7] = m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2];
    mat.m_elem[8] = m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0];

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Matrix<0, 0, t_type>::operator&(const Vector3<t_type, 3> &v) const
{ // matrix3 * [v]x, []x is skew-symmetric matrix
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "This method is only for 3 x 3 matrix");
    assert(m_col == 3 && "This method is only for 3 x 3 matrix");

    Matrix<0, 0, t_type> mat(3, 3);

    mat.m_elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    mat.m_elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    mat.m_elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    mat.m_elem[3] = m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1];
    mat.m_elem[4] = m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2];
    mat.m_elem[5] = m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0];

    mat.m_elem[6] = m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1];
    mat.m_elem[7] = m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2];
    mat.m_elem[8] = m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0];

    return mat;
}

/* Comparison operators */
template <typename t_type>
template <uint16_t row, uint16_t col>
inline bool Matrix<0, 0, t_type>::operator==(const Matrix<row, col, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row Dimensions do not matched");
    assert(m_col == col && "Col Dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - m.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 4] - m.m_elem[i + 4]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 6] - m.m_elem[i + 6]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - m.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - m.m_elem[i + 3]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 5] - m.m_elem[i + 5]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 7] - m.m_elem[i + 7]) > m_tolerance) return false;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <typename t_type>
template <uint16_t row, uint16_t col>
inline bool Matrix<0, 0, t_type>::operator!=(const Matrix<row, col, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row Dimensions do not matched");
    assert(m_col == col && "Col Dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - m.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 4] - m.m_elem[i + 4]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 6] - m.m_elem[i + 6]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - m.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - m.m_elem[i + 3]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 5] - m.m_elem[i + 5]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 7] - m.m_elem[i + 7]) > m_tolerance) return true;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <typename t_type>
inline bool Matrix<0, 0, t_type>::operator==(const Matrix<0, 0, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == m.m_row && "Row Dimensions do not matched");
    assert(m_col == m.m_col && "Col Dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - m.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 4] - m.m_elem[i + 4]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 6] - m.m_elem[i + 6]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - m.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - m.m_elem[i + 3]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 5] - m.m_elem[i + 5]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 7] - m.m_elem[i + 7]) > m_tolerance) return false;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <typename t_type>
inline bool Matrix<0, 0, t_type>::operator!=(const Matrix<0, 0, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == m.m_row && "Row Dimensions do not matched");
    assert(m_col == m.m_col && "Col Dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - m.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 4] - m.m_elem[i + 4]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 6] - m.m_elem[i + 6]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - m.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - m.m_elem[i + 3]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 5] - m.m_elem[i + 5]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 7] - m.m_elem[i + 7]) > m_tolerance) return true;
    }

    for (cnt = m_size % 8u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <typename t_type>
inline bool Matrix<0, 0, t_type>::operator==(const Matrix3<t_type, 3, 3> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row Dimensions do not matched");
    assert(m_col == 3 && "Col Dimensions do not matched");

    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance) return false;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance) return false;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance) return false;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance) return false;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance) return false;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance) return false;

    return true;
}

template <typename t_type>
inline bool Matrix<0, 0, t_type>::operator!=(const Matrix3<t_type, 3, 3> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row Dimensions do not matched");
    assert(m_col == 3 && "Col Dimensions do not matched");

    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance) return true;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance) return true;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance) return true;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance) return true;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance) return true;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance) return true;

    return false;
}

template <typename t_type>
inline bool Matrix<0, 0, t_type>::operator==(const Rotation<t_type, 3, 3> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row Dimensions do not matched");
    assert(m_col == 3 && "Col Dimensions do not matched");

    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance) return false;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance) return false;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance) return false;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance) return false;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance) return false;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance) return false;

    return true;
}

template <typename t_type>
inline bool Matrix<0, 0, t_type>::operator!=(const Rotation<t_type, 3, 3> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row Dimensions do not matched");
    assert(m_col == 3 && "Col Dimensions do not matched");

    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance) return true;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance) return true;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance) return true;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance) return true;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance) return true;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance) return true;

    return false;
}

template <typename t_type>
inline bool Matrix<0, 0, t_type>::operator==(const Transform<t_type, 4, 4> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 4 && "Row Dimensions do not matched");
    assert(m_col == 4 && "Col Dimensions do not matched");

    if (std::abs(m_elem[0] - m.m_R.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - m.m_R.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - m.m_R.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - m.m_p.m_elem[0]) > m_tolerance) return false;

    if (std::abs(m_elem[4] - m.m_R.m_elem[3]) > m_tolerance) return false;
    if (std::abs(m_elem[5] - m.m_R.m_elem[4]) > m_tolerance) return false;
    if (std::abs(m_elem[6] - m.m_R.m_elem[5]) > m_tolerance) return false;
    if (std::abs(m_elem[7] - m.m_p.m_elem[1]) > m_tolerance) return false;

    if (std::abs(m_elem[8] - m.m_R.m_elem[6]) > m_tolerance) return false;
    if (std::abs(m_elem[9] - m.m_R.m_elem[7]) > m_tolerance) return false;
    if (std::abs(m_elem[10] - m.m_R.m_elem[8]) > m_tolerance) return false;
    if (std::abs(m_elem[11] - m.m_p.m_elem[2]) > m_tolerance) return false;

    if (std::abs(m_elem[12]) > m_tolerance) return false;
    if (std::abs(m_elem[13]) > m_tolerance) return false;
    if (std::abs(m_elem[14]) > m_tolerance) return false;
    if (std::abs(1 - m_elem[15]) > m_tolerance) return false;

    return true;
}

template <typename t_type>
inline bool Matrix<0, 0, t_type>::operator!=(const Transform<t_type, 4, 4> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 4 && "Row Dimensions do not matched");
    assert(m_col == 4 && "Col Dimensions do not matched");

    if (std::abs(m_elem[0] - m.m_R.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - m.m_R.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - m.m_R.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - m.m_p.m_elem[0]) > m_tolerance) return true;

    if (std::abs(m_elem[4] - m.m_R.m_elem[3]) > m_tolerance) return true;
    if (std::abs(m_elem[5] - m.m_R.m_elem[4]) > m_tolerance) return true;
    if (std::abs(m_elem[6] - m.m_R.m_elem[5]) > m_tolerance) return true;
    if (std::abs(m_elem[7] - m.m_p.m_elem[1]) > m_tolerance) return true;

    if (std::abs(m_elem[8] - m.m_R.m_elem[6]) > m_tolerance) return true;
    if (std::abs(m_elem[9] - m.m_R.m_elem[7]) > m_tolerance) return true;
    if (std::abs(m_elem[10] - m.m_R.m_elem[8]) > m_tolerance) return true;
    if (std::abs(m_elem[11] - m.m_p.m_elem[2]) > m_tolerance) return true;

    if (std::abs(m_elem[12]) > m_tolerance) return true;
    if (std::abs(m_elem[13]) > m_tolerance) return true;
    if (std::abs(m_elem[14]) > m_tolerance) return true;
    if (std::abs(1 - m_elem[15]) > m_tolerance) return true;

    return false;
}

template <typename t_type>
inline void Matrix<0, 0, t_type>::Print(const char endChar)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
#if defined(ARDUINO)
    for (uint16_t irow = 0; irow < m_row; irow++)
    {
        for (uint16_t icol = 0; icol < m_col; icol++)
        {
            Serial.printf("%7.3f ", (t_type)(m_elem[irow * m_col + icol]));
        }
        Serial.write('\n');
    }
    Serial.write(endChar);
#else
    for (uint16_t irow = 0; irow < m_row; irow++)
    {
        for (uint16_t icol = 0; icol < m_col; icol++)
        {
            printf("%10.6f ", (t_type)(m_elem[irow * m_col + icol]));
        }
        printf("\n");
    }
    printf("%c", endChar);
#endif
}

//-- Template Function ------------------------------------------------------//
// scalar * matrix
template <typename type>
inline Matrix<0, 0, type> operator*(const type s, const Matrix<0, 0, type> &m)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    Matrix<0, 0, type> mat(m.m_row, m.m_col);

    for (cnt = m.m_size >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat.m_elem[i] = m.m_elem[i] * s;
        mat.m_elem[i + 2] = m.m_elem[i + 2] * s;
        mat.m_elem[i + 4] = m.m_elem[i + 4] * s;
        mat.m_elem[i + 6] = m.m_elem[i + 6] * s;
        mat.m_elem[i + 1] = m.m_elem[i + 1] * s;
        mat.m_elem[i + 3] = m.m_elem[i + 3] * s;
        mat.m_elem[i + 5] = m.m_elem[i + 5] * s;
        mat.m_elem[i + 7] = m.m_elem[i + 7] * s;
    }

    for (cnt = m.m_size % 8u; cnt > 0u; cnt--, i++)
    {
        mat.m_elem[i] = m.m_elem[i] * s;
    }

    return mat;
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTMATRIX0_TPP_
