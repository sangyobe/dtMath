/*!
\file       dtVector0.tpp
\brief      dtMath, Dynamic Memory Allocation General Vector(m x 1) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 04. 18
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTVECTOR0_TPP_
#define DTMATH_DTVECTOR0_TPP_

#include "dtMathMem.h"
#include "dtVector0.h"

#include <cassert>

namespace dt
{
namespace Math
{

template <typename t_type>
inline Vector<0, t_type>::Vector() : m_row(0), m_elem(nullptr)
{
}

template <typename t_type>
inline Vector<0, t_type>::Vector(const uint16_t row)
    : m_row(row)
{
    m_elem = MemAllocZeroInit<t_type>(m_row);
}

template <typename t_type>
inline Vector<0, t_type>::Vector(const uint16_t row, const t_type *element)
    : m_row(row)
{
    m_elem = MemAlloc<t_type>(m_row);
    memcpy(m_elem, element, sizeof(t_type) * m_row);
}

template <typename t_type>
inline Vector<0, t_type>::Vector(const Vector<0, t_type> &v)
    : m_row(v.m_row)
{
    m_elem = MemAlloc<t_type>(m_row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * m_row);
}

template <typename t_type>
inline Vector<0, t_type>::Vector(const Matrix<0, 0, t_type> &v)
    : m_row(v.m_row)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    m_elem = MemAlloc<t_type>(m_row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * m_row);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type>::Vector(const Vector<row, t_type> &v)
    : m_row(row)
{
    m_elem = MemAlloc<t_type>(row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * row);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type>::Vector(const Matrix<row, 1, t_type> &v)
    : m_row(row)
{
    m_elem = MemAlloc<t_type>(row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * row);
}

template <typename t_type>
inline Vector<0, t_type>::~Vector()
{
    if (m_elem)
    {
        MemFree<t_type>(m_elem);
        m_elem = nullptr;
    }
}

template <typename t_type>
inline void Vector<0, t_type>::NewSize(const uint16_t row)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    if (m_elem) Release();

    m_row = row;
    m_elem = MemAllocZeroInit<t_type>(m_row);
}

template <typename t_type>
inline void Vector<0, t_type>::NewSize(const uint16_t row, const t_type *element)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    if (m_elem) Release();

    m_row = row;
    m_elem = MemAlloc<t_type>(m_row);
    memcpy(m_elem, element, sizeof(t_type) * m_row);
}

template <typename t_type>
inline void Vector<0, t_type>::NewSize(const Vector<0, t_type> &v)
{
    assert(m_elem == nullptr && "Memory has been allocated");
    assert(v.m_elem != nullptr && "Argument Vector has not been allocated");

    if (m_elem) Release();

    m_row = v.m_row;
    m_elem = MemAlloc<t_type>(m_row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * m_row);
}

template <typename t_type>
inline void Vector<0, t_type>::NewSize(const Matrix<0, 0, t_type> &v)
{
    assert(m_elem == nullptr && "Memory has been allocated");
    assert(v.m_elem != nullptr && "Argument Vector has not been allocated");
    assert(v.m_col == 1 && "Column size is not 1");

    if (m_elem) Release();

    m_row = v.m_row;
    m_elem = MemAlloc<t_type>(m_row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * m_row);
}

template <typename t_type>
template <uint16_t row>
inline void Vector<0, t_type>::NewSize(const Vector<row, t_type> &v)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    if (m_elem) Release();

    m_row = row;
    m_elem = MemAlloc<t_type>(row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * row);
}

template <typename t_type>
template <uint16_t row>
inline void Vector<0, t_type>::NewSize(const Matrix<row, 1, t_type> &v)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    if (m_elem) Release();

    m_row = row;
    m_elem = MemAlloc<t_type>(row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * row);
}

template <typename t_type>
inline void Vector<0, t_type>::ReSize(const uint16_t row)
{
    assert(m_elem != nullptr && "Memory has been allocated");

    m_row = row;
    m_elem = MemReAlloc<t_type>(m_elem, m_row);
}

template <typename t_type>
inline void Vector<0, t_type>::ReSize(const uint16_t row, const t_type *element)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(element != nullptr && "Argument is empty");

    m_row = row;
    m_elem = MemReAlloc<t_type>(m_elem, m_row);
    memcpy(m_elem, element, sizeof(t_type) * m_row);
}

template <typename t_type>
inline void Vector<0, t_type>::ReSize(const Vector<0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has been allocated");
    assert(v.m_elem != nullptr && "Argument Vector has not been allocated");

    m_row = v.m_row;
    m_elem = MemReAlloc<t_type>(m_elem, m_row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * m_row);
}

template <typename t_type>
inline void Vector<0, t_type>::ReSize(const Matrix<0, 0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Argument Vector has not been allocated");
    assert(v.m_col == 1 && "Column size is not 1");

    m_row = v.m_row;
    m_elem = MemReAlloc<t_type>(m_elem, m_row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * m_row);
}

template <typename t_type>
template <uint16_t row>
inline void Vector<0, t_type>::ReSize(const Vector<row, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    m_row = row;
    m_elem = MemReAlloc<t_type>(m_elem, row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * row);
}

template <typename t_type>
template <uint16_t row>
inline void Vector<0, t_type>::ReSize(const Matrix<row, 1, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    m_row = row;
    m_elem = MemReAlloc<t_type>(m_elem, row);
    memcpy(m_elem, v.m_elem, sizeof(t_type) * row);
}

template <typename t_type>
inline void Vector<0, t_type>::Release()
{
    if (m_elem)
    {
        MemFree<t_type>(m_elem);
        m_elem = nullptr;
        m_row = 0;
    }
}

template <typename t_type>
inline void Vector<0, t_type>::SetZero()
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    memset(m_elem, 0, sizeof(t_type) * m_row);
}

template <typename t_type>
inline void Vector<0, t_type>::SetFill(const t_type value)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] = value;
        m_elem[irow + 1] = value;
        m_elem[irow + 2] = value;
        m_elem[irow + 3] = value;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
        m_elem[irow] = value;
}

template <typename t_type>
inline void Vector<0, t_type>::SetElement(const t_type *element, const size_t n_byte)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(element != nullptr && "element is nullptr");

    size_t vecSz = sizeof(t_type) * m_row;

    if (vecSz > n_byte) memcpy(m_elem, element, n_byte);
    else memcpy(m_elem, element, vecSz);
}

template <typename t_type>
inline void Vector<0, t_type>::SetElement(const t_type *element, const uint16_t row)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(element != nullptr && "element is nullptr");

    if (m_row > row) memcpy(m_elem, element, sizeof(t_type) * row);
    else memcpy(m_elem, element, sizeof(t_type) * m_row);
}

template <typename t_type>
inline void Vector<0, t_type>::SetBlock(const uint16_t idxRow, const Vector<0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");

    if (idxRow >= m_row) return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > v.m_row) rowSz = v.m_row;

    memcpy(&m_elem[idxRow], v.m_elem, sizeof(t_type) * rowSz);
}

template <typename t_type>
inline void Vector<0, t_type>::SetBlock(const uint16_t idxRow, const Vector3<t_type, 3> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");

    if (idxRow >= m_row) return;

    switch (m_row - idxRow)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[0];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        break;
    default:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    }
}

template <typename t_type>
inline void Vector<0, t_type>::SetBlock(const uint16_t idxRow, const Vector4<t_type, 4> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");

    if (idxRow >= m_row) return;

    switch (m_row - idxRow)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[0];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    default:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        break;
    }
}

template <typename t_type>
inline void Vector<0, t_type>::SetBlock(const uint16_t idxRow, const Vector6<t_type, 6> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");

    if (idxRow >= m_row) return;

    switch (m_row - idxRow)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[0];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    case 4:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        break;
    case 5:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        m_elem[idxRow + 4] = v.m_elem[4];
        break;
    default:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        m_elem[idxRow + 4] = v.m_elem[4];
        m_elem[idxRow + 5] = v.m_elem[5];
        break;
    }
}

template <typename t_type>
inline void Vector<0, t_type>::SetBlock(const uint16_t idxRow, const Matrix<0, 0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    if (idxRow >= m_row) return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > v.m_row) rowSz = v.m_row;

    memcpy(&m_elem[idxRow], v.m_elem, sizeof(t_type) * rowSz);
}

template <typename t_type>
template <uint16_t row>
inline void Vector<0, t_type>::SetBlock(const uint16_t idxRow, const Vector<row, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");

    if (idxRow >= m_row) return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > row) rowSz = row;

    memcpy(&m_elem[idxRow], v.m_elem, sizeof(t_type) * rowSz);
}

template <typename t_type>
template <uint16_t row>
inline void Vector<0, t_type>::SetBlock(const uint16_t idxRow, const Matrix<row, 1, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idxRow && "Index out of range");

    if (idxRow >= m_row) return;

    uint16_t rowSz = m_row - idxRow;
    if (rowSz > row) rowSz = row;

    memcpy(&m_elem[idxRow], v.m_elem, sizeof(t_type) * rowSz);
}

template <typename t_type>
inline void Vector<0, t_type>::SetSwap(const uint16_t i, const uint16_t j)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > i && "Index out of range");
    assert(m_row > j && "Index out of range");

    t_type elem = m_elem[i];
    m_elem[i] = m_elem[j];
    m_elem[j] = elem;
}

template <typename t_type>
inline void Vector<0, t_type>::SetNormalize()
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    uint16_t cnt = 0;
    uint16_t irow = 0;
    t_type norm = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        norm += m_elem[irow] * m_elem[irow];
        norm += m_elem[irow + 1] * m_elem[irow + 1];
        norm += m_elem[irow + 2] * m_elem[irow + 2];
        norm += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        norm += m_elem[irow] * m_elem[irow];
    }

    norm = std::sqrt(norm);

    if (norm < std::numeric_limits<t_type>::epsilon())
        norm = std::numeric_limits<t_type>::epsilon();

    for (cnt = m_row >> 2u, irow = 0; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] /= norm;
        m_elem[irow + 1] /= norm;
        m_elem[irow + 2] /= norm;
        m_elem[irow + 3] /= norm;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
        m_elem[irow] /= norm;
}

template <typename t_type>
inline const t_type *const Vector<0, t_type>::GetElementsAddr() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    return m_elem;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::GetBlock(const uint16_t idx, const uint16_t row)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idx && "Index out of range");

    Vector<0, t_type> v(row);
    uint16_t rowSize = m_row - idx;

    if (rowSize > row) rowSize = row;

    memcpy(v.m_elem, &m_elem[idx], sizeof(t_type) * rowSize);

    return v;
}

template <typename t_type>
inline Vector3<t_type, 3> Vector<0, t_type>::GetBlockVec3(const uint16_t idx)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idx && "Index out of range");

    t_type elem[3]{0};

    if (idx >= m_row) return Vector3<t_type, 3>(elem);

    switch (m_row - idx)
    {
    case 1:
        elem[0] = m_elem[idx];
        break;
    case 2:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        break;
    default:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
    };

    return Vector3<t_type, 3>(elem);
}

template <typename t_type>
inline Vector4<t_type, 4> Vector<0, t_type>::GetBlockVec4(const uint16_t idx)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idx && "Index out of range");

    t_type elem[4]{0};

    if (idx >= m_row) return Vector4<t_type, 4>(elem);

    switch (m_row - idx)
    {
    case 1:
        elem[0] = m_elem[idx];
        break;
    case 2:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        break;
    case 3:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        break;
    default:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        elem[3] = m_elem[idx + 3];
    };

    return Vector4<t_type, 4>(elem);
}

template <typename t_type>
inline Vector6<t_type, 6> Vector<0, t_type>::GetBlockVec6(const uint16_t idx)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idx && "Index out of range");

    t_type elem[6]{0};

    if (idx >= m_row) return Vector6<t_type, 6>(elem);

    switch (m_row - idx)
    {
    case 1:
        elem[0] = m_elem[idx];
        break;
    case 2:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        break;
    case 3:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        break;
    case 4:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        elem[3] = m_elem[idx + 3];
        break;
    case 5:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        elem[3] = m_elem[idx + 3];
        elem[4] = m_elem[idx + 4];
        break;
    default:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        elem[3] = m_elem[idx + 3];
        elem[4] = m_elem[idx + 4];
        elem[5] = m_elem[idx + 5];
    };

    return Vector6<t_type, 6>(elem);
}

template <typename t_type>
template <uint16_t row>
inline Vector<row, t_type> Vector<0, t_type>::GetBlock(const uint16_t idx)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idx && "Index out of range");

    t_type elem[row]{0};
    uint16_t rowSize = m_row - idx;

    if (idx >= m_row) return Vector<row, t_type>(elem); // Index out of range, return zero vec
    if (rowSize > row) rowSize = row;

    memcpy(elem, &m_elem[idx], sizeof(t_type) * rowSize);

    return Vector<row, t_type>(elem);
}

template <typename t_type>
inline int8_t Vector<0, t_type>::GetBlock(const uint16_t idx, Vector<0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idx && "Index out of range");

    uint16_t rowSize = m_row - idx;

    if (idx >= m_row) return -1;
    if (rowSize > v.m_row) rowSize = v.m_row;

    memcpy(v.m_elem, &m_elem[idx], sizeof(t_type) * rowSize);

    return 0;
}

template <typename t_type>
inline int8_t Vector<0, t_type>::GetBlockVec3(const uint16_t idx, Vector3<t_type, 3> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idx && "Index out of range");

    if (idx >= m_row) return -1;

    switch (m_row - idx)
    {
    case 1:
        v.m_elem[0] = m_elem[idx];
        break;
    case 2:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        break;
    default:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
    };

    return 0;
}

template <typename t_type>
inline int8_t Vector<0, t_type>::GetBlockVec4(const uint16_t idx, Vector4<t_type, 4> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idx && "Index out of range");

    if (idx >= m_row) return -1;

    switch (m_row - idx)
    {
    case 1:
        v.m_elem[0] = m_elem[idx];
        break;
    case 2:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        break;
    case 3:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        break;
    default:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        v.m_elem[3] = m_elem[idx + 3];
    };

    return 0;
}

template <typename t_type>
inline int8_t Vector<0, t_type>::GetBlockVec6(const uint16_t idx, Vector6<t_type, 6> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idx && "Index out of range");

    if (idx >= m_row) return -1;

    switch (m_row - idx)
    {
    case 1:
        v.m_elem[0] = m_elem[idx];
        break;
    case 2:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        break;
    case 3:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        break;
    case 4:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        v.m_elem[3] = m_elem[idx + 3];
        break;
    case 5:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        v.m_elem[3] = m_elem[idx + 3];
        v.m_elem[4] = m_elem[idx + 4];
        break;
    default:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        v.m_elem[3] = m_elem[idx + 3];
        v.m_elem[4] = m_elem[idx + 4];
        v.m_elem[5] = m_elem[idx + 5];
    };

    return 0;
}

template <typename t_type>
template <uint16_t row>
inline int8_t Vector<0, t_type>::GetBlock(const uint16_t idx, Vector<row, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > idx && "Index out of range");

    uint16_t rowSize = m_row - idx;

    if (idx >= m_row) return -1;
    if (rowSize > row) rowSize = row;

    memcpy(v.m_elem, &m_elem[idx], sizeof(t_type) * rowSize);

    return 0;
}

template <typename t_type>
inline t_type Vector<0, t_type>::GetNorm() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    t_type sqSum = 0;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        sqSum += m_elem[irow] * m_elem[irow];
        sqSum += m_elem[irow + 1] * m_elem[irow + 1];
        sqSum += m_elem[irow + 2] * m_elem[irow + 2];
        sqSum += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        sqSum += m_elem[irow] * m_elem[irow];
    }

    return std::sqrt(sqSum);
}

template <typename t_type>
inline t_type Vector<0, t_type>::GetSqNorm() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    t_type sqSum = 0;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        sqSum += m_elem[irow] * m_elem[irow];
        sqSum += m_elem[irow + 1] * m_elem[irow + 1];
        sqSum += m_elem[irow + 2] * m_elem[irow + 2];
        sqSum += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        sqSum += m_elem[irow] * m_elem[irow];
    }

    return sqSum;
}

template <typename t_type>
inline t_type Vector<0, t_type>::GetLpNorm(const int p) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, i = 0;
    t_type powSum = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        powSum += std::pow(std::abs(m_elem[i]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 1]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 2]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 3]), (t_type)p);
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, i++)
    {
        powSum += std::pow(std::abs(m_elem[i]), (t_type)p);
    }

    return std::pow(powSum, (t_type)1 / p);
}

template <typename t_type>
inline t_type Vector<0, t_type>::GetSum() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    t_type rtn = 0;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        rtn += m_elem[irow];
        rtn += m_elem[irow + 1];
        rtn += m_elem[irow + 2];
        rtn += m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        rtn += m_elem[irow];
    }

    return rtn;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::GetNormalized() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt = 0;
    uint16_t irow = 0;
    t_type norm = 0;
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        norm += m_elem[irow] * m_elem[irow];
        norm += m_elem[irow + 1] * m_elem[irow + 1];
        norm += m_elem[irow + 2] * m_elem[irow + 2];
        norm += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        norm += m_elem[irow] * m_elem[irow];
    }

    norm = std::sqrt(norm);

    if (norm < std::numeric_limits<t_type>::epsilon())
        norm = std::numeric_limits<t_type>::epsilon();

    for (cnt = m_row >> 2u, irow = 0; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] / norm;
        vec.m_elem[irow + 1] = m_elem[irow + 1] / norm;
        vec.m_elem[irow + 2] = m_elem[irow + 2] / norm;
        vec.m_elem[irow + 3] = m_elem[irow + 3] / norm;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
        vec.m_elem[irow] = m_elem[irow] / norm;

    return vec;
}

template <typename t_type>
inline Matrix<3, 3, t_type> Vector<0, t_type>::GetSkew() const
{
    static_assert(m_row == 3, "This method is only for 3 x 1 vector");
    assert(m_elem != nullptr && "Memory has not been allocated");

    return Matrix<3, 3, t_type>(
        0, -m_elem[2], m_elem[1],
        m_elem[2], 0, -m_elem[0],
        -m_elem[1], m_elem[0], 0);
}

template <typename t_type>
inline Matrix<0, 0, t_type> Vector<0, t_type>::Transpose() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    return Matrix<0, 0, t_type>(1, m_row, m_elem);
}

template <typename t_type>
inline void Vector<0, t_type>::Transpose(Matrix<0, 0, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == 1 && "Row dimensions do not matched");
    assert(m.m_col == m_row && "Col dimensions do not matched");

    memcpy(m.m_elem, m_elem, sizeof(t_type) * m_row);
}

template <typename t_type>
template <uint16_t col>
inline void Vector<0, t_type>::Transpose(Matrix<1, col, t_type> &m) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == col && "Col dimensions do not matched");
    memcpy(m.m_elem, m_elem, sizeof(t_type) * m_row);
}

/* Member access operators */
template <typename t_type>
inline t_type &Vector<0, t_type>::operator()(uint16_t irow)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > irow && "Index out of range");
    return m_elem[irow];
}

template <typename t_type>
inline const t_type &Vector<0, t_type>::operator()(uint16_t irow) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row > irow && "Index out of range");
    return m_elem[irow];
}

/* Assignment operators */
template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator=(const Vector<0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    memcpy(m_elem, v.m_elem, sizeof(t_type) * m_row);

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator+=(const Vector<0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += v.m_elem[irow];
        m_elem[irow + 1] += v.m_elem[irow + 1];
        m_elem[irow + 2] += v.m_elem[irow + 2];
        m_elem[irow + 3] += v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator-=(const Vector<0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= v.m_elem[irow];
        m_elem[irow + 1] -= v.m_elem[irow + 1];
        m_elem[irow + 2] -= v.m_elem[irow + 2];
        m_elem[irow + 3] -= v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator*=(const Vector<0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= v.m_elem[irow];
        m_elem[irow + 1] *= v.m_elem[irow + 1];
        m_elem[irow + 2] *= v.m_elem[irow + 2];
        m_elem[irow + 3] *= v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator/=(const Vector<0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    t_type den[4];

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        den[0] = v.m_elem[irow];
        den[1] = v.m_elem[irow + 1];
        den[2] = v.m_elem[irow + 2];
        den[3] = v.m_elem[irow + 3];

        if (std::abs(den[0]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[0] < 0) den[0] = -std::numeric_limits<t_type>::epsilon();
            else den[0] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[1]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[1] < 0) den[1] = -std::numeric_limits<t_type>::epsilon();
            else den[1] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[2]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[2] < 0) den[2] = -std::numeric_limits<t_type>::epsilon();
            else den[2] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[3]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[3] < 0) den[3] = -std::numeric_limits<t_type>::epsilon();
            else den[3] = std::numeric_limits<t_type>::epsilon();
        }
        m_elem[irow] /= den[0];
        m_elem[irow + 1] /= den[1];
        m_elem[irow + 2] /= den[2];
        m_elem[irow + 3] /= den[3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        den[0] = v.m_elem[irow];
        if (std::abs(den[0]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[0] < 0) den[0] = -std::numeric_limits<t_type>::epsilon();
            else den[0] = std::numeric_limits<t_type>::epsilon();
        }
        m_elem[irow] /= den[0];
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator=(const Vector3<t_type, 3> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator+=(const Vector3<t_type, 3> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator-=(const Vector3<t_type, 3> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator*=(const Vector3<t_type, 3> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator/=(const Vector3<t_type, 3> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    t_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[0] /= den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[1] /= den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[2] /= den;

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator=(const Vector4<t_type, 4> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator+=(const Vector4<t_type, 4> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator-=(const Vector4<t_type, 4> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator*=(const Vector4<t_type, 4> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator/=(const Vector4<t_type, 4> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    t_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[0] /= den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[1] /= den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[2] /= den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[3] /= den;

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator=(const Vector6<t_type, 6> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator+=(const Vector6<t_type, 6> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];
    m_elem[4] += v.m_elem[4];
    m_elem[5] += v.m_elem[5];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator-=(const Vector6<t_type, 6> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];
    m_elem[4] -= v.m_elem[4];
    m_elem[5] -= v.m_elem[5];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator*=(const Vector6<t_type, 6> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];
    m_elem[4] *= v.m_elem[4];
    m_elem[5] *= v.m_elem[5];

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator/=(const Vector6<t_type, 6> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    t_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[0] /= den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[1] /= den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[2] /= den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[3] /= den;

    den = v.m_elem[4];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[4] /= den;

    den = v.m_elem[5];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    m_elem[5] /= den;

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator=(const Matrix<0, 0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    memcpy(m_elem, v.m_elem, sizeof(t_type) * m_row);

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator+=(const Matrix<0, 0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += v.m_elem[irow];
        m_elem[irow + 1] += v.m_elem[irow + 1];
        m_elem[irow + 2] += v.m_elem[irow + 2];
        m_elem[irow + 3] += v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator-=(const Matrix<0, 0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= v.m_elem[irow];
        m_elem[irow + 1] -= v.m_elem[irow + 1];
        m_elem[irow + 2] -= v.m_elem[irow + 2];
        m_elem[irow + 3] -= v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator*=(const Matrix<0, 0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= v.m_elem[irow];
        m_elem[irow + 1] *= v.m_elem[irow + 1];
        m_elem[irow + 2] *= v.m_elem[irow + 2];
        m_elem[irow + 3] *= v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator/=(const Matrix<0, 0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    t_type den[4];
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        den[0] = v.m_elem[irow];
        den[1] = v.m_elem[irow + 1];
        den[2] = v.m_elem[irow + 2];
        den[3] = v.m_elem[irow + 3];

        if (std::abs(den[0]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[0] < 0) den[0] = -std::numeric_limits<t_type>::epsilon();
            else den[0] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[1]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[1] < 0) den[1] = -std::numeric_limits<t_type>::epsilon();
            else den[1] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[2]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[2] < 0) den[2] = -std::numeric_limits<t_type>::epsilon();
            else den[2] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[3]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[3] < 0) den[3] = -std::numeric_limits<t_type>::epsilon();
            else den[3] = std::numeric_limits<t_type>::epsilon();
        }
        m_elem[irow] /= den[0];
        m_elem[irow + 1] /= den[1];
        m_elem[irow + 2] /= den[2];
        m_elem[irow + 3] /= den[3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        den[0] = v.m_elem[irow];
        if (std::abs(den[0]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[0] < 0) den[0] = -std::numeric_limits<t_type>::epsilon();
            else den[0] = std::numeric_limits<t_type>::epsilon();
        }
        m_elem[irow] /= den[0];
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator=(const t_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] = s;
        m_elem[irow + 1] = s;
        m_elem[irow + 2] = s;
        m_elem[irow + 3] = s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] = s;
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator+=(const t_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += s;
        m_elem[irow + 1] += s;
        m_elem[irow + 2] += s;
        m_elem[irow + 3] += s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += s;
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator-=(const t_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= s;
        m_elem[irow + 1] -= s;
        m_elem[irow + 2] -= s;
        m_elem[irow + 3] -= s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= s;
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator*=(const t_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= s;
        m_elem[irow + 1] *= s;
        m_elem[irow + 2] *= s;
        m_elem[irow + 3] *= s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= s;
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator/=(const t_type s)
{
    t_type scalar = s;
    uint16_t cnt, irow = 0;

    if (std::abs(scalar) < std::numeric_limits<t_type>::epsilon())
    {
        if (scalar < 0) scalar = -std::numeric_limits<t_type>::epsilon();
        else scalar = std::numeric_limits<t_type>::epsilon();
    }

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] /= scalar;
        m_elem[irow + 1] /= scalar;
        m_elem[irow + 2] /= scalar;
        m_elem[irow + 3] /= scalar;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] /= scalar;
    }

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator=(const Vector<row, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    memcpy(m_elem, v.m_elem, sizeof(t_type) * row);

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator+=(const Vector<row, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += v.m_elem[irow];
        m_elem[irow + 1] += v.m_elem[irow + 1];
        m_elem[irow + 2] += v.m_elem[irow + 2];
        m_elem[irow + 3] += v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator-=(const Vector<row, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= v.m_elem[irow];
        m_elem[irow + 1] -= v.m_elem[irow + 1];
        m_elem[irow + 2] -= v.m_elem[irow + 2];
        m_elem[irow + 3] -= v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator*=(const Vector<row, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= v.m_elem[irow];
        m_elem[irow + 1] *= v.m_elem[irow + 1];
        m_elem[irow + 2] *= v.m_elem[irow + 2];
        m_elem[irow + 3] *= v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator/=(const Vector<row, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    t_type den[row]{0};
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 1]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 1] < 0) den[irow + 1] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 1] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 2]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 2] < 0) den[irow + 2] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 2] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 3]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 3] < 0) den[irow + 3] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 3] = std::numeric_limits<t_type>::epsilon();
        }
        m_elem[irow] /= den[irow];
        m_elem[irow + 1] /= den[irow + 1];
        m_elem[irow + 2] /= den[irow + 2];
        m_elem[irow + 3] /= den[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        m_elem[irow] /= den[irow];
    }

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator=(const Matrix<row, 1, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    memcpy(m_elem, v.m_elem, sizeof(t_type) * row);

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator+=(const Matrix<row, 1, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += v.m_elem[irow];
        m_elem[irow + 1] += v.m_elem[irow + 1];
        m_elem[irow + 2] += v.m_elem[irow + 2];
        m_elem[irow + 3] += v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator-=(const Matrix<row, 1, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= v.m_elem[irow];
        m_elem[irow + 1] -= v.m_elem[irow + 1];
        m_elem[irow + 2] -= v.m_elem[irow + 2];
        m_elem[irow + 3] -= v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator*=(const Matrix<row, 1, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= v.m_elem[irow];
        m_elem[irow + 1] *= v.m_elem[irow + 1];
        m_elem[irow + 2] *= v.m_elem[irow + 2];
        m_elem[irow + 3] *= v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= v.m_elem[irow];
    }

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator/=(const Matrix<row, 1, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    t_type den[row]{0};
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 1]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 1] < 0) den[irow + 1] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 1] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 2]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 2] < 0) den[irow + 2] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 2] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 3]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 3] < 0) den[irow + 3] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 3] = std::numeric_limits<t_type>::epsilon();
        }
        m_elem[irow] /= den[irow];
        m_elem[irow + 1] /= den[irow + 1];
        m_elem[irow + 2] /= den[irow + 2];
        m_elem[irow + 3] /= den[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        m_elem[irow] /= den[irow];
    }

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator&=(const Vector<0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == 3 && "Row dimensions do not matched");

    CrossProduct(v.m_elem);

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator&=(const Vector3<t_type, 3> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");

    CrossProduct(v.m_elem);

    return (*this);
}

template <typename t_type>
inline Vector<0, t_type> &Vector<0, t_type>::operator&=(const Matrix<0, 0, t_type> &v)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == 3 && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    CrossProduct(v.m_elem);

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator&=(const Vector<row, t_type> &v)
{
    static_assert(row == 3, "Row dimensions do not matched");
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");

    CrossProduct(v.m_elem);

    return (*this);
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> &Vector<0, t_type>::operator&=(const Matrix<row, 1, t_type> &v)
{
    static_assert(row == 3, "Row dimensions do not matched");
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");

    CrossProduct(v.m_elem);

    return (*this);
}

template <typename t_type>
inline CommaInit<0, t_type> Vector<0, t_type>::operator<<(const t_type s)
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    m_elem[0] = s;
    return CommaInit<0, t_type>(m_elem, m_row);
}

/* Arithmetic operators */
template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator-() const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = -m_elem[irow];
        vec.m_elem[irow + 1] = -m_elem[irow + 1];
        vec.m_elem[irow + 2] = -m_elem[irow + 2];
        vec.m_elem[irow + 3] = -m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = -m_elem[irow];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator+(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] + v.m_elem[irow];
        vec.m_elem[irow + 1] = m_elem[irow + 1] + v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] + v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] + v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] + v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator-(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] - v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] - v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] - v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator*(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] * v.m_elem[irow];
        vec.m_elem[irow + 1] = m_elem[irow + 1] * v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] * v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] * v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator/(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    t_type den[4];
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        den[0] = v.m_elem[irow];
        den[1] = v.m_elem[irow + 1];
        den[2] = v.m_elem[irow + 2];
        den[3] = v.m_elem[irow + 3];

        if (std::abs(den[0]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[0] < 0) den[0] = -std::numeric_limits<t_type>::epsilon();
            else den[0] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[1]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[1] < 0) den[1] = -std::numeric_limits<t_type>::epsilon();
            else den[1] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[2]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[2] < 0) den[2] = -std::numeric_limits<t_type>::epsilon();
            else den[2] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[3]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[3] < 0) den[3] = -std::numeric_limits<t_type>::epsilon();
            else den[3] = std::numeric_limits<t_type>::epsilon();
        }
        vec.m_elem[irow] = m_elem[irow] / den[0];
        vec.m_elem[irow + 1] = m_elem[irow + 1] / den[1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] / den[2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] / den[3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        den[0] = v.m_elem[irow];
        if (std::abs(den[0]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[0] < 0) den[0] = -std::numeric_limits<t_type>::epsilon();
            else den[0] = std::numeric_limits<t_type>::epsilon();
        }
        vec.m_elem[irow] = m_elem[irow] / den[0];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator+(const Vector3<t_type, 3> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    Vector<0, t_type> vec(m_row);

    vec.m_elem[0] = m_elem[0] + v.m_elem[0];
    vec.m_elem[1] = m_elem[1] + v.m_elem[1];
    vec.m_elem[2] = m_elem[2] + v.m_elem[2];

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator-(const Vector3<t_type, 3> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    Vector<0, t_type> vec(m_row);

    vec[0] = m_elem[0] - v.m_elem[0];
    vec[1] = m_elem[1] - v.m_elem[1];
    vec[2] = m_elem[2] - v.m_elem[2];

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator*(const Vector3<t_type, 3> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    Vector<0, t_type> vec(m_row);

    vec[0] = m_elem[0] * v.m_elem[0];
    vec[1] = m_elem[1] * v.m_elem[1];
    vec[2] = m_elem[2] * v.m_elem[2];

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator/(const Vector3<t_type, 3> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    t_type den;
    Vector<0, t_type> vec(m_row);

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[0] = m_elem[0] / den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[1] = m_elem[1] / den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[2] = m_elem[2] / den;

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator+(const Vector4<t_type, 4> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    Vector<0, t_type> vec(m_row);

    vec.m_elem[0] = m_elem[0] + v.m_elem[0];
    vec.m_elem[1] = m_elem[1] + v.m_elem[1];
    vec.m_elem[2] = m_elem[2] + v.m_elem[2];
    vec.m_elem[3] = m_elem[3] + v.m_elem[3];

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator-(const Vector4<t_type, 4> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    Vector<0, t_type> vec(m_row);

    vec.m_elem[0] = m_elem[0] - v.m_elem[0];
    vec.m_elem[1] = m_elem[1] - v.m_elem[1];
    vec.m_elem[2] = m_elem[2] - v.m_elem[2];
    vec.m_elem[3] = m_elem[3] - v.m_elem[3];

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator*(const Vector4<t_type, 4> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    Vector<0, t_type> vec(m_row);

    vec.m_elem[0] = m_elem[0] * v.m_elem[0];
    vec.m_elem[1] = m_elem[1] * v.m_elem[1];
    vec.m_elem[2] = m_elem[2] * v.m_elem[2];
    vec.m_elem[3] = m_elem[3] * v.m_elem[3];

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator/(const Vector4<t_type, 4> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    t_type den;
    Vector<0, t_type> vec(m_row);

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[0] = m_elem[0] / den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[1] = m_elem[1] / den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[2] = m_elem[2] / den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[3] = m_elem[3] / den;

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator+(const Vector6<t_type, 6> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    Vector<0, t_type> vec(m_row);

    vec.m_elem[0] = m_elem[0] + v.m_elem[0];
    vec.m_elem[1] = m_elem[1] + v.m_elem[1];
    vec.m_elem[2] = m_elem[2] + v.m_elem[2];
    vec.m_elem[3] = m_elem[3] + v.m_elem[3];
    vec.m_elem[4] = m_elem[4] + v.m_elem[4];
    vec.m_elem[5] = m_elem[5] + v.m_elem[5];

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator-(const Vector6<t_type, 6> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    Vector<0, t_type> vec(m_row);

    vec.m_elem[0] = m_elem[0] - v.m_elem[0];
    vec.m_elem[1] = m_elem[1] - v.m_elem[1];
    vec.m_elem[2] = m_elem[2] - v.m_elem[2];
    vec.m_elem[3] = m_elem[3] - v.m_elem[3];
    vec.m_elem[4] = m_elem[4] - v.m_elem[4];
    vec.m_elem[5] = m_elem[5] - v.m_elem[5];

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator*(const Vector6<t_type, 6> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    Vector<0, t_type> vec(m_row);

    vec.m_elem[0] = m_elem[0] * v.m_elem[0];
    vec.m_elem[1] = m_elem[1] * v.m_elem[1];
    vec.m_elem[2] = m_elem[2] * v.m_elem[2];
    vec.m_elem[3] = m_elem[3] * v.m_elem[3];
    vec.m_elem[4] = m_elem[4] * v.m_elem[4];
    vec.m_elem[5] = m_elem[5] * v.m_elem[5];

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator/(const Vector6<t_type, 6> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    t_type den;
    Vector<0, t_type> vec(m_row);

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[0] = m_elem[0] / den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[1] = m_elem[1] / den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[2] = m_elem[2] / den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[3] = m_elem[3] / den;

    den = v.m_elem[4];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[4] = m_elem[4] / den;

    den = v.m_elem[5];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec.m_elem[5] = m_elem[5] / den;

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator+(const Matrix<0, 0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] + v.m_elem[irow];
        vec.m_elem[irow + 1] = m_elem[irow + 1] + v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] + v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] + v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] + v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator-(const Matrix<0, 0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] - v.m_elem[irow];
        vec.m_elem[irow + 1] = m_elem[irow + 1] - v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] - v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] - v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] - v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator*(const Matrix<0, 0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] * v.m_elem[irow];
        vec.m_elem[irow + 1] = m_elem[irow + 1] * v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] * v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] * v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator/(const Matrix<0, 0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, irow = 0;
    t_type den[4];
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        den[0] = v.m_elem[irow];
        den[1] = v.m_elem[irow + 1];
        den[2] = v.m_elem[irow + 2];
        den[3] = v.m_elem[irow + 3];

        if (std::abs(den[0]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[0] < 0) den[0] = -std::numeric_limits<t_type>::epsilon();
            else den[0] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[1]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[1] < 0) den[1] = -std::numeric_limits<t_type>::epsilon();
            else den[1] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[2]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[2] < 0) den[2] = -std::numeric_limits<t_type>::epsilon();
            else den[2] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[3]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[3] < 0) den[3] = -std::numeric_limits<t_type>::epsilon();
            else den[3] = std::numeric_limits<t_type>::epsilon();
        }
        vec.m_elem[irow] = m_elem[irow] / den[0];
        vec.m_elem[irow + 1] = m_elem[irow + 1] / den[1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] / den[2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] / den[3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        den[0] = v.m_elem[irow];
        if (std::abs(den[0]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[0] < 0) den[0] = -std::numeric_limits<t_type>::epsilon();
            else den[0] = std::numeric_limits<t_type>::epsilon();
        }
        vec.m_elem[irow] = m_elem[irow] / den[0];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator+(const t_type s) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] + s;
        vec.m_elem[irow + 1] = m_elem[irow + 1] + s;
        vec.m_elem[irow + 2] = m_elem[irow + 2] + s;
        vec.m_elem[irow + 3] = m_elem[irow + 3] + s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] + s;
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator-(const t_type s) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] - s;
        vec.m_elem[irow + 1] = m_elem[irow + 1] - s;
        vec.m_elem[irow + 2] = m_elem[irow + 2] - s;
        vec.m_elem[irow + 3] = m_elem[irow + 3] - s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] - s;
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator*(const t_type s) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(m_row);

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] * s;
        vec.m_elem[irow + 1] = m_elem[irow + 1] * s;
        vec.m_elem[irow + 2] = m_elem[irow + 2] * s;
        vec.m_elem[irow + 3] = m_elem[irow + 3] * s;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] * s;
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator/(const t_type s) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    t_type scalar = s;
    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(m_row);

    if (std::abs(scalar) < std::numeric_limits<t_type>::epsilon())
    {
        if (scalar < 0) scalar = -std::numeric_limits<t_type>::epsilon();
        else scalar = std::numeric_limits<t_type>::epsilon();
    }

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] / scalar;
        vec.m_elem[irow + 1] = m_elem[irow + 1] / scalar;
        vec.m_elem[irow + 2] = m_elem[irow + 2] / scalar;
        vec.m_elem[irow + 3] = m_elem[irow + 3] / scalar;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] / scalar;
    }

    return vec;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> Vector<0, t_type>::operator+(const Vector<row, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(row);

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] + v.m_elem[irow];
        vec.m_elem[irow + 1] = m_elem[irow + 1] + v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] + v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] + v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] + v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> Vector<0, t_type>::operator-(const Vector<row, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(row);

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] - v.m_elem[irow];
        vec.m_elem[irow + 1] = m_elem[irow + 1] - v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] - v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] - v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] - v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> Vector<0, t_type>::operator*(const Vector<row, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(row);

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] * v.m_elem[irow];
        vec.m_elem[irow + 1] = m_elem[irow + 1] * v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] * v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] * v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> Vector<0, t_type>::operator/(const Vector<row, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(row);
    t_type den[row]{0};

    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 1]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 1] < 0) den[irow + 1] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 1] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 2]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 2] < 0) den[irow + 2] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 2] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 3]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 3] < 0) den[irow + 3] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 3] = std::numeric_limits<t_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
        vec[irow + 1] = m_elem[irow + 1] / den[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] / den[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] / den[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
    }

    return vec;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> Vector<0, t_type>::operator+(const Matrix<row, 1, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(row);

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] + v.m_elem[irow];
        vec.m_elem[irow + 1] = m_elem[irow + 1] + v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] + v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] + v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] + v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> Vector<0, t_type>::operator-(const Matrix<row, 1, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(row);

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] - v.m_elem[irow];
        vec.m_elem[irow + 1] = m_elem[irow + 1] - v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] - v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] - v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] - v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> Vector<0, t_type>::operator*(const Matrix<row, 1, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(row);

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = m_elem[irow] * v.m_elem[irow];
        vec.m_elem[irow + 1] = m_elem[irow + 1] * v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = m_elem[irow + 2] * v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = m_elem[irow] * v.m_elem[irow];
    }

    return vec;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> Vector<0, t_type>::operator/(const Matrix<row, 1, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    Vector<0, t_type> vec(row);
    t_type den[row]{0};

    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 1]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 1] < 0) den[irow + 1] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 1] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 2]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 2] < 0) den[irow + 2] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 2] = std::numeric_limits<t_type>::epsilon();
        }
        if (std::abs(den[irow + 3]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow + 3] < 0) den[irow + 3] = -std::numeric_limits<t_type>::epsilon();
            else den[irow + 3] = std::numeric_limits<t_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
        vec[irow + 1] = m_elem[irow + 1] / den[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] / den[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] / den[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
    }

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator&(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "This method is only for 3 x 1 vector");
    assert(v.m_row == 3 && "This method is only for 3 x 1 vector");

    Vector<0, t_type> vec(3);

    vec.m_elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    vec.m_elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    vec.m_elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator&(const Vector3<t_type, 3> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "This method is only for 3 x 1 vector");

    Vector<0, t_type> vec(3);

    vec.m_elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    vec.m_elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    vec.m_elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return vec;
}

template <typename t_type>
inline Vector<0, t_type> Vector<0, t_type>::operator&(const Matrix<0, 0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "This method is only for 3 x 1 vector");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == 3 && "This method is only for 3 x 1 vector");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    Vector<0, t_type> vec(3);

    vec.m_elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    vec.m_elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    vec.m_elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return vec;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> Vector<0, t_type>::operator&(const Vector<row, t_type> &v) const
{
    static_assert(row == 3, "This method is only for 3 x 1 vector");
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "This method is only for 3 x 1 vector");

    Vector<0, t_type> vec(3);

    vec.m_elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    vec.m_elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    vec.m_elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return vec;
}

template <typename t_type>
template <uint16_t row>
inline Vector<0, t_type> Vector<0, t_type>::operator&(const Matrix<row, 1, t_type> &v) const
{
    static_assert(row == 3, "This method is only for 3 x 1 vector");
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "This method is only for 3 x 1 vector");

    Vector<0, t_type> vec(3);

    vec.m_elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    vec.m_elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    vec.m_elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return vec;
}

template <typename t_type>
inline Matrix3<t_type, 3, 3> Vector<0, t_type>::operator&(const Matrix3<t_type, 3, 3> &m) const
{ // [v]x * Mat3, []x is skew-symmetric matrix
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "This method is only for 3 x 1 vector");

    return Matrix3<t_type, 3, 3>(
        m.m_elem[6] * m_elem[1] - m.m_elem[3] * m_elem[2], m.m_elem[7] * m_elem[1] - m.m_elem[4] * m_elem[2], m.m_elem[8] * m_elem[1] - m.m_elem[5] * m_elem[2],
        m.m_elem[0] * m_elem[2] - m.m_elem[6] * m_elem[0], m.m_elem[1] * m_elem[2] - m.m_elem[7] * m_elem[0], m.m_elem[2] * m_elem[2] - m.m_elem[8] * m_elem[0],
        m.m_elem[3] * m_elem[0] - m.m_elem[0] * m_elem[1], m.m_elem[4] * m_elem[0] - m.m_elem[1] * m_elem[1], m.m_elem[5] * m_elem[0] - m.m_elem[2] * m_elem[1]);
}

template <typename t_type>
inline Rotation<t_type, 3, 3> Vector<0, t_type>::operator&(const Rotation<t_type, 3, 3> &m) const
{ // [v]x * RotMat, []x is skew-symmetric matrix
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "This method is only for 3 x 1 vector");

    return Rotation<t_type, 3, 3>(
        m.m_elem[6] * m_elem[1] - m.m_elem[3] * m_elem[2], m.m_elem[7] * m_elem[1] - m.m_elem[4] * m_elem[2], m.m_elem[8] * m_elem[1] - m.m_elem[5] * m_elem[2],
        m.m_elem[0] * m_elem[2] - m.m_elem[6] * m_elem[0], m.m_elem[1] * m_elem[2] - m.m_elem[7] * m_elem[0], m.m_elem[2] * m_elem[2] - m.m_elem[8] * m_elem[0],
        m.m_elem[3] * m_elem[0] - m.m_elem[0] * m_elem[1], m.m_elem[4] * m_elem[0] - m.m_elem[1] * m_elem[1], m.m_elem[5] * m_elem[0] - m.m_elem[2] * m_elem[1]);
}

template <typename t_type>
inline Matrix<0, 0, t_type> Vector<0, t_type>::outer(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt;
    uint16_t irow, icol;
    Matrix<0, 0, t_type> mat(m_row, v.m_row);

    for (irow = 0; irow < m_row; irow++)
    {
        for (cnt = v.m_row >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            mat.m_elem[irow * v.m_row + icol] = m_elem[irow] * v.m_elem[icol];
            mat.m_elem[irow * v.m_row + icol + 1] = m_elem[irow] * v.m_elem[icol + 1];
            mat.m_elem[irow * v.m_row + icol + 2] = m_elem[irow] * v.m_elem[icol + 2];
            mat.m_elem[irow * v.m_row + icol + 3] = m_elem[irow] * v.m_elem[icol + 3];
        }

        for (cnt = v.m_row % 4u; cnt > 0u; cnt--, icol++)
            mat.m_elem[irow * v.m_row + icol] = m_elem[irow] * v.m_elem[icol];
    }

    return mat;
}

template <typename t_type>
inline Matrix<0, 0, t_type> Vector<0, t_type>::outer(const Matrix<0, 0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_col == 1 && "Row dimensions do not matched");

    uint16_t cnt;
    uint16_t irow, icol;
    Matrix<0, 0, t_type> mat(m_row, v.m_row);

    for (irow = 0; irow < m_row; irow++)
    {
        for (cnt = v.m_row >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            mat.m_elem[irow * v.m_row + icol] = m_elem[irow] * v.m_elem[icol];
            mat.m_elem[irow * v.m_row + icol + 1] = m_elem[irow] * v.m_elem[icol + 1];
            mat.m_elem[irow * v.m_row + icol + 2] = m_elem[irow] * v.m_elem[icol + 2];
            mat.m_elem[irow * v.m_row + icol + 3] = m_elem[irow] * v.m_elem[icol + 3];
        }

        for (cnt = v.m_row % 4u; cnt > 0u; cnt--, icol++)
            mat.m_elem[irow * v.m_row + icol] = m_elem[irow] * v.m_elem[icol];
    }

    return mat;
}

template <typename t_type>
template <uint16_t row>
inline Matrix<0, 0, t_type> Vector<0, t_type>::outer(const Vector<row, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt;
    uint16_t irow, icol;
    Matrix<0, 0, t_type> mat(m_row, row);

    for (irow = 0; irow < m_row; irow++)
    {
        for (cnt = row >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            mat.m_elem[irow * row + icol] = m_elem[irow] * v.m_elem[icol];
            mat.m_elem[irow * row + icol + 1] = m_elem[irow] * v.m_elem[icol + 1];
            mat.m_elem[irow * row + icol + 2] = m_elem[irow] * v.m_elem[icol + 2];
            mat.m_elem[irow * row + icol + 3] = m_elem[irow] * v.m_elem[icol + 3];
        }

        for (cnt = row % 4u; cnt > 0u; cnt--, icol++)
            mat.m_elem[irow * row + icol] = m_elem[irow] * v.m_elem[icol];
    }

    return mat;
}

template <typename t_type>
template <uint16_t row>
inline Matrix<0, 0, t_type> Vector<0, t_type>::outer(const Matrix<row, 1, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt;
    uint16_t irow, icol;
    Matrix<0, 0, t_type> mat(m_row, row);

    for (irow = 0; irow < m_row; irow++)
    {
        for (cnt = row >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            mat.m_elem[irow * row + icol] = m_elem[irow] * v.m_elem[icol];
            mat.m_elem[irow * row + icol + 1] = m_elem[irow] * v.m_elem[icol + 1];
            mat.m_elem[irow * row + icol + 2] = m_elem[irow] * v.m_elem[icol + 2];
            mat.m_elem[irow * row + icol + 3] = m_elem[irow] * v.m_elem[icol + 3];
        }

        for (cnt = row % 4u; cnt > 0u; cnt--, icol++)
            mat.m_elem[irow * row + icol] = m_elem[irow] * v.m_elem[icol];
    }

    return mat;
}

template <typename t_type>
inline t_type Vector<0, t_type>::inner(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    t_type result = 0;
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        result += m_elem[irow] * v.m_elem[irow];
        result += m_elem[irow + 1] * v.m_elem[irow + 1];
        result += m_elem[irow + 2] * v.m_elem[irow + 2];
        result += m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        result += m_elem[irow] * v.m_elem[irow];
    }

    return result;
}

template <typename t_type>
inline t_type Vector<0, t_type>::inner(const Vector3<t_type, 3> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");

    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2]);
}

template <typename t_type>
inline t_type Vector<0, t_type>::inner(const Vector4<t_type, 4> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 4 && "Row dimensions do not matched");

    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3]);
}

template <typename t_type>
inline t_type Vector<0, t_type>::inner(const Vector6<t_type, 6> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 6 && "Row dimensions do not matched");

    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3] +
        m_elem[4] * v.m_elem[4] +
        m_elem[5] * v.m_elem[5]);
}

template <typename t_type>
inline t_type Vector<0, t_type>::inner(const Matrix<0, 0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    t_type result = 0;
    uint16_t cnt, irow = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        result += m_elem[irow] * v.m_elem[irow];
        result += m_elem[irow + 1] * v.m_elem[irow + 1];
        result += m_elem[irow + 2] * v.m_elem[irow + 2];
        result += m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        result += m_elem[irow] * v.m_elem[irow];
    }

    return result;
}

template <typename t_type>
template <uint16_t row>
inline t_type Vector<0, t_type>::inner(const Vector<row, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    t_type result = 0;
    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        result += m_elem[irow] * v.m_elem[irow];
        result += m_elem[irow + 1] * v.m_elem[irow + 1];
        result += m_elem[irow + 2] * v.m_elem[irow + 2];
        result += m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        result += m_elem[irow] * v.m_elem[irow];
    }

    return result;
}

template <typename t_type>
template <uint16_t row>
inline t_type Vector<0, t_type>::inner(const Matrix<row, 1, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    t_type result = 0;
    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        result += m_elem[irow] * v.m_elem[irow];
        result += m_elem[irow + 1] * v.m_elem[irow + 1];
        result += m_elem[irow + 2] * v.m_elem[irow + 2];
        result += m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        result += m_elem[irow] * v.m_elem[irow];
    }

    return result;
}

/* Comparison operators */
template <typename t_type>
inline bool Vector<0, t_type>::operator==(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return false;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <typename t_type>
inline bool Vector<0, t_type>::operator!=(const Vector<0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == v.m_row && "Row dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return true;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <typename t_type>
inline bool Vector<0, t_type>::operator==(const Vector3<t_type, 3> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");

    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;

    return true;
}

template <typename t_type>
inline bool Vector<0, t_type>::operator!=(const Vector3<t_type, 3> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 3 && "Row dimensions do not matched");

    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;

    return false;
}

template <typename t_type>
inline bool Vector<0, t_type>::operator==(const Vector4<t_type, 4> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 4 && "Row dimensions do not matched");

    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;

    return true;
}

template <typename t_type>
inline bool Vector<0, t_type>::operator!=(const Vector4<t_type, 4> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 4 && "Row dimensions do not matched");

    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;

    return false;
}

template <typename t_type>
inline bool Vector<0, t_type>::operator==(const Vector6<t_type, 6> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 6 && "Row dimensions do not matched");

    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;
    if (std::abs(m_elem[4] - v.m_elem[4]) > m_tolerance) return false;
    if (std::abs(m_elem[5] - v.m_elem[5]) > m_tolerance) return false;

    return true;
}

template <typename t_type>
inline bool Vector<0, t_type>::operator!=(const Vector6<t_type, 6> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == 6 && "Row dimensions do not matched");

    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;
    if (std::abs(m_elem[4] - v.m_elem[4]) > m_tolerance) return true;
    if (std::abs(m_elem[5] - v.m_elem[5]) > m_tolerance) return true;

    return false;
}

template <typename t_type>
inline bool Vector<0, t_type>::operator==(const Matrix<0, 0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return false;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <typename t_type>
inline bool Vector<0, t_type>::operator!=(const Matrix<0, 0, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == m_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = m_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return true;
    }

    for (cnt = m_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <typename t_type>
template <uint16_t row>
inline bool Vector<0, t_type>::operator==(const Vector<row, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return false;
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <typename t_type>
template <uint16_t row>
inline bool Vector<0, t_type>::operator!=(const Vector<row, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return true;
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <typename t_type>
template <uint16_t row>
inline bool Vector<0, t_type>::operator==(const Matrix<row, 1, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return false;
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <typename t_type>
template <uint16_t row>
inline bool Vector<0, t_type>::operator!=(const Matrix<row, 1, t_type> &v) const
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    assert(m_row == row && "Row dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return true;
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <typename t_type>
inline void Vector<0, t_type>::Print(const char endChar)
{
#if defined(ARDUINO)
    for (uint16_t irow = 0; irow < t_row; irow++)
    {
        Serial.printf("%7.3f\n", (t_type)m_elem[irow]);
    }
    Serial.write(endChar);
#else
    for (uint16_t irow = 0; irow < m_row; irow++)
    {
        printf("%7.3f\n", (t_type)m_elem[irow]);
    }
    printf("%c", endChar);
#endif
}

//-- Private Member Function ------------------------------------------------//
template <typename t_type>
inline void Vector<0, t_type>::CrossProduct(const t_type *v)
{
    t_type elem[3];

    elem[0] = m_elem[1] * v[2] - m_elem[2] * v[1];
    elem[1] = m_elem[2] * v[0] - m_elem[0] * v[2];
    elem[2] = m_elem[0] * v[1] - m_elem[1] * v[0];

    m_elem[0] = elem[0];
    m_elem[1] = elem[1];
    m_elem[2] = elem[2];
}

//-- Template Function ------------------------------------------------------//
// scalar + vector
template <typename type>
inline Vector<0, type> operator+(const type s, const Vector<0, type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    Vector<0, type> vec(v.m_row);

    for (cnt = v.m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = v.m_elem[irow] + s;
        vec.m_elem[irow + 1] = v.m_elem[irow + 1] + s;
        vec.m_elem[irow + 2] = v.m_elem[irow + 2] + s;
        vec.m_elem[irow + 3] = v.m_elem[irow + 3] + s;
    }

    for (cnt = v.m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = v.m_elem[irow] + s;
    }

    return vec;
}

// scalar - vector
template <uint16_t row, typename type>
inline Vector<0, type> operator-(const type s, const Vector<0, type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    Vector<0, type> vec(v.m_row);

    for (cnt = v.m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[irow] = s - v.m_elem[irow];
        vec.m_elem[irow + 1] = s - v.m_elem[irow + 1];
        vec.m_elem[irow + 2] = s - v.m_elem[irow + 2];
        vec.m_elem[irow + 3] = s - v.m_elem[irow + 3];
    }

    for (cnt = v.m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[irow] = s - v.m_elem[irow];
    }

    return vec;
}

// scalar * vector
template <uint16_t row, typename type>
inline Vector<0, type> operator*(const type s, const Vector<0, type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    Vector<0, type> vec(v.m_row);

    for (cnt = v.m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = v.m_elem[irow] * s;
        vec[irow + 1] = v.m_elem[irow + 1] * s;
        vec[irow + 2] = v.m_elem[irow + 2] * s;
        vec[irow + 3] = v.m_elem[irow + 3] * s;
    }

    for (cnt = v.m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = v.m_elem[irow] * s;
    }

    return vec;
}

// scalar / vector
template <uint16_t row, typename type>
inline Vector<0, type> operator/(const type s, const Vector<0, type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt, irow = 0;
    type den[4];
    Vector<0, type> vec(v.m_row);

    for (cnt = v.m_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        den[0] = v.m_elem[irow];
        den[1] = v.m_elem[irow + 1];
        den[2] = v.m_elem[irow + 2];
        den[3] = v.m_elem[irow + 3];

        if (std::abs(den[0]) < std::numeric_limits<type>::epsilon())
        {
            if (den[0] < 0) den[0] = -std::numeric_limits<type>::epsilon();
            else den[0] = std::numeric_limits<type>::epsilon();
        }
        if (std::abs(den[1]) < std::numeric_limits<type>::epsilon())
        {
            if (den[1] < 0) den[1] = -std::numeric_limits<type>::epsilon();
            else den[1] = std::numeric_limits<type>::epsilon();
        }
        if (std::abs(den[2]) < std::numeric_limits<type>::epsilon())
        {
            if (den[2] < 0) den[2] = -std::numeric_limits<type>::epsilon();
            else den[2] = std::numeric_limits<type>::epsilon();
        }
        if (std::abs(den[3]) < std::numeric_limits<type>::epsilon())
        {
            if (den[3] < 0) den[3] = -std::numeric_limits<type>::epsilon();
            else den[3] = std::numeric_limits<type>::epsilon();
        }
        vec.m_elem[irow] = s / den[0];
        vec.m_elem[irow + 1] = s / den[1];
        vec.m_elem[irow + 2] = s / den[2];
        vec.m_elem[irow + 3] = s / den[3];
    }

    for (cnt = v.m_row % 4u; cnt > 0u; cnt--, irow++)
    {
        den[0] = v.m_elem[irow];
        if (std::abs(den[0]) < std::numeric_limits<type>::epsilon())
        {
            if (den[0] < 0) den[0] = -std::numeric_limits<type>::epsilon();
            else den[0] = std::numeric_limits<type>::epsilon();
        }
        vec.m_elem[irow] = s / den[0];
    }

    return vec;
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTVECTOR0_TPP_