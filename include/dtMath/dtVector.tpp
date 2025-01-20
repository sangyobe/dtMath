/*!
\file       dtVector.tpp
\brief      dtMath, General Vector(m x 1) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTVECTOR_TPP_
#define DTMATH_DTVECTOR_TPP_

#include "dtVector.h"

#include <cassert>

namespace dt
{
namespace Math
{

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type>::Vector() : m_elem()
{
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type>::Vector(const t_type *element)
{
    memcpy(m_elem, element, sizeof(t_type) * t_row);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type>::Vector(const t_type *element, const size_t n_byte)
{
    size_t vecSz = sizeof(t_type) * t_row;

    if (vecSz > n_byte)
    {
        memset(m_elem, 0, sizeof(t_type) * t_row);
        memcpy(m_elem, element, n_byte);
    }
    else
        memcpy(m_elem, element, vecSz);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type>::Vector(const Vector &v)
{
    memcpy(m_elem, v.m_elem, sizeof(t_type) * t_row);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type>::Vector(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    memcpy(m_elem, v.m_elem, sizeof(t_type) * t_row);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type>::Vector(const Matrix<t_row, 1, t_type> &v)
{
    memcpy(m_elem, v.m_elem, sizeof(t_type) * t_row);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type>::Vector(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    memcpy(m_elem, v.m_elem, sizeof(t_type) * t_row);
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetZero()
{
    memset(m_elem, 0, sizeof(t_type) * t_row);
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetFill(const t_type value)
{
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] = value;
        m_elem[irow + 1] = value;
        m_elem[irow + 2] = value;
        m_elem[irow + 3] = value;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
        m_elem[irow] = value;
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetElement(const t_type *element, const size_t n_byte)
{
    size_t vecSz = sizeof(t_type) * t_row;

    if (vecSz > n_byte) memcpy(m_elem, element, n_byte);
    else memcpy(m_elem, element, vecSz);
}

// template <uint16_t t_row, typename t_type>
// inline void Vector<t_row, t_type>::SetElement(const Vector &v)
// {
//     memcpy(m_elem, v.m_elem, sizeof(t_type) * t_row);
// }

// template <uint16_t t_row, typename t_type>
// inline void Vector<t_row, t_type>::SetElement(const Matrix<t_row, 1, t_type> &v)
// {
//     memcpy(m_elem, v.m_elem, sizeof(t_type) * t_row);
// }

// template <uint16_t t_row, typename t_type>
// inline void Vector<t_row, t_type>::SetElement(const Matrix<0, 0, t_type> &v)
// {
//     assert(v.m_elem != nullptr && "Memory has not been allocated");
//     assert(v.m_row == t_row && "Row dimensions do not matched");
//     assert(v.m_col == 1 && "Col dimensions do not matched");

//     memcpy(m_elem, v.m_elem, sizeof(t_type) * t_row);
// }

// template<uint16_t t_row, typename t_type>
// inline void Vector<t_row, t_type>::SetElement(const t_type elem0, ...)
//{
//     va_list ap;
//     va_start(ap, elem0);
//
//     m_elem[0] = elem0;
//
//     for (uint16_t irow = 1; irow < t_row; ++irow)
//     {
//         m_elem[irow] = (t_type)(va_arg(ap, double));
//     }
//
//     va_end(ap);
// }

template <uint16_t t_row, typename t_type>
template <uint16_t row>
inline void Vector<t_row, t_type>::SetBlock(const uint16_t idxRow, const Vector<row, t_type> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(row > jdx && "Index out of range");
    assert(((int16_t)row - jdx) >= size && "over size");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (rowSz > size) rowSz = size;

    memcpy(&m_elem[idxRow], &v.m_elem[jdx], sizeof(t_type) * rowSz);
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetBlock(const uint16_t idxRow, const t_type *v, const size_t n_byte)
{
    assert(t_row > idxRow && "Index out of range");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    uint16_t row = (uint16_t)(n_byte / sizeof(t_type));
    if (rowSz > row) rowSz = row;

    memcpy(&m_elem[idxRow], v, sizeof(t_type) * rowSz);
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetBlock(const uint16_t idxRow, const Vector3<t_type, 3> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(3 > jdx && "Index out of range");
    assert((3 - jdx) >= size && "over size");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (rowSz > size) rowSz = size;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[jdx];
        m_elem[idxRow + 1] = v.m_elem[jdx + 1];
        break;
    default:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        break;
    }
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetBlock(const uint16_t idxRow, const Vector4<t_type, 4> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "over size");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (rowSz > size) rowSz = size;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[jdx];
        m_elem[idxRow + 1] = v.m_elem[jdx + 1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[jdx];
        m_elem[idxRow + 1] = v.m_elem[jdx + 1];
        m_elem[idxRow + 2] = v.m_elem[jdx + 2];
        break;
    default:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        break;
    }
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetBlock(const uint16_t idxRow, const Quaternion<t_type, 4> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "over size");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (rowSz > size) rowSz = size;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[jdx];
        m_elem[idxRow + 1] = v.m_elem[jdx + 1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[jdx];
        m_elem[idxRow + 1] = v.m_elem[jdx + 1];
        m_elem[idxRow + 2] = v.m_elem[jdx + 2];
        break;
    default:
        m_elem[idxRow] = v.m_elem[0];
        m_elem[idxRow + 1] = v.m_elem[1];
        m_elem[idxRow + 2] = v.m_elem[2];
        m_elem[idxRow + 3] = v.m_elem[3];
        break;
    }
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetBlock(const uint16_t idxRow, const Vector6<t_type, 6> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(6 > jdx && "Index out of range");
    assert((6 - jdx) >= size && "over size");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (rowSz > size) rowSz = size;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow] = v.m_elem[jdx];
        m_elem[idxRow + 1] = v.m_elem[jdx + 1];
        break;
    case 3:
        m_elem[idxRow] = v.m_elem[jdx];
        m_elem[idxRow + 1] = v.m_elem[jdx + 1];
        m_elem[idxRow + 2] = v.m_elem[jdx + 2];
        break;
    case 4:
        m_elem[idxRow] = v.m_elem[jdx];
        m_elem[idxRow + 1] = v.m_elem[jdx + 1];
        m_elem[idxRow + 2] = v.m_elem[jdx + 2];
        m_elem[idxRow + 3] = v.m_elem[jdx + 3];
        break;
    case 5:
        m_elem[idxRow] = v.m_elem[jdx];
        m_elem[idxRow + 1] = v.m_elem[jdx + 1];
        m_elem[idxRow + 2] = v.m_elem[jdx + 2];
        m_elem[idxRow + 3] = v.m_elem[jdx + 3];
        m_elem[idxRow + 4] = v.m_elem[jdx + 4];
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

template <uint16_t t_row, typename t_type>
template <uint16_t row>
inline void Vector<t_row, t_type>::SetBlock(const uint16_t idxRow, const Matrix<row, 1, t_type> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(row > jdx && "Index out of range");
    assert(((int16_t)row - jdx) >= size && "over size");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (rowSz > size) rowSz = size;

    memcpy(&m_elem[idxRow], &v.m_elem[jdx], sizeof(t_type) * rowSz);
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetBlock(const uint16_t idxRow, const Matrix<0, 0, t_type> &v, const uint16_t jdx, int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(v.m_row > jdx && "Index out of range");
    assert(((int16_t)v.m_row - jdx) >= size && "over size");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (size < 0) size = v.m_row;
    if (rowSz > size) rowSz = size;

    memcpy(&m_elem[idxRow], &v.m_elem[jdx], sizeof(t_type) * rowSz);
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetBlock(const uint16_t idxRow, const Vector<0, t_type> &v, const uint16_t jdx, int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(v.m_row > jdx && "Index out of range");
    assert(((int16_t)v.m_row - jdx) >= size && "over size");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (size < 0) size = v.m_row;
    if (rowSz > size) rowSz = size;

    memcpy(&m_elem[idxRow], &v.m_elem[jdx], sizeof(t_type) * rowSz);
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetSwap(const uint16_t i, const uint16_t j)
{
    assert(t_row > i && "Index out of range");
    assert(t_row > j && "Index out of range");

    t_type elem = m_elem[i];
    m_elem[i] = m_elem[j];
    m_elem[j] = elem;
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::SetNormalize()
{
    uint16_t cnt = 0;
    uint16_t irow = 0;
    t_type norm = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        norm += m_elem[irow] * m_elem[irow];
        norm += m_elem[irow + 1] * m_elem[irow + 1];
        norm += m_elem[irow + 2] * m_elem[irow + 2];
        norm += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        norm += m_elem[irow] * m_elem[irow];
    }

    norm = std::sqrt(norm);

    if (norm < std::numeric_limits<t_type>::epsilon())
        norm = std::numeric_limits<t_type>::epsilon();

    for (cnt = t_row >> 2u, irow = 0; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] /= norm;
        m_elem[irow + 1] /= norm;
        m_elem[irow + 2] /= norm;
        m_elem[irow + 3] /= norm;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
        m_elem[irow] /= norm;
}

template <uint16_t t_row, typename t_type>
inline const t_type *const Vector<t_row, t_type>::GetElementsAddr() const
{
    return m_elem;
}

template <uint16_t t_row, typename t_type>
template <uint16_t row>
inline int8_t Vector<t_row, t_type>::GetBlock(const uint16_t idx, Vector<row, t_type> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(row > jdx && "Index out of range");
    assert(((int16_t)row - jdx) >= size && "over size");

    uint16_t rowSize = t_row - idx;

    if (idx >= t_row) return -1;
    if (rowSize > size) rowSize = size;

    memcpy(&v.m_elem[jdx], &m_elem[idx], sizeof(t_type) * rowSize);

    return 0;
}

template <uint16_t t_row, typename t_type>
inline int8_t Vector<t_row, t_type>::GetBlock(const uint16_t idx, Vector<0, t_type> &v, const uint16_t jdx, int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(v.m_row > jdx && "Index out of range");
    assert(((int16_t)v.m_row - jdx) >= size && "over size");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t rowSize = t_row - idx;

    if (idx >= t_row) return -1;
    if (size < 0) size = v.m_row;
    if (rowSize > size) rowSize = size;

    memcpy(&v.m_elem[jdx], &m_elem[idx], sizeof(t_type) * rowSize);

    return 0;
}

template <uint16_t t_row, typename t_type>
inline int8_t Vector<t_row, t_type>::GetBlock(const uint16_t idx, Vector3<t_type, 3> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(3 > jdx && "Index out of range");
    assert((3 - jdx) >= size && "over size");

    if (idx >= t_row) return -1;

    uint16_t rowSize = t_row - idx;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idx];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idx];
        v.m_elem[jdx + 1] = m_elem[idx + 1];
        break;
    default:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
    };

    return 0;
}

template <uint16_t t_row, typename t_type>
inline int8_t Vector<t_row, t_type>::GetBlock(const uint16_t idx, Vector4<t_type, 4> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "over size");

    if (idx >= t_row) return -1;

    uint16_t rowSize = t_row - idx;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idx];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idx];
        v.m_elem[jdx + 1] = m_elem[idx + 1];
        break;
    case 3:
        v.m_elem[jdx] = m_elem[idx];
        v.m_elem[jdx + 1] = m_elem[idx + 1];
        v.m_elem[jdx + 2] = m_elem[idx + 2];
        break;
    default:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        v.m_elem[3] = m_elem[idx + 3];
    };

    return 0;
}

template <uint16_t t_row, typename t_type>
inline int8_t Vector<t_row, t_type>::GetBlock(const uint16_t idx, Vector6<t_type, 6> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(6 > jdx && "Index out of range");
    assert((6 - jdx) >= size && "over size");

    if (idx >= t_row) return -1;

    uint16_t rowSize = t_row - idx;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idx];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idx];
        v.m_elem[jdx + 1] = m_elem[idx + 1];
        break;
    case 3:
        v.m_elem[jdx] = m_elem[idx];
        v.m_elem[jdx + 1] = m_elem[idx + 1];
        v.m_elem[jdx + 2] = m_elem[idx + 2];
        break;
    case 4:
        v.m_elem[jdx] = m_elem[idx];
        v.m_elem[jdx + 1] = m_elem[idx + 1];
        v.m_elem[jdx + 2] = m_elem[idx + 2];
        v.m_elem[jdx + 3] = m_elem[idx + 3];
        break;
    case 5:
        v.m_elem[jdx] = m_elem[idx];
        v.m_elem[jdx + 1] = m_elem[idx + 1];
        v.m_elem[jdx + 2] = m_elem[idx + 2];
        v.m_elem[jdx + 3] = m_elem[idx + 3];
        v.m_elem[jdx + 4] = m_elem[idx + 4];
        break;
    default:
        v.m_elem[jdx] = m_elem[idx];
        v.m_elem[jdx + 1] = m_elem[idx + 1];
        v.m_elem[jdx + 2] = m_elem[idx + 2];
        v.m_elem[jdx + 3] = m_elem[idx + 3];
        v.m_elem[jdx + 4] = m_elem[idx + 4];
        v.m_elem[jdx + 5] = m_elem[idx + 5];
    };

    return 0;
}

template <uint16_t t_row, typename t_type>
inline int8_t Vector<t_row, t_type>::GetBlock(const uint16_t idx, Quaternion<t_type, 4> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "over size");

    if (idx >= t_row) return -1;

    uint16_t rowSize = t_row - idx;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idx];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idx];
        v.m_elem[jdx + 1] = m_elem[idx + 1];
        break;
    case 3:
        v.m_elem[jdx] = m_elem[idx];
        v.m_elem[jdx + 1] = m_elem[idx + 1];
        v.m_elem[jdx + 2] = m_elem[idx + 2];
        break;
    default:
        v.m_elem[0] = m_elem[idx];
        v.m_elem[1] = m_elem[idx + 1];
        v.m_elem[2] = m_elem[idx + 2];
        v.m_elem[3] = m_elem[idx + 3];
    };

    return 0;
}

template <uint16_t t_row, typename t_type>
template <uint16_t row>
inline Vector<row, t_type> Vector<t_row, t_type>::GetBlockVec(const uint16_t idx, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(row > jdx && "Index out of range");
    assert(((int16_t)row - jdx) >= size && "over size");

    t_type elem[row]{0};
    uint16_t rowSize = t_row - idx;

    if (idx >= t_row) return Vector<row, t_type>(elem);
    if (rowSize > size) rowSize = size;

    memcpy(&elem[jdx], &m_elem[idx], sizeof(t_type) * rowSize);

    return Vector<row, t_type>(elem);
}

template <uint16_t t_row, typename t_type>
inline Vector<0, t_type> Vector<t_row, t_type>::GetBlockVec0(const uint16_t idx, const uint16_t row, const uint16_t jdx, int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(row > jdx && "Index out of range");
    assert(((int16_t)row - jdx) >= size && "over size");

    Vector<0, t_type> v(row);
    uint16_t rowSize = t_row - idx;

    if (size < 0) size = row;
    if (rowSize > size) rowSize = size;

    memcpy(&v.m_elem[jdx], &m_elem[idx], sizeof(t_type) * rowSize);

    return v;
}

template <uint16_t t_row, typename t_type>
inline Vector3<t_type, 3> Vector<t_row, t_type>::GetBlockVec3(const uint16_t idx, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(3 > jdx && "Index out of range");
    assert((3 - jdx) >= size && "over size");

    t_type elem[3]{0};

    if (idx >= t_row) return Vector3<t_type, 3>(elem);

    uint16_t rowSize = t_row - idx;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        elem[jdx] = m_elem[idx];
        break;
    case 2:
        elem[jdx] = m_elem[idx];
        elem[jdx + 1] = m_elem[idx + 1];
        break;
    default:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
    };

    return Vector3<t_type, 3>(elem);
}

template <uint16_t t_row, typename t_type>
inline Vector4<t_type, 4> Vector<t_row, t_type>::GetBlockVec4(const uint16_t idx, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "over size");

    t_type elem[4]{0};

    if (idx >= t_row) return Vector4<t_type, 4>(elem);

    uint16_t rowSize = t_row - idx;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        elem[jdx] = m_elem[idx];
        break;
    case 2:
        elem[jdx] = m_elem[idx];
        elem[jdx + 1] = m_elem[idx + 1];
        break;
    case 3:
        elem[jdx] = m_elem[idx];
        elem[jdx + 1] = m_elem[idx + 1];
        elem[jdx + 2] = m_elem[idx + 2];
        break;
    default:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        elem[3] = m_elem[idx + 3];
    };

    return Vector4<t_type, 4>(elem);
}

template <uint16_t t_row, typename t_type>
inline Vector6<t_type, 6> Vector<t_row, t_type>::GetBlockVec6(const uint16_t idx, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(6 > jdx && "Index out of range");
    assert((6 - jdx) >= size && "over size");

    t_type elem[6]{0};

    if (idx >= t_row) return Vector6<t_type, 6>(elem);

    uint16_t rowSize = t_row - idx;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        elem[jdx] = m_elem[idx];
        break;
    case 2:
        elem[jdx] = m_elem[idx];
        elem[jdx + 1] = m_elem[idx + 1];
        break;
    case 3:
        elem[jdx] = m_elem[idx];
        elem[jdx + 1] = m_elem[idx + 1];
        elem[jdx + 2] = m_elem[idx + 2];
        break;
    case 4:
        elem[jdx] = m_elem[idx];
        elem[jdx + 1] = m_elem[idx + 1];
        elem[jdx + 2] = m_elem[idx + 2];
        elem[jdx + 3] = m_elem[idx + 3];
        break;
    case 5:
        elem[jdx] = m_elem[idx];
        elem[jdx + 1] = m_elem[idx + 1];
        elem[jdx + 2] = m_elem[idx + 2];
        elem[jdx + 3] = m_elem[idx + 3];
        elem[jdx + 4] = m_elem[idx + 4];
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

template <uint16_t t_row, typename t_type>
inline Quaternion<t_type, 4> Vector<t_row, t_type>::GetBlockQuat(const uint16_t idx, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idx && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "over size");

    t_type elem[4]{0};

    if (idx >= t_row) return Vector4<t_type, 4>(elem);

    uint16_t rowSize = t_row - idx;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        elem[jdx] = m_elem[idx];
        break;
    case 2:
        elem[jdx] = m_elem[idx];
        elem[jdx + 1] = m_elem[idx + 1];
        break;
    case 3:
        elem[jdx] = m_elem[idx];
        elem[jdx + 1] = m_elem[idx + 1];
        elem[jdx + 2] = m_elem[idx + 2];
        break;
    default:
        elem[0] = m_elem[idx];
        elem[1] = m_elem[idx + 1];
        elem[2] = m_elem[idx + 2];
        elem[3] = m_elem[idx + 3];
    };

    return Quaternion<t_type, 4>(elem);
}

template <uint16_t t_row, typename t_type>
inline t_type Vector<t_row, t_type>::GetNorm() const
{
    t_type sqSum = 0;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        sqSum += m_elem[irow] * m_elem[irow];
        sqSum += m_elem[irow + 1] * m_elem[irow + 1];
        sqSum += m_elem[irow + 2] * m_elem[irow + 2];
        sqSum += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        sqSum += m_elem[irow] * m_elem[irow];
    }

    return std::sqrt(sqSum);
}

template <uint16_t t_row, typename t_type>
inline t_type Vector<t_row, t_type>::GetSqNorm() const
{
    t_type sqSum = 0;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        sqSum += m_elem[irow] * m_elem[irow];
        sqSum += m_elem[irow + 1] * m_elem[irow + 1];
        sqSum += m_elem[irow + 2] * m_elem[irow + 2];
        sqSum += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        sqSum += m_elem[irow] * m_elem[irow];
    }

    return sqSum;
}

template <uint16_t t_row, typename t_type>
inline t_type Vector<t_row, t_type>::GetLpNorm(const int p) const
{
    uint16_t cnt, i = 0;
    t_type powSum = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        powSum += std::pow(std::abs(m_elem[i]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 1]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 2]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 3]), (t_type)p);
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, i++)
    {
        powSum += std::pow(std::abs(m_elem[i]), (t_type)p);
    }

    return std::pow(powSum, (t_type)1 / p);
}

template <uint16_t t_row, typename t_type>
inline t_type Vector<t_row, t_type>::GetSum() const
{
    t_type rtn = 0;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        rtn += m_elem[irow];
        rtn += m_elem[irow + 1];
        rtn += m_elem[irow + 2];
        rtn += m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        rtn += m_elem[irow];
    }

    return rtn;
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::GetNormalized() const
{
    t_type vec[t_row];
    uint16_t cnt = 0;
    uint16_t irow = 0;
    t_type norm = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        norm += m_elem[irow] * m_elem[irow];
        norm += m_elem[irow + 1] * m_elem[irow + 1];
        norm += m_elem[irow + 2] * m_elem[irow + 2];
        norm += m_elem[irow + 3] * m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        norm += m_elem[irow] * m_elem[irow];
    }

    norm = std::sqrt(norm);

    if (norm < std::numeric_limits<t_type>::epsilon())
        norm = std::numeric_limits<t_type>::epsilon();

    for (cnt = t_row >> 2u, irow = 0; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] / norm;
        vec[irow + 1] = m_elem[irow + 1] / norm;
        vec[irow + 2] = m_elem[irow + 2] / norm;
        vec[irow + 3] = m_elem[irow + 3] / norm;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
        vec[irow] = m_elem[irow] / norm;

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, typename t_type>
inline Matrix<3, 3, t_type> Vector<t_row, t_type>::GetSkew() const
{
    static_assert(t_row == 3, "This method is only for 3 x 1 vector");

    return Matrix<3, 3, t_type>(
        0, -m_elem[2], m_elem[1],
        m_elem[2], 0, -m_elem[0],
        -m_elem[1], m_elem[0], 0);
}

template <uint16_t t_row, typename t_type>
inline Matrix<1, t_row, t_type> Vector<t_row, t_type>::Transpose() const
{
    return Matrix<1, t_row, t_type>(m_elem);
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::Transpose(Matrix<1, t_row, t_type> &m) const
{
    memcpy(m.m_elem, m_elem, sizeof(t_type) * t_row);
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::Transpose(Matrix<0, 0, t_type> &m) const
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == 1 && "Row dimensions do not matched");
    assert(m.m_col == t_row && "Col dimensions do not matched");

    memcpy(m.m_elem, m_elem, sizeof(t_type) * t_row);
}

/* Member access operators */
template <uint16_t t_row, typename t_type>
inline t_type &Vector<t_row, t_type>::operator()(uint16_t irow)
{
    assert(irow < t_row && "Index out of range");
    return m_elem[irow];
}

template <uint16_t t_row, typename t_type>
inline const t_type &Vector<t_row, t_type>::operator()(uint16_t irow) const
{
    assert(irow < t_row && "Index out of range");
    return m_elem[irow];
}

/* Assignment operators */
template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator=(const Vector &v)
{
    memcpy(m_elem, v.m_elem, sizeof(t_type) * t_row);

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator+=(const Vector &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += v.m_elem[irow];
        m_elem[irow + 1] += v.m_elem[irow + 1];
        m_elem[irow + 2] += v.m_elem[irow + 2];
        m_elem[irow + 3] += v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator-=(const Vector &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= v.m_elem[irow];
        m_elem[irow + 1] -= v.m_elem[irow + 1];
        m_elem[irow + 2] -= v.m_elem[irow + 2];
        m_elem[irow + 3] -= v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseMulEq(const Vector &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= v.m_elem[irow];
        m_elem[irow + 1] *= v.m_elem[irow + 1];
        m_elem[irow + 2] *= v.m_elem[irow + 2];
        m_elem[irow + 3] *= v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseDivEq(const Vector &v)
{
    uint16_t cnt, irow = 0;
    t_type den[t_row];
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
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

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
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

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator=(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == v.m_row && "Row dimensions do not matched");

    memcpy(m_elem, v.m_elem, sizeof(t_type) * t_row);

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator+=(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == v.m_row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += v.m_elem[irow];
        m_elem[irow + 1] += v.m_elem[irow + 1];
        m_elem[irow + 2] += v.m_elem[irow + 2];
        m_elem[irow + 3] += v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator-=(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == v.m_row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= v.m_elem[irow];
        m_elem[irow + 1] -= v.m_elem[irow + 1];
        m_elem[irow + 2] -= v.m_elem[irow + 2];
        m_elem[irow + 3] -= v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseMulEq(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == v.m_row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= v.m_elem[irow];
        m_elem[irow + 1] *= v.m_elem[irow + 1];
        m_elem[irow + 2] *= v.m_elem[irow + 2];
        m_elem[irow + 3] *= v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseDivEq(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == v.m_row && "Row dimensions do not matched");

    uint16_t cnt, irow = 0;
    t_type den[t_row];
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
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

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
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

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator=(const Vector3<t_type, t_row> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator+=(const Vector3<t_type, t_row> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator-=(const Vector3<t_type, t_row> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseMulEq(const Vector3<t_type, t_row> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseDivEq(const Vector3<t_type, t_row> &v)
{
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

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator=(const Vector4<t_type, t_row> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator+=(const Vector4<t_type, t_row> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator-=(const Vector4<t_type, t_row> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseMulEq(const Vector4<t_type, t_row> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseDivEq(const Vector4<t_type, t_row> &v)
{
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

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator=(const Vector6<t_type, t_row> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
    m_elem[4] = v.m_elem[4];
    m_elem[5] = v.m_elem[5];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator+=(const Vector6<t_type, t_row> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];
    m_elem[4] += v.m_elem[4];
    m_elem[5] += v.m_elem[5];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator-=(const Vector6<t_type, t_row> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];
    m_elem[4] -= v.m_elem[4];
    m_elem[5] -= v.m_elem[5];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseMulEq(const Vector6<t_type, t_row> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];
    m_elem[4] *= v.m_elem[4];
    m_elem[5] *= v.m_elem[5];

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseDivEq(const Vector6<t_type, t_row> &v)
{
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

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator=(const Matrix<t_row, 1, t_type> &v)
{
    memcpy(m_elem, v.m_elem, sizeof(t_type) * t_row);

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator+=(const Matrix<t_row, 1, t_type> &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += v.m_elem[irow];
        m_elem[irow + 1] += v.m_elem[irow + 1];
        m_elem[irow + 2] += v.m_elem[irow + 2];
        m_elem[irow + 3] += v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator-=(const Matrix<t_row, 1, t_type> &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= v.m_elem[irow];
        m_elem[irow + 1] -= v.m_elem[irow + 1];
        m_elem[irow + 2] -= v.m_elem[irow + 2];
        m_elem[irow + 3] -= v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseMulEq(const Matrix<t_row, 1, t_type> &v)
{
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= v.m_elem[irow];
        m_elem[irow + 1] *= v.m_elem[irow + 1];
        m_elem[irow + 2] *= v.m_elem[irow + 2];
        m_elem[irow + 3] *= v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseDivEq(const Matrix<t_row, 1, t_type> &v)
{
    uint16_t cnt, irow = 0;
    t_type den[t_row];
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
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

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
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

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator=(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    memcpy(m_elem, v.m_elem, sizeof(t_type) * t_row);

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator+=(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += v.m_elem[irow];
        m_elem[irow + 1] += v.m_elem[irow + 1];
        m_elem[irow + 2] += v.m_elem[irow + 2];
        m_elem[irow + 3] += v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator-=(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= v.m_elem[irow];
        m_elem[irow + 1] -= v.m_elem[irow + 1];
        m_elem[irow + 2] -= v.m_elem[irow + 2];
        m_elem[irow + 3] -= v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseMulEq(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= v.m_elem[irow];
        m_elem[irow + 1] *= v.m_elem[irow + 1];
        m_elem[irow + 2] *= v.m_elem[irow + 2];
        m_elem[irow + 3] *= v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= v.m_elem[irow];
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::CWiseDivEq(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, irow = 0;
    t_type den[t_row];
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
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

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
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

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator=(const t_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] = s;
        m_elem[irow + 1] = s;
        m_elem[irow + 2] = s;
        m_elem[irow + 3] = s;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] = s;
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator+=(const t_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] += s;
        m_elem[irow + 1] += s;
        m_elem[irow + 2] += s;
        m_elem[irow + 3] += s;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] += s;
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator-=(const t_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] -= s;
        m_elem[irow + 1] -= s;
        m_elem[irow + 2] -= s;
        m_elem[irow + 3] -= s;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] -= s;
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator*=(const t_type s)
{
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] *= s;
        m_elem[irow + 1] *= s;
        m_elem[irow + 2] *= s;
        m_elem[irow + 3] *= s;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] *= s;
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator/=(const t_type s)
{
    t_type scalar = s;
    uint16_t cnt, irow = 0;

    if (std::abs(scalar) < std::numeric_limits<t_type>::epsilon())
    {
        if (scalar < 0) scalar = -std::numeric_limits<t_type>::epsilon();
        else scalar = std::numeric_limits<t_type>::epsilon();
    }

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow] /= scalar;
        m_elem[irow + 1] /= scalar;
        m_elem[irow + 2] /= scalar;
        m_elem[irow + 3] /= scalar;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow] /= scalar;
    }

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator&=(const Vector3<t_type, t_row> &v)
{
    CrossProduct(v.m_elem);

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator&=(const Vector<3, t_type> &v)
{
    CrossProduct(v.m_elem);

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator&=(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == 3 && "Row dimensions do not matched");

    CrossProduct(v.m_elem);

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator&=(const Matrix<3, 1, t_type> &v)
{
    CrossProduct(v.m_elem);

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> &Vector<t_row, t_type>::operator&=(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == 3 && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    CrossProduct(v.m_elem);

    return (*this);
}

template <uint16_t t_row, typename t_type>
inline CommaInit<t_row, t_type> Vector<t_row, t_type>::operator<<(const t_type s)
{
    m_elem[0] = s;
    return CommaInit<t_row, t_type>(m_elem);
}

/* Arithmetic operators */
template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator-() const
{
    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = -m_elem[irow];
        vec[irow + 1] = -m_elem[irow + 1];
        vec[irow + 2] = -m_elem[irow + 2];
        vec[irow + 3] = -m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = -m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator+(const Vector &v) const
{
    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] + v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] + v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] + v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator-(const Vector &v) const
{
    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] - v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] - v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] - v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseMul(const Vector &v) const
{
    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] * v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] * v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseDiv(const Vector &v) const
{
    t_type vec[t_row];
    t_type den[t_row];
    uint16_t cnt, irow = 0;
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
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

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator+(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] + v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] + v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] + v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator-(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] - v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] - v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] - v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseMul(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] * v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] * v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseDiv(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    t_type vec[t_row];
    t_type den[t_row];
    uint16_t cnt, irow = 0;
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
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

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator+(const Vector3<t_type, t_row> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] + v.m_elem[0];
    vec[1] = m_elem[1] + v.m_elem[1];
    vec[2] = m_elem[2] + v.m_elem[2];

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator-(const Vector3<t_type, t_row> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] - v.m_elem[0];
    vec[1] = m_elem[1] - v.m_elem[1];
    vec[2] = m_elem[2] - v.m_elem[2];

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseMul(const Vector3<t_type, t_row> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] * v.m_elem[0];
    vec[1] = m_elem[1] * v.m_elem[1];
    vec[2] = m_elem[2] * v.m_elem[2];

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseDiv(const Vector3<t_type, t_row> &v) const
{
    t_type vec[t_row];
    t_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[0] = m_elem[0] / den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[1] = m_elem[1] / den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[2] = m_elem[2] / den;

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator+(const Vector4<t_type, t_row> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] + v.m_elem[0];
    vec[1] = m_elem[1] + v.m_elem[1];
    vec[2] = m_elem[2] + v.m_elem[2];
    vec[3] = m_elem[3] + v.m_elem[3];

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator-(const Vector4<t_type, t_row> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] - v.m_elem[0];
    vec[1] = m_elem[1] - v.m_elem[1];
    vec[2] = m_elem[2] - v.m_elem[2];
    vec[3] = m_elem[3] - v.m_elem[3];

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseMul(const Vector4<t_type, t_row> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] * v.m_elem[0];
    vec[1] = m_elem[1] * v.m_elem[1];
    vec[2] = m_elem[2] * v.m_elem[2];
    vec[3] = m_elem[3] * v.m_elem[3];

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseDiv(const Vector4<t_type, t_row> &v) const
{
    t_type vec[t_row];
    t_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[0] = m_elem[0] / den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[1] = m_elem[1] / den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[2] = m_elem[2] / den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[3] = m_elem[3] / den;

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator+(const Vector6<t_type, t_row> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] + v.m_elem[0];
    vec[1] = m_elem[1] + v.m_elem[1];
    vec[2] = m_elem[2] + v.m_elem[2];
    vec[3] = m_elem[3] + v.m_elem[3];
    vec[4] = m_elem[4] + v.m_elem[4];
    vec[5] = m_elem[5] + v.m_elem[5];

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator-(const Vector6<t_type, t_row> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] - v.m_elem[0];
    vec[1] = m_elem[1] - v.m_elem[1];
    vec[2] = m_elem[2] - v.m_elem[2];
    vec[3] = m_elem[3] - v.m_elem[3];
    vec[4] = m_elem[4] - v.m_elem[4];
    vec[5] = m_elem[5] - v.m_elem[5];

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseMul(const Vector6<t_type, t_row> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] * v.m_elem[0];
    vec[1] = m_elem[1] * v.m_elem[1];
    vec[2] = m_elem[2] * v.m_elem[2];
    vec[3] = m_elem[3] * v.m_elem[3];
    vec[4] = m_elem[4] * v.m_elem[4];
    vec[5] = m_elem[5] * v.m_elem[5];

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseDiv(const Vector6<t_type, t_row> &v) const
{
    t_type vec[t_row];
    t_type den;

    den = v.m_elem[0];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[0] = m_elem[0] / den;

    den = v.m_elem[1];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[1] = m_elem[1] / den;

    den = v.m_elem[2];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[2] = m_elem[2] / den;

    den = v.m_elem[3];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[3] = m_elem[3] / den;

    den = v.m_elem[4];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[4] = m_elem[4] / den;

    den = v.m_elem[5];
    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }
    vec[5] = m_elem[5] / den;

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator+(const Matrix<t_row, 1, t_type> &v) const
{
    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] + v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] + v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] + v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator-(const Matrix<t_row, 1, t_type> &v) const
{
    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] - v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] - v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] - v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseMul(const Matrix<t_row, 1, t_type> &v) const
{
    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] * v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] * v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseDiv(const Matrix<t_row, 1, t_type> &v) const
{
    t_type vec[t_row];
    t_type den[t_row];
    uint16_t cnt, irow = 0;
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
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

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator+(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] + v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] + v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] + v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] + v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator-(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] - v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] - v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] - v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] - v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseMul(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
        vec[irow + 1] = m_elem[irow + 1] * v.m_elem[irow + 1];
        vec[irow + 2] = m_elem[irow + 2] * v.m_elem[irow + 2];
        vec[irow + 3] = m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] * v.m_elem[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::CWiseDiv(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    t_type vec[t_row];
    t_type den[t_row];
    uint16_t cnt, irow = 0;
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
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

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<t_type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<t_type>::epsilon();
            else den[irow] = std::numeric_limits<t_type>::epsilon();
        }
        vec[irow] = m_elem[irow] / den[irow];
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator+(const t_type s) const
{
    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] + s;
        vec[irow + 1] = m_elem[irow + 1] + s;
        vec[irow + 2] = m_elem[irow + 2] + s;
        vec[irow + 3] = m_elem[irow + 3] + s;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] + s;
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator-(const t_type s) const
{
    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] - s;
        vec[irow + 1] = m_elem[irow + 1] - s;
        vec[irow + 2] = m_elem[irow + 2] - s;
        vec[irow + 3] = m_elem[irow + 3] - s;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] - s;
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator*(const t_type s) const
{
    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] * s;
        vec[irow + 1] = m_elem[irow + 1] * s;
        vec[irow + 2] = m_elem[irow + 2] * s;
        vec[irow + 3] = m_elem[irow + 3] * s;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] * s;
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator/(const t_type s) const
{
    t_type den = s;
    t_type vec[t_row];
    uint16_t cnt, irow = 0;

    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0) den = -std::numeric_limits<t_type>::epsilon();
        else den = std::numeric_limits<t_type>::epsilon();
    }

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = m_elem[irow] / den;
        vec[irow + 1] = m_elem[irow + 1] / den;
        vec[irow + 2] = m_elem[irow + 2] / den;
        vec[irow + 3] = m_elem[irow + 3] / den;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = m_elem[irow] / den;
    }

    return Vector(vec);
}

template <uint16_t t_row, typename t_type>
template <uint16_t col>
inline Matrix<t_row, col, t_type> Vector<t_row, t_type>::operator*(const Matrix<1, col, t_type> &m) const
{
    t_type mat[t_row * col];
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; irow++)
    {
        for (cnt = col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            mat[irow * col + icol] = m_elem[irow] * m.m_elem[icol];
            mat[irow * col + icol + 1] = m_elem[irow] * m.m_elem[icol + 1];
            mat[irow * col + icol + 2] = m_elem[irow] * m.m_elem[icol + 2];
            mat[irow * col + icol + 3] = m_elem[irow] * m.m_elem[icol + 3];
        }

        for (cnt = col % 4u; cnt > 0u; cnt--, icol++)
            mat[irow * col + icol] = m_elem[irow] * m.m_elem[icol];
    }

    return Matrix<t_row, col, t_type>(mat);
}

template <uint16_t t_row, typename t_type>
inline Matrix<0, 0, t_type> Vector<t_row, t_type>::operator*(const Matrix<0, 0, t_type> &m) const
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == 1 && "Row dimensions do not matched");

    Matrix<0, 0, t_type> mat(t_row, m.m_col);
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; irow++)
    {
        for (cnt = m.m_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            mat.m_elem[irow * m.m_col + icol] = m_elem[irow] * m.m_elem[icol];
            mat.m_elem[irow * m.m_col + icol + 1] = m_elem[irow] * m.m_elem[icol + 1];
            mat.m_elem[irow * m.m_col + icol + 2] = m_elem[irow] * m.m_elem[icol + 2];
            mat.m_elem[irow * m.m_col + icol + 3] = m_elem[irow] * m.m_elem[icol + 3];
        }

        for (cnt = m.m_col % 4u; cnt > 0u; cnt--, icol++)
            mat.m_elem[irow * m.m_col + icol] = m_elem[irow] * m.m_elem[icol];
    }

    return mat;
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator&(const Vector3<t_type, t_row> &v) const
{
    static_assert(t_row == 3, "This method is only for 3 x 1 vector");

    t_type elem[t_row];

    elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return Vector<t_row, t_type>(elem);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator&(const Vector<t_row, t_type> &v) const
{
    static_assert(t_row == 3, "This method is only for 3 x 1 vector");

    t_type elem[t_row];

    elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return Vector<t_row, t_type>(elem);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator&(const Vector<0, t_type> &v) const
{
    static_assert(t_row == 3, "This method is only for 3 x 1 vector");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    t_type elem[t_row];

    elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return Vector<t_row, t_type>(elem);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator&(const Matrix<t_row, 1, t_type> &v) const
{
    static_assert(t_row == 3, "This method is only for 3 x 1 vector");

    t_type elem[t_row];

    elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return Vector<t_row, t_type>(elem);
}

template <uint16_t t_row, typename t_type>
inline Vector<t_row, t_type> Vector<t_row, t_type>::operator&(const Matrix<0, 0, t_type> &v) const
{
    static_assert(t_row == 3, "This method is only for 3 x 1 vector");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    t_type elem[t_row];

    elem[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    elem[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    elem[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    return Vector<t_row, t_type>(elem);
}

template <uint16_t t_row, typename t_type>
inline Matrix3<t_type, 3, 3> Vector<t_row, t_type>::operator&(const Matrix3<t_type, 3, 3> &m) const
{ // [v]x * Mat3, []x is skew-symmetric matrix
    static_assert(t_row == 3, "This method is only for 3 x 1 vector");
    return Matrix3<t_type, 3, 3>(
        m.m_elem[6] * m_elem[1] - m.m_elem[3] * m_elem[2], m.m_elem[7] * m_elem[1] - m.m_elem[4] * m_elem[2], m.m_elem[8] * m_elem[1] - m.m_elem[5] * m_elem[2],
        m.m_elem[0] * m_elem[2] - m.m_elem[6] * m_elem[0], m.m_elem[1] * m_elem[2] - m.m_elem[7] * m_elem[0], m.m_elem[2] * m_elem[2] - m.m_elem[8] * m_elem[0],
        m.m_elem[3] * m_elem[0] - m.m_elem[0] * m_elem[1], m.m_elem[4] * m_elem[0] - m.m_elem[1] * m_elem[1], m.m_elem[5] * m_elem[0] - m.m_elem[2] * m_elem[1]);
}

template <uint16_t t_row, typename t_type>
inline Rotation<t_type, 3, 3> Vector<t_row, t_type>::operator&(const Rotation<t_type, 3, 3> &m) const
{ // [v]x * RotMat, []x is skew-symmetric matrix
    static_assert(t_row == 3, "This method is only for 3 x 1 vector");
    return Rotation<t_type, 3, 3>(
        m.m_elem[6] * m_elem[1] - m.m_elem[3] * m_elem[2], m.m_elem[7] * m_elem[1] - m.m_elem[4] * m_elem[2], m.m_elem[8] * m_elem[1] - m.m_elem[5] * m_elem[2],
        m.m_elem[0] * m_elem[2] - m.m_elem[6] * m_elem[0], m.m_elem[1] * m_elem[2] - m.m_elem[7] * m_elem[0], m.m_elem[2] * m_elem[2] - m.m_elem[8] * m_elem[0],
        m.m_elem[3] * m_elem[0] - m.m_elem[0] * m_elem[1], m.m_elem[4] * m_elem[0] - m.m_elem[1] * m_elem[1], m.m_elem[5] * m_elem[0] - m.m_elem[2] * m_elem[1]);
}

template <uint16_t t_row, typename t_type>
inline Matrix<t_row, 3, t_type> Vector<t_row, t_type>::Outer(const Vector3<t_type, 3> &v) const
{
    t_type mat[t_row * 3];
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; irow++)
    {
        mat[irow * 3] = m_elem[irow] * v.m_elem[0];
        mat[irow * 3 + 1] = m_elem[irow] * v.m_elem[1];
        mat[irow * 3 + 2] = m_elem[irow] * v.m_elem[2];
    }

    return Matrix<t_row, 3, t_type>(mat);
}

template <uint16_t t_row, typename t_type>
inline Matrix<t_row, 4, t_type> Vector<t_row, t_type>::Outer(const Vector4<t_type, 4> &v) const
{
    t_type mat[t_row * 4];
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; irow++)
    {
        mat[irow * 4] = m_elem[irow] * v.m_elem[0];
        mat[irow * 4 + 1] = m_elem[irow] * v.m_elem[1];
        mat[irow * 4 + 2] = m_elem[irow] * v.m_elem[2];
        mat[irow * 4 + 3] = m_elem[irow] * v.m_elem[3];
    }

    return Matrix<t_row, 4, t_type>(mat);
}

template <uint16_t t_row, typename t_type>
inline Matrix<t_row, 6, t_type> Vector<t_row, t_type>::Outer(const Vector6<t_type, 6> &v) const
{
    t_type mat[t_row * 6];
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; irow++)
    {
        mat[irow * 6] = m_elem[irow] * v.m_elem[0];
        mat[irow * 6 + 1] = m_elem[irow] * v.m_elem[1];
        mat[irow * 6 + 2] = m_elem[irow] * v.m_elem[2];
        mat[irow * 6 + 3] = m_elem[irow] * v.m_elem[3];
        mat[irow * 6 + 4] = m_elem[irow] * v.m_elem[4];
        mat[irow * 6 + 5] = m_elem[irow] * v.m_elem[5];
    }

    return Matrix<t_row, 6, t_type>(mat);
}

template <uint16_t t_row, typename t_type>
inline Matrix<0, 0, t_type> Vector<t_row, t_type>::Outer(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    Matrix<0, 0, t_type> mat(t_row, v.m_row);
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; irow++)
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

template <uint16_t t_row, typename t_type>
inline Matrix<0, 0, t_type> Vector<t_row, t_type>::Outer(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    Matrix<0, 0, t_type> mat(t_row, v.m_row);
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; irow++)
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

template <uint16_t t_row, typename t_type>
template <uint16_t row>
inline Matrix<t_row, row, t_type> Vector<t_row, t_type>::Outer(const Vector<row, t_type> &v) const
{
    t_type mat[t_row * row];
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; irow++)
    {
        for (cnt = row >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            mat[irow * row + icol] = m_elem[irow] * v.m_elem[icol];
            mat[irow * row + icol + 1] = m_elem[irow] * v.m_elem[icol + 1];
            mat[irow * row + icol + 2] = m_elem[irow] * v.m_elem[icol + 2];
            mat[irow * row + icol + 3] = m_elem[irow] * v.m_elem[icol + 3];
        }

        for (cnt = row % 4u; cnt > 0u; cnt--, icol++)
            mat[irow * row + icol] = m_elem[irow] * v.m_elem[icol];
    }

    return Matrix<t_row, row, t_type>(mat);
}

template <uint16_t t_row, typename t_type>
template <uint16_t row>
inline Matrix<t_row, row, t_type> Vector<t_row, t_type>::Outer(const Matrix<row, 1, t_type> &v) const
{
    t_type mat[t_row * row];
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; irow++)
    {
        for (cnt = row >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            mat[irow * row + icol] = m_elem[irow] * v.m_elem[icol];
            mat[irow * row + icol + 1] = m_elem[irow] * v.m_elem[icol + 1];
            mat[irow * row + icol + 2] = m_elem[irow] * v.m_elem[icol + 2];
            mat[irow * row + icol + 3] = m_elem[irow] * v.m_elem[icol + 3];
        }

        for (cnt = row % 4u; cnt > 0u; cnt--, icol++)
            mat[irow * row + icol] = m_elem[irow] * v.m_elem[icol];
    }

    return Matrix<t_row, row, t_type>(mat);
}

template <uint16_t t_row, typename t_type>
inline t_type Vector<t_row, t_type>::Inner(const Vector &v) const
{
    t_type result = 0;
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        result += m_elem[irow] * v.m_elem[irow];
        result += m_elem[irow + 1] * v.m_elem[irow + 1];
        result += m_elem[irow + 2] * v.m_elem[irow + 2];
        result += m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        result += m_elem[irow] * v.m_elem[irow];
    }

    return result;
}

template <uint16_t t_row, typename t_type>
inline t_type Vector<t_row, t_type>::Inner(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    t_type result = 0;
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        result += m_elem[irow] * v.m_elem[irow];
        result += m_elem[irow + 1] * v.m_elem[irow + 1];
        result += m_elem[irow + 2] * v.m_elem[irow + 2];
        result += m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        result += m_elem[irow] * v.m_elem[irow];
    }

    return result;
}

template <uint16_t t_row, typename t_type>
inline t_type Vector<t_row, t_type>::Inner(const Vector3<t_type, t_row> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2]);
}

template <uint16_t t_row, typename t_type>
inline t_type Vector<t_row, t_type>::Inner(const Vector4<t_type, t_row> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3]);
}

template <uint16_t t_row, typename t_type>
inline t_type Vector<t_row, t_type>::Inner(const Vector6<t_type, t_row> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3] +
        m_elem[4] * v.m_elem[4] +
        m_elem[5] * v.m_elem[5]);
}

template <uint16_t t_row, typename t_type>
inline t_type Vector<t_row, t_type>::Inner(const Matrix<t_row, 1, t_type> &v) const
{
    t_type result = 0;
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        result += m_elem[irow] * v.m_elem[irow];
        result += m_elem[irow + 1] * v.m_elem[irow + 1];
        result += m_elem[irow + 2] * v.m_elem[irow + 2];
        result += m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        result += m_elem[irow] * v.m_elem[irow];
    }

    return result;
}

template <uint16_t t_row, typename t_type>
inline t_type Vector<t_row, t_type>::Inner(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    t_type result = 0;
    uint16_t cnt, irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        result += m_elem[irow] * v.m_elem[irow];
        result += m_elem[irow + 1] * v.m_elem[irow + 1];
        result += m_elem[irow + 2] * v.m_elem[irow + 2];
        result += m_elem[irow + 3] * v.m_elem[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        result += m_elem[irow] * v.m_elem[irow];
    }

    return result;
}

/* Comparison operators */
template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator==(const Vector &v) const
{
    uint16_t cnt, i = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return false;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator!=(const Vector &v) const
{
    uint16_t cnt, i = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return true;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator==(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return false;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator!=(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return true;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator==(const Vector3<t_type, t_row> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;

    return true;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator!=(const Vector3<t_type, t_row> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;

    return false;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator==(const Vector4<t_type, t_row> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;

    return true;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator!=(const Vector4<t_type, t_row> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;

    return false;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator==(const Vector6<t_type, t_row> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;
    if (std::abs(m_elem[4] - v.m_elem[4]) > m_tolerance) return false;
    if (std::abs(m_elem[5] - v.m_elem[5]) > m_tolerance) return false;

    return true;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator!=(const Vector6<t_type, t_row> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;
    if (std::abs(m_elem[4] - v.m_elem[4]) > m_tolerance) return true;
    if (std::abs(m_elem[5] - v.m_elem[5]) > m_tolerance) return true;

    return false;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator==(const Matrix<t_row, 1, t_type> &v) const
{
    uint16_t cnt, i = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return false;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator!=(const Matrix<t_row, 1, t_type> &v) const
{
    uint16_t cnt, i = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return true;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator==(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return false;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <uint16_t t_row, typename t_type>
inline bool Vector<t_row, t_type>::operator!=(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, i += 4)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - v.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - v.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - v.m_elem[i + 3]) > m_tolerance) return true;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - v.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::Print(const char endChar)
{
#if defined(ARDUINO)
    for (uint16_t irow = 0; irow < t_row; irow++)
    {
        Serial.printf("%7.3f\n", (t_type)m_elem[irow]);
    }
    Serial.write(endChar);
#else
    for (uint16_t irow = 0; irow < t_row; irow++)
    {
        printf("%10.6f\n", (t_type)m_elem[irow]);
    }
    printf("%c", endChar);
#endif
}

//-- Private Member Function ------------------------------------------------//
template <uint16_t t_row, typename t_type>
inline void Vector<t_row, t_type>::CrossProduct(const t_type *v)
{
    static_assert(t_row == 3, "This method is only for 3 x 1 vector");

    t_type elem[t_row];

    elem[0] = m_elem[1] * v[2] - m_elem[2] * v[1];
    elem[1] = m_elem[2] * v[0] - m_elem[0] * v[2];
    elem[2] = m_elem[0] * v[1] - m_elem[1] * v[0];

    m_elem[0] = elem[0];
    m_elem[1] = elem[1];
    m_elem[2] = elem[2];
}

//-- Template Function ------------------------------------------------------//
// scalar + vector
template <uint16_t row, typename type>
inline Vector<row, type> operator+(const type s, const Vector<row, type> &v)
{
    type vec[row];
    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = v.m_elem[irow] + s;
        vec[irow + 1] = v.m_elem[irow + 1] + s;
        vec[irow + 2] = v.m_elem[irow + 2] + s;
        vec[irow + 3] = v.m_elem[irow + 3] + s;
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = v.m_elem[irow] + s;
    }

    return Vector<row, type>(vec);
}

// scalar - vector
template <uint16_t row, typename type>
inline Vector<row, type> operator-(const type s, const Vector<row, type> &v)
{
    type vec[row];
    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = s - v.m_elem[irow];
        vec[irow + 1] = s - v.m_elem[irow + 1];
        vec[irow + 2] = s - v.m_elem[irow + 2];
        vec[irow + 3] = s - v.m_elem[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = s - v.m_elem[irow];
    }

    return Vector<row, type>(vec);
}

// scalar * vector
template <uint16_t row, typename type>
inline Vector<row, type> operator*(const type s, const Vector<row, type> &v)
{
    type vec[row];
    uint16_t cnt, irow = 0;

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[irow] = v.m_elem[irow] * s;
        vec[irow + 1] = v.m_elem[irow + 1] * s;
        vec[irow + 2] = v.m_elem[irow + 2] * s;
        vec[irow + 3] = v.m_elem[irow + 3] * s;
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[irow] = v.m_elem[irow] * s;
    }

    return Vector<row, type>(vec);
}

// scalar / vector
template <uint16_t row, typename type>
inline Vector<row, type> operator/(const type s, const Vector<row, type> &v)
{
    type vec[row];
    type den[row];
    uint16_t cnt, irow = 0;
    memcpy(den, v.m_elem, sizeof(den));

    for (cnt = row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        if (std::abs(den[irow]) < std::numeric_limits<type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<type>::epsilon();
            else den[irow] = std::numeric_limits<type>::epsilon();
        }
        if (std::abs(den[irow + 1]) < std::numeric_limits<type>::epsilon())
        {
            if (den[irow + 1] < 0) den[irow + 1] = -std::numeric_limits<type>::epsilon();
            else den[irow + 1] = std::numeric_limits<type>::epsilon();
        }
        if (std::abs(den[irow + 2]) < std::numeric_limits<type>::epsilon())
        {
            if (den[irow + 2] < 0) den[irow + 2] = -std::numeric_limits<type>::epsilon();
            else den[irow + 2] = std::numeric_limits<type>::epsilon();
        }
        if (std::abs(den[irow + 3]) < std::numeric_limits<type>::epsilon())
        {
            if (den[irow + 3] < 0) den[irow + 3] = -std::numeric_limits<type>::epsilon();
            else den[irow + 3] = std::numeric_limits<type>::epsilon();
        }
        vec[irow] = s / den[irow];
        vec[irow + 1] = s / den[irow + 1];
        vec[irow + 2] = s / den[irow + 2];
        vec[irow + 3] = s / den[irow + 3];
    }

    for (cnt = row % 4u; cnt > 0u; cnt--, irow++)
    {
        if (std::abs(den[irow]) < std::numeric_limits<type>::epsilon())
        {
            if (den[irow] < 0) den[irow] = -std::numeric_limits<type>::epsilon();
            else den[irow] = std::numeric_limits<type>::epsilon();
        }
        vec[irow] = s / den[irow];
    }

    return Vector<row, type>(vec);
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTVECTOR_TPP_