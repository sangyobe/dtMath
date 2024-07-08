/*!
\file       dtVector4.h
\brief      dtMath, 4x1 Vector class, lighter and faster than general vector class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTVECTOR4_TPP_
#define DTMATH_DTVECTOR4_TPP_

#include "dtVector4.h"

#include <cassert>

namespace dt
{
namespace Math
{

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row>::Vector4()
{
    m_elem[0] = 0;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row>::Vector4(const t_type *element)
{
    m_elem[0] = element[0];
    m_elem[1] = element[1];
    m_elem[2] = element[2];
    m_elem[3] = element[3];
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row>::Vector4(const t_type *element, const size_t n_byte)
{
    switch (n_byte / sizeof(t_type))
    {
    case 1:
        m_elem[0] = element[0];
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = 0;
        break;
    case 2:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = 0;
        m_elem[3] = 0;
        break;
    case 3:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = 0;
        break;
    default:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        break;
    }
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row>::Vector4(const t_type v0, const t_type v1, const t_type v2, const t_type v3)
{
    m_elem[0] = v0;
    m_elem[1] = v1;
    m_elem[2] = v2;
    m_elem[3] = v3;
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row>::Vector4(const Vector4 &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row>::Vector4(const Vector<t_row, t_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row>::Vector4(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row>::Vector4(const Matrix<t_row, 1, t_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row>::Vector4(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetZero()
{
    m_elem[0] = 0;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetFill(const t_type value)
{
    m_elem[0] = value;
    m_elem[1] = value;
    m_elem[2] = value;
    m_elem[3] = value;
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetElement(const t_type *element, const size_t n_byte)
{
    switch (n_byte / sizeof(t_type))
    {
    case 1:
        m_elem[0] = element[0];
        break;
    case 2:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        break;
    case 3:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        break;
    default:
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        break;
    }
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetElement(const t_type v0, const t_type v1, const t_type v2, const t_type v3)
{
    m_elem[0] = v0;
    m_elem[1] = v1;
    m_elem[2] = v2;
    m_elem[3] = v3;
}

// template <typename t_type, uint16_t t_row>
// inline void Vector4<t_type, t_row>::SetElement(const Vector4 &v)
// {
//     m_elem[0] = v.m_elem[0];
//     m_elem[1] = v.m_elem[1];
//     m_elem[2] = v.m_elem[2];
//     m_elem[3] = v.m_elem[3];
// }

// template <typename t_type, uint16_t t_row>
// inline void Vector4<t_type, t_row>::SetElement(const Vector<t_row, t_type> &v)
// {
//     m_elem[0] = v.m_elem[0];
//     m_elem[1] = v.m_elem[1];
//     m_elem[2] = v.m_elem[2];
//     m_elem[3] = v.m_elem[3];
// }

// template <typename t_type, uint16_t t_row>
// inline void Vector4<t_type, t_row>::SetElement(const Matrix<t_row, 1, t_type> &v)
// {
//     m_elem[0] = v.m_elem[0];
//     m_elem[1] = v.m_elem[1];
//     m_elem[2] = v.m_elem[2];
//     m_elem[3] = v.m_elem[3];
// }

// template <typename t_type, uint16_t t_row>
// inline void Vector4<t_type, t_row>::SetElement(const Matrix<0, 0, t_type> &v)
// {
//     assert(v.m_elem != nullptr && "Memory has not been allocated");
//     assert(v.m_row == t_row && "Row dimensions do not matched");
//     assert(v.m_col == 1 && "Col dimensions do not matched");

//     m_elem[0] = v.m_elem[0];
//     m_elem[1] = v.m_elem[1];
//     m_elem[2] = v.m_elem[2];
//     m_elem[3] = v.m_elem[3];
// }

template <typename t_type, uint16_t t_row>
template <uint16_t row>
inline void Vector4<t_type, t_row>::SetBlock(const uint16_t idxRow, const Vector<row, t_type> &v)
{
    assert(t_row > idxRow && "Index out of range");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (rowSz > row) rowSz = row;

    switch (rowSz)
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
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        break;
    }
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetBlock(const uint16_t idxRow, const t_type *v, const size_t n_byte)
{
    assert(t_row > idxRow && "Index out of range");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    uint16_t row = n_byte / sizeof(t_type);
    if (rowSz > row) rowSz = row;

    switch (rowSz)
    {
    case 1:
        m_elem[idxRow] = v[0];
        break;
    case 2:
        m_elem[idxRow] = v[0];
        m_elem[idxRow + 1] = v[1];
        break;
    case 3:
        m_elem[idxRow] = v[0];
        m_elem[idxRow + 1] = v[1];
        m_elem[idxRow + 2] = v[2];
        break;
    default:
        m_elem[0] = v[0];
        m_elem[1] = v[1];
        m_elem[2] = v[2];
        m_elem[3] = v[3];
        break;
    }
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetBlock(const uint16_t idxRow, const Vector3<t_type, 3> &v)
{
    assert(t_row > idxRow && "Index out of range");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (rowSz > 3) rowSz = 3;

    switch (rowSz)
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
    }
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetBlock(const uint16_t idxRow, const Vector4<t_type, 4> &v)
{
    assert(t_row > idxRow && "Index out of range");

    if (idxRow >= t_row) return;

    switch (idxRow)
    {
    case 3:
        m_elem[3] = v.m_elem[0];
        break;
    case 2:
        m_elem[2] = v.m_elem[0];
        m_elem[3] = v.m_elem[1];
        break;
    case 1:
        m_elem[1] = v.m_elem[0];
        m_elem[2] = v.m_elem[1];
        m_elem[3] = v.m_elem[2];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        break;
    }
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetBlock(const uint16_t idxRow, const Vector6<t_type, 6> &v)
{
    assert(t_row > idxRow && "Index out of range");

    if (idxRow >= t_row) return;

    switch (idxRow)
    {
    case 3:
        m_elem[3] = v.m_elem[0];
        break;
    case 2:
        m_elem[2] = v.m_elem[0];
        m_elem[3] = v.m_elem[1];
        break;
    case 1:
        m_elem[1] = v.m_elem[0];
        m_elem[2] = v.m_elem[1];
        m_elem[3] = v.m_elem[2];
        break;
    default:
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        break;
    }
}

template <typename t_type, uint16_t t_row>
template <uint16_t row>
inline void Vector4<t_type, t_row>::SetBlock(const uint16_t idxRow, const Matrix<row, 1, t_type> &v)
{
    assert(t_row > idxRow && "Index out of range");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (rowSz > row) rowSz = row;

    switch (rowSz)
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
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        break;
    }
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetBlock(const uint16_t idxRow, const Matrix<0, 0, t_type> &v)
{
    assert(t_row > idxRow && "Index out of range");
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (rowSz > v.m_row) rowSz = v.m_row;

    switch (rowSz)
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
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        break;
    }
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetBlock(const uint16_t idxRow, const Vector<0, t_type> &v)
{
    assert(t_row > idxRow && "Index out of range");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    if (idxRow >= t_row) return;

    uint16_t rowSz = t_row - idxRow;
    if (rowSz > v.m_row) rowSz = v.m_row;

    switch (rowSz)
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
        m_elem[0] = v.m_elem[0];
        m_elem[1] = v.m_elem[1];
        m_elem[2] = v.m_elem[2];
        m_elem[3] = v.m_elem[3];
        break;
    }
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetSwap(const uint16_t i, const uint16_t j)
{
    assert(t_row > i && "Index out of range");
    assert(t_row > j && "Index out of range");

    t_type elem = m_elem[i];
    m_elem[i] = m_elem[j];
    m_elem[j] = elem;
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::SetNormalize()
{
    t_type norm = std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);

    if (norm < std::numeric_limits<t_type>::epsilon())
        norm = std::numeric_limits<t_type>::epsilon();

    m_elem[0] /= norm;
    m_elem[1] /= norm;
    m_elem[2] /= norm;
    m_elem[3] /= norm;
}

template <typename t_type, uint16_t t_row>
inline const t_type *const Vector4<t_type, t_row>::GetElementsAddr() const
{
    return m_elem;
}

template <typename t_type, uint16_t t_row>
template <uint16_t row>
inline Vector<row, t_type> Vector4<t_type, t_row>::GetBlock(const uint16_t idx)
{
    assert(t_row > idx && "Index out of range");

    t_type elem[row]{0};
    uint16_t rowSize = t_row - idx;

    if (idx >= t_row) return Vector<row, t_type>(elem);
    if (rowSize > row) rowSize = row;

    memcpy(elem, &m_elem[idx], sizeof(t_type) * rowSize);

    return Vector<row, t_type>(elem);
}

template <typename t_type, uint16_t t_row>
inline Vector<0, t_type> Vector4<t_type, t_row>::GetBlock(const uint16_t idx, const uint16_t row)
{
    assert(t_row > idx && "Index out of range");

    Vector<0, t_type> v(row);
    uint16_t rowSize = t_row - idx;

    if (rowSize > row) rowSize = row;

    memcpy(v.m_elem, &m_elem[idx], sizeof(t_type) * rowSize);

    return v;
}

template <typename t_type, uint16_t t_row>
inline Vector3<t_type, 3> Vector4<t_type, t_row>::GetBlockVec3(const uint16_t idx)
{
    assert(t_row > idx && "Index out of range");

    t_type elem[3]{0};

    if (idx >= t_row) return Vector3<t_type, 3>(elem);

    switch (t_row - idx)
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

template <typename t_type, uint16_t t_row>
template <uint16_t row>
inline int8_t Vector4<t_type, t_row>::GetBlock(const uint16_t idx, Vector<row, t_type> &v)
{
    assert(t_row > idx && "Index out of range");

    uint16_t rowSize = t_row - idx;

    if (idx >= t_row) return -1;
    if (rowSize > row) rowSize = row;

    memcpy(v.m_elem, &m_elem[idx], sizeof(t_type) * rowSize);

    return 0;
}

template <typename t_type, uint16_t t_row>
inline int8_t Vector4<t_type, t_row>::GetBlock(const uint16_t idx, Vector<0, t_type> &v)
{
    assert(t_row > idx && "Index out of range");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t rowSize = t_row - idx;

    if (idx >= t_row) return -1;
    if (rowSize > v.m_row) rowSize = v.m_row;

    memcpy(v.m_elem, &m_elem[idx], sizeof(t_type) * rowSize);

    return 0;
}

template <typename t_type, uint16_t t_row>
inline int8_t Vector4<t_type, t_row>::GetBlockVec3(const uint16_t idx, Vector3<t_type, 3> &v)
{
    assert(t_row > idx && "Index out of range");

    if (idx >= t_row) return -1;

    switch (t_row - idx)
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

template <typename t_type, uint16_t t_row>
inline t_type Vector4<t_type, t_row>::GetNorm() const
{
    return std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline t_type Vector4<t_type, t_row>::GetSqNorm() const
{
    return (
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline t_type Vector4<t_type, t_row>::GetLpNorm(const int p) const
{
    t_type powSum =
        std::pow(std::abs(m_elem[0]), (t_type)p) +
        std::pow(std::abs(m_elem[1]), (t_type)p) +
        std::pow(std::abs(m_elem[2]), (t_type)p) +
        std::pow(std::abs(m_elem[3]), (t_type)p);

    return std::pow(powSum, (t_type)1 / p);
}

template <typename t_type, uint16_t t_row>
inline t_type Vector4<t_type, t_row>::GetSum() const
{
    return (
        m_elem[0] +
        m_elem[1] +
        m_elem[2] +
        m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::GetNormalized() const
{
    t_type norm = std::sqrt(
        m_elem[0] * m_elem[0] +
        m_elem[1] * m_elem[1] +
        m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3]);

    if (norm < std::numeric_limits<t_type>::epsilon())
        norm = std::numeric_limits<t_type>::epsilon();

    return Vector4(
        m_elem[0] / norm,
        m_elem[1] / norm,
        m_elem[2] / norm,
        m_elem[3] / norm);
}

template <typename t_type, uint16_t t_row>
inline Matrix<1, t_row, t_type> Vector4<t_type, t_row>::Transpose() const
{
    return Matrix<1, t_row, t_type>(m_elem);
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::Transpose(Matrix<1, t_row, t_type> &m) const
{
    memcpy(m.m_elem, m_elem, sizeof(t_type) * t_row);
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::Transpose(Matrix<0, 0, t_type> &m) const
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == 1 && "Row dimensions do not matched");
    assert(m.m_col == t_row && "Col dimensions do not matched");

    memcpy(m.m_elem, m_elem, sizeof(t_type) * t_row);
}

/* Member access operators */
template <typename t_type, uint16_t t_row>
inline t_type &Vector4<t_type, t_row>::operator()(uint16_t irow)
{
    assert(irow < t_row && "Index out of range");
    return m_elem[irow];
}

template <typename t_type, uint16_t t_row>
inline const t_type &Vector4<t_type, t_row>::operator()(uint16_t irow) const
{
    assert(irow < t_row && "Index out of range");
    return m_elem[irow];
}

/* Assignment operators */
template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator=(const Vector4 &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator+=(const Vector4 &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator-=(const Vector4 &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::CWiseMulEq(const Vector4 &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::CWiseDivEq(const Vector4 &v)
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

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator=(const Vector<t_row, t_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator+=(const Vector<t_row, t_type> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator-=(const Vector<t_row, t_type> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::CWiseMulEq(const Vector<t_row, t_type> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::CWiseDivEq(const Vector<t_row, t_type> &v)
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

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator=(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator+=(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator-=(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == v.m_row && "Row dimensions do not matched");

    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator=(const Matrix<t_row, 1, t_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator+=(const Matrix<t_row, 1, t_type> &v)
{
    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator-=(const Matrix<t_row, 1, t_type> &v)
{
    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::CWiseMulEq(const Matrix<t_row, 1, t_type> &v)
{
    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::CWiseDivEq(const Matrix<t_row, 1, t_type> &v)
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

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator=(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    m_elem[0] = v.m_elem[0];
    m_elem[1] = v.m_elem[1];
    m_elem[2] = v.m_elem[2];
    m_elem[3] = v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator+=(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    m_elem[0] += v.m_elem[0];
    m_elem[1] += v.m_elem[1];
    m_elem[2] += v.m_elem[2];
    m_elem[3] += v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator-=(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    m_elem[0] -= v.m_elem[0];
    m_elem[1] -= v.m_elem[1];
    m_elem[2] -= v.m_elem[2];
    m_elem[3] -= v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::CWiseMulEq(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    m_elem[0] *= v.m_elem[0];
    m_elem[1] *= v.m_elem[1];
    m_elem[2] *= v.m_elem[2];
    m_elem[3] *= v.m_elem[3];

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::CWiseDivEq(const Matrix<0, 0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

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

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator=(const t_type s)
{
    m_elem[0] = s;
    m_elem[1] = s;
    m_elem[2] = s;
    m_elem[3] = s;

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator+=(const t_type s)
{
    m_elem[0] += s;
    m_elem[1] += s;
    m_elem[2] += s;
    m_elem[3] += s;

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator-=(const t_type s)
{
    m_elem[0] -= s;
    m_elem[1] -= s;
    m_elem[2] -= s;
    m_elem[3] -= s;

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator*=(const t_type s)
{
    m_elem[0] *= s;
    m_elem[1] *= s;
    m_elem[2] *= s;
    m_elem[3] *= s;

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> &Vector4<t_type, t_row>::operator/=(const t_type s)
{
    t_type den = s;

    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<t_type>::epsilon();
        else
            den = std::numeric_limits<t_type>::epsilon();
    }

    m_elem[0] /= den;
    m_elem[1] /= den;
    m_elem[2] /= den;
    m_elem[3] /= den;

    return (*this);
}

template <typename t_type, uint16_t t_row>
inline CommaInit<t_row, t_type> Vector4<t_type, t_row>::operator<<(const t_type s)
{
    m_elem[0] = s;
    return CommaInit<t_row, t_type>(m_elem);
}

/* Arithmetic operators */
template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator-() const
{
    return Vector4(
        -m_elem[0],
        -m_elem[1],
        -m_elem[2],
        -m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator+(const Vector4 &v) const
{
    return Vector4(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2],
        m_elem[3] + v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator-(const Vector4 &v) const
{
    return Vector4(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2],
        m_elem[3] - v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::CWiseMul(const Vector4 &v) const
{
    return Vector4(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2],
        m_elem[3] * v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::CWiseDiv(const Vector4 &v) const
{
    t_type den[4];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];

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

    return Vector4(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2],
        m_elem[3] / den[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator+(const Vector<t_row, t_type> &v) const
{
    return Vector4(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2],
        m_elem[3] + v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator-(const Vector<t_row, t_type> &v) const
{
    return Vector4(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2],
        m_elem[3] - v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::CWiseMul(const Vector<t_row, t_type> &v) const
{
    return Vector4(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2],
        m_elem[3] * v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::CWiseDiv(const Vector<t_row, t_type> &v) const
{
    t_type den[4];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];

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

    return Vector4(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2],
        m_elem[3] / den[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator+(const Matrix<t_row, 1, t_type> &v) const
{
    return Vector4(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2],
        m_elem[3] + v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator-(const Matrix<t_row, 1, t_type> &v) const
{
    return Vector4(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2],
        m_elem[3] - v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::CWiseMul(const Matrix<t_row, 1, t_type> &v) const
{
    return Vector4(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2],
        m_elem[3] * v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::CWiseDiv(const Matrix<t_row, 1, t_type> &v) const
{
    t_type den[4];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];

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

    return Vector4(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2],
        m_elem[3] / den[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator+(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    return Vector4(
        m_elem[0] + v.m_elem[0],
        m_elem[1] + v.m_elem[1],
        m_elem[2] + v.m_elem[2],
        m_elem[3] + v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator-(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    return Vector4(
        m_elem[0] - v.m_elem[0],
        m_elem[1] - v.m_elem[1],
        m_elem[2] - v.m_elem[2],
        m_elem[3] - v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::CWiseMul(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    return Vector4(
        m_elem[0] * v.m_elem[0],
        m_elem[1] * v.m_elem[1],
        m_elem[2] * v.m_elem[2],
        m_elem[3] * v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::CWiseDiv(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    t_type den[4];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];

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

    return Vector4(
        m_elem[0] / den[0],
        m_elem[1] / den[1],
        m_elem[2] / den[2],
        m_elem[3] / den[3]);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator+(const t_type s) const
{
    return Vector4(
        m_elem[0] + s,
        m_elem[1] + s,
        m_elem[2] + s,
        m_elem[3] + s);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator-(const t_type s) const
{
    return Vector4(
        m_elem[0] - s,
        m_elem[1] - s,
        m_elem[2] - s,
        m_elem[3] - s);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator*(const t_type s) const
{
    return Vector4(
        m_elem[0] * s,
        m_elem[1] * s,
        m_elem[2] * s,
        m_elem[3] * s);
}

template <typename t_type, uint16_t t_row>
inline Vector4<t_type, t_row> Vector4<t_type, t_row>::operator/(const t_type s) const
{
    t_type den = s;

    if (std::abs(den) < std::numeric_limits<t_type>::epsilon())
    {
        if (den < 0)
            den = -std::numeric_limits<t_type>::epsilon();
        else
            den = std::numeric_limits<t_type>::epsilon();
    }

    return Vector4(
        m_elem[0] / den,
        m_elem[1] / den,
        m_elem[2] / den,
        m_elem[3] / den);
}

template <typename t_type, uint16_t t_row>
template <uint16_t col>
inline Matrix<t_row, col, t_type> Vector4<t_type, t_row>::operator*(const Matrix<1, col, t_type> &m) const
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

template <typename t_type, uint16_t t_row>
inline Matrix<0, 0, t_type> Vector4<t_type, t_row>::operator*(const Matrix<0, 0, t_type> &m) const
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

template <typename t_type, uint16_t t_row>
inline Matrix<t_row, t_row, t_type> Vector4<t_type, t_row>::Outer(const Vector4 &v) const
{
    t_type mat[t_row * t_row];
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; irow++)
    {
        mat[irow * t_row] = m_elem[irow] * v.m_elem[0];
        mat[irow * t_row + 1] = m_elem[irow] * v.m_elem[1];
        mat[irow * t_row + 2] = m_elem[irow] * v.m_elem[2];
        mat[irow * t_row + 3] = m_elem[irow] * v.m_elem[3];
    }

    return Matrix<t_row, t_row, t_type>(mat);
}

template <typename t_type, uint16_t t_row>
inline Matrix<t_row, 3, t_type> Vector4<t_type, t_row>::Outer(const Vector3<t_type, 3> &v) const
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

template <typename t_type, uint16_t t_row>
inline Matrix<t_row, 6, t_type> Vector4<t_type, t_row>::Outer(const Vector6<t_type, 6> &v) const
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

template <typename t_type, uint16_t t_row>
inline Matrix<0, 0, t_type> Vector4<t_type, t_row>::Outer(const Vector<0, t_type> &v) const
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

template <typename t_type, uint16_t t_row>
inline Matrix<0, 0, t_type> Vector4<t_type, t_row>::Outer(const Matrix<0, 0, t_type> &v) const
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

template <typename t_type, uint16_t t_row>
template <uint16_t row>
inline Matrix<t_row, row, t_type> Vector4<t_type, t_row>::Outer(const Vector<row, t_type> &v) const
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

template <typename t_type, uint16_t t_row>
template <uint16_t row>
inline Matrix<t_row, row, t_type> Vector4<t_type, t_row>::Outer(const Matrix<row, 1, t_type> &v) const
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

template <typename t_type, uint16_t t_row>
inline t_type Vector4<t_type, t_row>::Inner(const Vector4 &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline t_type Vector4<t_type, t_row>::Inner(const Vector<t_row, t_type> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline t_type Vector4<t_type, t_row>::Inner(const Matrix<t_row, 1, t_type> &v) const
{
    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline t_type Vector4<t_type, t_row>::Inner(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3]);
}

template <typename t_type, uint16_t t_row>
inline t_type Vector4<t_type, t_row>::Inner(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    return (
        m_elem[0] * v.m_elem[0] +
        m_elem[1] * v.m_elem[1] +
        m_elem[2] * v.m_elem[2] +
        m_elem[3] * v.m_elem[3]);
}

/* Comparison operators */
template <typename t_type, uint16_t t_row>
inline bool Vector4<t_type, t_row>::operator==(const Vector4 &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;

    return true;
}

template <typename t_type, uint16_t t_row>
inline bool Vector4<t_type, t_row>::operator!=(const Vector4 &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;

    return false;
}

template <typename t_type, uint16_t t_row>
inline bool Vector4<t_type, t_row>::operator==(const Vector<t_row, t_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;

    return true;
}

template <typename t_type, uint16_t t_row>
inline bool Vector4<t_type, t_row>::operator!=(const Vector<t_row, t_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;

    return false;
}

template <typename t_type, uint16_t t_row>
inline bool Vector4<t_type, t_row>::operator==(const Matrix<t_row, 1, t_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;

    return true;
}

template <typename t_type, uint16_t t_row>
inline bool Vector4<t_type, t_row>::operator!=(const Matrix<t_row, 1, t_type> &v) const
{
    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;

    return false;
}

template <typename t_type, uint16_t t_row>
inline bool Vector4<t_type, t_row>::operator==(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;

    return true;
}

template <typename t_type, uint16_t t_row>
inline bool Vector4<t_type, t_row>::operator!=(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");

    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;

    return false;
}

template <typename t_type, uint16_t t_row>
inline bool Vector4<t_type, t_row>::operator==(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return false;

    return true;
}

template <typename t_type, uint16_t t_row>
inline bool Vector4<t_type, t_row>::operator!=(const Matrix<0, 0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == t_row && "Row dimensions do not matched");
    assert(v.m_col == 1 && "Col dimensions do not matched");

    if (std::abs(m_elem[0] - v.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - v.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - v.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - v.m_elem[3]) > m_tolerance) return true;

    return false;
}

template <typename t_type, uint16_t t_row>
inline void Vector4<t_type, t_row>::Print(const char endChar)
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
        printf("%7.3f\n", (t_type)m_elem[irow]);
    }
    printf("%c", endChar);
#endif
}

//-- Template Function ------------------------------------------------------//
// scalar + vector
template <typename type, uint16_t row>
inline Vector4<type, row> operator+(const type s, const Vector4<type, row> &v)
{
    return Vector4<type, row>(
        v.m_elem[0] + s,
        v.m_elem[1] + s,
        v.m_elem[2] + s,
        v.m_elem[3] + s);
}

// scalar - vector
template <typename type, uint16_t row>
inline Vector4<type, row> operator-(const type s, const Vector4<type, row> &v)
{
    return Vector4<type, row>(
        s - v.m_elem[0],
        s - v.m_elem[1],
        s - v.m_elem[2],
        s - v.m_elem[3]);
}

// scalar * vector
template <typename type, uint16_t row>
inline Vector4<type, row> operator*(const type s, const Vector4<type, row> &v)
{
    return Vector4<type, row>(
        v.m_elem[0] * s,
        v.m_elem[1] * s,
        v.m_elem[2] * s,
        v.m_elem[3] * s);
}

// scalar / vector
template <typename type, uint16_t row>
inline Vector4<type, row> operator/(const type s, const Vector4<type, row> &v)
{
    type den[4];

    den[0] = v.m_elem[0];
    den[1] = v.m_elem[1];
    den[2] = v.m_elem[2];
    den[3] = v.m_elem[3];

    if (std::abs(den[0]) < std::numeric_limits<type>::epsilon())
    {
        if (den[0] < 0)
            den[0] = -std::numeric_limits<type>::epsilon();
        else
            den[0] = std::numeric_limits<type>::epsilon();
    }
    if (std::abs(den[1]) < std::numeric_limits<type>::epsilon())
    {
        if (den[1] < 0)
            den[1] = -std::numeric_limits<type>::epsilon();
        else
            den[1] = std::numeric_limits<type>::epsilon();
    }
    if (std::abs(den[2]) < std::numeric_limits<type>::epsilon())
    {
        if (den[2] < 0)
            den[2] = -std::numeric_limits<type>::epsilon();
        else
            den[2] = std::numeric_limits<type>::epsilon();
    }
    if (std::abs(den[3]) < std::numeric_limits<type>::epsilon())
    {
        if (den[3] < 0)
            den[3] = -std::numeric_limits<type>::epsilon();
        else
            den[3] = std::numeric_limits<type>::epsilon();
    }

    return Vector4<type, row>(
        s / den[0],
        s / den[1],
        s / den[2],
        s / den[3]);
}

typedef Vector4<> dtVec4;

} // namespace Math
} // namespace dt

#endif // DTMATH_DTVECTOR4_TPP_