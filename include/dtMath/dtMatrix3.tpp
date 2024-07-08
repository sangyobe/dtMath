/*!
\file       dtMatrix3.h
\brief      dtMath, 3x3 Matrix class, lighter and faster than general matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATRIX3_TPP_
#define DTMATH_DTMATRIX3_TPP_

#include "dtMatrix3.h"

namespace dt
{
namespace Math
{

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col>::Matrix3()
{
    m_elem[0] = 0;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
    m_elem[4] = 0;
    m_elem[5] = 0;
    m_elem[6] = 0;
    m_elem[7] = 0;
    m_elem[8] = 0;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col>::Matrix3(const t_type *element)
{
    m_elem[0] = element[0];
    m_elem[1] = element[1];
    m_elem[2] = element[2];
    m_elem[3] = element[3];
    m_elem[4] = element[4];
    m_elem[5] = element[5];
    m_elem[6] = element[6];
    m_elem[7] = element[7];
    m_elem[8] = element[8];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col>::Matrix3(const t_type *element, const size_t n_byte)
{
    size_t matSz = sizeof(t_type) * t_row * t_col;

    if (matSz <= n_byte)
    {
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        m_elem[4] = element[4];
        m_elem[5] = element[5];
        m_elem[6] = element[6];
        m_elem[7] = element[7];
        m_elem[8] = element[8];
    }
    else
    {
        memset(m_elem, 0, matSz);
        memcpy(m_elem, element, n_byte);
    }
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col>::Matrix3(const char c, const t_type *element, const size_t n_byte)
{
    if (c == 'a')
    {
        size_t matSz = sizeof(t_type) * t_row * t_col;

        if (matSz <= n_byte)
        {
            m_elem[0] = element[0];
            m_elem[1] = element[1];
            m_elem[2] = element[2];
            m_elem[3] = element[3];
            m_elem[4] = element[4];
            m_elem[5] = element[5];
            m_elem[6] = element[6];
            m_elem[7] = element[7];
            m_elem[8] = element[8];
        }
        else
        {
            memset(m_elem, 0, matSz);
            memcpy(m_elem, element, n_byte);
        }
    }
    else if (c == 'd')
    {
        switch (n_byte / sizeof(t_type))
        {
        case 1:
            m_elem[0] = element[0];
            m_elem[1] = 0;
            m_elem[2] = 0;
            m_elem[3] = 0;
            m_elem[4] = 0;
            m_elem[5] = 0;
            m_elem[6] = 0;
            m_elem[7] = 0;
            m_elem[8] = 0;
            break;
        case 2:
            m_elem[0] = element[0];
            m_elem[1] = 0;
            m_elem[2] = 0;
            m_elem[3] = 0;
            m_elem[4] = element[1];
            m_elem[5] = 0;
            m_elem[6] = 0;
            m_elem[7] = 0;
            m_elem[8] = 0;
            break;
        default:
            m_elem[0] = element[0];
            m_elem[1] = 0;
            m_elem[2] = 0;
            m_elem[3] = 0;
            m_elem[4] = element[1];
            m_elem[5] = 0;
            m_elem[6] = 0;
            m_elem[7] = 0;
            m_elem[8] = element[2];
            break;
        }
    }
    else
    {
        m_elem[0] = 0;
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = 0;
        m_elem[4] = 0;
        m_elem[5] = 0;
        m_elem[6] = 0;
        m_elem[7] = 0;
        m_elem[8] = 0;
    }
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col>::Matrix3(
    const t_type m00, const t_type m01, const t_type m02,
    const t_type m10, const t_type m11, const t_type m12,
    const t_type m20, const t_type m21, const t_type m22)
{
    m_elem[0] = m00;
    m_elem[1] = m01;
    m_elem[2] = m02;
    m_elem[3] = m10;
    m_elem[4] = m11;
    m_elem[5] = m12;
    m_elem[6] = m20;
    m_elem[7] = m21;
    m_elem[8] = m22;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col>::Matrix3(const Matrix3 &m)
{
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col>::Matrix3(const Rotation<t_type, t_row, t_col> &m)
{
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col>::Matrix3(const Matrix<t_row, t_col, t_type> &m)
{
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col>::Matrix3(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == t_row && "Row dimensions do not matched");
    assert(m.m_col == t_col && "Col dimensions do not matched");

    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetZero()
{
    m_elem[0] = 0;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
    m_elem[4] = 0;
    m_elem[5] = 0;
    m_elem[6] = 0;
    m_elem[7] = 0;
    m_elem[8] = 0;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetIdentity()
{
    m_elem[0] = 1;
    m_elem[1] = 0;
    m_elem[2] = 0;
    m_elem[3] = 0;
    m_elem[4] = 1;
    m_elem[5] = 0;
    m_elem[6] = 0;
    m_elem[7] = 0;
    m_elem[8] = 1;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetDiagonal(const t_type d1, const t_type d2, const t_type d3)
{
    m_elem[0] = d1;
    m_elem[4] = d2;
    m_elem[8] = d3;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetDiagonal(const t_type *element, const size_t n_byte)
{
    switch (n_byte / sizeof(t_type))
    {
    case 1:
        m_elem[0] = element[0];
        break;
    case 2:
        m_elem[0] = element[0];
        m_elem[4] = element[1];
        break;
    default:
        m_elem[0] = element[0];
        m_elem[4] = element[1];
        m_elem[8] = element[2];
        break;
    }
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetDiagonal(const Vector<t_row, t_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[4] = v.m_elem[1];
    m_elem[8] = v.m_elem[2];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetDiagonal(const Vector<0, t_type> &v)
{
    assert(t_row == v.m_row && "Check dimensions");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    m_elem[0] = v.m_elem[0];
    m_elem[4] = v.m_elem[1];
    m_elem[8] = v.m_elem[2];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetDiagonal(const Vector3<t_type, t_row> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[4] = v.m_elem[1];
    m_elem[8] = v.m_elem[2];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetFill(const t_type value)
{
    m_elem[0] = value;
    m_elem[1] = value;
    m_elem[2] = value;
    m_elem[3] = value;
    m_elem[4] = value;
    m_elem[5] = value;
    m_elem[6] = value;
    m_elem[7] = value;
    m_elem[8] = value;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetElement(const t_type *element, const size_t n_byte)
{
    size_t matSz = sizeof(t_type) * t_row * t_col;

    if (matSz <= n_byte)
    {
        m_elem[0] = element[0];
        m_elem[1] = element[1];
        m_elem[2] = element[2];
        m_elem[3] = element[3];
        m_elem[4] = element[4];
        m_elem[5] = element[5];
        m_elem[6] = element[6];
        m_elem[7] = element[7];
        m_elem[8] = element[8];
    }
    else
        memcpy(m_elem, element, n_byte);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetElement(
    const t_type m00, const t_type m01, const t_type m02,
    const t_type m10, const t_type m11, const t_type m12,
    const t_type m20, const t_type m21, const t_type m22)
{
    m_elem[0] = m00;
    m_elem[1] = m01;
    m_elem[2] = m02;
    m_elem[3] = m10;
    m_elem[4] = m11;
    m_elem[5] = m12;
    m_elem[6] = m20;
    m_elem[7] = m21;
    m_elem[8] = m22;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetElement(const Matrix3 &m)
{
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetElement(const Rotation<t_type, t_row, t_col> &m)
{
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetElement(const Matrix<t_row, t_col, t_type> &m)
{
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetElement(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == t_row && "Row dimensions do not matched");
    assert(m.m_col == t_col && "Col dimensions do not matched");

    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
template <uint16_t col>
inline void Matrix3<t_type, t_row, t_col>::SetRowVec(const uint16_t idxRow, const Vector<col, t_type> &v)
{
    uint16_t maxCol = (t_col < col) ? t_col : col;

    for (uint16_t icol = 0; icol < maxCol; icol++)
    {
        m_elem[idxRow * t_col + icol] = v.m_elem[icol];
    }
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetRowVec(const uint16_t idxRow, const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t maxCol = (t_col < v.m_row) ? t_col : v.m_row;

    for (uint16_t icol = 0; icol < maxCol; icol++)
    {
        m_elem[idxRow * t_col + icol] = v.m_elem[icol];
    }
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetRowVec(const uint16_t idxRow, const Vector3<t_type, t_col> &v)
{
    m_elem[idxRow * t_col] = v.m_elem[0];
    m_elem[idxRow * t_col + 1] = v.m_elem[1];
    m_elem[idxRow * t_col + 2] = v.m_elem[2];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetRowVec(const uint16_t idxRow, const t_type *v, const size_t n_byte)
{
    uint16_t col = n_byte / sizeof(t_type);
    uint16_t maxCol = (t_col < col) ? t_col : col;

    for (uint16_t icol = 0; icol < maxCol; icol++)
    {
        m_elem[idxRow * t_col + icol] = v[icol];
    }
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
template <uint16_t row>
inline void Matrix3<t_type, t_row, t_col>::SetColVec(const uint16_t idxCol, const Vector<row, t_type> &v)
{
    uint16_t maxRow = (t_row < row) ? t_row : row;

    for (uint16_t irow = 0; irow < maxRow; irow++)
    {
        m_elem[irow * t_col + idxCol] = v.m_elem[irow];
    }
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetColVec(const uint16_t idxCol, const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t maxRow = (t_row < v.m_row) ? t_row : v.m_row;

    for (uint16_t irow = 0; irow < maxRow; irow++)
    {
        m_elem[irow * t_col + idxCol] = v.m_elem[irow];
    }
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetColVec(const uint16_t idxCol, const Vector3<t_type, t_row> &v)
{
    m_elem[idxCol] = v.m_elem[0];
    m_elem[t_col + idxCol] = v.m_elem[1];
    m_elem[t_col * 2 + idxCol] = v.m_elem[2];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetColVec(const uint16_t idxCol, const t_type *v, const size_t n_byte)
{
    uint16_t row = n_byte / sizeof(t_type);
    uint16_t maxRow = (t_row < row) ? t_row : row;

    for (uint16_t irow = 0; irow < maxRow; irow++)
    {
        m_elem[irow * t_col + idxCol] = v[irow];
    }
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2)
{
    t_type tmpElem[t_col];

    tmpElem[0] = m_elem[idxRow1 * t_col];
    m_elem[idxRow1 * t_col] = m_elem[idxRow2 * t_col];
    m_elem[idxRow2 * t_col] = tmpElem[0];

    tmpElem[1] = m_elem[idxRow1 * t_col + 1];
    m_elem[idxRow1 * t_col + 1] = m_elem[idxRow2 * t_col + 1];
    m_elem[idxRow2 * t_col + 1] = tmpElem[1];

    tmpElem[2] = m_elem[idxRow1 * t_col + 2];
    m_elem[idxRow1 * t_col + 2] = m_elem[idxRow2 * t_col + 2];
    m_elem[idxRow2 * t_col + 2] = tmpElem[2];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2)
{
    t_type tmpElem[t_row];

    tmpElem[0] = m_elem[idxCol1];
    m_elem[idxCol1] = m_elem[idxCol2];
    m_elem[idxCol2] = tmpElem[0];

    tmpElem[1] = m_elem[t_col + idxCol1];
    m_elem[t_col + idxCol1] = m_elem[t_col + idxCol2];
    m_elem[t_col + idxCol2] = tmpElem[1];

    tmpElem[2] = m_elem[2 * t_col + idxCol1];
    m_elem[2 * t_col + idxCol1] = m_elem[2 * t_col + idxCol2];
    m_elem[2 * t_col + idxCol2] = tmpElem[2];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline const t_type *const Matrix3<t_type, t_row, t_col>::GetElementsAddr() const
{
    return m_elem;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector3<t_type, t_col> Matrix3<t_type, t_row, t_col>::GetRowVec(const uint16_t idxRow) const
{
    return Vector3<t_type, t_row>(
        m_elem[idxRow * t_col],
        m_elem[idxRow * t_col + 1],
        m_elem[idxRow * t_col + 2]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector3<t_type, t_row> Matrix3<t_type, t_row, t_col>::GetColVec(const uint16_t idxCol) const
{
    return Vector3<t_type, t_row>(
        m_elem[idxCol],
        m_elem[1 * t_col + idxCol],
        m_elem[2 * t_col + idxCol]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline int8_t Matrix3<t_type, t_row, t_col>::GetRowVec(const uint16_t idxRow, Vector3<t_type, t_col> &v) const
{
    v.m_elem[0] = m_elem[idxRow * t_col];
    v.m_elem[1] = m_elem[idxRow * t_col + 1];
    v.m_elem[2] = m_elem[idxRow * t_col + 2];

    return 0;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline int8_t Matrix3<t_type, t_row, t_col>::GetColVec(const uint16_t idxCol, Vector3<t_type, t_row> &v) const
{
    v.m_elem[0] = m_elem[idxCol];
    v.m_elem[1] = m_elem[1 * t_col + idxCol];
    v.m_elem[2] = m_elem[2 * t_col + idxCol];

    return 0;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::Transpose() const
{
    return Matrix3(
        m_elem[0], m_elem[3], m_elem[6],
        m_elem[1], m_elem[4], m_elem[7],
        m_elem[2], m_elem[5], m_elem[8]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline t_type Matrix3<t_type, t_row, t_col>::Trace() const
{
    return (m_elem[0] + m_elem[4] + m_elem[8]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline t_type Matrix3<t_type, t_row, t_col>::GetNorm() const
{
    t_type sqSum =
        m_elem[0] * m_elem[0] + m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3] + m_elem[4] * m_elem[4] + m_elem[5] * m_elem[5] +
        m_elem[6] * m_elem[6] + m_elem[7] * m_elem[7] + m_elem[8] * m_elem[8];

    return std::sqrt(sqSum);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline t_type Matrix3<t_type, t_row, t_col>::GetSqNorm() const
{
    t_type sqSum =
        m_elem[0] * m_elem[0] + m_elem[1] * m_elem[1] + m_elem[2] * m_elem[2] +
        m_elem[3] * m_elem[3] + m_elem[4] * m_elem[4] + m_elem[5] * m_elem[5] +
        m_elem[6] * m_elem[6] + m_elem[7] * m_elem[7] + m_elem[8] * m_elem[8];

    return sqSum;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline NoPivLU<t_row, t_col, t_type> Matrix3<t_type, t_row, t_col>::GetNoPivLU() const
{
    return NoPivLU<t_row, t_col, t_type>(*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline PartialPivLU<t_row, t_col, t_type> Matrix3<t_type, t_row, t_col>::GetPartialPivLU() const
{
    return PartialPivLU<t_row, t_col, t_type>(*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline FullPivLU<t_row, t_col, t_type> Matrix3<t_type, t_row, t_col>::GetFullPivLU() const
{
    return FullPivLU<t_row, t_col, t_type>(*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline LLT<t_row, t_col, t_type> Matrix3<t_type, t_row, t_col>::GetLLT() const
{
    return LLT<t_row, t_col, t_type>(*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline LDLT<t_row, t_col, t_type> Matrix3<t_type, t_row, t_col>::GetLDLT() const
{
    return LDLT<t_row, t_col, t_type>(*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline QR<t_row, t_col, t_type> Matrix3<t_type, t_row, t_col>::GetQR() const
{
    return QR<t_row, t_col, t_type>(*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline SVD<t_row, t_col, t_type> Matrix3<t_type, t_row, t_col>::GetSVD() const
{
    return SVD<t_row, t_col, t_type>(*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::Inv(int8_t *isOk) const
{
    return Matrix3<t_type, t_row, t_col>(PartialPivLU<t_row, t_col, t_type>(*this).InverseArray(isOk));
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::PInv(int8_t *isOk, t_type tolerance) const
{
    return Matrix3<t_type, t_row, t_col>(SVD<t_row, t_col, t_type>(*this).InverseArray(isOk, tolerance));
}

/* Member access operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline t_type &Matrix3<t_type, t_row, t_col>::operator()(uint16_t irow, uint16_t icol)
{
    assert(irow < t_row && "Index out of range");
    assert(icol < t_col && "Index out of range");

    return m_elem[irow * t_col + icol];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline const t_type &Matrix3<t_type, t_row, t_col>::operator()(uint16_t irow, uint16_t icol) const
{
    assert(irow < t_row && "Index out of range");
    assert(icol < t_col && "Index out of range");

    return m_elem[irow * t_col + icol];
}

/* Assignment operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator=(const Matrix3 &m)
{
    // memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];

    return (*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator+=(const Matrix3 &m)
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator-=(const Matrix3 &m)
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator=(const Rotation<t_type, t_row, t_col> &m)
{
    // memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];

    return (*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator+=(const Rotation<t_type, t_row, t_col> &m)
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator-=(const Rotation<t_type, t_row, t_col> &m)
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator=(const Matrix<t_row, t_col, t_type> &m)
{
    // memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];

    return (*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator+=(const Matrix<t_row, t_col, t_type> &m)
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator-=(const Matrix<t_row, t_col, t_type> &m)
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator=(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == t_row && "Row dimensions do not matched");
    assert(m.m_col == t_col && "Col dimensions do not matched");

    m_elem[0] = m.m_elem[0];
    m_elem[1] = m.m_elem[1];
    m_elem[2] = m.m_elem[2];
    m_elem[3] = m.m_elem[3];
    m_elem[4] = m.m_elem[4];
    m_elem[5] = m.m_elem[5];
    m_elem[6] = m.m_elem[6];
    m_elem[7] = m.m_elem[7];
    m_elem[8] = m.m_elem[8];

    return (*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator+=(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == t_row && "Row dimensions do not matched");
    assert(m.m_col == t_col && "Col dimensions do not matched");

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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator-=(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == t_row && "Row dimensions do not matched");
    assert(m.m_col == t_col && "Col dimensions do not matched");

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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator=(const t_type s)
{
    m_elem[0] = s;
    m_elem[1] = s;
    m_elem[2] = s;
    m_elem[3] = s;
    m_elem[4] = s;
    m_elem[5] = s;
    m_elem[6] = s;
    m_elem[7] = s;
    m_elem[8] = s;

    return (*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator+=(const t_type s)
{
    m_elem[0] += s;
    m_elem[1] += s;
    m_elem[2] += s;
    m_elem[3] += s;
    m_elem[4] += s;
    m_elem[5] += s;
    m_elem[6] += s;
    m_elem[7] += s;
    m_elem[8] += s;

    return (*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator-=(const t_type s)
{
    m_elem[0] -= s;
    m_elem[1] -= s;
    m_elem[2] -= s;
    m_elem[3] -= s;
    m_elem[4] -= s;
    m_elem[5] -= s;
    m_elem[6] -= s;
    m_elem[7] -= s;
    m_elem[8] -= s;

    return (*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator*=(const t_type s)
{
    m_elem[0] *= s;
    m_elem[1] *= s;
    m_elem[2] *= s;
    m_elem[3] *= s;
    m_elem[4] *= s;
    m_elem[5] *= s;
    m_elem[6] *= s;
    m_elem[7] *= s;
    m_elem[8] *= s;

    return (*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> &Matrix3<t_type, t_row, t_col>::operator/=(const t_type s)
{
    t_type scalar = s;

    if (std::abs(scalar) < std::numeric_limits<t_type>::epsilon())
    {
        if (scalar < 0)
            scalar = -std::numeric_limits<t_type>::epsilon();
        else
            scalar = std::numeric_limits<t_type>::epsilon();
    }

    m_elem[0] /= scalar;
    m_elem[1] /= scalar;
    m_elem[2] /= scalar;
    m_elem[3] /= scalar;
    m_elem[4] /= scalar;
    m_elem[5] /= scalar;
    m_elem[6] /= scalar;
    m_elem[7] /= scalar;
    m_elem[8] /= scalar;

    return (*this);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline CommaInit<t_row * t_col, t_type> Matrix3<t_type, t_row, t_col>::operator<<(const t_type s)
{
    m_elem[0] = s;
    return CommaInit<t_row * t_col, t_type>(m_elem);
}

/* Arithmetic operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator-() const
{
    return Matrix3(*this) *= -1;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator+(const Matrix3 &m) const
{
    return Matrix3(*this) += m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator-(const Matrix3 &m) const
{
    return Matrix3(*this) -= m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator+(const Rotation<t_type, t_row, t_col> &m) const
{
    return Matrix3(*this) += m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator-(const Rotation<t_type, t_row, t_col> &m) const
{
    return Matrix3(*this) -= m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator+(const Matrix<t_row, t_col, t_type> &m) const
{
    return Matrix3(*this) += m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator-(const Matrix<t_row, t_col, t_type> &m) const
{
    return Matrix3(*this) -= m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator+(const Matrix<0, 0, t_type> &m) const
{
    return Matrix3(*this) += m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator-(const Matrix<0, 0, t_type> &m) const
{
    return Matrix3(*this) -= m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator+(const t_type s) const
{
    return Matrix3(*this) += s;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator-(const t_type s) const
{
    return Matrix3(*this) -= s;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator*(const t_type s) const
{
    return Matrix3(*this) *= s;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator/(const t_type s) const
{
    return Matrix3(*this) /= s;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
template <uint16_t col>
inline Matrix<t_row, col, t_type> Matrix3<t_type, t_row, t_col>::operator*(const Matrix<t_col, col, t_type> &m) const
{
    t_type mat[t_row * col];

    for (uint16_t irow = 0; irow < t_row; ++irow)
    {
        for (uint16_t icol = 0; icol < col; ++icol)
        {
            mat[irow * col + icol] = m_elem[irow * t_col] * m.m_elem[icol];
            mat[irow * col + icol] += m_elem[irow * t_col + 1] * m.m_elem[col + icol];
            mat[irow * col + icol] += m_elem[irow * t_col + 2] * m.m_elem[2 * col + icol];
        }
    }

    return Matrix<t_row, col, t_type>(mat);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix<0, 0, t_type> Matrix3<t_type, t_row, t_col>::operator*(const Matrix<0, 0, t_type> &m) const
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == 3 && "Dimensions do not matched");

    Matrix<0, 0, t_type> mat(t_row, m.m_col);

    for (uint16_t irow = 0; irow < t_row; ++irow)
    {
        for (uint16_t icol = 0; icol < m.m_col; ++icol)
        {
            mat[irow * m.m_col + icol] = m_elem[irow * t_col] * m.m_elem[icol];
            mat[irow * m.m_col + icol] += m_elem[irow * t_col + 1] * m.m_elem[m.m_col + icol];
            mat[irow * m.m_col + icol] += m_elem[irow * t_col + 2] * m.m_elem[2 * m.m_col + icol];
        }
    }

    return mat;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator*(const Matrix3 &m) const
{
    t_type mat[t_row * t_col];

    mat[0] = m_elem[0] * m.m_elem[0] + m_elem[1] * m.m_elem[3] + m_elem[2] * m.m_elem[6];
    mat[1] = m_elem[0] * m.m_elem[1] + m_elem[1] * m.m_elem[4] + m_elem[2] * m.m_elem[7];
    mat[2] = m_elem[0] * m.m_elem[2] + m_elem[1] * m.m_elem[5] + m_elem[2] * m.m_elem[8];

    mat[3] = m_elem[3] * m.m_elem[0] + m_elem[4] * m.m_elem[3] + m_elem[5] * m.m_elem[6];
    mat[4] = m_elem[3] * m.m_elem[1] + m_elem[4] * m.m_elem[4] + m_elem[5] * m.m_elem[7];
    mat[5] = m_elem[3] * m.m_elem[2] + m_elem[4] * m.m_elem[5] + m_elem[5] * m.m_elem[8];

    mat[6] = m_elem[6] * m.m_elem[0] + m_elem[7] * m.m_elem[3] + m_elem[8] * m.m_elem[6];
    mat[7] = m_elem[6] * m.m_elem[1] + m_elem[7] * m.m_elem[4] + m_elem[8] * m.m_elem[7];
    mat[8] = m_elem[6] * m.m_elem[2] + m_elem[7] * m.m_elem[5] + m_elem[8] * m.m_elem[8];

    return Matrix3<t_type, t_row, t_col>(mat);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator*(const Rotation<t_type, t_row, t_col> &m) const
{
    t_type mat[t_row * t_col];

    mat[0] = m_elem[0] * m.m_elem[0] + m_elem[1] * m.m_elem[3] + m_elem[2] * m.m_elem[6];
    mat[1] = m_elem[0] * m.m_elem[1] + m_elem[1] * m.m_elem[4] + m_elem[2] * m.m_elem[7];
    mat[2] = m_elem[0] * m.m_elem[2] + m_elem[1] * m.m_elem[5] + m_elem[2] * m.m_elem[8];

    mat[3] = m_elem[3] * m.m_elem[0] + m_elem[4] * m.m_elem[3] + m_elem[5] * m.m_elem[6];
    mat[4] = m_elem[3] * m.m_elem[1] + m_elem[4] * m.m_elem[4] + m_elem[5] * m.m_elem[7];
    mat[5] = m_elem[3] * m.m_elem[2] + m_elem[4] * m.m_elem[5] + m_elem[5] * m.m_elem[8];

    mat[6] = m_elem[6] * m.m_elem[0] + m_elem[7] * m.m_elem[3] + m_elem[8] * m.m_elem[6];
    mat[7] = m_elem[6] * m.m_elem[1] + m_elem[7] * m.m_elem[4] + m_elem[8] * m.m_elem[7];
    mat[8] = m_elem[6] * m.m_elem[2] + m_elem[7] * m.m_elem[5] + m_elem[8] * m.m_elem[8];

    return Matrix3<t_type, t_row, t_col>(mat);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector<t_row, t_type> Matrix3<t_type, t_row, t_col>::operator*(const Vector<t_col, t_type> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] * v.m_elem[0] + m_elem[1] * v.m_elem[1] + m_elem[2] * v.m_elem[2];
    vec[1] = m_elem[3] * v.m_elem[0] + m_elem[4] * v.m_elem[1] + m_elem[5] * v.m_elem[2];
    vec[2] = m_elem[6] * v.m_elem[0] + m_elem[7] * v.m_elem[1] + m_elem[8] * v.m_elem[2];

    return Vector<t_row, t_type>(vec);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector<t_row, t_type> Matrix3<t_type, t_row, t_col>::operator*(const Vector<0, t_type> &v) const
{
    assert(t_col == v.m_row && "Check dimensions");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    t_type vec[t_row];

    vec[0] = m_elem[0] * v.m_elem[0] + m_elem[1] * v.m_elem[1] + m_elem[2] * v.m_elem[2];
    vec[1] = m_elem[3] * v.m_elem[0] + m_elem[4] * v.m_elem[1] + m_elem[5] * v.m_elem[2];
    vec[2] = m_elem[6] * v.m_elem[0] + m_elem[7] * v.m_elem[1] + m_elem[8] * v.m_elem[2];

    return Vector<t_row, t_type>(vec);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector3<t_type, t_row> Matrix3<t_type, t_row, t_col>::operator*(const Vector3<t_type, t_col> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] * v.m_elem[0] + m_elem[1] * v.m_elem[1] + m_elem[2] * v.m_elem[2];
    vec[1] = m_elem[3] * v.m_elem[0] + m_elem[4] * v.m_elem[1] + m_elem[5] * v.m_elem[2];
    vec[2] = m_elem[6] * v.m_elem[0] + m_elem[7] * v.m_elem[1] + m_elem[8] * v.m_elem[2];

    return Vector3<t_type, t_row>(vec);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator&(const Vector<t_col, t_type> &v) const
{ // Mat3 * [v]x, []x is skew-symmetric matrix
    return Matrix3(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1], m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2], m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0],
        m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1], m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2], m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0],
        m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1], m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2], m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator&(const Vector<0, t_type> &v) const
{ // Mat3 * [v]x, []x is skew-symmetric matrix
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == v.m_row && "This method is only for 3 x 1 vector");

    return Matrix3(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1], m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2], m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0],
        m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1], m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2], m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0],
        m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1], m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2], m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Matrix3<t_type, t_row, t_col>::operator&(const Vector3<t_type, t_col> &v) const
{ // Mat3 * [v]x, []x is skew-symmetric matrix
    return Matrix3(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1], m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2], m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0],
        m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1], m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2], m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0],
        m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1], m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2], m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0]);
}

/* Comparison operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline bool Matrix3<t_type, t_row, t_col>::operator==(const Matrix3 &m) const
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline bool Matrix3<t_type, t_row, t_col>::operator!=(const Matrix3 &m) const
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline bool Matrix3<t_type, t_row, t_col>::operator==(const Rotation<t_type, t_row, t_col> &m) const
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline bool Matrix3<t_type, t_row, t_col>::operator!=(const Rotation<t_type, t_row, t_col> &m) const
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline bool Matrix3<t_type, t_row, t_col>::operator==(const Matrix<t_row, t_col, t_type> &m) const
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline bool Matrix3<t_type, t_row, t_col>::operator!=(const Matrix<t_row, t_col, t_type> &m) const
{
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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline bool Matrix3<t_type, t_row, t_col>::operator==(const Matrix<0, 0, t_type> &m) const
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == t_row && "Row dimensions do not matched");
    assert(m.m_col == t_col && "Col dimensions do not matched");

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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline bool Matrix3<t_type, t_row, t_col>::operator!=(const Matrix<0, 0, t_type> &m) const
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == t_row && "Row dimensions do not matched");
    assert(m.m_col == t_col && "Col dimensions do not matched");

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

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Matrix3<t_type, t_row, t_col>::Print(const char endChar)
{
#if defined(ARDUINO)
    for (uint16_t irow = 0; irow < t_row; irow++)
    {
        for (uint16_t icol = 0; icol < t_col; icol++)
        {
            Serial.printf("%7.3f ", (t_type)(m_elem[irow * t_col + icol]));
        }
        Serial.write('\n');
    }
    Serial.write(endChar);
#else
    for (uint16_t irow = 0; irow < t_row; irow++)
    {
        for (uint16_t icol = 0; icol < t_col; icol++)
        {
            printf("%7.3f ", (t_type)(m_elem[irow * t_col + icol]));
        }
        printf("\n");
    }
    printf("%c", endChar);
#endif
}

//-- Template Function ------------------------------------------------------//
// scalar * matrix
template <typename type, uint16_t row>
inline Matrix3<type, row> operator*(const type s, const Matrix3<type, row> &m)
{
    return Matrix3<type, row>(m) *= s;
}

typedef Matrix3<double> Mat3d;
typedef Matrix3<float> Mat3;

} // namespace Math
} // namespace dt

#endif // DTMATH_DTMATRIX3_TPP_
