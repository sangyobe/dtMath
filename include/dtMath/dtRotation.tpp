/*!
\file       dtRotation.h
\brief      dtMath, Rotation matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTROTATION_TPP_
#define DTMATH_DTROTATION_TPP_

#include "dtRotation.h"

namespace dt
{
namespace Math
{

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col>::Rotation()
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
inline Rotation<t_type, t_row, t_col>::Rotation(const t_type *element)
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
inline Rotation<t_type, t_row, t_col>::Rotation(const t_type *element, const size_t n_byte)
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
inline Rotation<t_type, t_row, t_col>::Rotation(char c, const t_type *element, const size_t n_byte)
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
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col>::Rotation(
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
inline Rotation<t_type, t_row, t_col>::Rotation(const uint16_t order, const t_type angle)
{
    t_type c = std::cos(angle);
    t_type s = std::sin(angle);

    switch (order)
    {
    case 0x0: // X-axis
        m_elem[0] = 1;
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = 0;
        m_elem[4] = c;
        m_elem[5] = -s;
        m_elem[6] = 0;
        m_elem[7] = s;
        m_elem[8] = c;
        break;
    case 0x1: // Y-axis
        m_elem[0] = c;
        m_elem[1] = 0;
        m_elem[2] = s;
        m_elem[3] = 0;
        m_elem[4] = 1;
        m_elem[5] = 0;
        m_elem[6] = -s;
        m_elem[7] = 0;
        m_elem[8] = c;
        break;
    case 0x2: // Z-axis
        m_elem[0] = c;
        m_elem[1] = -s;
        m_elem[2] = 0;
        m_elem[3] = s;
        m_elem[4] = c;
        m_elem[5] = 0;
        m_elem[6] = 0;
        m_elem[7] = 0;
        m_elem[8] = 1;
        break;
    default:
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
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col>::Rotation(const uint16_t order, const t_type angle1, const t_type angle2)
{
    SetElement(order & 0xF, angle1); // R1
    Rotation<t_type> R2((order >> 4) & 0xF, angle2);
    (*this) = (*this) * R2;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col>::Rotation(const uint16_t order, const t_type angle1, const t_type angle2, const t_type angle3)
{
    SetElement(order & 0xF, angle1); // R1
    Rotation<t_type> R2((order >> 4) & 0xF, angle2);
    Rotation<t_type> R3((order >> 8) & 0xF, angle3);
    (*this) = (*this) * R2 * R3;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col>::Rotation(const Rotation &m)
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
inline Rotation<t_type, t_row, t_col>::Rotation(const Matrix3<t_type, t_row, t_col> &m)
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
inline Rotation<t_type, t_row, t_col>::Rotation(const Matrix<t_row, t_col, t_type> &m)
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
inline Rotation<t_type, t_row, t_col>::Rotation(const Matrix<0, 0, t_type> &m)
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
inline Rotation<t_type, t_row, t_col>::Rotation(const uint16_t order, const Vector3<t_type, 3> &e)
{
    Euler2RotMat(order, e.m_elem);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col>::Rotation(const uint16_t order, const Vector<3, t_type> &e)
{
    Euler2RotMat(order, e.m_elem);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col>::Rotation(const uint16_t order, const Vector<0, t_type> &e)
{
    assert(e.m_elem != nullptr && "Memory has not been allocated");
    assert(e.m_row == 3 && "Row dimensions do not matched");

    Euler2RotMat(order, e.m_elem);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col>::Rotation(const Quaternion<t_type, 4> &q)
{
    Quat2RotMat(q.m_elem);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetZero()
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
inline void Rotation<t_type, t_row, t_col>::SetIdentity()
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
inline void Rotation<t_type, t_row, t_col>::SetDiagonal(const t_type d1, const t_type d2, const t_type d3)
{
    m_elem[0] = d1;
    m_elem[4] = d2;
    m_elem[8] = d3;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetDiagonal(const t_type *element, const size_t n_byte)
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
inline void Rotation<t_type, t_row, t_col>::SetDiagonal(const Vector<t_row, t_type> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[4] = v.m_elem[1];
    m_elem[8] = v.m_elem[2];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetDiagonal(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == 3 && "Row dimensions do not matched");

    m_elem[0] = v.m_elem[0];
    m_elem[4] = v.m_elem[1];
    m_elem[8] = v.m_elem[2];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetDiagonal(const Vector3<t_type, t_row> &v)
{
    m_elem[0] = v.m_elem[0];
    m_elem[4] = v.m_elem[1];
    m_elem[8] = v.m_elem[2];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetFill(const t_type value)
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
inline void Rotation<t_type, t_row, t_col>::SetElement(const t_type *element, const size_t n_byte)
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
inline void Rotation<t_type, t_row, t_col>::SetElement(
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
inline void Rotation<t_type, t_row, t_col>::SetElement(const uint16_t order, const t_type angle)
{
    t_type c = std::cos(angle);
    t_type s = std::sin(angle);

    switch (order)
    {
    case 0x0: // X-axis
        m_elem[0] = 1;
        m_elem[1] = 0;
        m_elem[2] = 0;
        m_elem[3] = 0;
        m_elem[4] = c;
        m_elem[5] = -s;
        m_elem[6] = 0;
        m_elem[7] = s;
        m_elem[8] = c;
        break;
    case 0x1: // Y-axis
        m_elem[0] = c;
        m_elem[1] = 0;
        m_elem[2] = s;
        m_elem[3] = 0;
        m_elem[4] = 1;
        m_elem[5] = 0;
        m_elem[6] = -s;
        m_elem[7] = 0;
        m_elem[8] = c;
        break;
    case 0x2: // Z-axis
        m_elem[0] = c;
        m_elem[1] = -s;
        m_elem[2] = 0;
        m_elem[3] = s;
        m_elem[4] = c;
        m_elem[5] = 0;
        m_elem[6] = 0;
        m_elem[7] = 0;
        m_elem[8] = 1;
        break;
    }
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const uint16_t order, const t_type angle1, const t_type angle2)
{
    SetElement(order & 0xF, angle1); // R1
    Rotation<t_type> R2((order >> 4) & 0xF, angle2);
    (*this) = (*this) * R2;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const uint16_t order, const t_type angle1, const t_type angle2, const t_type angle3)
{
    SetElement(order & 0xF, angle1);                 // R1
    Rotation<t_type> R2((order >> 4) & 0xF, angle2); // R2
    Rotation<t_type> R3((order >> 8) & 0xF, angle3); // R3
    (*this) = (*this) * R2 * R3;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const Rotation &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(m.m_row == t_row && "Row dimensions do not matched");
    assert(m.m_col == t_col && "Col dimensions do not matched");

    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const uint16_t order, const Vector3<t_type, 3> &e)
{
    Euler2RotMat(order, e.m_elem);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const uint16_t order, const Vector<3, t_type> &e)
{
    Euler2RotMat(order, e.m_elem);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const uint16_t order, const Vector<0, t_type> &e)
{
    assert(e.m_elem != nullptr && "Memory has not been allocated");
    assert(e.m_row == 3 && "Row dimensions do not matched");

    Euler2RotMat(order, e.m_elem);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const uint16_t order, const t_type *e)
{
    Euler2RotMat(order, e);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const Quaternion<t_type, 4> &q)
{
    Quat2RotMat(q.m_elem);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const t_type *q)
{
    Quat2RotMat(q);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetElement(const t_type w, const t_type x, const t_type y, const t_type z)
{
    t_type q[4] = {w, x, y, z};
    Quat2RotMat(q);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2)
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
inline void Rotation<t_type, t_row, t_col>::SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2)
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
inline const t_type *const Rotation<t_type, t_row, t_col>::GetElementsAddr() const
{
    return m_elem;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector3<t_type, 3> Rotation<t_type, t_row, t_col>::GetRowVec(const uint16_t idxRow) const
{
    return Vector3<t_type, 3>(
        m_elem[idxRow * t_col],
        m_elem[idxRow * t_col + 1],
        m_elem[idxRow * t_col + 2]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector3<t_type, 3> Rotation<t_type, t_row, t_col>::GetColVec(const uint16_t idxCol) const
{
    return Vector3<t_type, 3>(
        m_elem[0 * t_col + idxCol],
        m_elem[1 * t_col + idxCol],
        m_elem[2 * t_col + idxCol]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline int8_t Rotation<t_type, t_row, t_col>::GetRowVec(const uint16_t idxRow, Vector3<t_type, 3> &v) const
{
    v.m_elem[0] = m_elem[idxRow * t_col];
    v.m_elem[1] = m_elem[idxRow * t_col + 1];
    v.m_elem[2] = m_elem[idxRow * t_col + 2];

    return 0;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline int8_t Rotation<t_type, t_row, t_col>::GetColVec(const uint16_t idxCol, Vector3<t_type, 3> &v) const
{
    v.m_elem[0] = m_elem[0 * t_col + idxCol];
    v.m_elem[1] = m_elem[1 * t_col + idxCol];
    v.m_elem[2] = m_elem[2 * t_col + idxCol];

    return 0;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector3<t_type, 3> Rotation<t_type, t_row, t_col>::GetEulerAngles(uint16_t order) const
{
    /* Tait?Bryan angles */
    t_type vec[3];

    // 0:x, 1:y, 2:z
    // order = 0x012 -> zyx, 0x210 -> xyz, 0x102 -> zxy, inverse order!
    uint16_t o1 = order & 0xF;
    uint16_t o2 = (order >> 4) & 0xF;
    uint16_t o3 = (order >> 8) & 0xF;
    int sign = ((o1 + 1) == o2) ? 1 : -1;

    vec[1] = std::asin(sign * m_elem[o1 * t_col + o3]);

    if ((1 - std::fabs(m_elem[o1 * t_col + o3])) <= std::numeric_limits<t_type>::epsilon())
    {
        // case s2 = +-1, c2 = 0 s3 = 0 c3 = +-1
        vec[0] = std::atan2(sign * m_elem[o3 * t_col + o2], m_elem[o2 * t_col + o2]);
        vec[2] = 0;
    }
    else
    {
        // case s2 != 1 or s2 != -1
        vec[0] = std::atan2(-sign * m_elem[o2 * t_col + o3], m_elem[o3 * t_col + o3]);
        vec[2] = std::atan2(-sign * m_elem[o1 * t_col + o2], m_elem[o1 * t_col + o1]);
    }

    return Vector3<t_type, 3>(vec);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::Transpose() const
{
    return Rotation(
        m_elem[0], m_elem[3], m_elem[6],
        m_elem[1], m_elem[4], m_elem[7],
        m_elem[2], m_elem[5], m_elem[8]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::log() const
{
    /* Logarithmic map of rotation matrix */
    // log:SO(3) -> so(3); R -> log(R) = [u*phi]x
    // phi * u is Axis-angle representation
    // phi = arccos((trace(R)-1)/2)
    // u = (R-RT)V / 2sin(phi)
    // ()V is the inverse of []x, matrix([]x) -> vector(()V)
    // if phi is 0, singular

    t_type phi = std::acos((m_elem[0] + m_elem[4] + m_elem[8] - 1) * (t_type)(0.5));
    t_type alpha;
    t_type uphi[3]{};

    if (std::abs(phi) > std::numeric_limits<t_type>::epsilon())
    {
        alpha = phi / (2 * std::sin(phi));
        uphi[0] = (m_elem[7] - m_elem[5]) * alpha;
        uphi[1] = (m_elem[2] - m_elem[6]) * alpha;
        uphi[2] = (m_elem[3] - m_elem[1]) * alpha;
    }

    return Matrix3<t_type, t_row, t_col>(
        0, -uphi[2], uphi[1],
        uphi[2], 0, -uphi[0],
        -uphi[1], uphi[0], 0);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::expMap(Matrix3<t_type, t_row, t_col> &skew)
{
    /* Exponential map of rotation matirx */
    // exp:so(3) -> SO(3); skew = [PHI]x -> exp([PHI]x) = e^([PHI]x)
    // PHI = norm(PHI) * (PHI / norm(PHI)) = phi*uv
    // R{PHI} = exp([PHI]x) = exp(phi*[uv]x)
    // = I + sin(norm)[uv]x + (1 − cos(norm))*([uv]x)^2
    // = cos(norm)*I + sin(norm)*[uv]x + (1 - cos(norm))*uv*uvT

    t_type norm = std::sqrt(skew.m_elem[7] * skew.m_elem[7] + skew.m_elem[2] * skew.m_elem[2] + skew.m_elem[3] * skew.m_elem[3]);
    t_type uv[3] = {skew.m_elem[7] / norm, skew.m_elem[2] / norm, skew.m_elem[3] / norm};
    t_type cv = std::cos(norm);
    t_type sv = std::sin(norm);
    t_type x = (sv * uv[0]) - ((cv - 1) * uv[1] * uv[2]);
    t_type y = (sv * uv[1]) - ((cv - 1) * uv[2] * uv[0]);
    t_type z = (sv * uv[2]) - ((cv - 1) * uv[0] * uv[1]);

    return Rotation<t_type, t_row, t_col>(
        cv - (cv - 1) * uv[0] * uv[0], -z, y,
        z, cv - (cv - 1) * uv[1] * uv[1], -x,
        -y, x, cv - (cv - 1) * uv[2] * uv[2]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::ExpMap(Vector3<t_type, 3> &v)
{
    /* Exponential map of rotation matirx */
    // Exp:R3 -> SO(3); v = PHI -> Exp([PHI]x) = e^([PHI]x)
    // v = PHI = norm(PHI) * (PHI / norm(PHI)) = phi*uv
    // R{PHI} = exp([PHI]x) = exp(phi*[uv]x)
    // = I + sin(norm)[uv]x + (1 − cos(norm))*([uv]x)^2
    // = cos(norm)*I + sin(norm)*[uv]x + (1 - cos(norm))*uv*uvT

    t_type norm = std::sqrt(v.m_elem[0] * v.m_elem[0] + v.m_elem[1] * v.m_elem[1] + v.m_elem[2] * v.m_elem[2]);
    t_type uv[3] = {v.m_elem[0] / norm, v.m_elem[1] / norm, v.m_elem[2] / norm};
    t_type cv = std::cos(norm);
    t_type sv = std::sin(norm);
    t_type x = (sv * uv[0]) - ((cv - 1) * uv[1] * uv[2]);
    t_type y = (sv * uv[1]) - ((cv - 1) * uv[2] * uv[0]);
    t_type z = (sv * uv[2]) - ((cv - 1) * uv[0] * uv[1]);

    return Rotation<t_type, t_row, t_col>(
        cv - (cv - 1) * uv[0] * uv[0], -z, y,
        z, cv - (cv - 1) * uv[1] * uv[1], -x,
        -y, x, cv - (cv - 1) * uv[2] * uv[2]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::logMap() const
{
    /* Logarithmic map of rotation matrix */
    // log:SO(3) -> so(3); R -> log(R) = [u*phi]x
    // phi * u is Axis-angle representation
    // phi = arccos((trace(R)-1)/2)
    // u = (R-RT)V / 2sin(phi)
    // ()V is the inverse of []x, matrix([]x) -> vector(()V)
    // if phi is 0, singular

    t_type phi = std::acos((m_elem[0] + m_elem[4] + m_elem[8] - 1) * (t_type)(0.5));
    t_type alpha;
    t_type uphi[3]{};

    if (std::abs(phi) > std::numeric_limits<t_type>::epsilon())
    {
        alpha = phi / (2 * std::sin(phi));
        uphi[0] = (m_elem[7] - m_elem[5]) * alpha;
        uphi[1] = (m_elem[2] - m_elem[6]) * alpha;
        uphi[2] = (m_elem[3] - m_elem[1]) * alpha;
    }

    return Matrix3<t_type, t_row, t_col>(
        0, -uphi[2], uphi[1],
        uphi[2], 0, -uphi[0],
        -uphi[1], uphi[0], 0);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector3<t_type, 3> Rotation<t_type, t_row, t_col>::LogMap() const
{
    /* Logarithmic map of rotation matrix */
    // Log:SO(3) -> R3; R -> Log(R) = u*phi
    // phi * u is Axis-angle representation
    // phi = arccos((trace(R)-1)/2)
    // u = (R-RT)V / 2sin(phi)
    // ()V is the inverse of []x, matrix([]x) -> vector(()V)
    // if phi is 0, singular

    t_type phi = std::acos((m_elem[0] + m_elem[4] + m_elem[8] - 1) * (t_type)(0.5));
    t_type alpha;

    if (std::abs(phi) > std::numeric_limits<t_type>::epsilon())
        alpha = phi / (2 * std::sin(phi));
    else
        alpha = 0;

    return Vector3<t_type, 3>(
        (m_elem[7] - m_elem[5]) * alpha,
        (m_elem[2] - m_elem[6]) * alpha,
        (m_elem[3] - m_elem[1]) * alpha);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::ode(t_type wx, t_type wy, t_type wz) const
{
    /* Ordinary Differential Equation (ODE) */
    // d(R^{i-1}_{i})/dt = [w^{i-1}_{i}]x * R^{i-1}_{i}
    //                   = R^{i-1}_{i} * [w^{i}_{i-1}]x
    //
    // R is rotation matrix wrt frame i-1, R^{i-1}_{i}
    // omega(w) is angular velocity wrt frame i
    // dR/dt = R[w]x
    // where [w]x = [0 -wz wy]
    //              [wz 0 -wx]
    //              [-wy wx 0]

    return Rotation(
        m_elem[1] * wz - m_elem[2] * wy, -m_elem[0] * wz + m_elem[2] * wx, m_elem[0] * wy - m_elem[1] * wx,
        m_elem[4] * wz - m_elem[5] * wy, -m_elem[3] * wz + m_elem[5] * wx, m_elem[3] * wy - m_elem[4] * wx,
        m_elem[7] * wz - m_elem[8] * wy, -m_elem[6] * wz + m_elem[8] * wx, m_elem[6] * wy - m_elem[7] * wx);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::ode(t_type *w) const
{
    /* Ordinary Differential Equation (ODE) */
    // d(R^{i-1}_{i})/dt = [w^{i-1}_{i}]x * R^{i-1}_{i}
    //                   = R^{i-1}_{i} * [w^{i}_{i-1}]x
    //
    // R is rotation matrix wrt frame i-1, R^{i-1}_{i}
    // omega(w) is angular velocity wrt frame i
    // dR/dt = R[w]x
    // where [w]x = [0 -wz wy]
    //              [wz 0 -wx]
    //              [-wy wx 0]

    return Rotation(
        m_elem[1] * w[2] - m_elem[2] * w[1], -m_elem[0] * w[2] + m_elem[2] * w[0], m_elem[0] * w[1] - m_elem[1] * w[0],
        m_elem[4] * w[2] - m_elem[5] * w[1], -m_elem[3] * w[2] + m_elem[5] * w[0], m_elem[3] * w[1] - m_elem[4] * w[0],
        m_elem[7] * w[2] - m_elem[8] * w[1], -m_elem[6] * w[2] + m_elem[8] * w[0], m_elem[6] * w[1] - m_elem[7] * w[0]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::ode(Vector3<t_type, 3> w) const
{
    /* Ordinary Differential Equation (ODE) */
    // d(R^{i-1}_{i})/dt = [w^{i-1}_{i}]x * R^{i-1}_{i}
    //                   = R^{i-1}_{i} * [w^{i}_{i-1}]x
    //
    // R is rotation matrix wrt frame i-1, R^{i-1}_{i}
    // omega(w) is angular velocity wrt frame i
    // dR/dt = R[w]x
    // where [w]x = [0 -wz wy]
    //              [wz 0 -wx]
    //              [-wy wx 0]

    return Rotation(
        m_elem[1] * w.m_elem[2] - m_elem[2] * w.m_elem[1], -m_elem[0] * w.m_elem[2] + m_elem[2] * w.m_elem[0], m_elem[0] * w.m_elem[1] - m_elem[1] * w.m_elem[0],
        m_elem[4] * w.m_elem[2] - m_elem[5] * w.m_elem[1], -m_elem[3] * w.m_elem[2] + m_elem[5] * w.m_elem[0], m_elem[3] * w.m_elem[1] - m_elem[4] * w.m_elem[0],
        m_elem[7] * w.m_elem[2] - m_elem[8] * w.m_elem[1], -m_elem[6] * w.m_elem[2] + m_elem[8] * w.m_elem[0], m_elem[6] * w.m_elem[1] - m_elem[7] * w.m_elem[0]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::ode(Vector<3, t_type> w) const
{
    /* Ordinary Differential Equation (ODE) */
    // d(R^{i-1}_{i})/dt = [w^{i-1}_{i}]x * R^{i-1}_{i}
    //                   = R^{i-1}_{i} * [w^{i}_{i-1}]x
    //
    // R is rotation matrix wrt frame i-1, R^{i-1}_{i}
    // omega(w) is angular velocity wrt frame i
    // dR/dt = R[w]x
    // where [w]x = [0 -wz wy]
    //              [wz 0 -wx]
    //              [-wy wx 0]

    return Rotation(
        m_elem[1] * w.m_elem[2] - m_elem[2] * w.m_elem[1], -m_elem[0] * w.m_elem[2] + m_elem[2] * w.m_elem[0], m_elem[0] * w.m_elem[1] - m_elem[1] * w.m_elem[0],
        m_elem[4] * w.m_elem[2] - m_elem[5] * w.m_elem[1], -m_elem[3] * w.m_elem[2] + m_elem[5] * w.m_elem[0], m_elem[3] * w.m_elem[1] - m_elem[4] * w.m_elem[0],
        m_elem[7] * w.m_elem[2] - m_elem[8] * w.m_elem[1], -m_elem[6] * w.m_elem[2] + m_elem[8] * w.m_elem[0], m_elem[6] * w.m_elem[1] - m_elem[7] * w.m_elem[0]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::ode(Vector<0, t_type> w) const
{
    assert(w.m_elem != nullptr && "Memory has not been allocated");
    assert(w.m_row == 3 && "Row dimensions do not matched");

    /* Ordinary Differential Equation (ODE) */
    // d(R^{i-1}_{i})/dt = [w^{i-1}_{i}]x * R^{i-1}_{i}
    //                   = R^{i-1}_{i} * [w^{i}_{i-1}]x
    //
    // R is rotation matrix wrt frame i-1, R^{i-1}_{i}
    // omega(w) is angular velocity wrt frame i
    // dR/dt = R[w]x
    // where [w]x = [0 -wz wy]
    //              [wz 0 -wx]
    //              [-wy wx 0]

    return Rotation(
        m_elem[1] * w.m_elem[2] - m_elem[2] * w.m_elem[1], -m_elem[0] * w.m_elem[2] + m_elem[2] * w.m_elem[0], m_elem[0] * w.m_elem[1] - m_elem[1] * w.m_elem[0],
        m_elem[4] * w.m_elem[2] - m_elem[5] * w.m_elem[1], -m_elem[3] * w.m_elem[2] + m_elem[5] * w.m_elem[0], m_elem[3] * w.m_elem[1] - m_elem[4] * w.m_elem[0],
        m_elem[7] * w.m_elem[2] - m_elem[8] * w.m_elem[1], -m_elem[6] * w.m_elem[2] + m_elem[8] * w.m_elem[0], m_elem[6] * w.m_elem[1] - m_elem[7] * w.m_elem[0]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::Inv() const
{
    return Rotation(
        m_elem[0], m_elem[3], m_elem[6],
        m_elem[1], m_elem[4], m_elem[7],
        m_elem[2], m_elem[5], m_elem[8]);
}

/* Member access operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline t_type &Rotation<t_type, t_row, t_col>::operator()(uint16_t irow, uint16_t icol)
{
    assert(irow < t_row && "Index out of range");
    assert(icol < t_col && "Index out of range");

    return m_elem[irow * t_col + icol];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline const t_type &Rotation<t_type, t_row, t_col>::operator()(uint16_t irow, uint16_t icol) const
{
    assert(irow < t_row && "Index out of range");
    assert(icol < t_col && "Index out of range");

    return m_elem[irow * t_col + icol];
}

/* Assignment operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> &Rotation<t_type, t_row, t_col>::operator=(const Rotation &m)
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
inline CommaInit<t_row * t_col, t_type> Rotation<t_type, t_row, t_col>::operator<<(const t_type s)
{
    m_elem[0] = s;
    return CommaInit<t_row * t_col, t_type>(m_elem);
}

/* Arithmetic operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator-() const
{
    return Rotation(
        -m_elem[0], -m_elem[1], -m_elem[2],
        -m_elem[3], -m_elem[4], -m_elem[5],
        -m_elem[6], -m_elem[7], -m_elem[8]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator+(const Rotation &m) const
{
    return Matrix3<t_type, t_row, t_col>(*this) += m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator-(const Rotation &m) const
{
    return Matrix3<t_type, t_row, t_col>(*this) -= m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator+(const Matrix3<t_type, t_row, t_col> &m) const
{
    return Matrix3<t_type, t_row, t_col>(*this) += m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator-(const Matrix3<t_type, t_row, t_col> &m) const
{
    return Matrix3<t_type, t_row, t_col>(*this) -= m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator+(const Matrix<t_row, t_col, t_type> &m) const
{
    return Matrix3<t_type, t_row, t_col>(*this) += m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator-(const Matrix<t_row, t_col, t_type> &m) const
{
    return Matrix3<t_type, t_row, t_col>(*this) -= m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator+(const Matrix<0, 0, t_type> &m) const
{
    return Matrix3<t_type, t_row, t_col>(*this) += m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator-(const Matrix<0, 0, t_type> &m) const
{
    return Matrix3<t_type, t_row, t_col>(*this) -= m;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator*(const t_type s) const
{
    return Matrix3<t_type, t_row, t_col>(*this) *= s;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator/(const t_type s) const
{
    return Matrix3<t_type, t_row, t_col>(*this) /= s;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
template <uint16_t col>
inline Matrix<t_row, col, t_type> Rotation<t_type, t_row, t_col>::operator*(const Matrix<t_col, col, t_type> &m) const
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
inline Matrix<0, 0, t_type> Rotation<t_type, t_row, t_col>::operator*(const Matrix<0, 0, t_type> &m) const
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
inline Matrix3<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator*(const Matrix3<t_type, t_row, t_col> &m) const
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
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator*(const Rotation &m) const
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

    return Rotation(mat);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector<t_row, t_type> Rotation<t_type, t_row, t_col>::operator*(const Vector<t_col, t_type> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] * v.m_elem[0] + m_elem[1] * v.m_elem[1] + m_elem[2] * v.m_elem[2];
    vec[1] = m_elem[3] * v.m_elem[0] + m_elem[4] * v.m_elem[1] + m_elem[5] * v.m_elem[2];
    vec[2] = m_elem[6] * v.m_elem[0] + m_elem[7] * v.m_elem[1] + m_elem[8] * v.m_elem[2];

    return Vector<t_row, t_type>(vec);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector<t_row, t_type> Rotation<t_type, t_row, t_col>::operator*(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == 3 && "Dimensions do not matched");

    t_type vec[t_row];

    vec[0] = m_elem[0] * v.m_elem[0] + m_elem[1] * v.m_elem[1] + m_elem[2] * v.m_elem[2];
    vec[1] = m_elem[3] * v.m_elem[0] + m_elem[4] * v.m_elem[1] + m_elem[5] * v.m_elem[2];
    vec[2] = m_elem[6] * v.m_elem[0] + m_elem[7] * v.m_elem[1] + m_elem[8] * v.m_elem[2];

    return Vector<t_row, t_type>(vec);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector3<t_type, t_row> Rotation<t_type, t_row, t_col>::operator*(const Vector3<t_type, t_col> &v) const
{
    t_type vec[t_row];

    vec[0] = m_elem[0] * v.m_elem[0] + m_elem[1] * v.m_elem[1] + m_elem[2] * v.m_elem[2];
    vec[1] = m_elem[3] * v.m_elem[0] + m_elem[4] * v.m_elem[1] + m_elem[5] * v.m_elem[2];
    vec[2] = m_elem[6] * v.m_elem[0] + m_elem[7] * v.m_elem[1] + m_elem[8] * v.m_elem[2];

    return Vector3<t_type, t_row>(vec);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator&(const Vector<t_col, t_type> &v) const
{ // RotMat * [v]x, []x is skew-symmetric matrix
    return Rotation(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1], m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2], m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0],
        m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1], m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2], m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0],
        m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1], m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2], m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator&(const Vector<0, t_type> &v) const
{ // RotMat * [v]x, []x is skew-symmetric matrix
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(v.m_row == 3 && "Dimensions do not matched");

    return Rotation(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1], m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2], m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0],
        m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1], m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2], m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0],
        m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1], m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2], m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0]);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Rotation<t_type, t_row, t_col> Rotation<t_type, t_row, t_col>::operator&(const Vector3<t_type, t_col> &v) const
{ // RotMat * [v]x, []x is skew-symmetric matrix
    return Rotation(
        m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1], m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2], m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0],
        m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1], m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2], m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0],
        m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1], m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2], m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0]);
}

/* Comparison operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline bool Rotation<t_type, t_row, t_col>::operator==(const Rotation &m) const
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
inline bool Rotation<t_type, t_row, t_col>::operator!=(const Rotation &m) const
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
inline bool Rotation<t_type, t_row, t_col>::operator==(const Matrix3<t_type, t_row, t_col> &m) const
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
inline bool Rotation<t_type, t_row, t_col>::operator!=(const Matrix3<t_type, t_row, t_col> &m) const
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
inline bool Rotation<t_type, t_row, t_col>::operator==(const Matrix<t_row, t_col, t_type> &m) const
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
inline bool Rotation<t_type, t_row, t_col>::operator!=(const Matrix<t_row, t_col, t_type> &m) const
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
inline bool Rotation<t_type, t_row, t_col>::operator==(const Matrix<0, 0, t_type> &m) const
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
inline bool Rotation<t_type, t_row, t_col>::operator!=(const Matrix<0, 0, t_type> &m) const
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
inline void Rotation<t_type, t_row, t_col>::Print(const char endChar)
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

//-- Private Member Function ------------------------------------------------//
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::Euler2RotMat(const uint16_t order, const t_type *e)
{
    t_type s1 = std::sin(e[0]);
    t_type c1 = std::cos(e[0]);
    t_type s2 = std::sin(e[1]);
    t_type c2 = std::cos(e[1]);
    t_type s3 = std::sin(e[2]);
    t_type c3 = std::cos(e[2]);

    /* Only Tait?Bryan angles */
    switch (order)
    {
    case 0x120: // xzy
        m_elem[0] = c2 * c3;
        m_elem[1] = -s2;
        m_elem[2] = c2 * s3;
        m_elem[3] = s1 * s3 + c1 * c3 * s2;
        m_elem[4] = c1 * c2;
        m_elem[5] = c1 * s2 * s3 - c3 * s1;
        m_elem[6] = c3 * s1 * s2 - c1 * s3;
        m_elem[7] = c2 * s1;
        m_elem[8] = c1 * c3 + s1 * s2 * s3;
        break;
    case 0x210: // xyz
        m_elem[0] = c2 * c3;
        m_elem[1] = -c2 * s3;
        m_elem[2] = s2;
        m_elem[3] = c1 * s3 + c3 * s1 * s2;
        m_elem[4] = c1 * c3 - s1 * s2 * s3;
        m_elem[5] = -c2 * s1;
        m_elem[6] = s1 * s3 - c1 * c3 * s2;
        m_elem[7] = c3 * s1 + c1 * s2 * s3;
        m_elem[8] = c1 * c2;
        break;
    case 0x201: // yxz
        m_elem[0] = c1 * c3 + s1 * s2 * s3;
        m_elem[1] = c3 * s1 * s2 - c1 * s3;
        m_elem[2] = c2 * s1;
        m_elem[3] = c2 * s3;
        m_elem[4] = c2 * c3;
        m_elem[5] = -s2;
        m_elem[6] = c1 * s2 * s3 - c3 * s1;
        m_elem[7] = c1 * c3 * s2 + s1 * s3;
        m_elem[8] = c1 * c2;
        break;
    case 0x021: // yzx
        m_elem[0] = c1 * c2;
        m_elem[1] = s1 * s3 - c1 * c3 * s2;
        m_elem[2] = c3 * s1 + c1 * s2 * s3;
        m_elem[3] = s2;
        m_elem[4] = c2 * c3;
        m_elem[5] = -c2 * s3;
        m_elem[6] = -c2 * s1;
        m_elem[7] = c1 * s3 + c3 * s1 * s2;
        m_elem[8] = c1 * c3 - s1 * s2 * s3;
        break;
    case 0x012: // zyx
        m_elem[0] = c1 * c2;
        m_elem[1] = c1 * s2 * s3 - c3 * s1;
        m_elem[2] = s1 * s3 + c1 * c3 * s2;
        m_elem[3] = c2 * s1;
        m_elem[4] = c1 * c3 + s1 * s2 * s3;
        m_elem[5] = c3 * s1 * s2 - c1 * s3;
        m_elem[6] = -s2;
        m_elem[7] = c2 * s3;
        m_elem[8] = c2 * c3;
        break;
    case 0x102: // zxy
        m_elem[0] = c1 * c3 - s1 * s2 * s3;
        m_elem[1] = -c2 * s1;
        m_elem[2] = c1 * s3 + c3 * s1 * s2;
        m_elem[3] = c3 * s1 + c1 * s2 * s3;
        m_elem[4] = c1 * c2;
        m_elem[5] = s1 * s3 - c1 * c3 * s2;
        m_elem[6] = -c2 * s3;
        m_elem[7] = s2;
        m_elem[8] = c2 * c3;
        break;
    }
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Rotation<t_type, t_row, t_col>::Quat2RotMat(const t_type *q)
{
    // m_elem[0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
    m_elem[0] = 1 - 2 * (q[2] * q[2] + q[3] * q[3]);
    m_elem[1] = 2 * (q[1] * q[2] - q[0] * q[3]);
    m_elem[2] = 2 * (q[1] * q[3] + q[0] * q[2]);

    m_elem[3] = 2 * (q[1] * q[2] + q[0] * q[3]);
    // m_elem[4] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
    m_elem[4] = 1 - 2 * (q[1] * q[1] + q[3] * q[3]);
    m_elem[5] = 2 * (q[2] * q[3] - q[0] * q[1]);

    m_elem[6] = 2 * (q[1] * q[3] - q[0] * q[2]);
    m_elem[7] = 2 * (q[2] * q[3] + q[0] * q[1]);
    // m_elem[8] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
    m_elem[8] = 1 - 2 * (q[1] * q[1] + q[2] * q[2]);
}

//-- Template Function ------------------------------------------------------//
// scalar * matrix
template <typename type, uint16_t row, uint16_t col>
inline Matrix3<type, row, col> operator*(const type s, const Rotation<type, row, col> &m)
{
    return Matrix3<type, row, col>(m) *= s;
}

typedef Rotation<> dtRotMat;

} // namespace Math
} // namespace dt

#endif // DTMATH_DTROTATION_TPP_