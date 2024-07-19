/*!
\file       dtTransform.h
\brief      dtMath, Homogeneous transformation matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTTRANSFORM_TPP_
#define DTMATH_DTTRANSFORM_TPP_

#include "dtTransform.h"

namespace dt
{
namespace Math
{

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Transform<t_type, t_row, t_col>::Transform()
{
    m_dummy[0] = 0;
    m_dummy[1] = 0;
    m_dummy[2] = 0;
    m_dummy[3] = 1;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Transform<t_type, t_row, t_col>::Transform(const Rotation<t_type, 3, 3> &R, const Vector3<t_type, 3> &p)
{
    m_R = R;
    m_p = p;
    m_dummy[0] = 0;
    m_dummy[1] = 0;
    m_dummy[2] = 0;
    m_dummy[3] = 1;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Transform<t_type, t_row, t_col>::Transform(const Quaternion<t_type, 4> &q, const Vector3<t_type, 3> &p)
{
    m_R = Rotation<t_type, 3, 3>(q);
    m_p = p;
    m_dummy[0] = 0;
    m_dummy[1] = 0;
    m_dummy[2] = 0;
    m_dummy[3] = 1;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Transform<t_type, t_row, t_col>::Transform(const uint16_t order, const Vector3<t_type, 3> &e, const Vector3<t_type, 3> &p)
{
    m_R = Rotation<t_type, 3, 3>(order, e);
    m_p = p;
    m_dummy[0] = 0;
    m_dummy[1] = 0;
    m_dummy[2] = 0;
    m_dummy[3] = 1;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Transform<t_type, t_row, t_col>::Transform(const Transform &m)
{
    m_R = m.m_R;
    m_p = m.m_p;
    m_dummy[0] = 0;
    m_dummy[1] = 0;
    m_dummy[2] = 0;
    m_dummy[3] = 1;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Transform<t_type, t_row, t_col>::SetZero()
{
    m_R.SetZero();
    m_p.SetZero();
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Transform<t_type, t_row, t_col>::SetIdentity()
{
    m_R.SetIdentity();
    m_p.SetZero();
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Transform<t_type, t_row, t_col>::SetElement(const Vector3<t_type, 3> &p)
{
    m_p = p;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Transform<t_type, t_row, t_col>::SetElement(const Rotation<t_type, 3, 3> &R)
{
    m_R = R;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Transform<t_type, t_row, t_col>::SetElement(const Quaternion<t_type, 4> &q)
{
    m_R = Rotation<t_type, 3, 3>(q);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Transform<t_type, t_row, t_col>::SetElement(const uint16_t order, const Vector3<t_type, 3> &e)
{
    m_R = Rotation<t_type, 3, 3>(order, e);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Transform<t_type, t_row, t_col>::SetElement(const Rotation<t_type, 3, 3> &R, const Vector3<t_type, 3> &p)
{
    m_R = R;
    m_p = p;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Transform<t_type, t_row, t_col>::SetElement(const Quaternion<t_type, 4> &q, const Vector3<t_type, 3> &p)
{
    m_R = Rotation<t_type, 3, 3>(q);
    m_p = p;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Transform<t_type, t_row, t_col>::SetElement(const uint16_t order, const Vector3<t_type, 3> &e, const Vector3<t_type, 3> &p)
{
    m_R = Rotation<t_type, 3, 3>(order, e);
    m_p = p;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Transform<t_type, t_row, t_col>::SetElement(const Transform &m)
{
    m_R = m.m_R;
    m_p = m.m_p;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector6<t_type, 6> Transform<t_type, t_row, t_col>::GetError(const Transform &m) const
{
    // Vector3<t_type, 3> pos_error = m_p - m.m_p;
    // Quaternion<t_type, 4> qd(m_R);
    // Vector3<t_type, 3> ori_error = qd.GetOriErr(Quaternion<t_type, 4>(m.m_R));
    //
    // return Vector6<t_type, 6>(pos_error, ori_error);

    return Vector6<t_type, 6>(
        m_p - m.m_p,
        Quaternion<t_type, 4>(m_R).GetOriErr(Quaternion<t_type, 4>(m.m_R)));
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix<t_col, t_row, t_type> Transform<t_type, t_row, t_col>::Transpose() const
{
    t_type mat[t_col * t_row] = {
        m_R.m_elem[0], m_R.m_elem[3], m_R.m_elem[6], 0,
        m_R.m_elem[1], m_R.m_elem[4], m_R.m_elem[7], 0,
        m_R.m_elem[2], m_R.m_elem[5], m_R.m_elem[8], 0,
        m_p.m_elem[0], m_p.m_elem[1], m_p.m_elem[2], (t_type)(1)};

    return Matrix<t_col, t_row, t_type>(mat);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Transform<t_type, t_row, t_col> Transform<t_type, t_row, t_col>::Inv() const
{
    return Transform(m_R.Inv(), -(m_R.Inv() * m_p));
}

/* Member access operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline t_type &Transform<t_type, t_row, t_col>::operator()(uint16_t irow, uint16_t icol)
{
    assert(irow < t_row && "Index out of range");
    assert(icol < t_col && "Index out of range");

    if (irow < 3 && icol == 3)
        return m_p.m_elem[irow];
    else if (irow == 3)
        return m_dummy[icol];

    return m_R.m_elem[irow * 3 + icol];
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline const t_type &Transform<t_type, t_row, t_col>::operator()(uint16_t irow, uint16_t icol) const
{
    assert(irow < t_row && "Index out of range");
    assert(icol < t_col && "Index out of range");

    if (irow < 3 && icol == 3)
        return m_p.m_elem[irow];
    else if (irow == 3)
        return m_dummy[icol];

    return m_R.m_elem[irow * 3 + icol];
}

/* Assignment operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Transform<t_type, t_row, t_col> &Transform<t_type, t_row, t_col>::operator=(const Transform &m)
{
    m_R = m.m_R;
    m_p = m.m_p;

    return (*this);
}

/* Arithmetic operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix<t_row, t_col, t_type> Transform<t_type, t_row, t_col>::operator+(const Matrix<t_row, t_col, t_type> &m) const
{
    t_type mat[t_row * t_col] = {
        m_R.m_elem[0] + m.m_elem[0], m_R.m_elem[1] + m.m_elem[1], m_R.m_elem[2] + m.m_elem[2], m_p.m_elem[0] + m.m_elem[3],
        m_R.m_elem[3] + m.m_elem[4], m_R.m_elem[4] + m.m_elem[5], m_R.m_elem[5] + m.m_elem[6], m_p.m_elem[1] + m.m_elem[7],
        m_R.m_elem[6] + m.m_elem[8], m_R.m_elem[7] + m.m_elem[9], m_R.m_elem[8] + m.m_elem[10], m_p.m_elem[2] + m.m_elem[11],
        m.m_elem[12], m.m_elem[13], m.m_elem[14], (t_type)(1) + m.m_elem[15]};

    return Matrix<t_row, t_col, t_type>(mat);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Matrix<t_row, t_col, t_type> Transform<t_type, t_row, t_col>::operator-(const Matrix<t_row, t_col, t_type> &m) const
{
    t_type mat[t_row * t_col] = {
        m_R.m_elem[0] - m.m_elem[0], m_R.m_elem[1] - m.m_elem[1], m_R.m_elem[2] - m.m_elem[2], m_p.m_elem[0] - m.m_elem[3],
        m_R.m_elem[3] - m.m_elem[4], m_R.m_elem[4] - m.m_elem[5], m_R.m_elem[5] - m.m_elem[6], m_p.m_elem[1] - m.m_elem[7],
        m_R.m_elem[6] - m.m_elem[8], m_R.m_elem[7] - m.m_elem[9], m_R.m_elem[8] - m.m_elem[10], m_p.m_elem[2] - m.m_elem[11],
        -m.m_elem[12], -m.m_elem[13], -m.m_elem[14], (t_type)(1) - m.m_elem[15]};

    return Matrix<t_row, t_col, t_type>(mat);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
template <uint16_t col>
inline Matrix<t_row, col, t_type> Transform<t_type, t_row, t_col>::operator*(const Matrix<t_row, col, t_type> &m) const
{
    t_type mat[t_row * col];

    for (uint16_t irow = 0; irow < 3; ++irow)
    {
        for (uint16_t icol = 0; icol < col; ++icol)
        {
            mat[irow * col + icol] =
                m_R.m_elem[3 * irow] * m.m_elem[icol] +
                m_R.m_elem[3 * irow + 1] * m.m_elem[col + icol] +
                m_R.m_elem[3 * irow + 2] * m.m_elem[2 * col + icol] +
                m_p.m_elem[irow] * m.m_elem[3 * col + icol];
        }
    }

    memcpy(&mat[3 * col], &m.m_elem[3 * col], sizeof(t_type) * col);

    return Matrix<t_row, col, t_type>(mat);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Transform<t_type, t_row, t_col> Transform<t_type, t_row, t_col>::operator*(const Transform &m) const
{
    return Transform(m_R * m.m_R, m_R * m.m_p + m_p);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector3<t_type, 3> Transform<t_type, t_row, t_col>::operator*(const Vector<3, t_type> &v) const
{
    return Vector3<t_type, 3>(m_R * v + m_p);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector3<t_type, 3> Transform<t_type, t_row, t_col>::operator*(const Vector<0, t_type> &v) const
{
    return Vector3<t_type, 3>(m_R * v + m_p);
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline Vector3<t_type, 3> Transform<t_type, t_row, t_col>::operator*(const Vector3<t_type, 3> &v) const
{
    return Vector3<t_type, 3>(m_R * v + m_p);
}

/* Comparison operators */
template <typename t_type, uint16_t t_row, uint16_t t_col>
inline bool Transform<t_type, t_row, t_col>::operator==(const Transform &m) const
{
    if (m_R != m.m_R)
        return false;
    if (m_p != m.m_p)
        return false;

    return true;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline bool Transform<t_type, t_row, t_col>::operator!=(const Transform &m) const
{
    if (m_R != m.m_R)
        return true;
    if (m_p != m.m_p)
        return true;

    return false;
}

template <typename t_type, uint16_t t_row, uint16_t t_col>
inline void Transform<t_type, t_row, t_col>::Print(const char endChar)
{
#if defined(ARDUINO)
    for (uint16_t irow = 0; irow < 3; irow++)
    {
        Serial.printf("%7.3f %7.3f %7.3f %7.3f\n",
                      m_R.m_elem[irow * 3], m_R.m_elem[irow * 3 + 1], m_R.m_elem[irow * 3 + 2], m_p.m_elem[irow]);
    }
    Serial.printf("%7.3f %7.3f %7.3f %7.3f\n", m_dummy[0], m_dummy[1], m_dummy[2], m_dummy[3]);
    Serial.write(endChar);
#else
    for (uint16_t irow = 0; irow < 3; irow++)
    {
        printf("%7.3f %7.3f %7.3f %7.3f\n",
               m_R.m_elem[irow * 3], m_R.m_elem[irow * 3 + 1], m_R.m_elem[irow * 3 + 2], m_p.m_elem[irow]);
    }
    printf("%7.3f %7.3f %7.3f %7.3f\n", m_dummy[0], m_dummy[1], m_dummy[2], m_dummy[3]);
    printf("%c", endChar);
#endif
}

typedef Transform<> dtTMat;

} // namespace Math
} // namespace dt

#endif // DTMATH_DTTRANSFORM_TPP_
