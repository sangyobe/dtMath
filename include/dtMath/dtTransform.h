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

#ifndef DTMATH_DTTRANSFORM_H_
#define DTMATH_DTTRANSFORM_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>

namespace dt
{
namespace Math
{

template <uint16_t m_size, typename t_type> class CommaInit;
template <uint16_t t_row, typename t_type> class Vector;
template <typename t_type, uint16_t t_row> class Vector3;
template <typename t_type, uint16_t t_row> class Vector6;
template <typename t_type, uint16_t t_row> class Quaternion;
template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Rotation;

template <typename t_type = float, uint16_t t_row = 4, uint16_t t_col = 4>
class Transform
{
private:
    t_type m_tolerance = std::numeric_limits<t_type>::epsilon();
    Rotation<t_type, 3, 3> m_R;
    Vector3<t_type, 3> m_p;
    t_type m_dummy[4];

public:
    Transform();
    Transform(const Rotation<t_type, 3, 3> &R, const Vector3<t_type, 3> &p);
    Transform(const Quaternion<t_type, 4> &q, const Vector3<t_type, 3> &p);
    Transform(const uint16_t order, const Vector3<t_type, 3> &e, const Vector3<t_type, 3> &p);
    Transform(const Transform &m);
    ~Transform() {}

    void SetZero();
    void SetIdentity();
    void SetElement(const Vector3<t_type, 3> &p);
    void SetElement(const Rotation<t_type, 3, 3> &R);
    void SetElement(const Quaternion<t_type, 4> &q);
    void SetElement(const uint16_t order, const Vector3<t_type, 3> &e);
    void SetElement(const Rotation<t_type, 3, 3> &R, const Vector3<t_type, 3> &p);
    void SetElement(const Quaternion<t_type, 4> &q, const Vector3<t_type, 3> &p);
    void SetElement(const uint16_t order, const Vector3<t_type, 3> &e, const Vector3<t_type, 3> &p);
    void SetElement(const Transform &m);

    Quaternion<t_type, 4> q() const { return Quaternion<t_type, 4>(m_R); }
    Rotation<t_type, 3, 3> R() const { return m_R; }
    Vector3<t_type, 3> e(uint16_t order) const { return m_R.GetEulerAngles(order); }
    Vector3<t_type, 3> p() const { return m_p; }
    Vector6<t_type, 6> GetError(const Transform &m) const;
    Matrix<t_col, t_row, t_type> Transpose() const;
    Transform Inv() const;

    /* Member access operators */
    t_type &operator()(uint16_t irow, uint16_t icol);             // returns a row of modifiable elements
    const t_type &operator()(uint16_t irow, uint16_t icol) const; // returns a row of non-modifiable elements

    /* Assignment operators */
    Transform &operator=(const Transform &m); // transform matrix  = transform matrix

    /* Arithmetic operators */
    Matrix<t_row, t_col, t_type> operator+(const Matrix<t_row, t_col, t_type> &m) const; // transform matrix + matrix
    Matrix<t_row, t_col, t_type> operator-(const Matrix<t_row, t_col, t_type> &m) const; // transform matrix - matrix

    template <uint16_t col>
    Matrix<t_row, col, t_type> operator*(const Matrix<t_row, col, t_type> &m) const; // transform matrix * matrix
    Transform operator*(const Transform &m) const;                                   // transform matrix * transform matrix
    Vector3<t_type, 3> operator*(const Vector<3, t_type> &v) const;                  // matrix * vector
    Vector3<t_type, 3> operator*(const Vector<0, t_type> &v) const;                  // matrix * vector
    Vector3<t_type, 3> operator*(const Vector3<t_type, 3> &v) const;                 // matrix * vector3

    /* Comparison operators */
    bool operator==(const Transform &m) const; // (true or false) matrix1 == matrix
    bool operator!=(const Transform &m) const; // (true or false) matrix1 == matrix

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row, uint16_t col> friend class Transform;
    template <uint16_t row, uint16_t col, typename type> friend class Matrix;
};

} // namespace Math
} // namespace dt

#include "dtTransform.tpp"

#endif // DTMATH_DTTRANSFORM_H_
