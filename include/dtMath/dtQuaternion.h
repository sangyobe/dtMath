/*!
\file       dtQuaternion.h
\brief      dtMath, Quaternion class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTQUATERNION_H_
#define DTMATH_DTQUATERNION_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

namespace dt
{
namespace Math
{

template <uint16_t m_size, typename t_type> class CommaInit;
template <uint16_t t_row, typename t_type> class Vector;
template <typename t_type, uint16_t t_row> class Vector3;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Matrix3;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Rotation;

template <typename t_type = float, uint16_t t_row = 4>
class Quaternion
{
private:
    t_type m_tolerance = std::numeric_limits<t_type>::epsilon();
    t_type m_elem[t_row]; // [w, x, y, z]
    inline int SGN(const t_type x) const { return (x >= t_type(0)) ? 1 : -1; }
    inline void Euler2Quat(const uint16_t order, const t_type *e);
    inline void RotMat2Quat(const t_type *rm);

public:
    Quaternion();
    Quaternion(const t_type *element);
    Quaternion(const t_type w, const t_type x, const t_type y, const t_type z);
    Quaternion(const uint16_t order, const t_type angle);
    Quaternion(const uint16_t order, const t_type angle1, const t_type angle2);
    Quaternion(const uint16_t order, const t_type angle1, const t_type angle2, const t_type angle3);
    Quaternion(const Quaternion &q);
    Quaternion(const uint16_t order, const Vector3<t_type, 3> &e);
    Quaternion(const uint16_t order, const Vector<3, t_type> &e);
    Quaternion(const Rotation<t_type, 3, 3> &rm);
    ~Quaternion() {}

    void SetZero();
    void SetFill(const t_type value);
    void SetElement(const t_type *element);
    void SetElement(const t_type w, const t_type x, const t_type y, const t_type z);
    void SetElement(const uint16_t order, const t_type angle);
    void SetElement(const uint16_t order, const t_type angle1, const t_type angle2);
    void SetElement(const uint16_t order, const t_type angle1, const t_type angle2, const t_type angle3);
    void SetElement(const Quaternion &q);
    void SetElement(const uint16_t order, const Vector3<t_type, 3> &e);
    void SetElement(const uint16_t order, const Vector<3, t_type> &e);
    void SetElement(const Rotation<t_type, 3, 3> &rm);
    void SetSwap(const uint16_t i, const uint16_t j);
    void SetNormalize();

    const t_type *const GetElementsAddr() const;
    t_type GetNorm() const;
    t_type GetSqNorm() const;
    t_type GetSum() const;
    Quaternion GetNormalized() const;
    Quaternion GetConj() const;
    Vector3<t_type, 3> GetEulerAngles(const uint16_t order) const;
    Vector3<t_type, 3> GetOriErr(const Quaternion &q) const;
    Quaternion exp() const;                                // exp(q) = e^(q), Expoential of general quaternions
    Quaternion log() const;                                // log(q) = ln(q), Logarithm of general quaternions
    Vector3<t_type, 3> Log() const;                        // Log(q) = u*phi : S3 -> R3, here S3 is the 3-dimensional surface of the unit sphere of R4 = unit quaternions
    Quaternion ode(t_type wx, t_type wy, t_type wz) const; // dq/dt
    Quaternion ode(t_type *w) const;                       // dq/dt
    Quaternion ode(Vector3<t_type, 3> w) const;            // dq/dt
    Quaternion ode(Vector<3, t_type> w) const;             // dq/dt
    Quaternion Inv() const;

    /* Member access operators */
    t_type &operator()(uint16_t irow);             // returns a row of modifiable elements
    const t_type &operator()(uint16_t irow) const; // returns a row of non-modifiable elements

    /* Assignment operators */
    Quaternion &operator=(const Quaternion &q);          // quaternion1  = quaternion2
    Quaternion &operator+=(const Quaternion &q);         // quaternion1 += quaternion2
    Quaternion &operator-=(const Quaternion &q);         // quaternion1 -= quaternion2
    Quaternion &operator*=(const t_type s);              // quaternion1 *= scalar
    Quaternion &operator/=(const t_type s);              // quaternion1 /= scalar
    CommaInit<t_row, t_type> operator<<(const t_type s); // Init first matrix elements

    /* Arithmetic operators */
    Quaternion operator-() const;                    // -quaternion : minus sign
    Quaternion operator+(const Quaternion &q) const; // quaternion + quaternion
    Quaternion operator-(const Quaternion &q) const; // quaternion - quaternion
    Quaternion operator*(const t_type s) const;      // quaternion * scalar
    Quaternion operator/(const t_type s) const;      // quaternion / scalar
    Quaternion operator*(const Quaternion &q) const; // quaternion1 * quaternion2 : Quaternion Multiplication

    /* Comparison operators */
    bool operator==(const Quaternion &q) const; // (true or false) quaternion1 == quaternion2
    bool operator!=(const Quaternion &q) const; // (true or false) quaternion1 != quaternion2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row> friend class Quaternion;
    template <typename type, uint16_t row, uint16_t col> friend class Rotation;

    /* Friend template function */
    template <typename type, uint16_t row>
    friend Quaternion<type, row> operator*(const type s, const Quaternion<type, row> &v); // scalar * quaternion
};

} // namespace Math
} // namespace dt

#include "dtQuaternion.tpp"

#endif // DTMATH_DTQUATERNION_H_
