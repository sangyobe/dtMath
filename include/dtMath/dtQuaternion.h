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
template <typename t_type, uint16_t t_row> class Vector4;
template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Matrix3;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Rotation;

template <typename t_type = float, uint16_t t_row = 4>
class Quaternion
{
private:
    t_type m_tolerance = std::numeric_limits<t_type>::epsilon();
    t_type m_elem[t_row]; // [w, x, y, z]
    inline int SGN(const t_type x) const { return (x >= t_type(0)) ? 1 : -1; }
    inline void Euler2Quat(const uint16_t order, const t_type *e); // adh : private to public
    inline void RotMat2Quat(const t_type *rm);

public:
    Quaternion();
    Quaternion(const t_type *element);                                                               // input: quaternion elements, w->x->y->z
    Quaternion(const t_type w, const t_type x, const t_type y, const t_type z);                      // input: quaternion elements
    Quaternion(const uint16_t order, const t_type angle);                                            // input: angle and relative frame(body frame)
    Quaternion(const uint16_t order, const t_type angle1, const t_type angle2);                      // input: angle and relative frame(body frame)
    Quaternion(const uint16_t order, const t_type angle1, const t_type angle2, const t_type angle3); // input: angle and relative frame(body frame)
    Quaternion(const Quaternion &q);                                                                 // input: quaternion
    Quaternion(const Vector4<t_type, t_row> &q);                                                     // input: quaternion, w->x->y->z
    Quaternion(const uint16_t order, const Vector3<t_type, 3> &e);                                   // input: euler angles
    Quaternion(const uint16_t order, const Vector<3, t_type> &e);                                    // input: euler angles
    Quaternion(const Rotation<t_type, 3, 3> &rm);                                                    // input: rotation matrix
    Quaternion(const Vector3<t_type, 3> &v);                                                         // input: rotation vector(v in R3), v -> ExpMap(v) = quat(in S3)
    Quaternion(const Vector<3, t_type> &v);                                                          // input: rotation vector(v in R3), v -> ExpMap(v) = quat(in S3)
    ~Quaternion() {}

    void SetZero();
    void SetFill(const t_type value);
    void SetElement(const t_type *element);                                                               // input: quaternion elements, w->x->y->z
    void SetElement(const t_type w, const t_type x, const t_type y, const t_type z);                      // input: quaternion elements
    void SetElement(const uint16_t order, const t_type angle);                                            // input: angle and relative frame(body frame)
    void SetElement(const uint16_t order, const t_type angle1, const t_type angle2);                      // input: angle and relative frame(body frame)
    void SetElement(const uint16_t order, const t_type angle1, const t_type angle2, const t_type angle3); // input: angle and relative frame(body frame)
    void SetElement(const uint16_t order, const Vector3<t_type, 3> &e);                                   // input: euler angles
    void SetElement(const uint16_t order, const Vector<3, t_type> &e);                                    // input: euler angles
    void SetElement(const Rotation<t_type, 3, 3> &rm);                                                    // input: rotation matrix
    void SetElement(const Vector3<t_type, 3> &v);                                                         // input: rotation vector(v in R3), v -> ExpMap(v) = quat(in S3)
    void SetElement(const Vector<3, t_type> &v);                                                          // input: rotation vector(v in R3), v -> ExpMap(v) = quat(in S3)
    void SetSwap(const uint16_t i, const uint16_t j);
    void SetNormalize();

    const t_type *const GetElementsAddr() const;
    t_type GetNorm() const;
    t_type GetSqNorm() const;
    t_type GetSum() const;
    Quaternion GetNormalized() const;
    Quaternion GetConj() const;
    Vector3<t_type, 3> GetEulerAngles(const uint16_t order) const;
    Rotation<t_type, 3, 3> GetRotMat() const;
    Vector3<t_type, 3> GetOriErr(const Quaternion &q) const;
    Matrix<t_row, t_row, t_type> GetLmat() const;                             // Left-quaternion-product matrices, q1 * q2 = L(q1)q2
    Matrix<t_row, t_row, t_type> GetRmat() const;                             // Right-quaternion-product matrices, q1 * q2 = L(q2)q1
    Matrix<4, 3, t_type> GetGmat() const;                                     // the attitude Jacobian, utility for Jacobian of Vector-Valued Functions
    Matrix<3, 4, t_type> GetGTmat() const;                                    // transpose of G matrix
    static Matrix<t_row, t_row, t_type> GetLmat(const Quaternion<t_type> &q); // Left-quaternion-product matrices, q1 * q2 = L(q1)q2
    static Matrix<t_row, t_row, t_type> GetRmat(const Quaternion<t_type> &q); // Right-quaternion-product matrices, q1 * q2 = L(q2)q1
    static Matrix<4, 3, t_type> GetGmat(const Quaternion<t_type> &q);         // the attitude Jacobian, utility for Jacobian of Vector-Valued Functions
    // inline void Euler2Quat(const uint16_t order, const t_type *e);

    Matrix<1, t_row, t_type> Transpose() const;
    void Transpose(Matrix<1, t_row, t_type> &m) const;
    void Transpose(Matrix<0, 0, t_type> &m) const;

    Quaternion exp() const;                                // Exponential of general quaternions, e^(q)
    static Quaternion expMap(Quaternion<t_type> &V);       // Exponential map, Hp -> S3, exp(V) = e^(V) = q, V is pure quaternion(Hp)
    static Quaternion ExpMap(Vector3<t_type, 3> &v);       // Exponential map, R3 -> S3, Exp(v) = exp(v/2) e^(v/2) = q, v is rotateion vector(R3)
    Quaternion log() const;                                // Logarithm of general quaternions, log(q)
    Quaternion logMap() const;                             // Logarithmic map, S3 -> Hp, log(q) = ln(q) = u*th
    Vector3<t_type, 3> LogMap() const;                     // Logarithmic map, S3 -> R3, Log(q) = u*phi, S3 is the 3-dimensional surface of the unit sphere of R4 = unit quaternions
    Vector3<t_type, 3> InvCayleyMap() const;               // Inverse Cayley map, S3 -> R3, inv(Phi(q)) = qv / qw
    Quaternion ode(t_type wx, t_type wy, t_type wz) const; // dq/dt = 0.5*q*w, w is local angular velocity
    Quaternion ode(t_type *w) const;                       // dq/dt = 0.5*q*w, w is local angular velocity
    Quaternion ode(Vector3<t_type, 3> w) const;            // dq/dt = 0.5*q*w, w is local angular velocity
    Quaternion ode(Vector<3, t_type> w) const;             // dq/dt = 0.5*q*w, w is local angular velocity
    Quaternion Inv() const;

    /* Member access operators */
    t_type &operator()(uint16_t irow);             // returns a row of modifiable elements
    const t_type &operator()(uint16_t irow) const; // returns a row of non-modifiable elements

    /* Assignment operators */
    Quaternion &operator=(const Quaternion &q);             // quaternion1  = quaternion2
    Quaternion &operator+=(const Quaternion &q);            // quaternion1 += quaternion2
    Quaternion &operator-=(const Quaternion &q);            // quaternion1 -= quaternion2
    Quaternion &operator=(const Vector4<t_type, 4> &q);     // quaternion1  = quaternion2
    Quaternion &operator+=(const Vector4<t_type, 4> &q);    // quaternion1 += quaternion2
    Quaternion &operator-=(const Vector4<t_type, 4> &q);    // quaternion1 -= quaternion2
    Quaternion &operator=(const Vector<t_row, t_type> &q);  // quaternion1  = quaternion2
    Quaternion &operator+=(const Vector<t_row, t_type> &q); // quaternion1 += quaternion2
    Quaternion &operator-=(const Vector<t_row, t_type> &q); // quaternion1 -= quaternion2
    Quaternion &operator=(const Vector<0, t_type> &q);      // quaternion1  = quaternion2
    Quaternion &operator+=(const Vector<0, t_type> &q);     // quaternion1 += quaternion2
    Quaternion &operator-=(const Vector<0, t_type> &q);     // quaternion1 -= quaternion2
    Quaternion &operator*=(const t_type s);                 // quaternion1 *= scalar
    Quaternion &operator/=(const t_type s);                 // quaternion1 /= scalar
    CommaInit<t_row, t_type> operator<<(const t_type s);    // Init first matrix elements

    /* Arithmetic operators */
    Quaternion operator-() const;                    // -quaternion : minus sign
    Quaternion operator+(const Quaternion &q) const; // quaternion + quaternion
    Quaternion operator-(const Quaternion &q) const; // quaternion - quaternion
    Quaternion operator*(const t_type s) const;      // quaternion * scalar
    Quaternion operator/(const t_type s) const;      // quaternion / scalar
    Quaternion operator*(const Quaternion &q) const; // quaternion1 * quaternion2 : Quaternion Multiplication

    t_type Inner(const Quaternion &v) const;               // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2
    t_type Inner(const Vector4<t_type, t_row> &v) const;   // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2
    t_type Inner(const Vector<t_row, t_type> &v) const;    // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2
    t_type Inner(const Matrix<t_row, 1, t_type> &v) const; // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2
    t_type Inner(const Vector<0, t_type> &v) const;        // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2
    t_type Inner(const Matrix<0, 0, t_type> &v) const;     // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2

    /* Comparison operators */
    bool operator==(const Quaternion &q) const; // (true or false) quaternion1 == quaternion2
    bool operator!=(const Quaternion &q) const; // (true or false) quaternion1 != quaternion2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <uint16_t row, typename type> friend class Vector;
    template <typename type, uint16_t row> friend class Vector4;
    template <typename type, uint16_t row> friend class Quaternion;
    template <uint16_t row, uint16_t col, typename type> friend class Matrix;
    template <typename type, uint16_t row, uint16_t col> friend class Rotation;

    /* Friend template function */
    template <typename type, uint16_t row>
    friend Quaternion<type, row> operator*(const type s, const Quaternion<type, row> &v); // scalar * quaternion
};

} // namespace Math
} // namespace dt

#include "dtQuaternion.tpp"

#endif // DTMATH_DTQUATERNION_H_
