/*!
\file       dtQuaternion.h
\brief      dtMath, Quaternion class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
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

template <uint16_t m_size, typename m_type> class CdtCommaInit;
template <uint16_t m_row, typename m_type> class CdtVector;
template <typename m_type, uint16_t m_row> class CdtVector3;
template <typename m_type, uint16_t m_row, uint16_t m_col> class CdtMatrix3;
template <typename m_type, uint16_t m_row, uint16_t m_col> class CdtRotation;

template <typename m_type = float, uint16_t m_row = 4>
class CdtQuaternion
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row]; // [w, x, y, z]
    inline int SGN(const m_type x) const { return (x >= m_type(0)) ? 1 : -1; }
    inline void Euler2Quat(const uint16_t order, const m_type* e);
    inline void RotMat2Quat(const m_type* rm);

public:
    CdtQuaternion();
    CdtQuaternion(const m_type* element);
    CdtQuaternion(const m_type w, const m_type x, const m_type y, const m_type z);
    CdtQuaternion(const uint16_t order, const m_type angle);
    CdtQuaternion(const uint16_t order, const m_type angle1, const m_type angle2);
    CdtQuaternion(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3);
    CdtQuaternion(const CdtQuaternion& q);
    CdtQuaternion(const uint16_t order, const CdtVector3<m_type, 3>& e);
    CdtQuaternion(const uint16_t order, const CdtVector<3, m_type>& e);
    CdtQuaternion(const CdtRotation<m_type, 3, 3>& rm);
    ~CdtQuaternion() {}

    void SetZero();
    void SetFill(const m_type value);
    void SetElement(const m_type* element);
    void SetElement(const m_type w, const m_type x, const m_type y, const m_type z);
    void SetElement(const uint16_t order, const m_type angle);
    void SetElement(const uint16_t order, const m_type angle1, const m_type angle2);
    void SetElement(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3);
    void SetElement(const CdtQuaternion& q);
    void SetElement(const uint16_t order, const CdtVector3<m_type, 3>& e);
    void SetElement(const uint16_t order, const CdtVector<3, m_type>& e);
    void SetElement(const CdtRotation<m_type, 3, 3>& rm);
    void SetSwap(const uint16_t i, const uint16_t j);
    void SetNormalize();

    const m_type* const GetElementsAddr() const;
    m_type GetNorm() const;
    m_type GetSqNorm() const;
    m_type GetSum() const;
    CdtQuaternion GetNormalized() const;
    CdtQuaternion GetConj() const;
    CdtVector3<m_type, 3> GetEulerAngles(const uint16_t order) const;
    CdtVector3<m_type, 3> GetOriErr(const CdtQuaternion& q) const;
    CdtQuaternion exp() const;                                  // exp(q) = e^(q), Expoential of general quaternions
    CdtQuaternion log() const;                                  // log(q) = ln(q), Logarithm of general quaternions
    CdtVector3<m_type, 3> Log() const;                          // Log(q) = u*phi : S3 -> R3, here S3 is the 3-dimensional surface of the unit sphere of R4 = unit quaternions
    CdtQuaternion ode(m_type wx, m_type wy, m_type wz) const;   // dq/dt
    CdtQuaternion ode(m_type *w) const;                         // dq/dt
    CdtQuaternion ode(CdtVector3<m_type, 3> w) const;           // dq/dt
    CdtQuaternion ode(CdtVector<3, m_type> w) const;            // dq/dt
    CdtQuaternion Inv() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type& operator ()(uint16_t irow) { return m_elem[irow]; }
    // returns a row of non-modifiable elements
    const m_type& operator ()(uint16_t irow) const { return m_elem[irow]; }

    /* Assignment operators */
    CdtQuaternion& operator  =(const CdtQuaternion& q);     // quaternion1  = quaternion2
    CdtQuaternion& operator +=(const CdtQuaternion& q);     // quaternion1 += quaternion2
    CdtQuaternion& operator -=(const CdtQuaternion& q);     // quaternion1 -= quaternion2
    CdtQuaternion& operator *=(const m_type s);             // quaternion1 *= scalar
    CdtQuaternion& operator /=(const m_type s);             // quaternion1 /= scalar
    CdtCommaInit<m_row, m_type> operator <<(const m_type s);// Init first matrix elements

    /* Arithmetic operators */
    CdtQuaternion operator -() const;                       // -quaternion : minus sign
    CdtQuaternion operator +(const CdtQuaternion& q) const; // quaternion + quaternion
    CdtQuaternion operator -(const CdtQuaternion& q) const; // quaternion - quaternion
    CdtQuaternion operator *(const m_type s) const;         // quaternion * scalar
    CdtQuaternion operator /(const m_type s) const;         // quaternion / scalar
    CdtQuaternion operator *(const CdtQuaternion& q) const; // quaternion1 * quaternion2 : Quaternion Multiplication

    /* Comparison operators */
    bool operator ==(const CdtQuaternion& q) const;         // (true or false) quaternion1 == quaternion2
    bool operator !=(const CdtQuaternion& q) const;         // (true or false) quaternion1 != quaternion2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row> friend class CdtQuaternion;
    template <typename type, uint16_t row, uint16_t col> friend class CdtRotation;

    /* Friend template function */
    template <typename type, uint16_t row>
    friend CdtQuaternion<type, row> operator *(const type s, const CdtQuaternion<type, row>& v);  // scalar * quaternion
};

#include "dtQuaternion.tpp"

#endif // DTMATH_DTQUATERNION_H_
