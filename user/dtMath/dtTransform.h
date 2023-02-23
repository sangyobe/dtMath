/*!
\file       dtTransform.h
\brief      dtMath, Homogeneous transformation matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
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

template <uint16_t m_size, typename m_type> class CdtCommaInit;
template <uint16_t m_row, typename m_type> class CdtVector;
template <typename m_type, uint16_t m_row> class CdtVector3;
template <typename m_type, uint16_t m_row> class CdtVector6;
template <typename m_type, uint16_t m_row> class CdtQuaternion;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtMatrix;

template <typename m_type, uint16_t m_row, uint16_t m_col>
class CdtRotation;

template <typename m_type = float, uint16_t m_row = 4, uint16_t m_col = 4>
class CdtTransform
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    CdtRotation<m_type, 3, 3> m_R;
    CdtVector3<m_type, 3> m_p;
    m_type m_dummy[4];

public:
    CdtTransform();
    CdtTransform(const CdtRotation<m_type, 3, 3>& R, const CdtVector3<m_type, 3>& p);
    CdtTransform(const CdtQuaternion<m_type, 4>& q, const CdtVector3<m_type, 3>& p);
    CdtTransform(const uint16_t order, const CdtVector3<m_type, 3>& e, const CdtVector3<m_type, 3>& p);
    CdtTransform(const CdtTransform& m);
    ~CdtTransform() {}

    void SetZero();
    void SetIdentity();
    void SetElement(const CdtVector3<m_type, 3>& p);
    void SetElement(const CdtRotation<m_type, 3, 3>& R);
    void SetElement(const CdtQuaternion<m_type, 4>& q);
    void SetElement(const uint16_t order, const CdtVector3<m_type, 3>& e);
    void SetElement(const CdtRotation<m_type, 3, 3>& R, const CdtVector3<m_type, 3>& p);
    void SetElement(const CdtQuaternion<m_type, 4>& q, const CdtVector3<m_type, 3>& p);
    void SetElement(const uint16_t order, const CdtVector3<m_type, 3>& e, const CdtVector3<m_type, 3>& p);
    void SetElement(const CdtTransform& m);

    CdtQuaternion<m_type, 4> q() const { return CdtQuaternion<m_type, 4>(m_R); }
    CdtRotation<m_type, 3, 3> R() const { return m_R; }
    CdtVector3<m_type, 3> e(uint16_t order) const { return m_R.GetEulerAngles(order); }
    CdtVector3<m_type, 3> p() const { return m_p; }
    CdtVector6<m_type, 6> GetError(const CdtTransform& m) const;
    CdtMatrix<m_col, m_row, m_type> Transpose() const;
    CdtTransform Inv() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type& operator ()(uint16_t irow, uint16_t icol);
    // returns a row of non-modifiable elements
    const m_type& operator ()(uint16_t irow, uint16_t icol) const;

    /* Assignment operators */
    CdtTransform& operator  =(const CdtTransform& m);   // transform matrix  = transform matrix

    /* Arithmetic operators */
    CdtMatrix<m_row, m_col, m_type> operator +(const CdtMatrix<m_row, m_col, m_type>& m) const; // transform matrix + matrix
    CdtMatrix<m_row, m_col, m_type> operator -(const CdtMatrix<m_row, m_col, m_type>& m) const; // transform matrix - matrix

    template<uint16_t col>
    CdtMatrix<m_row, col, m_type> operator *(const CdtMatrix<m_row, col, m_type>& m) const; // transform matrix * matrix
    CdtTransform operator *(const CdtTransform& m) const;                                   // transform matrix * transform matrix
    CdtVector3<m_type, 3> operator *(const CdtVector<3, m_type>& v) const;                  // matrix * vector
    CdtVector3<m_type, 3> operator *(const CdtVector3<m_type, 3>& v) const;                 // matrix * vector3

    /* Comparison operators */
    bool operator == (const CdtTransform& m) const; // (true or false) matrix1 == matrix
    bool operator != (const CdtTransform& m) const; // (true or false) matrix1 == matrix

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row, uint16_t col> friend class CdtTransform;
    template <uint16_t row, uint16_t col, typename type> friend class CdtMatrix;
};

#include "dtTransform.tpp"

#endif // DTMATH_DTTRANSFORM_H_

