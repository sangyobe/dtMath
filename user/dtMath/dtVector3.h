/*!
\file       dtVector3.h
\brief      dtMath, 3x1 Vector class, lighter and faster than general vector class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTVECTOR3_H_
#define DTMATH_DTVECTOR3_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

template <uint16_t m_size, typename m_type> class CdtCommaInit;
template <uint16_t m_row, typename m_type> class CdtVector;
template <typename m_type, uint16_t m_row> class CdtVector4;
template <typename m_type, uint16_t m_row> class CdtVector6;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtMatrix;
template <typename m_type, uint16_t m_row, uint16_t m_col> class CdtMatrix3;
template <typename m_type, uint16_t m_row, uint16_t m_col> class CdtRotation;

template <typename m_type = float, uint16_t m_row = 3>
class CdtVector3
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row];
    CdtVector3(const m_type* element);
    inline void CrossProduct(const m_type* v);

public:
    CdtVector3();
    CdtVector3(const m_type* element, const size_t n_byte);
    CdtVector3(const m_type i, const m_type j, const m_type k);
    CdtVector3(const CdtVector3& v);
    CdtVector3(const CdtVector<m_row, m_type>& v);
    CdtVector3(const CdtMatrix<m_row, 1, m_type>& v);
    ~CdtVector3() {}

    void SetZero();
    void SetFill(const m_type value);
    void SetElement(const m_type* element, const size_t n_byte);
    void SetElement(const m_type i, const m_type j, const m_type k);
    void SetElement(const CdtVector3& v);
    void SetElement(const CdtVector<m_row, m_type>& v);
    void SetElement(const CdtMatrix<m_row, 1, m_type>& v);
    template <uint16_t row>
    void SetBlock(const uint16_t idxRow, const CdtVector<row, m_type>& v);
    void SetBlock(const uint16_t idxRow, const m_type *v, const size_t n_byte);
    void SetBlock(const uint16_t idxRow, const CdtVector3<m_type, 3>& v);
    void SetBlock(const uint16_t idxRow, const CdtVector4<m_type, 4>& v);
    void SetBlock(const uint16_t idxRow, const CdtVector6<m_type, 6>& v);
    template <uint16_t row>
    void SetBlock(const uint16_t idxRow, const CdtMatrix<row, 1, m_type>& v);
    void SetSwap(const uint16_t i, const uint16_t j);
    void SetNormalize();

    const m_type* const GetElementsAddr() const;
    template <uint16_t row>
    CdtVector<row, m_type> GetBlock(const uint16_t idx);
    template <uint16_t row>
    int8_t GetBlock(const uint16_t idx, CdtVector<row, m_type>& v);
    m_type GetNorm() const;
    m_type GetSqNorm() const;
    m_type GetSum() const;
    CdtVector3 GetNormalized() const;
    CdtMatrix3<m_type, m_row, m_row> GetSkew() const;
    CdtMatrix<1, m_row, m_type> Transpose() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type& operator ()(uint16_t irow) { return m_elem[irow]; }
    // returns a row of non-modifiable elements
    const m_type& operator ()(uint16_t irow) const { return m_elem[irow]; }

    /* Assignment operators */
    CdtVector3& operator  =(const CdtVector3& v);                   // vector1  = vector2
    CdtVector3& operator +=(const CdtVector3& v);                   // vector1 += vector2
    CdtVector3& operator -=(const CdtVector3& v);                   // vector1 -= vector2
    CdtVector3& operator *=(const CdtVector3& v);                   // vector1 *= vector2, vector1(i) *= vector2(i)
    CdtVector3& operator /=(const CdtVector3& v);                   // vector1 /= vector2, vector1(i) /= vector2(i)
    CdtVector3& operator  =(const CdtVector<m_row, m_type>& v);     // vector1  = vector2
    CdtVector3& operator +=(const CdtVector<m_row, m_type>& v);     // vector1 += vector2
    CdtVector3& operator -=(const CdtVector<m_row, m_type>& v);     // vector1 -= vector2
    CdtVector3& operator *=(const CdtVector<m_row, m_type>& v);     // vector1 *= vector2, vector1(i) *= vector2(i)
    CdtVector3& operator /=(const CdtVector<m_row, m_type>& v);     // vector1 /= vector2, vector1(i) /= vector2(i)
    CdtVector3& operator  =(const CdtMatrix<m_row, 1, m_type>& v);  // vector1  = vector2
    CdtVector3& operator +=(const CdtMatrix<m_row, 1, m_type>& v);  // vector1 += vector2
    CdtVector3& operator -=(const CdtMatrix<m_row, 1, m_type>& v);  // vector1 -= vector2
    CdtVector3& operator *=(const CdtMatrix<m_row, 1, m_type>& v);  // vector1 *= vector2, vector1(i) *= vector2(i)
    CdtVector3& operator /=(const CdtMatrix<m_row, 1, m_type>& v);  // vector1 /= vector2, vector1(i) /= vector2(i)
    CdtVector3& operator  =(const m_type s);                        // vector1  = scalar, all elements set scalar
    CdtVector3& operator +=(const m_type s);                        // vector1 += scalar, vector1(i) += scalar
    CdtVector3& operator -=(const m_type s);                        // vector1 -= scalar, vector1(i) -= scalar
    CdtVector3& operator *=(const m_type s);                        // vector1 *= scalar
    CdtVector3& operator /=(const m_type s);                        // vector1 /= scalar
    CdtVector3& operator &=(const CdtVector3& v);                   // vector1 x vector2 cross product
    CdtVector3& operator &=(const CdtVector<m_row, m_type>& v);     // vector1 x vector2 cross product
    CdtVector3& operator &=(const CdtMatrix<m_row, 1, m_type>& v);  // vector1 x vector2 cross product
    CdtCommaInit<m_row, m_type> operator <<(const m_type s);        // Init first matrix elements

    /* Arithmetic operators */
    CdtVector3 operator -() const;                                      // minus sign
    CdtVector3 operator +(const CdtVector3& v) const;                   // vector + vector
    CdtVector3 operator -(const CdtVector3& v) const;                   // vector - vector
    CdtVector3 operator *(const CdtVector3& v) const;                   // vector * vector, vector(i) = vector1(i) * vector2(i)
    CdtVector3 operator /(const CdtVector3& v) const;                   // vector / vector, vector(i) = vector1(i) / vector2(i)
    CdtVector3 operator +(const CdtVector<m_row, m_type>& v) const;     // vector + vector
    CdtVector3 operator -(const CdtVector<m_row, m_type>& v) const;     // vector - vector
    CdtVector3 operator *(const CdtVector<m_row, m_type>& v) const;     // vector * vector, vector(i) = vector1(i) * vector2(i)
    CdtVector3 operator /(const CdtVector<m_row, m_type>& v) const;     // vector / vector, vector(i) = vector1(i) / vector2(i)
    CdtVector3 operator +(const CdtMatrix<m_row, 1, m_type>& v) const;  // vector + vector
    CdtVector3 operator -(const CdtMatrix<m_row, 1, m_type>& v) const;  // vector - vector
    CdtVector3 operator *(const CdtMatrix<m_row, 1, m_type>& v) const;  // vector * vector, vector(i) = vector1(i) * vector2(i)
    CdtVector3 operator /(const CdtMatrix<m_row, 1, m_type>& v) const;  // vector / vector, vector(i) = vector1(i) / vector2(i)
    CdtVector3 operator +(const m_type s) const;                        // vector + scalar, vector(i) = vector1(i) + scalar
    CdtVector3 operator -(const m_type s) const;                        // vector - scalar, vector(i) = vector1(i) - scalar
    CdtVector3 operator *(const m_type s) const;                        // vector * scalar
    CdtVector3 operator /(const m_type s) const;                        // vector / scalar
    CdtVector3 operator &(const CdtVector3& v) const;                   // vector x vector cross product
    CdtVector3 operator &(const CdtVector<m_row, m_type>& v) const;     // vector x vector cross product
    CdtVector3 operator &(const CdtMatrix<m_row, 1, m_type>& v) const;  // vector x vector cross product
    CdtMatrix3<m_type, m_row, m_row> operator &(const CdtMatrix3<m_type, m_row, m_row> &m) const;   // [v]x * RotMat, []x is skew-symmetric matrix
    CdtRotation<m_type, m_row, m_row> operator &(const CdtRotation<m_type, m_row, m_row> &m) const; // [v]x * RotMat, []x is skew-symmetric matrix

    template <uint16_t col>
    CdtMatrix<m_row, col, m_type> operator *(const CdtMatrix<1, col, m_type>& m) const; // vector1 * matrix(1xcol) outer product

    m_type dot(const CdtVector3& v) const;                  // vector1 * vector2 dot(inner) product
    m_type dot(const CdtVector<m_row, m_type>& v) const;    // vector1 * vector2 dot(inner) product
    m_type dot(const CdtMatrix<m_row, 1, m_type>& v) const; // vector1 * vector2 dot(inner) product

    /* Comparison operators */
    bool operator ==(const CdtVector3& v) const;                    // (true or false) vector1 == vector2
    bool operator !=(const CdtVector3& v) const;                    // (true or false) vector1 != vector2
    bool operator ==(const CdtVector<m_row, m_type>& v) const;      // (true or false) vector1 == vector2
    bool operator !=(const CdtVector<m_row, m_type>& v) const;      // (true or false) vector1 != vector2
    bool operator ==(const CdtMatrix<m_row, 1, m_type>& v) const;   // (true or false) vector1 == vector2
    bool operator !=(const CdtMatrix<m_row, 1, m_type>& v) const;   // (true or false) vector1 != vector2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row> friend class CdtVector3;
    template <uint16_t row, uint16_t col, typename type> friend class CdtMatrix;
    template <typename type, uint16_t row, uint16_t col> friend class CdtMatrix3;
    template <typename type, uint16_t row, uint16_t col> friend class CdtRotation;
    template <typename type, uint16_t row, uint16_t col> friend class CdtTransform;

    template <uint16_t row, typename type> friend class CdtVector;
    template <typename type, uint16_t row> friend class CdtVector4;
    template <typename type, uint16_t row> friend class CdtVector6;
    template <typename type, uint16_t row> friend class CdtQuaternion;

    /* Friend template function */
    template <typename type, uint16_t row>
    friend CdtVector3<type, row> operator+(const type s, const CdtVector3<type, row>& v);   // scalar + vector, scalar + vector(i)
    template <typename type, uint16_t row>
    friend CdtVector3<type, row> operator-(const type s, const CdtVector3<type, row>& v);   // scalar - vector, scalar - vector(i)
    template <typename type, uint16_t row>
    friend CdtVector3<type, row> operator*(const type s, const CdtVector3<type, row>& v);   // scalar * vector, scalar * vector(i)
    template <typename type, uint16_t row>
    friend CdtVector3<type, row> operator/(const type s, const CdtVector3<type, row>& v);   // scalar / vector, scalar / vector(i)
};

#include "dtVector3.tpp"

#endif // DTMATH_DTVECTOR3_H_
