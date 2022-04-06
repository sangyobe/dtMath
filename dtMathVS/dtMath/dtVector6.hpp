/*!
\file       dtVector6.hpp
\brief      dtMath, 6x1 Vector class, lighter and faster than general vector class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#pragma once

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
template <typename m_type, uint16_t m_row> class CdtVector3;
template <typename m_type, uint16_t m_row> class CdtVector4;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtMatrix;

template <typename m_type = float, uint16_t m_row = 6>
class CdtVector6
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row];
    CdtVector6(const m_type* element);

public:
    CdtVector6();
    CdtVector6(const m_type* element, const size_t n_byte);
    CdtVector6(const m_type px, const m_type py, const m_type pz, const m_type ox, const m_type oy, const m_type oz);
    CdtVector6(const CdtVector6& v);
    CdtVector6(const CdtVector3<m_type, 3>& p, const CdtVector3<m_type, 3>& o);
    CdtVector6(const CdtVector<m_row, m_type>& v);
    CdtVector6(const CdtMatrix<m_row, 1, m_type>& v);
    ~CdtVector6() {}

    void SetZero();
    void SetFill(const m_type value);
    void SetElement(const m_type* element, const size_t n_byte);
    void SetElement(const m_type px, const m_type py, const m_type pz, const m_type ox, const m_type oy, const m_type oz);
    void SetElement(const CdtVector6& v);
    void SetElement(const CdtVector3<m_type, 3>& p, const CdtVector3<m_type, 3>& o);
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
    CdtVector3<m_type, 3> GetBlockVec3(const uint16_t idx);
    CdtVector4<m_type, 4> GetBlockVec4(const uint16_t idx);
    template <uint16_t row>
    int8_t GetBlock(const uint16_t idx, CdtVector<row, m_type>& v);
    int8_t GetBlockVec3(const uint16_t idx, CdtVector3<m_type, 3>& v);
    int8_t GetBlockVec4(const uint16_t idx, CdtVector4<m_type, 4>& v);
    CdtVector3<m_type, 3> GetPos() const;
    CdtVector3<m_type, 3> GetOri() const;
    m_type GetNorm() const;
    m_type GetSqNorm() const;
    m_type GetSum() const;
    CdtVector6 GetNormalized() const;
    CdtMatrix<1, m_row, m_type> Transpose() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type& operator ()(uint16_t irow) { return m_elem[irow]; }
    // returns a row of non-modifiable elements
    const m_type& operator ()(uint16_t irow) const { return m_elem[irow]; }

    /* Assignment operators */
    CdtVector6& operator  =(const CdtVector6& v);                   // vector1  = vector2
    CdtVector6& operator +=(const CdtVector6& v);                   // vector1 += vector2
    CdtVector6& operator -=(const CdtVector6& v);                   // vector1 -= vector2
    CdtVector6& operator *=(const CdtVector6& v);                   // vector1 *= vector2, vector1(i) *= vector2(i)
    CdtVector6& operator /=(const CdtVector6& v);                   // vector1 /= vector2, vector1(i) /= vector2(i)
    CdtVector6& operator  =(const CdtVector<m_row, m_type>& v);     // vector1  = vector2
    CdtVector6& operator +=(const CdtVector<m_row, m_type>& v);     // vector1 += vector2
    CdtVector6& operator -=(const CdtVector<m_row, m_type>& v);     // vector1 -= vector2
    CdtVector6& operator *=(const CdtVector<m_row, m_type>& v);     // vector1 *= vector2, vector1(i) *= vector2(i)
    CdtVector6& operator /=(const CdtVector<m_row, m_type>& v);     // vector1 /= vector2, vector1(i) /= vector2(i)
    CdtVector6& operator  =(const CdtMatrix<m_row, 1, m_type>& v);  // vector1  = vector2
    CdtVector6& operator +=(const CdtMatrix<m_row, 1, m_type>& v);  // vector1 += vector2
    CdtVector6& operator -=(const CdtMatrix<m_row, 1, m_type>& v);  // vector1 -= vector2
    CdtVector6& operator *=(const CdtMatrix<m_row, 1, m_type>& v);  // vector1 *= vector2, vector1(i) *= vector2(i)
    CdtVector6& operator /=(const CdtMatrix<m_row, 1, m_type>& v);  // vector1 /= vector2, vector1(i) /= vector2(i)
    CdtVector6& operator  =(const m_type s);                        // vector1  = scalar, all elements set scalar
    CdtVector6& operator +=(const m_type s);                        // vector1 += scalar, vector1(i) += scalar
    CdtVector6& operator -=(const m_type s);                        // vector1 -= scalar, vector1(i) -= scalar
    CdtVector6& operator *=(const m_type s);                        // vector1 *= scalar
    CdtVector6& operator /=(const m_type s);                        // vector1 /= scalar
    CdtCommaInit<m_row, m_type> operator <<(const m_type s);        // Init first matrix elements

    /* Arithmetic operators */
    CdtVector6 operator -() const;                                      // minus sign
    CdtVector6 operator +(const CdtVector6& v) const;                   // vector + vector
    CdtVector6 operator -(const CdtVector6& v) const;                   // vector - vector
    CdtVector6 operator *(const CdtVector6& v) const;                   // vector * vector, vector(i) = vector1(i) * vector2(i)
    CdtVector6 operator /(const CdtVector6& v) const;                   // vector / vector, vector(i) = vector1(i) / vector2(i)
    CdtVector6 operator +(const CdtVector<m_row, m_type>& v) const;     // vector + vector
    CdtVector6 operator -(const CdtVector<m_row, m_type>& v) const;     // vector - vector
    CdtVector6 operator *(const CdtVector<m_row, m_type>& v) const;     // vector * vector, vector(i) = vector1(i) * vector2(i)
    CdtVector6 operator /(const CdtVector<m_row, m_type>& v) const;     // vector / vector, vector(i) = vector1(i) / vector2(i)
    CdtVector6 operator +(const CdtMatrix<m_row, 1, m_type>& v) const;  // vector + vector
    CdtVector6 operator -(const CdtMatrix<m_row, 1, m_type>& v) const;  // vector - vector
    CdtVector6 operator *(const CdtMatrix<m_row, 1, m_type>& v) const;  // vector * vector, vector(i) = vector1(i) * vector2(i)
    CdtVector6 operator /(const CdtMatrix<m_row, 1, m_type>& v) const;  // vector / vector, vector(i) = vector1(i) / vector2(i)
    CdtVector6 operator +(const m_type s) const;                        // vector + scalar, vector(i) = vector1(i) + scalar
    CdtVector6 operator -(const m_type s) const;                        // vector - scalar, vector(i) = vector1(i) - scalar
    CdtVector6 operator *(const m_type s) const;                        // vector * scalar
    CdtVector6 operator /(const m_type s) const;                        // vector / scalar

    template <uint16_t col>
    CdtMatrix<m_row, col, m_type> operator *(const CdtMatrix<1, col, m_type>& m) const;

    m_type dot(const CdtVector6& v) const;                  // vector * vector dot(inner) product
    m_type dot(const CdtVector<m_row, m_type>& v) const;    // vector * vector dot(inner) product
    m_type dot(const CdtMatrix<m_row, 1, m_type>& v) const; // Vector * vector dot(inner) product

    /* Comparison operators */
    bool operator ==(const CdtVector6& v) const;                    // (true or false) vector1 == vector2
    bool operator !=(const CdtVector6& v) const;                    // (true or false) vector1 != vector2
    bool operator ==(const CdtVector<m_row, m_type>& v) const;      // (true or false) vector1 == vector2
    bool operator !=(const CdtVector<m_row, m_type>& v) const;      // (true or false) vector1 != vector2
    bool operator ==(const CdtMatrix<m_row, 1, m_type>& v) const;   // (true or false) vector1 == vector2
    bool operator !=(const CdtMatrix<m_row, 1, m_type>& v) const;   // (true or false) vector1 != vector2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row> friend class CdtVector6;
    template <uint16_t row, uint16_t col, typename type> friend class CdtMatrix;
    template <uint16_t row, typename type> friend class CdtVector;
    template <typename type, uint16_t row> friend class CdtVector3;
    template <typename type, uint16_t row> friend class CdtVector4;

    /* Friend template function */
    template <typename type, uint16_t row>
    friend CdtVector6<type, row> operator +(const type s, const CdtVector6<type, row>& v);  // scalar + vector, scalar + vector(i)
    template <typename type, uint16_t row>
    friend CdtVector6<type, row> operator -(const type s, const CdtVector6<type, row>& v);  // scalar - vector, scalar - vector(i)
    template <typename type, uint16_t row>
    friend CdtVector6<type, row> operator *(const type s, const CdtVector6<type, row>& v);  // scalar * vector, scalar * vector(i)
    template <typename type, uint16_t row>
    friend CdtVector6<type, row> operator /(const type s, const CdtVector6<type, row>& v);  // scalar / vector, scalar / vector(i)
};

#include "dtVector6.ipp"
