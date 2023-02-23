/*!
\file       dtVector.h
\brief      dtMath, General Vector(m x 1) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTVECTOR_H_
#define DTMATH_DTVECTOR_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#include <string.h>
//#include <cstdarg>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

template <uint16_t m_size, typename m_type> class CdtCommaInit;
template <typename m_type, uint16_t m_row> class CdtVector3;
template <typename m_type, uint16_t m_row> class CdtVector4;
template <typename m_type, uint16_t m_row> class CdtVector6;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtMatrix;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtCscMatrix;
template <typename m_type, uint16_t m_row, uint16_t m_col> class CdtRotation;

template <uint16_t m_row, typename m_type = float>
class CdtVector
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row];
    CdtVector(const m_type* element);
    inline void CrossProduct(const m_type* v);

public:
    CdtVector();
    CdtVector(const m_type* element, const size_t n_byte);
    CdtVector(const CdtVector& v);
    CdtVector(const CdtMatrix<m_row, 1, m_type>& v);
    ~CdtVector() {}

    void SetZero();
    void SetFill(const m_type value);
    void SetElement(const m_type* element, const size_t n_byte);
    void SetElement(const CdtVector& v);
    void SetElement(const CdtMatrix<m_row, 1, m_type>& v);
    //void SetElement(const m_type elem0, ...);
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
    uint16_t GetDim() const { return m_row; }
    template <uint16_t row>
    CdtVector<row, m_type> GetBlock(const uint16_t idx);
    CdtVector3<m_type, 3> GetBlockVec3(const uint16_t idx);
    CdtVector4<m_type, 4> GetBlockVec4(const uint16_t idx);
    CdtVector6<m_type, 6> GetBlockVec6(const uint16_t idx);
    template <uint16_t row>
    int8_t GetBlock(const uint16_t idx, CdtVector<row, m_type> &v);
    int8_t GetBlockVec3(const uint16_t idx, CdtVector3<m_type, 3> &v);
    int8_t GetBlockVec4(const uint16_t idx, CdtVector4<m_type, 4> &v);
    int8_t GetBlockVec6(const uint16_t idx, CdtVector6<m_type, 6> &v);
    m_type GetNorm() const;
    m_type GetSqNorm() const;
    m_type GetSum() const;
    CdtVector GetNormalized() const;
    CdtMatrix<3, 3, m_type> GetSkew() const;
    CdtMatrix<1, m_row, m_type> Transpose() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type& operator ()(uint16_t irow) { return m_elem[irow]; }
    // returns a row of non-modifiable elements
    const m_type& operator ()(uint16_t irow) const { return m_elem[irow]; }

    /* Assignment operators */
    CdtVector& operator  =(const CdtVector& v);                     // vector1  = vector2
    CdtVector& operator +=(const CdtVector& v);                     // vector1 += vector2
    CdtVector& operator -=(const CdtVector& v);                     // vector1 -= vector2
    CdtVector& operator *=(const CdtVector& v);                     // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    CdtVector& operator /=(const CdtVector& v);                     // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    CdtVector& operator  =(const CdtVector3<m_type, m_row>& v);     // vector1  = vector2
    CdtVector& operator +=(const CdtVector3<m_type, m_row>& v);     // vector1 += vector2
    CdtVector& operator -=(const CdtVector3<m_type, m_row>& v);     // vector1 -= vector2
    CdtVector& operator *=(const CdtVector3<m_type, m_row>& v);     // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    CdtVector& operator /=(const CdtVector3<m_type, m_row>& v);     // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    CdtVector& operator  =(const CdtVector4<m_type, m_row>& v);     // vector1  = vector2
    CdtVector& operator +=(const CdtVector4<m_type, m_row>& v);     // vector1 += vector2
    CdtVector& operator -=(const CdtVector4<m_type, m_row>& v);     // vector1 -= vector2
    CdtVector& operator *=(const CdtVector4<m_type, m_row>& v);     // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    CdtVector& operator /=(const CdtVector4<m_type, m_row>& v);     // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    CdtVector& operator  =(const CdtVector6<m_type, m_row>& v);     // vector1  = vector2
    CdtVector& operator +=(const CdtVector6<m_type, m_row>& v);     // vector1 += vector2
    CdtVector& operator -=(const CdtVector6<m_type, m_row>& v);     // vector1 -= vector2
    CdtVector& operator *=(const CdtVector6<m_type, m_row>& v);     // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    CdtVector& operator /=(const CdtVector6<m_type, m_row>& v);     // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    CdtVector& operator  =(const CdtMatrix<m_row, 1, m_type>& v);   // vector1  = vector2
    CdtVector& operator +=(const CdtMatrix<m_row, 1, m_type>& v);   // vector1 += vector2
    CdtVector& operator -=(const CdtMatrix<m_row, 1, m_type>& v);   // vector1 -= vector2
    CdtVector& operator *=(const CdtMatrix<m_row, 1, m_type>& v);   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    CdtVector& operator /=(const CdtMatrix<m_row, 1, m_type>& v);   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    CdtVector& operator  =(const m_type s);                         // vector1  = scalar, all elements set scalar
    CdtVector& operator +=(const m_type s);                         // vector1 += scalar, vector1(i) = vector1(i) + scalar
    CdtVector& operator -=(const m_type s);                         // vector1 -= scalar, vector1(i) = vector1(i) - scalar
    CdtVector& operator *=(const m_type s);                         // vector1 *= scalar
    CdtVector& operator /=(const m_type s);                         // vector1 /= scalar
    CdtVector& operator &=(const CdtVector3<m_type, m_row>& v);     // vector1 x vector2 cross product
    CdtVector& operator &=(const CdtVector<m_row, m_type>& v);      // vector1 x vector2 cross product
    CdtVector& operator &=(const CdtMatrix<m_row, 1, m_type>& v);   // vector1 x vector2 cross product
    CdtCommaInit<m_row, m_type> operator <<(const m_type s);        // Init first matrix elements

    /* Arithmetic operators */
    CdtVector operator -() const;                                       // minus sign
    CdtVector operator +(const CdtVector& v) const;                     // vector + vector
    CdtVector operator -(const CdtVector& v) const;                     // vector - vector
    CdtVector operator *(const CdtVector& v) const;                     // vector * vector = m_elem(i) * v.m_elem(i)
    CdtVector operator /(const CdtVector& v) const;                     // vector / vector = m_elem(i) / v.m_elem(i)
    CdtVector operator +(const CdtVector3<m_type, m_row>& v) const;     // vector + vector
    CdtVector operator -(const CdtVector3<m_type, m_row>& v) const;     // vector - vector
    CdtVector operator *(const CdtVector3<m_type, m_row>& v) const;     // vector * vector = m_elem(i) * v.m_elem(i)
    CdtVector operator /(const CdtVector3<m_type, m_row>& v) const;     // vector / vector = m_elem(i) / v.m_elem(i)
    CdtVector operator +(const CdtVector4<m_type, m_row>& v) const;     // vector + vector
    CdtVector operator -(const CdtVector4<m_type, m_row>& v) const;     // vector - vector
    CdtVector operator *(const CdtVector4<m_type, m_row>& v) const;     // vector * vector = m_elem(i) * v.m_elem(i)
    CdtVector operator /(const CdtVector4<m_type, m_row>& v) const;     // vector / vector = m_elem(i) / v.m_elem(i)
    CdtVector operator +(const CdtVector6<m_type, m_row>& v) const;     // vector + vector
    CdtVector operator -(const CdtVector6<m_type, m_row>& v) const;     // vector - vector
    CdtVector operator *(const CdtVector6<m_type, m_row>& v) const;     // vector * vector = m_elem(i) * v.m_elem(i)
    CdtVector operator /(const CdtVector6<m_type, m_row>& v) const;     // vector / vector = m_elem(i) / v.m_elem(i)
    CdtVector operator +(const CdtMatrix<m_row, 1, m_type>& v) const;   // vector + vector
    CdtVector operator -(const CdtMatrix<m_row, 1, m_type>& v) const;   // vector - vector
    CdtVector operator *(const CdtMatrix<m_row, 1, m_type>& v) const;   // vector * vector = m_elem(i) * v.m_elem(i)
    CdtVector operator /(const CdtMatrix<m_row, 1, m_type>& v) const;   // vector / vector = m_elem(i) / v.m_elem(i)
    CdtVector operator +(const m_type s) const;                         // vector + scalar = m_elem(i) + scalar
    CdtVector operator -(const m_type s) const;                         // vector - scalar = m_elem(i) - scalar
    CdtVector operator *(const m_type s) const;                         // vector * scalar
    CdtVector operator /(const m_type s) const;                         // vector / scalar
    CdtVector operator &(const CdtVector3<m_type, m_row>& v) const;     // vector x vector cross product
    CdtVector operator &(const CdtVector<m_row, m_type>& v) const;      // vector x vector cross product
    CdtVector operator &(const CdtMatrix<m_row, 1, m_type>& v) const;   // vector x vector cross product
    CdtMatrix3<m_type, m_row, m_row> operator &(const CdtMatrix3<m_type, m_row, m_row> &m) const;   // [v]x * RotMat, []x is skew-symmetric matrix
    CdtRotation<m_type, m_row, m_row> operator &(const CdtRotation<m_type, m_row, m_row> &m) const; // [v]x * RotMat, []x is skew-symmetric matrix

    template <uint16_t col>
    CdtMatrix<m_row, col, m_type> operator *(const CdtMatrix<1, col, m_type>& m) const; // vector1 * matrix(1xcol) outer product

    m_type dot(const CdtVector& v) const;                   // vector1 * vector2 dot(inner) product
    m_type dot(const CdtVector3<m_type, m_row>& v) const;   // vector1 * vector2 dot(inner) product
    m_type dot(const CdtVector4<m_type, m_row>& v) const;   // vector1 * vector2 dot(inner) product
    m_type dot(const CdtVector6<m_type, m_row>& v) const;   // vector1 * vector2 dot(inner) product
    m_type dot(const CdtMatrix<m_row, 1, m_type>& v) const; // vector1 * vector2 dot(inner) product

    /* Comparison operators */
    bool operator ==(const CdtVector& v) const;                     // (true or false) vector1 == vector2
    bool operator !=(const CdtVector& v) const;                     // (true or false) vector1 != vector2
    bool operator ==(const CdtVector3<m_type, m_row>& v) const;     // (true or false) vector1 == vector2
    bool operator !=(const CdtVector3<m_type, m_row>& v) const;     // (true or false) vector1 != vector2
    bool operator ==(const CdtVector4<m_type, m_row>& v) const;     // (true or false) vector1 == vector2
    bool operator !=(const CdtVector4<m_type, m_row>& v) const;     // (true or false) vector1 != vector2
    bool operator ==(const CdtVector6<m_type, m_row>& v) const;     // (true or false) vector1 == vector2
    bool operator !=(const CdtVector6<m_type, m_row>& v) const;     // (true or false) vector1 != vector2
    bool operator ==(const CdtMatrix<m_row, 1, m_type>& v) const;   // (true or false) vector1 == vector2
    bool operator !=(const CdtMatrix<m_row, 1, m_type>& v) const;   // (true or false) vector1 != vector2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <uint16_t row, typename type> friend class CdtVector;
    template <uint16_t row, uint16_t col, typename type> friend class CdtMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class CdtCscMatrix;
    template <typename type, uint16_t row, uint16_t col> friend class CdtMatrix3;
    template <typename type, uint16_t row, uint16_t col> friend class CdtRotation;
    template <typename type, uint16_t row> friend class CdtVector3;
    template <typename type, uint16_t row> friend class CdtVector4;
    template <typename type, uint16_t row> friend class CdtVector6;
    template <typename type, uint16_t row> friend class CdtQuaternion;

    template <uint16_t row, uint16_t col, typename type> friend class CdtLowerTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class CdtUpperTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class CdtNoPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class CdtPartialPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class CdtLLT;
    template <uint16_t row, uint16_t col, typename type> friend class CdtLDLT;
    template <uint16_t row, uint16_t col, typename type> friend class CdtSVD;

    /* Friend template function */
    template <uint16_t row, typename type>
    friend CdtVector<row, type> operator +(const type s, const CdtVector<row, type>& v);    // scalar + vector = scalar + m_elem(i)
    template <uint16_t row, typename type>
    friend CdtVector<row, type> operator -(const type s, const CdtVector<row, type>& v);    // scalar - vector = scalar - m_elem(i)
    template <uint16_t row, typename type>
    friend CdtVector<row, type> operator *(const type s, const CdtVector<row, type>& v);    // scalar * vector = scalar * m_elem(i)
    template <uint16_t row, typename type>
    friend CdtVector<row, type> operator /(const type s, const CdtVector<row, type>& v);    // scalar / vector = scalar / m_elem(i)
};

#include "dtVector.tpp"

#endif // DTMATH_DTVECTOR_H_
