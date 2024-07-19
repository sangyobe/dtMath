/*!
\file       dtVector.h
\brief      dtMath, General Vector(m x 1) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTVECTOR_H_
#define DTMATH_DTVECTOR_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#include <string.h>
// #include <cstdarg>
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
template <typename t_type, uint16_t t_row> class Vector3;
template <typename t_type, uint16_t t_row> class Vector4;
template <typename t_type, uint16_t t_row> class Vector6;
template <typename t_type, uint16_t t_row> class Quaternion;
template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;
template <uint16_t t_row, uint16_t t_col, typename t_type> class CscMatrix;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Matrix3;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Rotation;

template <uint16_t t_row, typename t_type = float>
class Vector
{
private:
    t_type m_tolerance = std::numeric_limits<t_type>::epsilon();
    t_type m_elem[t_row];
    Vector(const t_type *element);
    inline void CrossProduct(const t_type *v);

public:
    Vector();
    Vector(const t_type *element, const size_t n_byte);
    Vector(const Vector &v);
    Vector(const Vector<0, t_type> &v);
    Vector(const Matrix<t_row, 1, t_type> &v);
    Vector(const Matrix<0, 0, t_type> &v);
    ~Vector() {}

    void SetZero();
    void SetFill(const t_type value);
    void SetElement(const t_type *element, const size_t n_byte);
    // void SetElement(const Vector &v);
    // void SetElement(const Matrix<t_row, 1, t_type> &v);
    // void SetElement(const Matrix<0, 0, t_type> &v);
    // void SetElement(const t_type elem0, ...);
    template <uint16_t row>
    void SetBlock(const uint16_t idxRow, const Vector<row, t_type> &v, const uint16_t jdx = 0, const int16_t size = row);
    void SetBlock(const uint16_t idxRow, const t_type *v, const size_t n_byte);
    void SetBlock(const uint16_t idxRow, const Vector3<t_type, 3> &v, const uint16_t jdx = 0, const int16_t size = 3);
    void SetBlock(const uint16_t idxRow, const Vector4<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4);
    void SetBlock(const uint16_t idxRow, const Vector6<t_type, 6> &v, const uint16_t jdx = 0, const int16_t size = 6);
    void SetBlock(const uint16_t idxRow, const Quaternion<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4);
    template <uint16_t row>
    void SetBlock(const uint16_t idxRow, const Matrix<row, 1, t_type> &v, const uint16_t jdx = 0, const int16_t size = row);
    void SetBlock(const uint16_t idxRow, const Matrix<0, 0, t_type> &v, const uint16_t jdx = 0, int16_t size = -1);
    void SetBlock(const uint16_t idxRow, const Vector<0, t_type> &v, const uint16_t jdx = 0, int16_t size = -1);
    void SetSwap(const uint16_t i, const uint16_t j);
    void SetNormalize();

    const t_type *const GetElementsAddr() const;
    uint16_t GetDim() const { return t_row; }
    template <uint16_t row>
    int8_t GetBlock(const uint16_t idx, Vector<row, t_type> &v, const uint16_t jdx = 0, const int16_t size = row);
    int8_t GetBlock(const uint16_t idx, Vector<0, t_type> &v, const uint16_t jdx = 0, int16_t size = -1);
    int8_t GetBlock(const uint16_t idx, Vector3<t_type, 3> &v, const uint16_t jdx = 0, const int16_t size = 3);
    int8_t GetBlock(const uint16_t idx, Vector4<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4);
    int8_t GetBlock(const uint16_t idx, Vector6<t_type, 6> &v, const uint16_t jdx = 0, const int16_t size = 6);
    int8_t GetBlock(const uint16_t idx, Quaternion<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4);
    template <uint16_t row>
    Vector<row, t_type> GetBlockVec(const uint16_t idx, const uint16_t jdx = 0, const int16_t size = row);
    Vector<0, t_type> GetBlockVec0(const uint16_t idx, const uint16_t row, const uint16_t jdx = 0, int16_t size = -1);
    Vector3<t_type, 3> GetBlockVec3(const uint16_t idx, const uint16_t jdx = 0, const int16_t size = 3);
    Vector4<t_type, 4> GetBlockVec4(const uint16_t idx, const uint16_t jdx = 0, const int16_t size = 4);
    Vector6<t_type, 6> GetBlockVec6(const uint16_t idx, const uint16_t jdx = 0, const int16_t size = 6);
    Quaternion<t_type, 4> GetBlockQuat(const uint16_t idx, const uint16_t jdx = 0, const int16_t size = 4);

    t_type GetNorm() const;              // Frobenius Norm (Euclidean norm, L2 Norm)
    t_type GetSqNorm() const;            // Squared Frobenius Norm (Euclidean norm, Squared L2 Norm)
    t_type GetLpNorm(const int p) const; // Generalized Norm (Lp Norm)
    t_type GetSum() const;
    Vector GetNormalized() const;
    Matrix<3, 3, t_type> GetSkew() const;
    Matrix<1, t_row, t_type> Transpose() const;
    void Transpose(Matrix<1, t_row, t_type> &m) const;
    void Transpose(Matrix<0, 0, t_type> &m) const;

    /* Member access operators */
    t_type &operator()(uint16_t irow);             // returns a row of modifiable elements
    const t_type &operator()(uint16_t irow) const; // returns a row of non-modifiable elements

    /* Assignment operators */
    Vector &operator=(const Vector &v);                    // vector1  = vector2
    Vector &operator+=(const Vector &v);                   // vector1 += vector2
    Vector &operator-=(const Vector &v);                   // vector1 -= vector2
    Vector &CWiseMulEq(const Vector &v);                   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &CWiseDivEq(const Vector &v);                   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const Vector<0, t_type> &v);         // vector1  = vector2
    Vector &operator+=(const Vector<0, t_type> &v);        // vector1 += vector2
    Vector &operator-=(const Vector<0, t_type> &v);        // vector1 -= vector2
    Vector &CWiseMulEq(const Vector<0, t_type> &v);        // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &CWiseDivEq(const Vector<0, t_type> &v);        // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const Vector3<t_type, t_row> &v);    // vector1  = vector2
    Vector &operator+=(const Vector3<t_type, t_row> &v);   // vector1 += vector2
    Vector &operator-=(const Vector3<t_type, t_row> &v);   // vector1 -= vector2
    Vector &CWiseMulEq(const Vector3<t_type, t_row> &v);   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &CWiseDivEq(const Vector3<t_type, t_row> &v);   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const Vector4<t_type, t_row> &v);    // vector1  = vector2
    Vector &operator+=(const Vector4<t_type, t_row> &v);   // vector1 += vector2
    Vector &operator-=(const Vector4<t_type, t_row> &v);   // vector1 -= vector2
    Vector &CWiseMulEq(const Vector4<t_type, t_row> &v);   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &CWiseDivEq(const Vector4<t_type, t_row> &v);   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const Vector6<t_type, t_row> &v);    // vector1  = vector2
    Vector &operator+=(const Vector6<t_type, t_row> &v);   // vector1 += vector2
    Vector &operator-=(const Vector6<t_type, t_row> &v);   // vector1 -= vector2
    Vector &CWiseMulEq(const Vector6<t_type, t_row> &v);   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &CWiseDivEq(const Vector6<t_type, t_row> &v);   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const Matrix<t_row, 1, t_type> &v);  // vector1  = vector2
    Vector &operator+=(const Matrix<t_row, 1, t_type> &v); // vector1 += vector2
    Vector &operator-=(const Matrix<t_row, 1, t_type> &v); // vector1 -= vector2
    Vector &CWiseMulEq(const Matrix<t_row, 1, t_type> &v); // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &CWiseDivEq(const Matrix<t_row, 1, t_type> &v); // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const Matrix<0, 0, t_type> &v);      // vector1  = vector2
    Vector &operator+=(const Matrix<0, 0, t_type> &v);     // vector1 += vector2
    Vector &operator-=(const Matrix<0, 0, t_type> &v);     // vector1 -= vector2
    Vector &CWiseMulEq(const Matrix<0, 0, t_type> &v);     // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &CWiseDivEq(const Matrix<0, 0, t_type> &v);     // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const t_type s);                     // vector1  = scalar, all elements set scalar
    Vector &operator+=(const t_type s);                    // vector1 += scalar, vector1(i) = vector1(i) + scalar
    Vector &operator-=(const t_type s);                    // vector1 -= scalar, vector1(i) = vector1(i) - scalar
    Vector &operator*=(const t_type s);                    // vector1 *= scalar
    Vector &operator/=(const t_type s);                    // vector1 /= scalar
    Vector &operator&=(const Vector3<t_type, t_row> &v);   // vector1 x vector2 cross product
    Vector &operator&=(const Vector<3, t_type> &v);        // vector1 x vector2 cross product
    Vector &operator&=(const Vector<0, t_type> &v);        // vector1 x vector2 cross product
    Vector &operator&=(const Matrix<3, 1, t_type> &v);     // vector1 x vector2 cross product
    Vector &operator&=(const Matrix<0, 0, t_type> &v);     // vector1 x vector2 cross product
    CommaInit<t_row, t_type> operator<<(const t_type s);   // Init first matrix elements

    /* Arithmetic operators */
    Vector operator-() const;                                  // minus sign
    Vector operator+(const Vector &v) const;                   // vector + vector
    Vector operator-(const Vector &v) const;                   // vector - vector
    Vector CWiseMul(const Vector &v) const;                    // vector * vector = m_elem(i) * v.m_elem(i)
    Vector CWiseDiv(const Vector &v) const;                    // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const Vector<0, t_type> &v) const;        // vector + vector
    Vector operator-(const Vector<0, t_type> &v) const;        // vector - vector
    Vector CWiseMul(const Vector<0, t_type> &v) const;         // vector * vector = m_elem(i) * v.m_elem(i)
    Vector CWiseDiv(const Vector<0, t_type> &v) const;         // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const Vector3<t_type, t_row> &v) const;   // vector + vector
    Vector operator-(const Vector3<t_type, t_row> &v) const;   // vector - vector
    Vector CWiseMul(const Vector3<t_type, t_row> &v) const;    // vector * vector = m_elem(i) * v.m_elem(i)
    Vector CWiseDiv(const Vector3<t_type, t_row> &v) const;    // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const Vector4<t_type, t_row> &v) const;   // vector + vector
    Vector operator-(const Vector4<t_type, t_row> &v) const;   // vector - vector
    Vector CWiseMul(const Vector4<t_type, t_row> &v) const;    // vector * vector = m_elem(i) * v.m_elem(i)
    Vector CWiseDiv(const Vector4<t_type, t_row> &v) const;    // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const Vector6<t_type, t_row> &v) const;   // vector + vector
    Vector operator-(const Vector6<t_type, t_row> &v) const;   // vector - vector
    Vector CWiseMul(const Vector6<t_type, t_row> &v) const;    // vector * vector = m_elem(i) * v.m_elem(i)
    Vector CWiseDiv(const Vector6<t_type, t_row> &v) const;    // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const Matrix<t_row, 1, t_type> &v) const; // vector + vector
    Vector operator-(const Matrix<t_row, 1, t_type> &v) const; // vector - vector
    Vector CWiseMul(const Matrix<t_row, 1, t_type> &v) const;  // vector * vector = m_elem(i) * v.m_elem(i)
    Vector CWiseDiv(const Matrix<t_row, 1, t_type> &v) const;  // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const Matrix<0, 0, t_type> &v) const;     // vector + vector
    Vector operator-(const Matrix<0, 0, t_type> &v) const;     // vector - vector
    Vector CWiseMul(const Matrix<0, 0, t_type> &v) const;      // vector * vector = m_elem(i) * v.m_elem(i)
    Vector CWiseDiv(const Matrix<0, 0, t_type> &v) const;      // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const t_type s) const;                    // vector + scalar = m_elem(i) + scalar
    Vector operator-(const t_type s) const;                    // vector - scalar = m_elem(i) - scalar
    Vector operator*(const t_type s) const;                    // vector * scalar
    Vector operator/(const t_type s) const;                    // vector / scalar
    template <uint16_t col>
    Matrix<t_row, col, t_type> operator*(const Matrix<1, col, t_type> &m) const; // vector(mx1) * matrix(1xcol) = matrix(mxcol)
    Matrix<0, 0, t_type> operator*(const Matrix<0, 0, t_type> &m) const;         // vector(mx1) * matrix(1xn) = matrix(mxn)
    Vector operator&(const Vector3<t_type, t_row> &v) const;                     // vector x vector cross product
    Vector operator&(const Vector<t_row, t_type> &v) const;                      // vector x vector cross product
    Vector operator&(const Vector<0, t_type> &v) const;                          // vector x vector cross product
    Vector operator&(const Matrix<t_row, 1, t_type> &v) const;                   // vector x vector cross product
    Vector operator&(const Matrix<0, 0, t_type> &v) const;                       // vector x vector cross product
    Matrix3<t_type, 3, 3> operator&(const Matrix3<t_type, 3, 3> &m) const;       // [v]x * RotMat, []x is skew-symmetric matrix
    Rotation<t_type, 3, 3> operator&(const Rotation<t_type, 3, 3> &m) const;     // [v]x * RotMat, []x is skew-symmetric matrix

    Matrix<t_row, 3, t_type> Outer(const Vector3<t_type, 3> &v) const; // vector1 * vector2, outer product, vector1 * vector2.Transpos()
    Matrix<t_row, 4, t_type> Outer(const Vector4<t_type, 4> &v) const; // vector1 * vector2, outer product, vector1 * vector2.Transpos()
    Matrix<t_row, 6, t_type> Outer(const Vector6<t_type, 6> &v) const; // vector1 * vector2, outer product, vector1 * vector2.Transpos()
    Matrix<0, 0, t_type> Outer(const Vector<0, t_type> &v) const;      // vector1 * vector2 outer product, vector1 * vector2.Transpos()
    Matrix<0, 0, t_type> Outer(const Matrix<0, 0, t_type> &v) const;   // vector1 * vector2 outer product, vector1 * vector2.Transpos()
    template <uint16_t row>
    Matrix<t_row, row, t_type> Outer(const Vector<row, t_type> &v) const; // vector1 * vector2 outer product, vector1 * vector2.Transpos()
    template <uint16_t row>
    Matrix<t_row, row, t_type> Outer(const Matrix<row, 1, t_type> &v) const; // vector1 * vector2 outer product, vector1 * vector2.Transpose()

    t_type Inner(const Vector &v) const;                   // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2
    t_type Inner(const Vector<0, t_type> &v) const;        // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2
    t_type Inner(const Vector3<t_type, t_row> &v) const;   // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2
    t_type Inner(const Vector4<t_type, t_row> &v) const;   // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2
    t_type Inner(const Vector6<t_type, t_row> &v) const;   // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2
    t_type Inner(const Matrix<t_row, 1, t_type> &v) const; // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2
    t_type Inner(const Matrix<0, 0, t_type> &v) const;     // vector1 * vector2 inner(dot) product, vector1.Transpose() * vector2

    /* Comparison operators */
    bool operator==(const Vector &v) const;                   // (true or false) vector1 == vector2
    bool operator!=(const Vector &v) const;                   // (true or false) vector1 != vector2
    bool operator==(const Vector<0, t_type> &v) const;        // (true or false) vector1 == vector2
    bool operator!=(const Vector<0, t_type> &v) const;        // (true or false) vector1 != vector2
    bool operator==(const Vector3<t_type, t_row> &v) const;   // (true or false) vector1 == vector2
    bool operator!=(const Vector3<t_type, t_row> &v) const;   // (true or false) vector1 != vector2
    bool operator==(const Vector4<t_type, t_row> &v) const;   // (true or false) vector1 == vector2
    bool operator!=(const Vector4<t_type, t_row> &v) const;   // (true or false) vector1 != vector2
    bool operator==(const Vector6<t_type, t_row> &v) const;   // (true or false) vector1 == vector2
    bool operator!=(const Vector6<t_type, t_row> &v) const;   // (true or false) vector1 != vector2
    bool operator==(const Matrix<t_row, 1, t_type> &v) const; // (true or false) vector1 == vector2
    bool operator!=(const Matrix<t_row, 1, t_type> &v) const; // (true or false) vector1 != vector2
    bool operator==(const Matrix<0, 0, t_type> &v) const;     // (true or false) vector1 == vector2
    bool operator!=(const Matrix<0, 0, t_type> &v) const;     // (true or false) vector1 != vector2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <uint16_t row, typename type> friend class Vector;
    template <uint16_t row, uint16_t col, typename type> friend class Matrix;
    template <uint16_t row, uint16_t col, typename type> friend class CscMatrix;
    template <typename type, uint16_t row, uint16_t col> friend class Matrix3;
    template <typename type, uint16_t row, uint16_t col> friend class Rotation;
    template <typename type, uint16_t row> friend class Vector3;
    template <typename type, uint16_t row> friend class Vector4;
    template <typename type, uint16_t row> friend class Vector6;
    template <typename type, uint16_t row> friend class Quaternion;

    template <uint16_t row, uint16_t col, typename type> friend class LowerMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class UpperMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class NoPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class PartialPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class FullPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class LLT;
    template <uint16_t row, uint16_t col, typename type> friend class LDLT;
    template <uint16_t row, uint16_t col, typename type> friend class SVD;

    template <uint16_t row, uint16_t col, typename type> friend class CscAMD;
    template <uint16_t row, uint16_t col, typename type> friend class CscLLT;

    /* Friend template function */
    template <uint16_t row, typename type>
    friend Vector<row, type> operator+(const type s, const Vector<row, type> &v); // scalar + vector = scalar + m_elem(i)
    template <uint16_t row, typename type>
    friend Vector<row, type> operator-(const type s, const Vector<row, type> &v); // scalar - vector = scalar - m_elem(i)
    template <uint16_t row, typename type>
    friend Vector<row, type> operator*(const type s, const Vector<row, type> &v); // scalar * vector = scalar * m_elem(i)
    template <uint16_t row, typename type>
    friend Vector<row, type> operator/(const type s, const Vector<row, type> &v); // scalar / vector = scalar / m_elem(i)
};

} // namespace Math
} // namespace dt

#include "dtVector.tpp"

#include "dtVector0.h"

#endif // DTMATH_DTVECTOR_H_
