/*!
\file       dtVector3.h
\brief      dtMath, 3x1 Vector class, lighter and faster than general vector class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
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

namespace dt
{
namespace Math
{

template <uint16_t m_size, typename t_type> class CommaInit;
template <uint16_t t_row, typename t_type> class Vector;
template <typename t_type, uint16_t t_row> class Vector4;
template <typename t_type, uint16_t t_row> class Vector6;
template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Matrix3;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Rotation;

template <typename t_type = float, uint16_t t_row = 3>
class Vector3
{
private:
    t_type m_tolerance = std::numeric_limits<t_type>::epsilon();
    t_type m_elem[t_row];
    Vector3(const t_type *element);
    inline void CrossProduct(const t_type *v);

public:
    Vector3();
    Vector3(const t_type *element, const size_t n_byte);
    Vector3(const t_type i, const t_type j, const t_type k);
    Vector3(const Vector3 &v);
    Vector3(const Vector<t_row, t_type> &v);
    Vector3(const Vector<0, t_type> &v);
    Vector3(const Matrix<t_row, 1, t_type> &v);
    Vector3(const Matrix<0, 0, t_type> &v);
    ~Vector3() {}

    void SetZero();
    void SetFill(const t_type value);
    void SetElement(const t_type *element, const size_t n_byte);
    void SetElement(const t_type i, const t_type j, const t_type k);
    // void SetElement(const Vector3 &v);
    // void SetElement(const Vector<t_row, t_type> &v);
    // void SetElement(const Matrix<t_row, 1, t_type> &v);
    // void SetElement(const Matrix<0, 0, t_type> &v);
    template <uint16_t row>
    void SetBlock(const uint16_t idxRow, const Vector<row, t_type> &v);
    void SetBlock(const uint16_t idxRow, const t_type *v, const size_t n_byte);
    void SetBlock(const uint16_t idxRow, const Vector3<t_type, 3> &v);
    void SetBlock(const uint16_t idxRow, const Vector4<t_type, 4> &v);
    void SetBlock(const uint16_t idxRow, const Vector6<t_type, 6> &v);
    template <uint16_t row>
    void SetBlock(const uint16_t idxRow, const Matrix<row, 1, t_type> &v);
    void SetBlock(const uint16_t idxRow, const Matrix<0, 0, t_type> &v);
    void SetBlock(const uint16_t idxRow, const Vector<0, t_type> &v);
    void SetSwap(const uint16_t i, const uint16_t j);
    void SetNormalize();

    const t_type *const GetElementsAddr() const;
    template <uint16_t row>
    Vector<row, t_type> GetBlock(const uint16_t idx);
    Vector<0, t_type> GetBlock(const uint16_t idx, const uint16_t row);
    template <uint16_t row>
    int8_t GetBlock(const uint16_t idx, Vector<row, t_type> &v);
    int8_t GetBlock(const uint16_t idx, Vector<0, t_type> &v);

    t_type GetNorm() const;              // Frobenius Norm (Euclidean norm, L2 Norm)
    t_type GetSqNorm() const;            // Squared Frobenius Norm (Euclidean norm, Squared L2 Norm)
    t_type GetLpNorm(const int p) const; // Generalized Norm (Lp Norm)
    t_type GetSum() const;
    Vector3 GetNormalized() const;
    Matrix3<t_type, t_row, t_row> GetSkew() const;
    Matrix<1, t_row, t_type> Transpose() const;
    void Transpose(Matrix<1, t_row, t_type> &m) const;
    void Transpose(Matrix<0, 0, t_type> &m) const;

    /* Member access operators */
    t_type &operator()(uint16_t irow);             // returns a row of modifiable elements
    const t_type &operator()(uint16_t irow) const; // returns a row of non-modifiable elements

    /* Assignment operators */
    Vector3 &operator=(const Vector3 &v);                   // vector1  = vector2
    Vector3 &operator+=(const Vector3 &v);                  // vector1 += vector2
    Vector3 &operator-=(const Vector3 &v);                  // vector1 -= vector2
    Vector3 &CWiseMulEq(const Vector3 &v);                  // vector1 *= vector2, vector1(i) *= vector2(i)
    Vector3 &CWiseDivEq(const Vector3 &v);                  // vector1 /= vector2, vector1(i) /= vector2(i)
    Vector3 &operator=(const Vector<t_row, t_type> &v);     // vector1  = vector2
    Vector3 &operator+=(const Vector<t_row, t_type> &v);    // vector1 += vector2
    Vector3 &operator-=(const Vector<t_row, t_type> &v);    // vector1 -= vector2
    Vector3 &CWiseMulEq(const Vector<t_row, t_type> &v);    // vector1 *= vector2, vector1(i) *= vector2(i)
    Vector3 &CWiseDivEq(const Vector<t_row, t_type> &v);    // vector1 /= vector2, vector1(i) /= vector2(i)
    Vector3 &operator=(const Vector<0, t_type> &v);         // vector1  = vector2
    Vector3 &operator+=(const Vector<0, t_type> &v);        // vector1 += vector2
    Vector3 &operator-=(const Vector<0, t_type> &v);        // vector1 -= vector2
    Vector3 &CWiseMulEq(const Vector<0, t_type> &v);        // vector1 *= vector2, vector1(i) *= vector2(i)
    Vector3 &CWiseDivEq(const Vector<0, t_type> &v);        // vector1 /= vector2, vector1(i) /= vector2(i)
    Vector3 &operator=(const Matrix<t_row, 1, t_type> &v);  // vector1  = vector2
    Vector3 &operator+=(const Matrix<t_row, 1, t_type> &v); // vector1 += vector2
    Vector3 &operator-=(const Matrix<t_row, 1, t_type> &v); // vector1 -= vector2
    Vector3 &CWiseMulEq(const Matrix<t_row, 1, t_type> &v); // vector1 *= vector2, vector1(i) *= vector2(i)
    Vector3 &CWiseDivEq(const Matrix<t_row, 1, t_type> &v); // vector1 /= vector2, vector1(i) /= vector2(i)
    Vector3 &operator=(const Matrix<0, 0, t_type> &v);      // vector1  = vector2
    Vector3 &operator+=(const Matrix<0, 0, t_type> &v);     // vector1 += vector2
    Vector3 &operator-=(const Matrix<0, 0, t_type> &v);     // vector1 -= vector2
    Vector3 &CWiseMulEq(const Matrix<0, 0, t_type> &v);     // vector1 *= vector2, vector1(i) *= vector2(i)
    Vector3 &CWiseDivEq(const Matrix<0, 0, t_type> &v);     // vector1 /= vector2, vector1(i) /= vector2(i)
    Vector3 &operator=(const t_type s);                     // vector1  = scalar, all elements set scalar
    Vector3 &operator+=(const t_type s);                    // vector1 += scalar, vector1(i) += scalar
    Vector3 &operator-=(const t_type s);                    // vector1 -= scalar, vector1(i) -= scalar
    Vector3 &operator*=(const t_type s);                    // vector1 *= scalar
    Vector3 &operator/=(const t_type s);                    // vector1 /= scalar
    Vector3 &operator&=(const Vector3 &v);                  // vector1 x vector2 cross product
    Vector3 &operator&=(const Vector<t_row, t_type> &v);    // vector1 x vector2 cross product
    Vector3 &operator&=(const Vector<0, t_type> &v);        // vector1 x vector2 cross product
    Vector3 &operator&=(const Matrix<t_row, 1, t_type> &v); // vector1 x vector2 cross product
    Vector3 &operator&=(const Matrix<0, 0, t_type> &v);     // vector1 x vector2 cross product
    CommaInit<t_row, t_type> operator<<(const t_type s);    // Init first matrix elements

    /* Arithmetic operators */
    Vector3 operator-() const;                                  // minus sign
    Vector3 operator+(const Vector3 &v) const;                  // vector + vector
    Vector3 operator-(const Vector3 &v) const;                  // vector - vector
    Vector3 CWiseMul(const Vector3 &v) const;                   // vector * vector, vector(i) = vector1(i) * vector2(i)
    Vector3 CWiseDiv(const Vector3 &v) const;                   // vector / vector, vector(i) = vector1(i) / vector2(i)
    Vector3 operator+(const Vector<t_row, t_type> &v) const;    // vector + vector
    Vector3 operator-(const Vector<t_row, t_type> &v) const;    // vector - vector
    Vector3 CWiseMul(const Vector<t_row, t_type> &v) const;     // vector * vector, vector(i) = vector1(i) * vector2(i)
    Vector3 CWiseDiv(const Vector<t_row, t_type> &v) const;     // vector / vector, vector(i) = vector1(i) / vector2(i)
    Vector3 operator+(const Vector<0, t_type> &v) const;        // vector + vector
    Vector3 operator-(const Vector<0, t_type> &v) const;        // vector - vector
    Vector3 CWiseMul(const Vector<0, t_type> &v) const;         // vector * vector, vector(i) = vector1(i) * vector2(i)
    Vector3 CWiseDiv(const Vector<0, t_type> &v) const;         // vector / vector, vector(i) = vector1(i) / vector2(i)
    Vector3 operator+(const Matrix<t_row, 1, t_type> &v) const; // vector + vector
    Vector3 operator-(const Matrix<t_row, 1, t_type> &v) const; // vector - vector
    Vector3 CWiseMul(const Matrix<t_row, 1, t_type> &v) const;  // vector * vector, vector(i) = vector1(i) * vector2(i)
    Vector3 CWiseDiv(const Matrix<t_row, 1, t_type> &v) const;  // vector / vector, vector(i) = vector1(i) / vector2(i)
    Vector3 operator+(const Matrix<0, 0, t_type> &v) const;     // vector + vector
    Vector3 operator-(const Matrix<0, 0, t_type> &v) const;     // vector - vector
    Vector3 CWiseMul(const Matrix<0, 0, t_type> &v) const;      // vector * vector, vector(i) = vector1(i) * vector2(i)
    Vector3 CWiseDiv(const Matrix<0, 0, t_type> &v) const;      // vector / vector, vector(i) = vector1(i) / vector2(i)
    Vector3 operator+(const t_type s) const;                    // vector + scalar, vector(i) = vector1(i) + scalar
    Vector3 operator-(const t_type s) const;                    // vector - scalar, vector(i) = vector1(i) - scalar
    Vector3 operator*(const t_type s) const;                    // vector * scalar
    Vector3 operator/(const t_type s) const;                    // vector / scalar
    template <uint16_t col>
    Matrix<t_row, col, t_type> operator*(const Matrix<1, col, t_type> &m) const;             // vector(3x1) * matrix(1xcol) = matrix(3xcol)
    Matrix<0, 0, t_type> operator*(const Matrix<0, 0, t_type> &m) const;                     // vector(3x1) * matrix(1xn) = matrix(3xn)
    Vector3 operator&(const Vector3 &v) const;                                               // vector x vector cross product
    Vector3 operator&(const Vector<t_row, t_type> &v) const;                                 // vector x vector cross product
    Vector3 operator&(const Vector<0, t_type> &v) const;                                     // vector x vector cross product
    Vector3 operator&(const Matrix<t_row, 1, t_type> &v) const;                              // vector x vector cross product
    Vector3 operator&(const Matrix<0, 0, t_type> &v) const;                                  // vector x vector cross product
    Matrix3<t_type, t_row, t_row> operator&(const Matrix3<t_type, t_row, t_row> &m) const;   // [v]x * RotMat, []x is skew-symmetric matrix
    Rotation<t_type, t_row, t_row> operator&(const Rotation<t_type, t_row, t_row> &m) const; // [v]x * RotMat, []x is skew-symmetric matrix

    Matrix<t_row, t_row, t_type> Outer(const Vector3 &v) const;        // vector1 * vector2, outer product
    Matrix<t_row, 4, t_type> Outer(const Vector4<t_type, 4> &v) const; // vector1 * vector2, outer product
    Matrix<t_row, 6, t_type> Outer(const Vector6<t_type, 6> &v) const; // vector1 * vector2, outer product
    Matrix<0, 0, t_type> Outer(const Vector<0, t_type> &v) const;      // vector1 * vector2 outer product
    Matrix<0, 0, t_type> Outer(const Matrix<0, 0, t_type> &v) const;   // vector1 * vector2 outer product
    template <uint16_t row>
    Matrix<t_row, row, t_type> Outer(const Vector<row, t_type> &v) const; // vector1 * vector2, outer product
    template <uint16_t row>
    Matrix<t_row, row, t_type> Outer(const Matrix<row, 1, t_type> &v) const; // vector1 * vector2 outer product

    t_type Inner(const Vector3 &v) const;                  // vector1 * vector2 inner(dot) product
    t_type Inner(const Vector<t_row, t_type> &v) const;    // vector1 * vector2 inner(dot) product
    t_type Inner(const Matrix<t_row, 1, t_type> &v) const; // vector1 * vector2 inner(dot) product
    t_type Inner(const Vector<0, t_type> &v) const;        // vector1 * vector2 inner(dot) product
    t_type Inner(const Matrix<0, 0, t_type> &v) const;     // vector1 * vector2 inner(dot) product

    /* Comparison operators */
    bool operator==(const Vector3 &v) const;                  // (true or false) vector1 == vector2
    bool operator!=(const Vector3 &v) const;                  // (true or false) vector1 != vector2
    bool operator==(const Vector<t_row, t_type> &v) const;    // (true or false) vector1 == vector2
    bool operator!=(const Vector<t_row, t_type> &v) const;    // (true or false) vector1 != vector2
    bool operator==(const Matrix<t_row, 1, t_type> &v) const; // (true or false) vector1 == vector2
    bool operator!=(const Matrix<t_row, 1, t_type> &v) const; // (true or false) vector1 != vector2
    bool operator==(const Vector<0, t_type> &v) const;        // (true or false) vector1 == vector2
    bool operator!=(const Vector<0, t_type> &v) const;        // (true or false) vector1 != vector2
    bool operator==(const Matrix<0, 0, t_type> &v) const;     // (true or false) vector1 == vector2
    bool operator!=(const Matrix<0, 0, t_type> &v) const;     // (true or false) vector1 != vector2

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row> friend class Vector3;
    template <uint16_t row, uint16_t col, typename type> friend class Matrix;
    template <uint16_t row, uint16_t col, typename type> friend class CscMatrix;
    template <typename type, uint16_t row, uint16_t col> friend class Matrix3;
    template <typename type, uint16_t row, uint16_t col> friend class Rotation;
    template <typename type, uint16_t row, uint16_t col> friend class Transform;

    template <uint16_t row, typename type> friend class Vector;
    template <typename type, uint16_t row> friend class Vector4;
    template <typename type, uint16_t row> friend class Vector6;
    template <typename type, uint16_t row> friend class Quaternion;

    /* Friend template function */
    template <typename type, uint16_t row>
    friend Vector3<type, row> operator+(const type s, const Vector3<type, row> &v); // scalar + vector, scalar + vector(i)
    template <typename type, uint16_t row>
    friend Vector3<type, row> operator-(const type s, const Vector3<type, row> &v); // scalar - vector, scalar - vector(i)
    template <typename type, uint16_t row>
    friend Vector3<type, row> operator*(const type s, const Vector3<type, row> &v); // scalar * vector, scalar * vector(i)
    template <typename type, uint16_t row>
    friend Vector3<type, row> operator/(const type s, const Vector3<type, row> &v); // scalar / vector, scalar / vector(i)
};

} // namespace Math
} // namespace dt

#include "dtVector3.tpp"

#endif // DTMATH_DTVECTOR3_H_
