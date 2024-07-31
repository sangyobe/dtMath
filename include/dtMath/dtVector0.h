/*!
\file       dtVector.h
\brief      dtMath, Dynamic Memory Allocation General Vector(m x 1) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 04. 18
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTVECTOR0_H_
#define DTMATH_DTVECTOR0_H_

#include "dtDefine.h"
#include "dtVector.h"

#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
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

template <typename t_type>
class Vector<0, t_type>
{
private:
    t_type m_tolerance = std::numeric_limits<t_type>::epsilon();
    uint16_t m_row;
    t_type *m_elem;

    inline void CrossProduct(const t_type *v);

public:
    Vector();
    Vector(const uint16_t row);
    Vector(const uint16_t row, const t_type *element);
    Vector(const Vector &v);
    Vector(const Matrix<0, 0, t_type> &v);
    template <uint16_t row>
    Vector(const Vector<row, t_type> &v);
    template <uint16_t row>
    Vector(const Matrix<row, 1, t_type> &v);
    ~Vector();

    void NewSize(const uint16_t row);
    void NewSize(const uint16_t row, const t_type *element);
    void NewSize(const Vector &v);
    void NewSize(const Matrix<0, 0, t_type> &v);
    template <uint16_t row>
    void NewSize(const Vector<row, t_type> &v);
    template <uint16_t row>
    void NewSize(const Matrix<row, 1, t_type> &v);

    void ReSize(const uint16_t row);
    void ReSize(const uint16_t row, const t_type *element);
    void ReSize(const Vector &v);
    void ReSize(const Matrix<0, 0, t_type> &v);
    template <uint16_t row>
    void ReSize(const Vector<row, t_type> &v);
    template <uint16_t row>
    void ReSize(const Matrix<row, 1, t_type> &v);

    void Release(); // Free Memory(m_elem)

    void SetZero();
    void SetFill(const t_type value);
    void SetElement(const t_type *element, const size_t n_byte);
    void SetElement(const t_type *element, const uint16_t row);
    void SetBlock(const uint16_t idxRow, const Vector &v);
    void SetBlock(const uint16_t idxRow, const Vector3<t_type, 3> &v);
    void SetBlock(const uint16_t idxRow, const Vector4<t_type, 4> &v);
    void SetBlock(const uint16_t idxRow, const Vector6<t_type, 6> &v);
    void SetBlock(const uint16_t idxRow, const Matrix<0, 0, t_type> &v);
    template <uint16_t row>
    void SetBlock(const uint16_t idxRow, const Vector<row, t_type> &v);
    template <uint16_t row>
    void SetBlock(const uint16_t idxRow, const Matrix<row, 1, t_type> &v);
    void SetSwap(const uint16_t i, const uint16_t j);
    void SetNormalize();

    const t_type *const GetElementsAddr() const;
    uint16_t GetDim() const { return m_row; }
    Vector<0, t_type> GetBlock(const uint16_t idx, const uint16_t row);
    Vector3<t_type, 3> GetBlockVec3(const uint16_t idx);
    Vector4<t_type, 4> GetBlockVec4(const uint16_t idx);
    Vector6<t_type, 6> GetBlockVec6(const uint16_t idx);
    template <uint16_t row>
    Vector<row, t_type> GetBlock(const uint16_t idx);
    int8_t GetBlock(const uint16_t idx, Vector<0, t_type> &v);
    int8_t GetBlockVec3(const uint16_t idx, Vector3<t_type, 3> &v);
    int8_t GetBlockVec4(const uint16_t idx, Vector4<t_type, 4> &v);
    int8_t GetBlockVec6(const uint16_t idx, Vector6<t_type, 6> &v);
    template <uint16_t row>
    int8_t GetBlock(const uint16_t idx, Vector<row, t_type> &v);
    t_type GetNorm() const;              // Frobenius Norm (Euclidean norm, L2 Norm)
    t_type GetSqNorm() const;            // Squared Frobenius Norm (Euclidean norm, Squared L2 Norm)
    t_type GetLpNorm(const int p) const; // Generalized Norm (Lp Norm)
    t_type GetSum() const;
    Vector GetNormalized() const;
    Matrix<3, 3, t_type> GetSkew() const;
    Matrix<0, 0, t_type> Transpose() const;
    void Transpose(Matrix<0, 0, t_type> &m) const;
    template <uint16_t col>
    void Transpose(Matrix<1, col, t_type> &m) const;

    /* Member access operators */
    t_type &operator()(uint16_t irow);             // returns a row of modifiable elements
    const t_type &operator()(uint16_t irow) const; // returns a row of non-modifiable elements

    /* Assignment operators */
    Vector &operator=(const Vector &v);                // vector1  = vector2
    Vector &operator+=(const Vector &v);               // vector1 += vector2
    Vector &operator-=(const Vector &v);               // vector1 -= vector2
    Vector &operator*=(const Vector &v);               // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &operator/=(const Vector &v);               // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const Vector3<t_type, 3> &v);    // vector1  = vector2
    Vector &operator+=(const Vector3<t_type, 3> &v);   // vector1 += vector2
    Vector &operator-=(const Vector3<t_type, 3> &v);   // vector1 -= vector2
    Vector &operator*=(const Vector3<t_type, 3> &v);   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &operator/=(const Vector3<t_type, 3> &v);   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const Vector4<t_type, 4> &v);    // vector1  = vector2
    Vector &operator+=(const Vector4<t_type, 4> &v);   // vector1 += vector2
    Vector &operator-=(const Vector4<t_type, 4> &v);   // vector1 -= vector2
    Vector &operator*=(const Vector4<t_type, 4> &v);   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &operator/=(const Vector4<t_type, 4> &v);   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const Vector6<t_type, 6> &v);    // vector1  = vector2
    Vector &operator+=(const Vector6<t_type, 6> &v);   // vector1 += vector2
    Vector &operator-=(const Vector6<t_type, 6> &v);   // vector1 -= vector2
    Vector &operator*=(const Vector6<t_type, 6> &v);   // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &operator/=(const Vector6<t_type, 6> &v);   // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const Matrix<0, 0, t_type> &v);  // vector1  = vector2
    Vector &operator+=(const Matrix<0, 0, t_type> &v); // vector1 += vector2
    Vector &operator-=(const Matrix<0, 0, t_type> &v); // vector1 -= vector2
    Vector &operator*=(const Matrix<0, 0, t_type> &v); // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    Vector &operator/=(const Matrix<0, 0, t_type> &v); // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    Vector &operator=(const t_type s);                 // vector1  = scalar, all elements set scalar
    Vector &operator+=(const t_type s);                // vector1 += scalar, vector1(i) = vector1(i) + scalar
    Vector &operator-=(const t_type s);                // vector1 -= scalar, vector1(i) = vector1(i) - scalar
    Vector &operator*=(const t_type s);                // vector1 *= scalar
    Vector &operator/=(const t_type s);                // vector1 /= scalar
    template <uint16_t row>
    Vector &operator=(const Vector<row, t_type> &v); // vector1  = vector2
    template <uint16_t row>
    Vector &operator+=(const Vector<row, t_type> &v); // vector1 += vector2
    template <uint16_t row>
    Vector &operator-=(const Vector<row, t_type> &v); // vector1 -= vector2
    template <uint16_t row>
    Vector &operator*=(const Vector<row, t_type> &v); // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    template <uint16_t row>
    Vector &operator/=(const Vector<row, t_type> &v); // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)
    template <uint16_t row>
    Vector &operator=(const Matrix<row, 1, t_type> &v); // vector1  = vector2
    template <uint16_t row>
    Vector &operator+=(const Matrix<row, 1, t_type> &v); // vector1 += vector2
    template <uint16_t row>
    Vector &operator-=(const Matrix<row, 1, t_type> &v); // vector1 -= vector2
    template <uint16_t row>
    Vector &operator*=(const Matrix<row, 1, t_type> &v); // vector1 *= vector2, vector1(i) = vector1(i) * vector2(i)
    template <uint16_t row>
    Vector &operator/=(const Matrix<row, 1, t_type> &v); // vector1 /= vector2, vector1(i) = vector1(i) / vector2(i)

    Vector &operator&=(const Vector &v);               // vector1 x vector2 cross product
    Vector &operator&=(const Vector3<t_type, 3> &v);   // vector1 x vector2 cross product
    Vector &operator&=(const Matrix<0, 0, t_type> &v); // vector1 x vector2 cross product
    template <uint16_t row>
    Vector &operator&=(const Vector<row, t_type> &v); // vector1 x vector2 cross product
    template <uint16_t row>
    Vector &operator&=(const Matrix<row, 1, t_type> &v); // vector1 x vector2 cross product

    CommaInit<0, t_type> operator<<(const t_type s); // Init first matrix elements

    /* Arithmetic operators */
    Vector operator-() const;                              // minus sign
    Vector operator+(const Vector &v) const;               // vector + vector
    Vector operator-(const Vector &v) const;               // vector - vector
    Vector operator*(const Vector &v) const;               // vector * vector = m_elem(i) * v.m_elem(i)
    Vector operator/(const Vector &v) const;               // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const Vector3<t_type, 3> &v) const;   // vector + vector
    Vector operator-(const Vector3<t_type, 3> &v) const;   // vector - vector
    Vector operator*(const Vector3<t_type, 3> &v) const;   // vector * vector = m_elem(i) * v.m_elem(i)
    Vector operator/(const Vector3<t_type, 3> &v) const;   // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const Vector4<t_type, 4> &v) const;   // vector + vector
    Vector operator-(const Vector4<t_type, 4> &v) const;   // vector - vector
    Vector operator*(const Vector4<t_type, 4> &v) const;   // vector * vector = m_elem(i) * v.m_elem(i)
    Vector operator/(const Vector4<t_type, 4> &v) const;   // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const Vector6<t_type, 6> &v) const;   // vector + vector
    Vector operator-(const Vector6<t_type, 6> &v) const;   // vector - vector
    Vector operator*(const Vector6<t_type, 6> &v) const;   // vector * vector = m_elem(i) * v.m_elem(i)
    Vector operator/(const Vector6<t_type, 6> &v) const;   // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const Matrix<0, 0, t_type> &v) const; // vector + vector
    Vector operator-(const Matrix<0, 0, t_type> &v) const; // vector - vector
    Vector operator*(const Matrix<0, 0, t_type> &v) const; // vector * vector = m_elem(i) * v.m_elem(i)
    Vector operator/(const Matrix<0, 0, t_type> &v) const; // vector / vector = m_elem(i) / v.m_elem(i)
    Vector operator+(const t_type s) const;                // vector + scalar = m_elem(i) + scalar
    Vector operator-(const t_type s) const;                // vector - scalar = m_elem(i) - scalar
    Vector operator*(const t_type s) const;                // vector * scalar
    Vector operator/(const t_type s) const;                // vector / scalar
    template <uint16_t row>
    Vector operator+(const Vector<row, t_type> &v) const; // vector + vector
    template <uint16_t row>
    Vector operator-(const Vector<row, t_type> &v) const; // vector - vector
    template <uint16_t row>
    Vector operator*(const Vector<row, t_type> &v) const; // vector * vector = m_elem(i) * v.m_elem(i)
    template <uint16_t row>
    Vector operator/(const Vector<row, t_type> &v) const; // vector / vector = m_elem(i) / v.m_elem(i)
    template <uint16_t row>
    Vector operator+(const Matrix<row, 1, t_type> &v) const; // vector + vector
    template <uint16_t row>
    Vector operator-(const Matrix<row, 1, t_type> &v) const; // vector - vector
    template <uint16_t row>
    Vector operator*(const Matrix<row, 1, t_type> &v) const; // vector * vector = m_elem(i) * v.m_elem(i)
    template <uint16_t row>
    Vector operator/(const Matrix<row, 1, t_type> &v) const; // vector / vector = m_elem(i) / v.m_elem(i)

    Vector operator&(const Vector &v) const;               // vector x vector cross product
    Vector operator&(const Vector3<t_type, 3> &v) const;   // vector x vector cross product
    Vector operator&(const Matrix<0, 0, t_type> &v) const; // vector x vector cross product
    template <uint16_t row>
    Vector operator&(const Vector<row, t_type> &v) const; // vector x vector cross product
    template <uint16_t row>
    Vector operator&(const Matrix<row, 1, t_type> &v) const; // vector x vector cross product

    Matrix3<t_type, 3, 3> operator&(const Matrix3<t_type, 3, 3> &m) const;   // [v]x * RotMat, []x is skew-symmetric matrix
    Rotation<t_type, 3, 3> operator&(const Rotation<t_type, 3, 3> &m) const; // [v]x * RotMat, []x is skew-symmetric matrix

    Matrix<0, 0, t_type> outer(const Vector &v) const;               // vector1 * vector2 outer product
    Matrix<0, 0, t_type> outer(const Matrix<0, 0, t_type> &v) const; // vector1 * matrix(1xcol) outer product
    template <uint16_t row>
    Matrix<0, 0, t_type> outer(const Vector<row, t_type> &v) const; // vector1 * vector2 outer product
    template <uint16_t row>
    Matrix<0, 0, t_type> outer(const Matrix<row, 1, t_type> &v) const; // vector1 * matrix(1xcol) outer product

    t_type inner(const Vector &v) const;               // vector1 * vector2 inner(dot) product
    t_type inner(const Vector3<t_type, 3> &v) const;   // vector1 * vector2 inner(dot) product
    t_type inner(const Vector4<t_type, 4> &v) const;   // vector1 * vector2 inner(dot) product
    t_type inner(const Vector6<t_type, 6> &v) const;   // vector1 * vector2 inner(dot) product
    t_type inner(const Matrix<0, 0, t_type> &v) const; // vector1 * vector2 inner(dot) product
    template <uint16_t row>
    t_type inner(const Vector<row, t_type> &v) const; // vector1 * vector2 inner(dot) product
    template <uint16_t row>
    t_type inner(const Matrix<row, 1, t_type> &v) const; // vector1 * vector2 inner(dot) product

    /* Comparison operators */
    bool operator==(const Vector &v) const;               // (true or false) vector1 == vector2
    bool operator!=(const Vector &v) const;               // (true or false) vector1 != vector2
    bool operator==(const Vector3<t_type, 3> &v) const;   // (true or false) vector1 == vector2
    bool operator!=(const Vector3<t_type, 3> &v) const;   // (true or false) vector1 != vector2
    bool operator==(const Vector4<t_type, 4> &v) const;   // (true or false) vector1 == vector2
    bool operator!=(const Vector4<t_type, 4> &v) const;   // (true or false) vector1 != vector2
    bool operator==(const Vector6<t_type, 6> &v) const;   // (true or false) vector1 == vector2
    bool operator!=(const Vector6<t_type, 6> &v) const;   // (true or false) vector1 != vector2
    bool operator==(const Matrix<0, 0, t_type> &v) const; // (true or false) vector1 == vector2
    bool operator!=(const Matrix<0, 0, t_type> &v) const; // (true or false) vector1 != vector2
    template <uint16_t row>
    bool operator==(const Vector<row, t_type> &v) const; // (true or false) vector1 == vector2
    template <uint16_t row>
    bool operator!=(const Vector<row, t_type> &v) const; // (true or false) vector1 != vector2
    template <uint16_t row>
    bool operator==(const Matrix<row, 1, t_type> &v) const; // (true or false) vector1 == vector2
    template <uint16_t row>
    bool operator!=(const Matrix<row, 1, t_type> &v) const; // (true or false) vector1 != vector2

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

    /* Friend template function */
    template <typename type>
    friend Vector<0, type> operator+(const type s, const Vector<0, type> &v); // scalar + vector = scalar + m_elem(i)
    template <typename type>
    friend Vector<0, type> operator-(const type s, const Vector<0, type> &v); // scalar - vector = scalar - m_elem(i)
    template <typename type>
    friend Vector<0, type> operator*(const type s, const Vector<0, type> &v); // scalar * vector = scalar * m_elem(i)
    template <typename type>
    friend Vector<0, type> operator/(const type s, const Vector<0, type> &v); // scalar / vector = scalar / m_elem(i)
};

} // namespace Math
} // namespace dt

#include "dtVector0.tpp"

#endif // DTMATH_DTVECTOR_H_
