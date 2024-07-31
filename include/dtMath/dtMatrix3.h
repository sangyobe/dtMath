/*!
\file       dtMatrix3.h
\brief      dtMath, 3x3 Matrix class, lighter and faster than general matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATRIX3_H_
#define DTMATH_DTMATRIX3_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
#include <stdint.h>
#include <stdio.h>
#include <string.h>
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
template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;
template <uint16_t t_row, uint16_t t_col, typename t_type> class NoPivLU;
template <uint16_t t_row, uint16_t t_col, typename t_type> class PartialPivLU;
template <uint16_t t_row, uint16_t t_col, typename t_type> class FullPivLU;
template <uint16_t t_row, uint16_t t_col, typename t_type> class LLT;
template <uint16_t t_row, uint16_t t_col, typename t_type> class LDLT;
template <uint16_t t_row, uint16_t t_col, typename t_type> class QR;
template <uint16_t t_row, uint16_t t_col, typename t_type> class SVD;
template <typename type, uint16_t row, uint16_t col> class Rotation;

template <typename t_type = float, uint16_t t_row = 3, uint16_t t_col = 3>
class Matrix3
{
private:
    t_type m_tolerance = std::numeric_limits<t_type>::epsilon();
    t_type m_elem[t_row * t_col];
    Matrix3(const t_type *element);

public:
    Matrix3();
    Matrix3(const t_type *element, const size_t n_byte);
    Matrix3(const char c, const t_type *element, const size_t n_byte);
    Matrix3(
        const t_type m00, const t_type m01, const t_type m02,
        const t_type m10, const t_type m11, const t_type m12,
        const t_type m20, const t_type m21, const t_type m22);
    Matrix3(const Matrix3 &m);
    Matrix3(const Rotation<t_type, t_row, t_col> &m);
    Matrix3(const Matrix<t_row, t_col, t_type> &m);
    Matrix3(const Matrix<0, 0, t_type> &m);
    ~Matrix3() {}

    void SetZero();
    void SetIdentity();
    void SetDiagonal(const t_type d1, const t_type d2, const t_type d3);
    void SetDiagonal(const t_type *element, const size_t n_byte);
    void SetDiagonal(const t_type value);
    void SetDiagonal(const Vector<t_row, t_type> &v);
    void SetDiagonal(const Vector<0, t_type> &v);
    void SetDiagonal(const Vector3<t_type, t_row> &v);
    void SetFill(const t_type value);
    void SetElement(const t_type *element, const size_t n_byte);
    void SetElement(
        const t_type m00, const t_type m01, const t_type m02,
        const t_type m10, const t_type m11, const t_type m12,
        const t_type m20, const t_type m21, const t_type m22);
    void SetElement(const Matrix3 &m);
    void SetElement(const Rotation<t_type, t_row, t_col> &m);
    void SetElement(const Matrix<t_row, t_col, t_type> &m);
    void SetElement(const Matrix<0, 0, t_type> &m);
    template <uint16_t col>
    void SetRowVec(const uint16_t idxRow, const Vector<col, t_type> &v);
    void SetRowVec(const uint16_t idxRow, const Vector<0, t_type> &v);
    void SetRowVec(const uint16_t idxRow, const Vector3<t_type, t_col> &v);
    void SetRowVec(const uint16_t idxRow, const t_type *v, const size_t n_byte);
    template <uint16_t row>
    void SetColVec(const uint16_t idxCol, const Vector<row, t_type> &v);
    void SetColVec(const uint16_t idxCol, const Vector<0, t_type> &v);
    void SetColVec(const uint16_t idxCol, const Vector3<t_type, t_row> &v);
    void SetColVec(const uint16_t idxCol, const t_type *v, const size_t n_byte);
    void SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2);
    void SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2);

    const t_type *const GetElementsAddr() const;
    Vector3<t_type, t_col> GetRowVec(const uint16_t idxRow) const;
    Vector3<t_type, t_row> GetColVec(const uint16_t idxCol) const;
    int8_t GetRowVec(const uint16_t idxRow, Vector3<t_type, t_col> &v) const;
    int8_t GetColVec(const uint16_t idxCol, Vector3<t_type, t_row> &v) const;
    Matrix3 Transpose() const;

    t_type Trace() const;
    t_type GetNorm() const;   // Frobenius Norm (Euclidean norm)
    t_type GetSqNorm() const; // Squared Frobenius Norm (Euclidean norm)

    NoPivLU<t_row, t_col, t_type> GetNoPivLU() const;
    PartialPivLU<t_row, t_col, t_type> GetPartialPivLU() const;
    FullPivLU<t_row, t_col, t_type> GetFullPivLU() const;
    LLT<t_row, t_col, t_type> GetLLT() const;
    LDLT<t_row, t_col, t_type> GetLDLT() const;
    QR<t_row, t_col, t_type> GetQR() const;
    SVD<t_row, t_col, t_type> GetSVD() const;

    Matrix3 Inv(int8_t *isOk = nullptr) const;
    Matrix3 PInv(int8_t *isOk = nullptr, t_type tolerance = 0) const;

    /* Member access operators */
    t_type &operator()(uint16_t irow, uint16_t icol);             // returns a row of modifiable elements
    const t_type &operator()(uint16_t irow, uint16_t icol) const; // returns a row of non-modifiable elements

    /* Assignment operators */
    Matrix3 &operator=(const Matrix3 &m);                         // matrix  = matrix
    Matrix3 &operator+=(const Matrix3 &m);                        // matrix += matrix
    Matrix3 &operator-=(const Matrix3 &m);                        // matrix -= matrix
    Matrix3 &operator=(const Rotation<t_type, t_row, t_col> &m);  // matrix  = matrix
    Matrix3 &operator+=(const Rotation<t_type, t_row, t_col> &m); // matrix += matrix
    Matrix3 &operator-=(const Rotation<t_type, t_row, t_col> &m); // matrix -= matrix
    Matrix3 &operator=(const Matrix<t_row, t_col, t_type> &m);    // matrix  = matrix
    Matrix3 &operator+=(const Matrix<t_row, t_col, t_type> &m);   // matrix += matrix
    Matrix3 &operator-=(const Matrix<t_row, t_col, t_type> &m);   // matrix -= matrix
    Matrix3 &operator=(const Matrix<0, 0, t_type> &m);            // matrix  = matrix
    Matrix3 &operator+=(const Matrix<0, 0, t_type> &m);           // matrix += matrix
    Matrix3 &operator-=(const Matrix<0, 0, t_type> &m);           // matrix -= matrix
    Matrix3 &operator=(const t_type s);                           // matrix  = scalar, all elements set scalar
    Matrix3 &operator+=(const t_type s);                          // matrix += scalar, matrix(i) += scalar
    Matrix3 &operator-=(const t_type s);                          // matrix -= scalar, matrix(i) -= scalar
    Matrix3 &operator*=(const t_type s);                          // matrix *= scalar
    Matrix3 &operator/=(const t_type s);                          // matrix /= scalar
    CommaInit<t_row * t_col, t_type> operator<<(const t_type s);  // Init first matrix elements

    /* Arithmetic operators */
    Matrix3 operator-() const;                                        // minus sign
    Matrix3 operator+(const Matrix3 &m) const;                        // matrix + matrix
    Matrix3 operator-(const Matrix3 &m) const;                        // matrix - matrix
    Matrix3 operator+(const Rotation<t_type, t_row, t_col> &m) const; // matrix + matrix
    Matrix3 operator-(const Rotation<t_type, t_row, t_col> &m) const; // matrix - matrix
    Matrix3 operator+(const Matrix<t_row, t_col, t_type> &m) const;   // matrix + matrix
    Matrix3 operator-(const Matrix<t_row, t_col, t_type> &m) const;   // matrix - matrix
    Matrix3 operator+(const Matrix<0, 0, t_type> &m) const;           // matrix + matrix
    Matrix3 operator-(const Matrix<0, 0, t_type> &m) const;           // matrix - matrix
    Matrix3 operator+(const t_type s) const;                          // matrix + scalar, matrix(i) + scalar
    Matrix3 operator-(const t_type s) const;                          // matrix - scalar, matrix(i) - scalar
    Matrix3 operator*(const t_type s) const;                          // matrix * scalar
    Matrix3 operator/(const t_type s) const;                          // matrix / scalar

    template <uint16_t col>
    Matrix<t_row, col, t_type> operator*(const Matrix<t_col, col, t_type> &m) const; // matrix * matrix
    Matrix<0, 0, t_type> operator*(const Matrix<0, 0, t_type> &m) const;             // matrix * matrix
    Matrix3 operator*(const Matrix3 &m) const;                                       // matrix * matrix
    Matrix3 operator*(const Rotation<t_type, t_row, t_col> &m) const;                // matrix * rotation
    Vector<t_row, t_type> operator*(const Vector<t_col, t_type> &v) const;           // matrix * vector
    Vector<t_row, t_type> operator*(const Vector<0, t_type> &v) const;               // matrix * vector
    Vector3<t_type, t_row> operator*(const Vector3<t_type, t_col> &v) const;         // matrix * vector
    Matrix3 operator&(const Vector<t_col, t_type> &v) const;                         // matrix * [v]x, []x is skew-symmetric matrix
    Matrix3 operator&(const Vector<0, t_type> &v) const;                             // matrix * [v]x, []x is skew-symmetric matrix
    Matrix3 operator&(const Vector3<t_type, t_col> &v) const;                        // matrix * [v]x, []x is skew-symmetric matrix

    /* Comparison operators */
    bool operator==(const Matrix3 &m) const;                        // (true or false) matrix == matrix
    bool operator!=(const Matrix3 &m) const;                        // (true or false) matrix != matrix
    bool operator==(const Rotation<t_type, t_row, t_col> &m) const; // (true or false) matrix == matrix
    bool operator!=(const Rotation<t_type, t_row, t_col> &m) const; // (true or false) matrix != matrix
    bool operator==(const Matrix<t_row, t_col, t_type> &m) const;   // (true or false) matrix == matrix
    bool operator!=(const Matrix<t_row, t_col, t_type> &m) const;   // (true or false) matrix != matrix
    bool operator==(const Matrix<0, 0, t_type> &m) const;           // (true or false) matrix == matrix
    bool operator!=(const Matrix<0, 0, t_type> &m) const;           // (true or false) matrix != matrix

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row, uint16_t col> friend class Matrix3;
    template <uint16_t row, uint16_t col, typename type> friend class CscMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class Matrix;
    template <typename type, uint16_t row, uint16_t col> friend class Rotation;

    template <uint16_t row, uint16_t col, typename type> friend class LowerMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class UpperMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class NoPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class PartialPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class FullPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class LLT;
    template <uint16_t row, uint16_t col, typename type> friend class LDLT;
    template <uint16_t row, uint16_t col, typename type> friend class QR;
    template <uint16_t row, uint16_t col, typename type> friend class SVD;

    /* Friend template function */
    template <typename type, uint16_t row, uint16_t col>
    friend Matrix3<type, row, col> operator*(const type s, const Matrix3<type, row, col> &m); // scalar * matrix
};

} // namespace Math
} // namespace dt

#include "dtMatrix3.tpp"

#endif // DTMATH_DTMATRIX3_H_
