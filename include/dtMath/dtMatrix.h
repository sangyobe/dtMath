/*!
\file       dtMatrix.h
\brief      dtMath, General Matrix(m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATRIX_H_
#define DTMATH_DTMATRIX_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
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
template <typename t_type, uint16_t t_row> class Vector4;
template <typename t_type, uint16_t t_row> class Quaternion;
template <typename t_type, uint16_t t_row> class Vector6;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Matrix3;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Rotation;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Transform;
template <uint16_t t_row, uint16_t t_col, typename t_type> class NoPivLU;
template <uint16_t t_row, uint16_t t_col, typename t_type> class PartialPivLU;
template <uint16_t t_row, uint16_t t_col, typename t_type> class FullPivLU;
template <uint16_t t_row, uint16_t t_col, typename t_type> class LLT;
template <uint16_t t_row, uint16_t t_col, typename t_type> class LDLT;
template <uint16_t t_row, uint16_t t_col, typename t_type> class QR;
template <uint16_t t_row, uint16_t t_col, typename t_type> class SVD;
template <uint16_t t_row, uint16_t t_col, typename t_type> class CscMatrix;

template <uint16_t t_row, uint16_t t_col, typename t_type = float>
class Matrix
{
private:
    t_type m_tolerance = std::numeric_limits<t_type>::epsilon();
    t_type m_elem[t_row * t_col];
    Matrix(const t_type *element);

public:
    Matrix();
    Matrix(const t_type *element, const size_t n_byte);
    Matrix(const char c, const t_type *element, const size_t n_byte); // JhJo(230424) : modified to handle default n_byte as matrix size.
    Matrix(const Matrix &m);
    ~Matrix() {}

    void SetZero();
    void SetIdentity();
    void SetDiagonal(const t_type *element, const size_t n_byte);
    void SetDiagonal(const t_type value);
    template <uint16_t row>
    void SetDiagonal(const Vector<row, t_type> &v);
    void SetDiagonal(const Vector<0, t_type> &v);
    void SetDiagonal(const Vector3<t_type, 3> &v);
    void SetDiagonal(const Vector4<t_type, 4> &v);
    void SetDiagonal(const Vector6<t_type, 6> &v);
    void SetDiagonal(const Quaternion<t_type, 4> &v);
    void SetFill(const t_type value);
    void SetFillRow(const uint16_t idxRow, const t_type value);
    void SetFillCol(const uint16_t idxCol, const t_type value);
    void SetElement(const t_type *element, const size_t n_byte);
    template <uint16_t row, uint16_t col>
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix<row, col, t_type> &m, const uint16_t jdxRow = 0, const uint16_t jdxCol = 0, const int16_t rSize = row, const int16_t cSize = col);
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix<0, 0, t_type> &m, const uint16_t jdxRow = 0, const uint16_t jdxCol = 0, int16_t rSize = -1, int16_t cSize = -1);
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix3<t_type, 3, 3> &m);
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Rotation<t_type, 3, 3> &m);
    template <uint16_t col>
    void SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Vector<col, t_type> &v, const uint16_t jdx = 0, const int16_t size = col);
    void SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Vector<0, t_type> &v, const uint16_t jdx = 0, int16_t size = -1);
    void SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Vector3<t_type, 3> &v, const uint16_t jdx = 0, const int16_t size = 3);
    void SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Vector4<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4);
    void SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Quaternion<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4);
    void SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Vector6<t_type, 6> &v, const uint16_t jdx = 0, const int16_t size = 6);
    void SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const t_type *v, const size_t n_byte);
    template <uint16_t row>
    void SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Vector<row, t_type> &v, const uint16_t jdx = 0, const int16_t size = row);
    void SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Vector<0, t_type> &v, const uint16_t jdx = 0, int16_t size = -1);
    void SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Vector3<t_type, 3> &v, const uint16_t jdx = 0, const int16_t size = 3);
    void SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Vector4<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4);
    void SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Quaternion<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4);
    void SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Vector6<t_type, 6> &v, const uint16_t jdx = 0, const int16_t size = 6);
    void SetColVec(const uint16_t idxRow, const uint16_t idxCol, const t_type *v, const size_t n_byte);
    void SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2);
    void SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2);

    const t_type *const GetElementsAddr() const;
    uint16_t GetRowSize() const { return t_row; } // size of row
    uint16_t GetColSize() const { return t_col; } // size of colum
    template <uint16_t row>
    void GetDiagonal(Vector<row, t_type> &v);
    void GetDiagonal(Vector<0, t_type> &v);
    void GetDiagonal(Vector3<t_type, 3> &v);
    void GetDiagonal(Vector4<t_type, 4> &v);
    void GetDiagonal(Vector6<t_type, 6> &v);
    void GetDiagonal(Quaternion<t_type, 4> &v);
    template <uint16_t row>
    Vector<row, t_type> GetDiagonal();
    Vector<0, t_type> GetDiagonalVec0();
    Vector3<t_type, 3> GetDiagonalVec3();
    Vector4<t_type, 4> GetDiagonalVec4();
    Vector6<t_type, 6> GetDiagonalVec6();
    Quaternion<t_type, 4> GetDiagonalQuat();
    template <uint16_t row, uint16_t col>
    Matrix<row, col, t_type> GetBlock(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdxRow = 0, const uint16_t jdxCol = 0, const int16_t rSize = row, const int16_t cSize = col);
    Matrix<0, 0, t_type> GetBlock(const uint16_t idxRow, const uint16_t idxCol, const uint16_t row, const uint16_t col, const uint16_t jdxRow = 0, const uint16_t jdxCol = 0, int16_t rSize = -1, int16_t cSize = -1);
    template <uint16_t row, uint16_t col>
    int8_t GetBlock(const uint16_t idxRow, const uint16_t idxCol, Matrix<row, col, t_type> &m, const uint16_t jdxRow = 0, const uint16_t jdxCol = 0, const int16_t rSize = row, const int16_t cSize = col);
    int8_t GetBlock(const uint16_t idxRow, const uint16_t idxCol, Matrix3<t_type, 3, 3> &m, const uint16_t jdxRow = 0, const uint16_t jdxCol = 0, const int16_t rSize = 3, const int16_t cSize = 3);
    int8_t GetBlock(const uint16_t idxRow, const uint16_t idxCol, Matrix<0, 0, t_type> &m, const uint16_t jdxRow = 0, const uint16_t jdxCol = 0, int16_t rSize = -1, int16_t cSize = -1);
    template <uint16_t col>
    int8_t GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Vector<col, t_type> &v, const uint16_t jdx = 0, const int16_t size = col) const;
    int8_t GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Vector<0, t_type> &v, const uint16_t jdx = 0, int16_t size = -1) const;
    int8_t GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Vector3<t_type, 3> &v, const uint16_t jdx = 0, const int16_t size = 3) const;
    int8_t GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Vector4<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4) const;
    int8_t GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Vector6<t_type, 6> &v, const uint16_t jdx = 0, const int16_t size = 6) const;
    int8_t GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Quaternion<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4) const;
    template <uint16_t row>
    int8_t GetColVec(const uint16_t idxRow, const uint16_t idxCol, Vector<row, t_type> &v, const uint16_t jdx = 0, const int16_t size = row) const;
    int8_t GetColVec(const uint16_t idxRow, const uint16_t idxCol, Vector<0, t_type> &v, const uint16_t jdx = 0, int16_t size = -1) const;
    int8_t GetColVec(const uint16_t idxRow, const uint16_t idxCol, Vector3<t_type, 3> &v, const uint16_t jdx = 0, const int16_t size = 3) const;
    int8_t GetColVec(const uint16_t idxRow, const uint16_t idxCol, Vector4<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4) const;
    int8_t GetColVec(const uint16_t idxRow, const uint16_t idxCol, Vector6<t_type, 6> &v, const uint16_t jdx = 0, const int16_t size = 6) const;
    int8_t GetColVec(const uint16_t idxRow, const uint16_t idxCol, Quaternion<t_type, 4> &v, const uint16_t jdx = 0, const int16_t size = 4) const;
    template <uint16_t col>
    Vector<col, t_type> GetRowVec(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx = 0, const int16_t size = col) const;
    Vector<0, t_type> GetRowVec0(const uint16_t idxRow, const uint16_t idxCol, const uint16_t col, const uint16_t jdx = 0, int16_t size = -1) const;
    Vector3<t_type, 3> GetRowVec3(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx = 0, const int16_t size = 3) const;
    Vector4<t_type, 4> GetRowVec4(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx = 0, const int16_t size = 4) const;
    Vector6<t_type, 6> GetRowVec6(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx = 0, const int16_t size = 6) const;
    Quaternion<t_type, 4> GetRowQuat(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx = 0, const int16_t size = 4) const;
    template <uint16_t row>
    Vector<row, t_type> GetColVec(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx = 0, const int16_t size = row) const;
    Vector<0, t_type> GetColVec0(const uint16_t idxRow, const uint16_t idxCol, const uint16_t row, const uint16_t jdx = 0, const int16_t size = -1) const;
    Vector3<t_type, 3> GetColVec3(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx = 0, const int16_t size = 3) const;
    Vector4<t_type, 4> GetColVec4(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx = 0, const int16_t size = 4) const;
    Vector6<t_type, 6> GetColVec6(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx = 0, const int16_t size = 6) const;
    Quaternion<t_type, 4> GetColQuat(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx = 0, const int16_t size = 4) const;

    CscMatrix<t_row, t_col, t_type> GetCscMat() const; // Compressed Sparse Column Matrix
    Matrix<t_col, t_row, t_type> Transpose() const;

    t_type Trace() const;
    t_type GetNorm() const;              // Frobenius Norm (Euclidean norm, L2 Norm)
    t_type GetSqNorm() const;            // Squared Frobenius Norm (Euclidean norm, Squared L2 Norm)
    t_type GetLpNorm(const int p) const; // Generalized Norm (Lp Norm)
    t_type Determinant() const;          // From LU Decomposition

    NoPivLU<t_row, t_col, t_type> GetNoPivLU() const;
    PartialPivLU<t_row, t_col, t_type> GetPartialPivLU() const;
    FullPivLU<t_row, t_col, t_type> GetFullPivLU() const;
    LLT<t_row, t_col, t_type> GetLLT() const;
    LDLT<t_row, t_col, t_type> GetLDLT() const;
    QR<t_row, t_col, t_type> GetQR() const;
    SVD<t_row, t_col, t_type> GetSVD() const;

    Matrix<t_row, t_col, t_type> Inv(int8_t *isOk = nullptr) const;  // Inverse Using LU Partial Pivoting
    Matrix<t_row, t_col, t_type> FInv(int8_t *isOk = nullptr) const; // Inverse Using LU Full Pivoting
    Matrix<t_col, t_row, t_type> PInv(int8_t *isOk = nullptr, t_type tolerance = std::numeric_limits<t_type>::epsilon()) const;

    t_type Max(uint16_t *irow = nullptr, uint16_t *icol = nullptr);
    t_type Min(uint16_t *irow = nullptr, uint16_t *icol = nullptr);
    t_type AbsMax(uint16_t *irow = nullptr, uint16_t *icol = nullptr);
    t_type AbsMin(uint16_t *irow = nullptr, uint16_t *icol = nullptr);

    /* Member access operators */
    t_type &operator()(uint16_t irow, uint16_t icol);             // returns a row of modifiable elements
    const t_type &operator()(uint16_t irow, uint16_t icol) const; // returns a row of non-modifiable elements

    /* Assignment operators */
    Matrix &operator=(const Matrix &m);                           // matrix  = matrix
    Matrix &operator+=(const Matrix &m);                          // matrix += matrix
    Matrix &operator-=(const Matrix &m);                          // matrix -= matrix
    Matrix &operator=(const Matrix<0, 0, t_type> &m);             // matrix  = matrix
    Matrix &operator+=(const Matrix<0, 0, t_type> &m);            // matrix += matrix
    Matrix &operator-=(const Matrix<0, 0, t_type> &m);            // matrix -= matrix
    Matrix &operator=(const Matrix3<t_type, t_row, t_col> &m);    // matrix  = matrix3
    Matrix &operator+=(const Matrix3<t_type, t_row, t_col> &m);   // matrix += matrix3
    Matrix &operator-=(const Matrix3<t_type, t_row, t_col> &m);   // matrix -= matrix3
    Matrix &operator=(const Rotation<t_type, t_row, t_col> &m);   // matrix  = RotMat
    Matrix &operator+=(const Rotation<t_type, t_row, t_col> &m);  // matrix += RotMat
    Matrix &operator-=(const Rotation<t_type, t_row, t_col> &m);  // matrix -= RotMat
    Matrix &operator=(const Transform<t_type, t_row, t_col> &m);  // matrix  = Transform
    Matrix &operator+=(const Transform<t_type, t_row, t_col> &m); // matrix += Transform
    Matrix &operator-=(const Transform<t_type, t_row, t_col> &m); // matrix -= Transform
    Matrix &operator=(const t_type s);                            // matrix  = scalar, all elements set scalar
    Matrix &operator+=(const t_type s);                           // matrix += scalar, matrix(i) += scalar
    Matrix &operator-=(const t_type s);                           // matrix -= scalar, matrix(i) -= scalar
    Matrix &operator*=(const t_type s);                           // matrix *= scalar
    Matrix &operator/=(const t_type s);                           // matrix /= scalar
    CommaInit<t_row * t_col, t_type> operator<<(const t_type s);  // Init first matrix elements

    /* Arithmetic operators */
    Matrix operator-() const;                                         // minus sign
    Matrix operator+(const Matrix &m) const;                          // matrix + matrix
    Matrix operator-(const Matrix &m) const;                          // matrix - matrix
    Matrix operator+(const Matrix<0, 0, t_type> &m) const;            // matrix + matrix
    Matrix operator-(const Matrix<0, 0, t_type> &m) const;            // matrix - matrix
    Matrix operator+(const Matrix3<t_type, t_row, t_col> &m) const;   // matrix + matrix3
    Matrix operator-(const Matrix3<t_type, t_row, t_col> &m) const;   // matrix - matrix3
    Matrix operator+(const Rotation<t_type, t_row, t_col> &m) const;  // matrix + RotMat
    Matrix operator-(const Rotation<t_type, t_row, t_col> &m) const;  // matrix - RotMat
    Matrix operator+(const Transform<t_type, t_row, t_col> &m) const; // matrix + Transform
    Matrix operator-(const Transform<t_type, t_row, t_col> &m) const; // matrix - Transform
    Matrix operator+(const t_type s) const;                           // matrix + scalar, matrix(i) + scalar
    Matrix operator-(const t_type s) const;                           // matrix - scalar, matrix(i) - scalar
    Matrix operator*(const t_type s) const;                           // matrix * scalar
    Matrix operator/(const t_type s) const;                           // matrix / scalar

    template <uint16_t col>
    Matrix<t_row, col, t_type> operator*(const Matrix<t_col, col, t_type> &m) const;        // matrix * matrix
    Matrix<0, 0, t_type> operator*(const Matrix<0, 0, t_type> &m) const;                    // matrix * matrix
    Matrix<t_row, t_col, t_type> operator*(const Matrix3<t_type, t_col, t_col> &m) const;   // matrix * matrix
    Matrix<t_row, t_col, t_type> operator*(const Rotation<t_type, t_col, t_col> &m) const;  // matrix * RotMat
    Matrix<t_row, t_col, t_type> operator*(const Transform<t_type, t_col, t_col> &m) const; // matrix * Transform
    Vector<t_row, t_type> operator*(const Vector<t_col, t_type> &v) const;                  // matrix * vector
    Vector<t_row, t_type> operator*(const Vector<0, t_type> &v) const;                      // matrix * vector
    Vector<t_row, t_type> operator*(const Vector3<t_type, t_col> &v) const;                 // matrix * vector3
    Vector<t_row, t_type> operator*(const Vector4<t_type, t_col> &v) const;                 // matrix * vector4
    Vector<t_row, t_type> operator*(const Vector6<t_type, t_col> &v) const;                 // matrix * vector6
    Vector<t_row, t_type> operator*(const Quaternion<t_type, t_col> &v) const;              // matrix * quaternion
    Matrix<t_row, t_col, t_type> operator&(const Vector<t_col, t_type> &v) const;           // matrix3 * [v]x, []x is skew-symmetric matrix
    Matrix<t_row, t_col, t_type> operator&(const Vector<0, t_type> &v) const;               // matrix3 * [v]x, []x is skew-symmetric matrix
    Matrix<t_row, t_col, t_type> operator&(const Vector3<t_type, t_col> &v) const;          // matrix3 * [v]x, []x is skew-symmetric matrix

    /* Comparison operators */
    bool operator==(const Matrix &m) const;                          // (true or false) matrix1 == matrix
    bool operator!=(const Matrix &m) const;                          // (true or false) matrix1 != matrix
    bool operator==(const Matrix<0, 0, t_type> &m) const;            // (true or false) matrix1 == matrix
    bool operator!=(const Matrix<0, 0, t_type> &m) const;            // (true or false) matrix1 != matrix
    bool operator==(const Matrix3<t_type, t_row, t_col> &m) const;   // (true or false) matrix == matrix3
    bool operator!=(const Matrix3<t_type, t_row, t_col> &m) const;   // (true or false) matrix != matrix3
    bool operator==(const Rotation<t_type, t_row, t_col> &m) const;  // (true or false) matrix == RotMat
    bool operator!=(const Rotation<t_type, t_row, t_col> &m) const;  // (true or false) matrix != RotMat
    bool operator==(const Transform<t_type, t_row, t_col> &m) const; // (true or false) matrix == Transform
    bool operator!=(const Transform<t_type, t_row, t_col> &m) const; // (true or false) matrix != Transform

    void Print(const char endChar = 0);

    /* Friend classes */
    template <uint16_t row, uint16_t col, typename type> friend class Matrix;
    template <uint16_t row, uint16_t col, typename type> friend class CscMatrix;
    template <typename type, uint16_t row, uint16_t col> friend class Matrix3;
    template <typename type, uint16_t row, uint16_t col> friend class Rotation;
    template <typename type, uint16_t row, uint16_t col> friend class Transform;

    template <uint16_t row, typename type> friend class Vector;
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
    template <uint16_t row, uint16_t col, typename type> friend class QR;
    template <uint16_t row, uint16_t col, typename type> friend class SVD;

    /* Friend template function */
    template <uint16_t row, uint16_t col, typename type>
    friend Matrix<row, col, type> operator*(const type s, const Matrix<row, col, type> &m); // scalar * matrix
};

} // namespace Math
} // namespace dt

#include "dtMatrix.tpp"

#include "dtMatrix0.h"

#endif // DTMATH_DTMATRIX_H_
