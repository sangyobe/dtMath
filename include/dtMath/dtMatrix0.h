/*!
\file       dtMatrix0.h
\brief      dtMath, Dynamic Memory Alloction General Matrix(m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 04. 01
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATRIX0_H_
#define DTMATH_DTMATRIX0_H_

#include "dtDefine.h"
#include "dtMatrix.h"

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

template <typename t_type>
class Matrix<0, 0, t_type>
{
private:
    const t_type m_tolerance = std::numeric_limits<t_type>::epsilon();
    uint16_t m_row, m_col, m_size;
    t_type *m_elem;

public:
    Matrix();
    Matrix(const uint16_t row, const uint16_t col);
    Matrix(const uint16_t row, const uint16_t col, const t_type *element);
    Matrix(const Matrix &m);
    template <uint16_t row, uint16_t col>
    Matrix(const Matrix<row, col, t_type> &m);
    ~Matrix();

    void NewSize(const uint16_t row, const uint16_t col);                        // Resize Memory(m_elem)
    void NewSize(const uint16_t row, const uint16_t col, const t_type *element); // Resize Memory(m_elem) & set element
    void NewSize(const Matrix &m);
    template <uint16_t row, uint16_t col>
    void NewSize(const Matrix<row, col, t_type> &m);

    void ReSize(const uint16_t row, const uint16_t col);                        // Resize Memory(m_elem)
    void ReSize(const uint16_t row, const uint16_t col, const t_type *element); // Resize Memory(m_elem) & set element
    void ReSize(const Matrix &m);
    template <uint16_t row, uint16_t col>
    void ReSize(const Matrix<row, col, t_type> &m);

    void Release(); // Free Memory(m_elem)

    void SetZero();
    void SetIdentity();
    void SetDiagonal(const t_type *element, const size_t n_byte);
    void SetFill(const t_type value);
    void SetFillRow(const uint16_t idxRow, const t_type value);
    void SetFillCol(const uint16_t idxCol, const t_type value);
    void SetElement(const t_type *element, const size_t n_byte);
    void SetElement(const t_type *element, const uint16_t row, const uint16_t col);
    template <uint16_t row, uint16_t col>
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix<row, col, t_type> &m);
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix<0, 0, t_type> &m);
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix3<t_type, 3, 3> &m);
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Rotation<t_type, 3, 3> &m);
    template <uint16_t col>
    void SetRowVec(const uint16_t idxRow, const Vector<col, t_type> &v);
    void SetRowVec(const uint16_t idxRow, const Vector<0, t_type> &v);
    void SetRowVec(const uint16_t idxRow, const Vector3<t_type, 3> &v);
    void SetRowVec(const uint16_t idxRow, const Vector4<t_type, 4> &v);
    void SetRowVec(const uint16_t idxRow, const Vector6<t_type, 6> &v);
    void SetRowVec(const uint16_t idxRow, const t_type *v, const size_t n_byte);
    template <uint16_t row>
    void SetColVec(const uint16_t idxCol, const Vector<row, t_type> &v);
    void SetColVec(const uint16_t idxCol, const Vector<0, t_type> &v);
    void SetColVec(const uint16_t idxCol, const Vector3<t_type, 3> &v);
    void SetColVec(const uint16_t idxCol, const Vector4<t_type, 4> &v);
    void SetColVec(const uint16_t idxCol, const Vector6<t_type, 6> &v);
    void SetColVec(const uint16_t idxCol, const t_type *v, const size_t n_byte);
    void SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2);
    void SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2);

    const t_type *const GetElementsAddr() const;
    uint16_t GetRowSize() const { return m_row; } // size of row
    uint16_t GetColSize() const { return m_col; } // size of colum
    size_t GetMemSize() const;
    template <uint16_t row, uint16_t col>
    Matrix<row, col, t_type> GetBlock(const uint16_t idxRow, const uint16_t idxCol);
    Matrix<0, 0, t_type> GetBlock(const uint16_t idxRow, const uint16_t idxCol, const uint16_t row, const uint16_t col);
    template <uint16_t row, uint16_t col>
    int8_t GetBlock(const uint16_t idxRow, const uint16_t idxCol, Matrix<row, col, t_type> &m);
    int8_t GetBlock(const uint16_t idxRow, const uint16_t idxCol, Matrix<0, 0, t_type> &m);
    template <uint16_t col> Vector<col, t_type> GetRowVec(const uint16_t idxRow) const;
    template <uint16_t row> Vector<row, t_type> GetColVec(const uint16_t idxCol) const;
    Vector<0, t_type> GetRowVec(const uint16_t idxRow, const uint16_t col) const;
    Vector<0, t_type> GetColVec(const uint16_t idxCol, const uint16_t row) const;
    template <uint16_t col> int8_t GetRowVec(const uint16_t idxRow, Vector<col, t_type> &v) const;
    template <uint16_t row> int8_t GetColVec(const uint16_t idxCol, Vector<row, t_type> &v) const;
    int8_t GetRowVec(const uint16_t idxRow, Vector<0, t_type> &v) const;
    int8_t GetColVec(const uint16_t idxCol, Vector<0, t_type> &v) const;
    template <uint16_t row, uint16_t col>
    CscMatrix<row, col, t_type> GetCscMat() const; // Compressed Sparse Column MatrixX
    template <uint16_t row, uint16_t col>
    int8_t GetCscMat(CscMatrix<row, col, t_type> &m) const; // Compressed Sparse Column MatrixX

    Matrix<0, 0, t_type> Transpose() const;
    t_type Trace() const;
    t_type GetNorm() const;              // Frobenius Norm (Euclidean norm, L2 Norm)
    t_type GetSqNorm() const;            // Squared Frobenius Norm (Euclidean norm, Squared L2 Norm)
    t_type GetLpNorm(const int p) const; // Generalized Norm (Lp Norm)
    t_type Determinant() const;          // From LU Decomposition

    NoPivLU<0, 0, t_type> GetNoPivLU() const;
    PartialPivLU<0, 0, t_type> GetPartialPivLU() const;
    FullPivLU<0, 0, t_type> GetFullPivLU() const;
    LLT<0, 0, t_type> GetLLT() const;
    LDLT<0, 0, t_type> GetLDLT() const;
    // QR<t_row, t_col, t_type> GetQR() const;
    SVD<0, 0, t_type> GetSVD() const;

    Matrix<0, 0, t_type> Inv(int8_t *isOk = nullptr) const;
    Matrix<0, 0, t_type> FInv(int8_t *isOk = nullptr) const;
    Matrix<0, 0, t_type> PInv(int8_t *isOk = nullptr, t_type tolerance = std::numeric_limits<t_type>::epsilon()) const;

    /* Member access operators */
    t_type &operator()(uint16_t irow, uint16_t icol);             // returns a row of modifiable elements
    const t_type &operator()(uint16_t irow, uint16_t icol) const; // returns a row of non-modifiable elements

    /* Assignment operators */
    Matrix &operator=(const Matrix &m);  // matrix  = matrix
    Matrix &operator+=(const Matrix &m); // matrix += matrix
    Matrix &operator-=(const Matrix &m); // matrix -= matrix
    template <uint16_t row, uint16_t col>
    Matrix &operator=(const Matrix<row, col, t_type> &m); // matrix  = matrix
    template <uint16_t row, uint16_t col>
    Matrix &operator+=(const Matrix<row, col, t_type> &m); // matrix += matrix
    template <uint16_t row, uint16_t col>
    Matrix &operator-=(const Matrix<row, col, t_type> &m); // matrix -= matrix
    Matrix &operator=(const Matrix3<t_type, 3, 3> &m);     // matrix  = matrix3
    Matrix &operator+=(const Matrix3<t_type, 3, 3> &m);    // matrix += matrix3
    Matrix &operator-=(const Matrix3<t_type, 3, 3> &m);    // matrix -= matrix3
    Matrix &operator=(const Rotation<t_type, 3, 3> &m);    // matrix  = RotMat
    Matrix &operator+=(const Rotation<t_type, 3, 3> &m);   // matrix += RotMat
    Matrix &operator-=(const Rotation<t_type, 3, 3> &m);   // matrix -= RotMat
    Matrix &operator=(const Transform<t_type, 4, 4> &m);   // matrix  = Transform
    Matrix &operator+=(const Transform<t_type, 4, 4> &m);  // matrix += Transform
    Matrix &operator-=(const Transform<t_type, 4, 4> &m);  // matrix -= Transform
    Matrix &operator=(const t_type s);                     // matrix  = scalar, all elements set scalar
    Matrix &operator+=(const t_type s);                    // matrix += scalar, matrix(i) += scalar
    Matrix &operator-=(const t_type s);                    // matrix -= scalar, matrix(i) -= scalar
    Matrix &operator*=(const t_type s);                    // matrix *= scalar
    Matrix &operator/=(const t_type s);                    // matrix /= scalar
    CommaInit<0, t_type> operator<<(const t_type s);       // Init first matrix elements

    /* Arithmetic operators */
    Matrix operator-() const;                // minus sign
    Matrix operator+(const Matrix &m) const; // matrix + matrix
    Matrix operator-(const Matrix &m) const; // matrix - matrix
    template <uint16_t row, uint16_t col>
    Matrix operator+(const Matrix<row, col, t_type> &m) const; // matrix + matrix
    template <uint16_t row, uint16_t col>
    Matrix operator-(const Matrix<row, col, t_type> &m) const; // matrix - matrix
    Matrix operator+(const Matrix3<t_type, 3, 3> &m) const;    // matrix + matrix3
    Matrix operator-(const Matrix3<t_type, 3, 3> &m) const;    // matrix - matrix3
    Matrix operator+(const Rotation<t_type, 3, 3> &m) const;   // matrix + RotMat
    Matrix operator-(const Rotation<t_type, 3, 3> &m) const;   // matrix - RotMat
    Matrix operator+(const Transform<t_type, 4, 4> &m) const;  // matrix + Transform
    Matrix operator-(const Transform<t_type, 4, 4> &m) const;  // matrix - Transform
    Matrix operator+(const t_type s) const;                    // matrix + scalar, matrix(i) + scalar
    Matrix operator-(const t_type s) const;                    // matrix - scalar, matrix(i) - scalar
    Matrix operator*(const t_type s) const;                    // matrix * scalar
    Matrix operator/(const t_type s) const;                    // matrix / scalar

    Matrix<0, 0, t_type> operator*(const Matrix<0, 0, t_type> &m) const;    // matrix * matrix
    Matrix<0, 0, t_type> operator*(const Matrix3<t_type, 3, 3> &m) const;   // matrix * matrix
    Matrix<0, 0, t_type> operator*(const Rotation<t_type, 3, 3> &m) const;  // matrix * RotMat
    Matrix<0, 0, t_type> operator*(const Transform<t_type, 4, 4> &m) const; // matrix * Transform
    template <uint16_t row, uint16_t col>
    Matrix<0, 0, t_type> operator*(const Matrix<row, col, t_type> &m) const; // matrix * matrix
    Vector<0, t_type> operator*(const Vector<0, t_type> &v) const;           // matrix * vector
    Vector<0, t_type> operator*(const Vector3<t_type, 3> &v) const;          // matrix * vector3
    Vector<0, t_type> operator*(const Vector4<t_type, 4> &v) const;          // matrix * vector4
    Vector<0, t_type> operator*(const Vector6<t_type, 6> &v) const;          // matrix * vector6
    template <uint16_t col>
    Vector<0, t_type> operator*(const Vector<col, t_type> &v) const;   // matrix * vector
    Matrix<0, 0, t_type> operator&(const Vector<0, t_type> &v) const;  // matrix3 * [v]x, []x is skew-symmetric matrix
    Matrix<0, 0, t_type> operator&(const Vector<3, t_type> &v) const;  // matrix3 * [v]x, []x is skew-symmetric matrix
    Matrix<0, 0, t_type> operator&(const Vector3<t_type, 3> &v) const; // matrix3 * [v]x, []x is skew-symmetric matrix

    /* Comparison operators */
    template <uint16_t row, uint16_t col>
    bool operator==(const Matrix<row, col, t_type> &m) const; // (true or false) matrix1 == matrix
    template <uint16_t row, uint16_t col>
    bool operator!=(const Matrix<row, col, t_type> &m) const; // (true or false) matrix1 != matrix
    bool operator==(const Matrix &m) const;                   // (true or false) matrix1 == matrix
    bool operator!=(const Matrix &m) const;                   // (true or false) matrix1 != matrix
    bool operator==(const Matrix3<t_type, 3, 3> &m) const;    // (true or false) matrix == matrix3
    bool operator!=(const Matrix3<t_type, 3, 3> &m) const;    // (true or false) matrix != matrix3
    bool operator==(const Rotation<t_type, 3, 3> &m) const;   // (true or false) matrix == RotMat
    bool operator!=(const Rotation<t_type, 3, 3> &m) const;   // (true or false) matrix != RotMat
    bool operator==(const Transform<t_type, 4, 4> &m) const;  // (true or false) matrix == Transform
    bool operator!=(const Transform<t_type, 4, 4> &m) const;  // (true or false) matrix != Transform

    void Print(const char endChar = 0);

    // /* Friend classes */
    template <uint16_t row, uint16_t col, typename type> friend class Matrix;
    template <uint16_t row, uint16_t col, typename type> friend class CscMatrix;
    // template <typename type, uint16_t row, uint16_t col> friend class Matrix3;
    // template <typename type, uint16_t row, uint16_t col> friend class Rotation;
    // template <typename type, uint16_t row, uint16_t col> friend class Transform;

    template <uint16_t row, typename type> friend class Vector;
    template <typename type, uint16_t row> friend class Vector3;
    template <typename type, uint16_t row> friend class Vector4;
    template <typename type, uint16_t row> friend class Vector6;

    template <uint16_t row, uint16_t col, typename type> friend class LowerMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class UpperMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class NoPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class PartialPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class FullPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class LLT;
    template <uint16_t row, uint16_t col, typename type> friend class LDLT;
    // template <uint16_t row, uint16_t col, typename type> friend class QR;
    template <uint16_t row, uint16_t col, typename type> friend class SVD;

    /* Friend template function */
    template <typename type>
    friend Matrix<0, 0, type> operator*(const type s, const Matrix<0, 0, type> &m); // scalar * matrix
};

} // namespace Math
} // namespace dt

#include "dtMatrix0.tpp"

#endif // DTMATH_DTMATRIX0_H_
