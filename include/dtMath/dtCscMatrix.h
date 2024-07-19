/*!
\file       dtCscMatrix.h
\brief      dtMath, Compressed Sparse Column Matrix (m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 05. 29
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCSC_MATRIX_H_
#define DTMATH_DTCSC_MATRIX_H_

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

template <uint16_t t_row, typename t_type> class Vector;
template <typename t_type, uint16_t t_row> class Vector3;
template <typename t_type, uint16_t t_row> class Vector4;
template <typename t_type, uint16_t t_row> class Vector6;
template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Matrix3;

template <uint16_t t_row, uint16_t t_col, typename t_type = float>
class CscMatrix
{
private:
    int m_elemNum;
    t_type m_elem[t_row * t_col];
    int m_rowIdx[t_row * t_col]; // row indices
    int m_colPtr[t_col + 1];     // column index pointer
    CscMatrix(const t_type *element, const int elemNum, const int *rowIdx, const int *colPtr);

public:
    CscMatrix();
    CscMatrix(const t_type *element, const size_t n_byte);
    CscMatrix(const CscMatrix &m);
    CscMatrix(const CscMatrix<0, 0, t_type> &m);
    CscMatrix(const Matrix<t_row, t_col, t_type> &m);
    CscMatrix(const Matrix<0, 0, t_type> &m);
    ~CscMatrix() {}

    void SetElement(const t_type *element, const size_t n_byte);
    void SetElement(const Matrix<t_row, t_col, t_type> &m);
    void SetElement(const Matrix3<t_type, t_row, t_col> &m);
    void SetElement(const Matrix<0, 0, t_type> &m);

    /* Data access functions for sparse matrix */
    const t_type *const GetDataAddr() const;
    uint16_t GetNoOfElem() const { return m_elemNum; }
    const int *const GetRowIdx() const;
    const int *const GetColPtr() const;
    uint16_t GetRowSize() const { return t_row; } // size of row
    uint16_t GetColSize() const { return t_col; } // size of colum
    Vector<t_col, t_type> GetRowVec(const uint16_t idxRow) const;
    Vector<t_row, t_type> GetColVec(const uint16_t idxCol) const;
    Vector<0, t_type> GetRowVec(const uint16_t idxRow, const uint16_t col) const;
    Vector<0, t_type> GetColVec(const uint16_t idxCol, const uint16_t row) const;
    int8_t GetRowVec(const uint16_t idxRow, Vector<t_col, t_type> &v) const;
    int8_t GetColVec(const uint16_t idxCol, Vector<t_row, t_type> &v) const;
    int8_t GetRowVec(const uint16_t idxRow, Vector<0, t_type> &v) const;
    int8_t GetColVec(const uint16_t idxCol, Vector<0, t_type> &v) const;
    Matrix<t_row, t_col, t_type> GetDenseMat() const;
    CscMatrix<t_row, t_col, t_type> GetUpperTriangular() const;
    CscMatrix<t_row, t_col, t_type> GetLowerTriangular() const;
    int ComparePattern(const CscMatrix<t_row, t_col, t_type> &b) const;

    CscMatrix<t_col, t_row, t_type> Transpose() const;
    t_type GetNorm() const;
    t_type GetSqNorm() const;
    t_type GetLpNorm(const int p) const;

    /* Assignment operators */
    CscMatrix &operator=(const CscMatrix &m);               // matrix  = matrix
    CscMatrix &operator=(const CscMatrix<0, 0, t_type> &m); // matrix  = matrix
    CscMatrix &operator*=(const t_type s);                  // matrix *= scalar
    CscMatrix &operator/=(const t_type s);                  // matrix /= scalar

    /* Arithmetic operators */
    CscMatrix operator-() const;               // minus sign
    CscMatrix operator*(const t_type s) const; // matrix * scalar
    CscMatrix operator/(const t_type s) const; // matrix / scalar

    CscMatrix operator+(const CscMatrix<t_row, t_col, t_type> &m) const; // Sparse matrix + Sparse matrix
    CscMatrix operator-(const CscMatrix<t_row, t_col, t_type> &m) const; // Sparse matrix - Sparse matrix
    template <uint16_t col>
    CscMatrix<t_row, col, t_type> operator*(const CscMatrix<t_col, col, t_type> &m) const; // Sparse matrix * Sparse matrix
    Vector<t_row, t_type> operator*(const Vector<t_col, t_type> &v) const;                 // matrix * vector
    Vector<t_row, t_type> operator*(const Vector<0, t_type> &v) const;                     // matrix * vector
    Vector<t_row, t_type> operator*(const Vector3<t_type, t_col> &v) const;                // matrix * vector3
    Vector<t_row, t_type> operator*(const Vector4<t_type, t_col> &v) const;                // matrix * vector4
    Vector<t_row, t_type> operator*(const Vector6<t_type, t_col> &v) const;                // matrix * vector6
    Vector<t_col, t_type> TposeVec(const Vector<t_row, t_type> &v) const;                  // matrix^T * vector
    Vector<t_col, t_type> TposeVec(const Vector<0, t_type> &v) const;                      // matrix^T * vector
    Vector<t_col, t_type> TposeVec(const Vector3<t_type, t_row> &v) const;                 // matrix^T * vector3
    Vector<t_col, t_type> TposeVec(const Vector4<t_type, t_row> &v) const;                 // matrix^T * vector4
    Vector<t_col, t_type> TposeVec(const Vector6<t_type, t_row> &v) const;                 // matrix^T * vector6

    /* Basic Solve Operations */
    Vector<t_row, t_type> LowerTriangularSolve(const Vector<t_row, t_type> &v) const;
    Vector<t_row, t_type> LTTriangularSolve(const Vector<t_row, t_type> &v) const;

    Vector<t_row, t_type> UpperTriangularSolve(const Vector<t_row, t_type> &v) const;
    Vector<t_row, t_type> UTTriangularSolve(const Vector<t_row, t_type> &v) const;

    /* Permutations */

    CscMatrix<t_row, t_col, t_type> PermuteRow(const Vector<t_row, int> &p) const; // Row permutation
    CscMatrix<t_row, t_col, t_type> PermuteCol(const Vector<t_col, int> &p) const; // Column permutation

    /* Factorization */

    // CscLLT<t_row, t_col, t_type> CscLLT() const; // Cholesky factorization
    // CscLDLT<t_row, t_col, t_type> CscLDLT() const; // Cholesky factorization
    // CscAMD<t_row, t_col, t_type> CscAMD() const; // Simplified Cholesky factorization

    void Print(const char endChar = 0);

    /* Friend classes */
    template <uint16_t row, uint16_t col, typename type> friend class CscMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class CscAMD;
    template <uint16_t row, uint16_t col, typename type> friend class CscLLT;
};

} // namespace Math
} // namespace dt

#include "dtCscMatrix.tpp"

#include "dtCscMatrix0.h"

#endif // DTMATH_DTCSC_MATRIX_H_
