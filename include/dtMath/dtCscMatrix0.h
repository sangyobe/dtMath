/*!
\file       dtCscMatrix0.h
\brief      dtMath, Dynamic Memory Allocation Compressed Sparse Column Matrix (m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 06. 10
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCSC_MATRIX0_H_
#define DTMATH_DTCSC_MATRIX0_H_

#include "dtCscMatrix.h"
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

template <typename t_type>
class CscMatrix<0, 0, t_type>
{
private:
    uint16_t m_row, m_col, m_size;
    int m_elemNum;
    t_type *m_elem;
    int *m_rowIdx; // row indices
    int *m_colPtr; // column index pointer
    bool m_compressed;

public:
    CscMatrix();
    CscMatrix(const uint16_t row, const uint16_t col, const t_type *element, const int elemNum, const int *rowIdx, const int *colPtr);
    CscMatrix(const uint16_t row, const uint16_t col);
    CscMatrix(const uint16_t row, const uint16_t col, const t_type *element);
    CscMatrix(const CscMatrix &m);
    CscMatrix(const Matrix<0, 0, t_type> &m);
    template <uint16_t row, uint16_t col>
    CscMatrix(const CscMatrix<row, col, t_type> &m);
    template <uint16_t row, uint16_t col>
    CscMatrix(const Matrix<row, col, t_type> &m);
    ~CscMatrix();

    void NewSize(const uint16_t row, const uint16_t col);
    void NewSize(const uint16_t row, const uint16_t col, const t_type *element);
    void NewSize(const CscMatrix &m);
    void NewSize(const Matrix<0, 0, t_type> &m);
    template <uint16_t row, uint16_t col>
    void NewSize(const CscMatrix<row, col, t_type> &m);
    template <uint16_t row, uint16_t col>
    void NewSize(const Matrix<row, col, t_type> &m);

    void ReSize(const uint16_t row, const uint16_t col);
    void ReSize(const uint16_t row, const uint16_t col, const t_type *element);
    void ReSize(const CscMatrix &m);
    void ReSize(const Matrix<0, 0, t_type> &m);
    template <uint16_t row, uint16_t col>
    void ReSize(const CscMatrix<row, col, t_type> &m);
    template <uint16_t row, uint16_t col>
    void ReSize(const Matrix<row, col, t_type> &m);

    void Compress();
    void Release();

    void SetElement(const t_type *element, const size_t n_byte);
    void SetElement(const t_type *element, const uint16_t row, const uint16_t col);
    void SetElement(const Matrix<0, 0, t_type> &m);
    void SetElement(const Matrix3<t_type, 3, 3> &m);
    template <uint16_t row, uint16_t col>
    void SetElement(const Matrix<row, col, t_type> &m);

    const t_type *const GetDataAddr() const;
    uint16_t GetNoOfElem() const { return m_elemNum; }
    const int *const GetRowIdx() const;
    const int *const GetColPtr() const;
    uint16_t GetRowSize() const { return m_row; } // size of row
    uint16_t GetColSize() const { return m_col; } // size of colum
    size_t GetMemSize() const;
    template <uint16_t col> Vector<col, t_type> GetRowVec(const uint16_t idxRow) const;
    template <uint16_t row> Vector<row, t_type> GetColVec(const uint16_t idxCol) const;
    Vector<0, t_type> GetRowVec(const uint16_t idxRow, const uint16_t col) const;
    Vector<0, t_type> GetColVec(const uint16_t idxCol, const uint16_t row) const;
    template <uint16_t col> int8_t GetRowVec(const uint16_t idxRow, Vector<col, t_type> &v) const;
    template <uint16_t row> int8_t GetColVec(const uint16_t idxCol, Vector<row, t_type> &v) const;
    int8_t GetRowVec(const uint16_t idxRow, Vector<0, t_type> &v) const;
    int8_t GetColVec(const uint16_t idxCol, Vector<0, t_type> &v) const;
    Matrix<0, 0, t_type> GetDenseMat() const;
    CscMatrix<0, 0, t_type> GetUpperTriangular() const;
    CscMatrix<0, 0, t_type> GetLowerTriangular() const;
    int ComparePattern(const CscMatrix<0, 0, t_type> &b) const;

    CscMatrix<0, 0, t_type> Transpose() const;
    t_type GetNorm() const;
    t_type GetSqNorm() const;
    t_type GetLpNorm(const int p) const;

    /* Assignment operators */
    template <uint16_t row, uint16_t col>
    CscMatrix &operator=(const CscMatrix<row, col, t_type> &m); // matrix  = matrix
    CscMatrix &operator=(const CscMatrix &m);                   // matrix  = matrix
    CscMatrix &operator*=(const t_type s);                      // matrix *= scalar
    CscMatrix &operator/=(const t_type s);                      // matrix /= scalar

    /* Arithmetic operators */
    CscMatrix operator-() const;               // minus sign
    CscMatrix operator*(const t_type s) const; // matrix * scalar
    CscMatrix operator/(const t_type s) const; // matrix / scalar

    CscMatrix operator+(const CscMatrix<0, 0, t_type> &m) const; // Sparse matrix + Sparse matrix
    CscMatrix operator-(const CscMatrix<0, 0, t_type> &m) const; // Sparse matrix - Sparse matrix
    CscMatrix operator*(const CscMatrix<0, 0, t_type> &m) const; // Sparse matrix * Sparse matrix
    template <uint16_t col>
    Vector<0, t_type> operator*(const Vector<col, t_type> &v) const; // matrix * vector
    Vector<0, t_type> operator*(const Vector<0, t_type> &v) const;   // matrix * vector
    Vector<0, t_type> operator*(const Vector3<t_type, 3> &v) const;  // matrix * vector3
    Vector<0, t_type> operator*(const Vector4<t_type, 4> &v) const;  // matrix * vector4
    Vector<0, t_type> operator*(const Vector6<t_type, 6> &v) const;  // matrix * vector6
    template <uint16_t row>
    Vector<0, t_type> TposeVec(const Vector<row, t_type> &v) const; // matrix^T * vector
    Vector<0, t_type> TposeVec(const Vector<0, t_type> &v) const;   // matrix^T * vector
    Vector<0, t_type> TposeVec(const Vector3<t_type, 3> &v) const;  // matrix^T * vector3
    Vector<0, t_type> TposeVec(const Vector4<t_type, 4> &v) const;  // matrix^T * vector4
    Vector<0, t_type> TposeVec(const Vector6<t_type, 6> &v) const;  // matrix^T * vector6

    /* Basic Solve Operations */
    Vector<0, t_type> LowerTriangularSolve(const Vector<0, t_type> &v);
    Vector<0, t_type> LTTriangularSolve(const Vector<0, t_type> &v);

    Vector<0, t_type> UpperTriangularSolve(const Vector<0, t_type> &v);
    Vector<0, t_type> UTTriangularSolve(const Vector<0, t_type> &v);

    /* Permutations */

    template <uint16_t row>
    CscMatrix<0, 0, t_type> PermuteRow(const Vector<row, int> &p) const; // Row permutation
    template <uint16_t col>
    CscMatrix<0, 0, t_type> PermuteCol(const Vector<col, int> &p) const; // Column permutation

    /* Factorization */

    // CscLLT<t_row, t_col, t_type> CscLLT() const; // Cholesky factorization
    // CscLDLT<t_row, t_col, t_type> CscLDLT() const; // Cholesky factorization
    // CscAMD<t_row, t_col, t_type> CscAMD() const; // Simplified Cholesky factorization

    void Print(const char endChar = 0);

    /* Friend classes */
};

} // namespace Math
} // namespace dt

#include "dtCscMatrix0.tpp"

#endif // DTMATH_DTCSC_MATRIX0_H_
