/*!
\file       dtMatrix.h
\brief      dtMath, General Matrix(m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
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

template <uint16_t m_size, typename m_type> class CdtCommaInit;
template <uint16_t m_row, typename m_type> class CdtVector;
template <typename m_type, uint16_t m_row> class CdtVector3;
template <typename m_type, uint16_t m_row> class CdtVector4;
template <typename m_type, uint16_t m_row> class CdtVector6;
template <typename m_type, uint16_t m_row, uint16_t m_col> class CdtMatrix3;
template <typename m_type, uint16_t m_row, uint16_t m_col> class CdtRotation;
template <typename m_type, uint16_t m_row, uint16_t m_col> class CdtTransform;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtNoPivLU;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtPartialPivLU;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtLLT;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtLDLT;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtQR;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtSVD;

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class CdtMatrix
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row * m_col];
    CdtMatrix(const m_type* element);

public:
    CdtMatrix();
    CdtMatrix(const m_type* element, const size_t n_byte);
    CdtMatrix(const char c, const m_type* element, const size_t n_byte);
    CdtMatrix(const CdtMatrix& m);
    ~CdtMatrix() {}

    void SetZero();
    void SetIdentity();
    void SetDiagonal(const m_type* element, const size_t n_byte);
    void SetFill(const m_type value);
    void SetElement(const m_type* element, const size_t n_byte);
    template <uint16_t row, uint16_t col>
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const CdtMatrix<row, col, m_type>& m);
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const CdtMatrix3<m_type, 3, 3>& m);
    void SetBlock(const uint16_t idxRow, const uint16_t idxCol, const CdtRotation<m_type, 3, 3>& m);
    template <uint16_t col>
    void SetRowVec(const uint16_t idxRow, const CdtVector<col, m_type>& v);
    void SetRowVec(const uint16_t idxRow, const CdtVector3<m_type, 3>&v);
    void SetRowVec(const uint16_t idxRow, const CdtVector4<m_type, 4>& v);
    void SetRowVec(const uint16_t idxRow, const CdtVector6<m_type, 6>& v);
    void SetRowVec(const uint16_t idxRow, const m_type *v, const size_t n_byte);
    template <uint16_t row>
    void SetColVec(const uint16_t idxCol, const CdtVector<row, m_type>& v);
    void SetColVec(const uint16_t idxCol, const CdtVector3<m_type, 3>& v);
    void SetColVec(const uint16_t idxCol, const CdtVector4<m_type, 4>& v);
    void SetColVec(const uint16_t idxCol, const CdtVector6<m_type, 6>& v);
    void SetColVec(const uint16_t idxCol, const m_type *v, const size_t n_byte);
    void SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2);
    void SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2);

    const m_type* const GetElementsAddr() const;
    uint16_t GetRowSize() const { return m_row; }   // size of row
    uint16_t GetColSize() const { return m_col; }   // size of colum
    template <uint16_t row, uint16_t col>
    CdtMatrix<row, col, m_type> GetBlock(const uint16_t idxRow, const uint16_t idxCol);
    CdtVector<m_col, m_type> GetRowVec(const uint16_t idxRow) const;
    CdtVector<m_row, m_type> GetColVec(const uint16_t idxCol) const;
    template <uint16_t row, uint16_t col>
    int8_t GetBlock(const uint16_t idxRow, const uint16_t idxCol, CdtMatrix<row, col, m_type>& m);
    int8_t GetRowVec(const uint16_t idxRow, CdtVector<m_col, m_type>& v) const;
    int8_t GetColVec(const uint16_t idxCol, CdtVector<m_row, m_type>& v) const;
    CdtCscMatrix<m_row, m_col, m_type> GetCscMat() const; // Compressed Sparse Column Matrix
    CdtMatrix<m_col, m_row, m_type> Transpose() const;

    m_type Trace() const;
    m_type GetNorm() const;     // Frobenius Norm (Euclidean norm)
    m_type GetSqNorm() const;   // Squared Frobenius Norm (Euclidean norm)
    m_type Determinant() const; // From LU Decomposition

    CdtNoPivLU<m_row, m_col, m_type> NoPivLU() const;
    CdtPartialPivLU<m_row, m_col, m_type> PartialPivLU() const;
    CdtLLT<m_row, m_col, m_type> LLT() const;
    CdtLDLT<m_row, m_col, m_type> LDLT() const;
    CdtQR<m_row, m_col, m_type> QR() const;
    CdtSVD<m_row, m_col, m_type> SVD() const;

    CdtMatrix<m_row, m_col, m_type> Inv(int8_t *isOk = nullptr) const;
    CdtMatrix<m_col, m_row, m_type> PInv(int8_t *isOk = nullptr, m_type tolerance = std::numeric_limits<m_type>::epsilon()) const;

    /* Member access operators */
    // returns a row of modifiable elements
    inline m_type& operator ()(uint16_t irow, uint16_t icol) { return m_elem[irow * m_col + icol]; }
    // returns a row of non-modifiable elements
    inline const m_type& operator ()(uint16_t irow, uint16_t icol) const { return m_elem[irow * m_col + icol]; }

    /* Assignment operators */
    CdtMatrix& operator  =(const CdtMatrix& m);                             // matrix  = matrix
    CdtMatrix& operator +=(const CdtMatrix& m);                             // matrix += matrix
    CdtMatrix& operator -=(const CdtMatrix& m);                             // matrix -= matrix
    CdtMatrix& operator  =(const CdtMatrix3<m_type, m_row, m_col>& m);      // matrix  = matrix3
    CdtMatrix& operator +=(const CdtMatrix3<m_type, m_row, m_col>& m);      // matrix += matrix3
    CdtMatrix& operator -=(const CdtMatrix3<m_type, m_row, m_col>& m);      // matrix -= matrix3
    CdtMatrix& operator  =(const CdtRotation<m_type, m_row, m_col>& m);     // matrix  = RotMat
    CdtMatrix& operator +=(const CdtRotation<m_type, m_row, m_col>& m);     // matrix += RotMat
    CdtMatrix& operator -=(const CdtRotation<m_type, m_row, m_col>& m);     // matrix -= RotMat
    CdtMatrix& operator  =(const CdtTransform<m_type, m_row, m_col>& m);    // matrix + Transform
    CdtMatrix& operator +=(const CdtTransform<m_type, m_row, m_col>& m);    // matrix + Transform
    CdtMatrix& operator -=(const CdtTransform<m_type, m_row, m_col>& m);    // matrix - Transform
    CdtMatrix& operator  =(const m_type s);                                 // matrix  = scalar, all elements set scalar
    CdtMatrix& operator +=(const m_type s);                                 // matrix += scalar, matrix(i) += scalar
    CdtMatrix& operator -=(const m_type s);                                 // matrix -= scalar, matrix(i) -= scalar
    CdtMatrix& operator *=(const m_type s);                                 // matrix *= scalar
    CdtMatrix& operator /=(const m_type s);                                 // matrix /= scalar
    CdtCommaInit<m_row*m_col, m_type> operator <<(const m_type s);          // Init first matrix elements

    /* Arithmetic operators */
    CdtMatrix operator -() const;                                            // minus sign
    CdtMatrix operator +(const CdtMatrix& m) const;                          // matrix + matrix
    CdtMatrix operator -(const CdtMatrix& m) const;                          // matrix - matrix
    CdtMatrix operator +(const CdtMatrix3<m_type, m_row, m_col>& m) const;   // matrix + matrix3
    CdtMatrix operator -(const CdtMatrix3<m_type, m_row, m_col>& m) const;   // matrix - matrix3
    CdtMatrix operator +(const CdtRotation<m_type, m_row, m_col>& m) const;  // matrix + RotMat
    CdtMatrix operator -(const CdtRotation<m_type, m_row, m_col>& m) const;  // matrix - RotMat
    CdtMatrix operator +(const CdtTransform<m_type, m_row, m_col>& m) const; // matrix + Transform
    CdtMatrix operator -(const CdtTransform<m_type, m_row, m_col>& m) const; // matrix - Transform
    CdtMatrix operator +(const m_type s) const;                              // matrix + scalar, matrix(i) + scalar
    CdtMatrix operator -(const m_type s) const;                              // matrix - scalar, matrix(i) - scalar
    CdtMatrix operator *(const m_type s) const;                              // matrix * scalar
    CdtMatrix operator /(const m_type s) const;                              // matrix / scalar

    template <uint16_t col>
    CdtMatrix<m_row, col, m_type> operator *(const CdtMatrix<m_col, col, m_type>& m) const;       // matrix * matrix
    CdtMatrix<m_row, m_col, m_type> operator *(const CdtMatrix3<m_type, m_col, m_col>& m) const;  // matrix * matrix
    CdtMatrix<m_row, m_col, m_type> operator *(const CdtRotation<m_type, m_col, m_col>& m) const; // matrix * RotMat
    CdtMatrix<m_row, m_col, m_type> operator *(const CdtTransform<m_type, m_col, m_col>& m) const;// matrix * Transform
    CdtVector<m_row, m_type> operator *(const CdtVector<m_col, m_type>& v) const;                 // matrix * vector
    CdtVector<m_row, m_type> operator *(const CdtVector3<m_type, m_col>& v) const;                // matrix * vector3
    CdtVector<m_row, m_type> operator *(const CdtVector4<m_type, m_col>& v) const;                // matrix * vector4
    CdtVector<m_row, m_type> operator *(const CdtVector6<m_type, m_col>& v) const;                // matrix * vector6
    CdtMatrix<m_row, m_col, m_type> operator &(const CdtVector<m_col, m_type>& v) const;          // matrix3 * [v]x, []x is skew-symmetric matrix
    CdtMatrix<m_row, m_col, m_type> operator &(const CdtVector3<m_type, m_col>& v) const;         // matrix3 * [v]x, []x is skew-symmetric matrix

    /* Comparison operators */
    bool operator ==(const CdtMatrix& m) const;                             // (true or false) matrix1 == matrix
    bool operator !=(const CdtMatrix& m) const;                             // (true or false) matrix1 != matrix
    bool operator ==(const CdtMatrix3<m_type, m_row, m_col>& m) const;      // (true or false) matrix == matrix3
    bool operator !=(const CdtMatrix3<m_type, m_row, m_col>& m) const;      // (true or false) matrix != matrix3
    bool operator ==(const CdtRotation<m_type, m_row, m_col>& m) const;     // (true or false) matrix == RotMat
    bool operator !=(const CdtRotation<m_type, m_row, m_col>& m) const;     // (true or false) matrix != RotMat
    bool operator ==(const CdtTransform<m_type, m_row, m_col>& m) const;    // (true or false) matrix == Transform
    bool operator !=(const CdtTransform<m_type, m_row, m_col>& m) const;    // (true or false) matrix != Transform

    void Print(const char endChar = 0);

    /* Friend classes */
    template <uint16_t row, uint16_t col, typename type> friend class CdtMatrix;
    template <uint16_t row, uint16_t col, typename type> friend class CdtCscMatrix;
    template <typename type, uint16_t row, uint16_t col> friend class CdtMatrix3;
    template <typename type, uint16_t row, uint16_t col> friend class CdtRotation;
    template <typename type, uint16_t row, uint16_t col> friend class CdtTransform;

    template <uint16_t row, typename type> friend class CdtVector;
    template <typename type, uint16_t row> friend class CdtVector3;
    template <typename type, uint16_t row> friend class CdtVector4;
    template <typename type, uint16_t row> friend class CdtVector6;

    template <uint16_t row, uint16_t col, typename type> friend class CdtLowerTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class CdtUpperTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class CdtNoPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class CdtPartialPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class CdtLLT;
    template <uint16_t row, uint16_t col, typename type> friend class CdtLDLT;
    template <uint16_t row, uint16_t col, typename type> friend class CdtQR;
    template <uint16_t row, uint16_t col, typename type> friend class CdtSVD;

    /* Friend template function */
    template <uint16_t row, uint16_t col, typename type>
    friend CdtMatrix<row, col, type> operator*(const type s, const CdtMatrix<row, col, type>& m); // scalar * matrix
};

#include "dtMatrix.tpp"

#endif // DTMATH_DTMATRIX_H_
