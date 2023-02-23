/*!
\file       dtMatrix3.h
\brief      dtMath, 3x3 Matrix class, lighter and faster than general matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATRIX3_H_
#define DTMATH_DTMATRIX3_H_

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
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtMatrix;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtNoPivLU;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtPartialPivLU;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtLLT;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtLDLT;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtQR;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtSVD;

template <typename m_type = float, uint16_t m_row = 3, uint16_t m_col = 3>
class CdtMatrix3
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row * m_col];
    CdtMatrix3(const m_type* element);

public:
    CdtMatrix3();
    CdtMatrix3(const m_type* element, const size_t n_byte);
    CdtMatrix3(const char c, const m_type* element, const size_t n_byte);
    CdtMatrix3(
        const m_type m00, const m_type m01, const m_type m02,
        const m_type m10, const m_type m11, const m_type m12,
        const m_type m20, const m_type m21, const m_type m22);
    CdtMatrix3(const CdtMatrix3& m);
    CdtMatrix3(const CdtRotation<m_type, m_row, m_col>& m);
    CdtMatrix3(const CdtMatrix<m_row, m_col, m_type>& m);
    ~CdtMatrix3() {}

    void SetZero();
    void SetIdentity();
    void SetDiagonal(const m_type d1, const m_type d2, const m_type d3);
    void SetDiagonal(const m_type* element, const size_t n_byte);
    void SetDiagonal(const CdtVector<m_row, m_type>& v);
    void SetDiagonal(const CdtVector3<m_type, m_row>& v);
    void SetFill(const m_type value);
    void SetElement(const m_type* element, const size_t n_byte);
    void SetElement(
        const m_type m00, const m_type m01, const m_type m02,
        const m_type m10, const m_type m11, const m_type m12,
        const m_type m20, const m_type m21, const m_type m22);
    void SetElement(const CdtMatrix3& m);
    void SetElement(const CdtRotation<m_type, m_row, m_col>& m);
    void SetElement(const CdtMatrix<m_row, m_col, m_type>& m);
    template <uint16_t col>
    void SetRowVec(const uint16_t idxRow, const CdtVector<col, m_type>& v);
    void SetRowVec(const uint16_t idxRow, const CdtVector3<m_type, m_col>& v);
    void SetRowVec(const uint16_t idxRow, const m_type *v, const size_t n_byte);
    template <uint16_t row>
    void SetColVec(const uint16_t idxCol, const CdtVector<row, m_type>& v);
    void SetColVec(const uint16_t idxCol, const CdtVector3<m_type, m_row>& v);
    void SetColVec(const uint16_t idxCol, const m_type *v, const size_t n_byte);
    void SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2);
    void SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2);

    const m_type* const GetElementsAddr() const;
    CdtVector3<m_type, m_col> GetRowVec(const uint16_t idxRow) const;
    CdtVector3<m_type, m_row> GetColVec(const uint16_t idxCol) const;
    int8_t GetRowVec(const uint16_t idxRow, CdtVector3<m_type, m_col>& v) const;
    int8_t GetColVec(const uint16_t idxCol, CdtVector3<m_type, m_row>& v) const;
    CdtMatrix3 Transpose() const;

    m_type Trace() const;
    m_type GetNorm() const;     // Frobenius Norm (Euclidean norm)
    m_type GetSqNorm() const;   // Squared Frobenius Norm (Euclidean norm)

    CdtNoPivLU<m_row, m_col, m_type> NoPivLU() const;
    CdtPartialPivLU<m_row, m_col, m_type> PartialPivLU() const;
    CdtLLT<m_row, m_col, m_type> LLT() const;
    CdtLDLT<m_row, m_col, m_type> LDLT() const;
    CdtQR<m_row, m_col, m_type> QR() const;
    CdtSVD<m_row, m_col, m_type> SVD() const;

    CdtMatrix3 Inv(int8_t *isOk = nullptr) const;
    CdtMatrix3 PInv(int8_t *isOk = nullptr, m_type tolerance = 0) const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type& operator ()(uint16_t irow, uint16_t icol) { return m_elem[irow * m_col + icol]; }
    // returns a row of non-modifiable elements
    const m_type& operator ()(uint16_t irow, uint16_t icol) const { return m_elem[irow * m_col + icol]; }

    /* Assignment operators */
    CdtMatrix3& operator  =(const CdtMatrix3& m);                           // matrix  = matrix
    CdtMatrix3& operator +=(const CdtMatrix3& m);                           // matrix += matrix
    CdtMatrix3& operator -=(const CdtMatrix3& m);                           // matrix -= matrix
    CdtMatrix3& operator  =(const CdtRotation<m_type, m_row, m_col>& m);    // matrix  = matrix
    CdtMatrix3& operator +=(const CdtRotation<m_type, m_row, m_col>& m);    // matrix += matrix
    CdtMatrix3& operator -=(const CdtRotation<m_type, m_row, m_col>& m);    // matrix -= matrix
    CdtMatrix3& operator  =(const CdtMatrix<m_row, m_col, m_type>& m);      // matrix  = matrix
    CdtMatrix3& operator +=(const CdtMatrix<m_row, m_col, m_type>& m);      // matrix += matrix
    CdtMatrix3& operator -=(const CdtMatrix<m_row, m_col, m_type>& m);      // matrix -= matrix
    CdtMatrix3& operator  =(const m_type s);                                // matrix  = scalar, all elements set scalar
    CdtMatrix3& operator +=(const m_type s);                                // matrix += scalar, matrix(i) += scalar
    CdtMatrix3& operator -=(const m_type s);                                // matrix -= scalar, matrix(i) -= scalar
    CdtMatrix3& operator *=(const m_type s);                                // matrix *= scalar
    CdtMatrix3& operator /=(const m_type s);                                // matrix /= scalar
    CdtCommaInit<m_row*m_col, m_type> operator <<(const m_type s);          // Init first matrix elements

    /* Arithmetic operators */
    CdtMatrix3 operator -() const;                                          // minus sign
    CdtMatrix3 operator +(const CdtMatrix3& m) const;                       // matrix + matrix
    CdtMatrix3 operator -(const CdtMatrix3& m) const;                       // matrix - matrix
    CdtMatrix3 operator +(const CdtRotation<m_type, m_row, m_col>& m) const;// matrix + matrix
    CdtMatrix3 operator -(const CdtRotation<m_type, m_row, m_col>& m) const;// matrix - matrix
    CdtMatrix3 operator +(const CdtMatrix<m_row, m_col, m_type>& m) const;  // matrix + matrix
    CdtMatrix3 operator -(const CdtMatrix<m_row, m_col, m_type>& m) const;  // matrix - matrix
    CdtMatrix3 operator +(const m_type s) const;                            // matrix + scalar, matrix(i) + scalar
    CdtMatrix3 operator -(const m_type s) const;                            // matrix - scalar, matrix(i) - scalar
    CdtMatrix3 operator *(const m_type s) const;                            // matrix * scalar
    CdtMatrix3 operator /(const m_type s) const;                            // matrix / scalar

    template <uint16_t col>
    CdtMatrix<m_row, col, m_type> operator *(const CdtMatrix<m_col, col, m_type>& m) const; // matrix * matrix
    CdtMatrix3 operator *(const CdtMatrix3& m) const;                                       // matrix * matrix
    CdtMatrix3 operator *(const CdtRotation<m_type, m_row, m_col>& m) const;                // matrix * rotation
    CdtVector<m_row, m_type> operator *(const CdtVector<m_col, m_type>& v) const;           // matrix * vector
    CdtVector3<m_type, m_row> operator *(const CdtVector3<m_type, m_col>& v) const;         // matrix * vector
    CdtMatrix3 operator &(const CdtVector<m_col, m_type>& v) const;                         // matrix * [v]x, []x is skew-symmetric matrix
    CdtMatrix3 operator &(const CdtVector3<m_type, m_col>& v) const;                        // matrix * [v]x, []x is skew-symmetric matrix

    /* Comparison operators */
    bool operator ==(const CdtMatrix3& m) const;                            // (true or false) matrix == matrix
    bool operator !=(const CdtMatrix3& m) const;                            // (true or false) matrix != matrix
    bool operator ==(const CdtRotation<m_type, m_row, m_col>& m) const;     // (true or false) matrix == matrix
    bool operator !=(const CdtRotation<m_type, m_row, m_col>& m) const;     // (true or false) matrix != matrix
    bool operator ==(const CdtMatrix<m_row, m_col, m_type>& m) const;       // (true or false) matrix == matrix
    bool operator !=(const CdtMatrix<m_row, m_col, m_type>& m) const;       // (true or false) matrix != matrix

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row, uint16_t col> friend class CdtMatrix3;
    template <uint16_t row, uint16_t col, typename type> friend class CdtMatrix;
    template <typename type, uint16_t row, uint16_t col> friend class CdtRotation;

    template <uint16_t row, uint16_t col, typename type> friend class CdtLowerTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class CdtUpperTriangular;
    template <uint16_t row, uint16_t col, typename type> friend class CdtNoPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class CdtPartialPivLU;
    template <uint16_t row, uint16_t col, typename type> friend class CdtLLT;
    template <uint16_t row, uint16_t col, typename type> friend class CdtLDLT;
    template <uint16_t row, uint16_t col, typename type> friend class CdtQR;
    template <uint16_t row, uint16_t col, typename type> friend class CdtSVD;

    /* Friend template function */
    template <typename type, uint16_t row, uint16_t col>
    friend CdtMatrix3<type, row, col> operator*(const type s, const CdtMatrix3<type, row, col>& m); // scalar * matrix
};

#include "dtMatrix3.tpp"

#endif // DTMATH_DTMATRIX3_H_
