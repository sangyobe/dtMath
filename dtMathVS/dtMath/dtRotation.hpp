/*!
\file       dtRotation.hpp
\brief      dtMath, Rotation matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#pragma once

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
template <typename m_type, uint16_t m_row> class CdtQuaternion;
template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtMatrix;
template <typename m_type, uint16_t m_row, uint16_t m_col> class CdtMatrix3;

template <typename m_type = float, uint16_t m_row = 3, uint16_t m_col = 3>
class CdtRotation
{
private:
    m_type m_tolerance = std::numeric_limits<m_type>::epsilon();
    m_type m_elem[m_row * m_col];
    CdtRotation(const m_type* element);
    inline void Euler2RotMat(const uint16_t order, const m_type* e);
    inline void Quat2RotMat(const m_type* q);

public:
    CdtRotation();
    CdtRotation(const m_type* element, const size_t n_byte);
    CdtRotation(const char c, const m_type* element, const size_t n_byte);
    CdtRotation(
        const m_type m00, const m_type m01, const m_type m02,
        const m_type m10, const m_type m11, const m_type m12,
        const m_type m20, const m_type m21, const m_type m22);
    CdtRotation(const uint16_t order, const m_type angle);
    CdtRotation(const uint16_t order, const m_type angle1, const m_type angle2);
    CdtRotation(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3);
    CdtRotation(const CdtRotation& m);
    CdtRotation(const CdtMatrix3<m_type, m_row, m_col>& m);
    CdtRotation(const CdtMatrix<m_row, m_col, m_type>& m);
    CdtRotation(const uint16_t order, const CdtVector3<m_type, 3>& e);
    CdtRotation(const uint16_t order, const CdtVector<3, m_type>& e);
    CdtRotation(const CdtQuaternion<m_type, 4>& q);
    ~CdtRotation() {}

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

    void SetElement(const uint16_t order, const m_type angle);
    void SetElement(const uint16_t order, const m_type angle1, const m_type angle2);
    void SetElement(const uint16_t order, const m_type angle1, const m_type angle2, const m_type angle3);

    void SetElement(const CdtRotation& m);
    void SetElement(const CdtMatrix3<m_type, m_row, m_col>& m);
    void SetElement(const CdtMatrix<m_row, m_col, m_type>& m);

    void SetElement(const uint16_t order, const CdtVector3<m_type, 3>& e);
    void SetElement(const uint16_t order, const CdtVector<3, m_type>& e);
    void SetElement(const uint16_t order, const m_type *e);

    void SetElement(const CdtQuaternion<m_type, 4>& q);
    void SetElement(const m_type *q);
    void SetElement(const m_type w, const m_type x, const m_type y, const m_type z);

    void SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2);
    void SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2);

    const m_type* const GetElementsAddr() const;
    CdtVector3<m_type, 3> GetRowVec(const uint16_t idxRow) const;
    CdtVector3<m_type, 3> GetColVec(const uint16_t idxCol) const;
    int8_t GetRowVec(const uint16_t idxRow, CdtVector3<m_type, 3>& v) const;
    int8_t GetColVec(const uint16_t idxCol, CdtVector3<m_type, 3>& v) const;
    CdtVector3<m_type, 3> GetEulerAngles(uint16_t order) const;
    CdtRotation Transpose() const;
    CdtRotation log() const;                                // log(R) = ln(R) : SO(3) -> so(3), SO(3) is Special Orthogonal Group, so(3) is the set of skew-symmetric 3x3 matrices
    CdtVector3<m_type, 3> Log() const;                      // Log(R) = u*phi : SO(3) -> R3
    CdtRotation ode(m_type wx, m_type wy, m_type wz) const; // dR/dt = R*[w]x
    CdtRotation ode(m_type *w) const;                       // dR/dt = R*[w]x
    CdtRotation ode(CdtVector3<m_type, 3> w) const;         // dR/dt = R*[w]x
    CdtRotation ode(CdtVector<3, m_type> w) const;          // dR/dt = R*[w]x
    CdtRotation Inv() const;

    /* Member access operators */
    // returns a row of modifiable elements
    m_type& operator ()(uint16_t irow, uint16_t icol) { return m_elem[irow * m_col + icol]; }
    // returns a row of non-modifiable elements
    const m_type& operator ()(uint16_t irow, uint16_t icol) const { return m_elem[irow * m_col + icol]; }

    /* Assignment operators */
    CdtRotation& operator =(const CdtRotation& m);                  // matrix = matrix
    CdtCommaInit<m_row*m_col, m_type> operator <<(const m_type s);  // Init first matrix elements

    /* Arithmetic operators */
    CdtRotation operator -() const;  // minus sign
    CdtMatrix3<m_type, m_row, m_col> operator +(const CdtRotation& m) const;
    CdtMatrix3<m_type, m_row, m_col> operator -(const CdtRotation& m) const;
    CdtMatrix3<m_type, m_row, m_col> operator +(const CdtMatrix3<m_type, m_row, m_col>& m) const;
    CdtMatrix3<m_type, m_row, m_col> operator -(const CdtMatrix3<m_type, m_row, m_col>& m) const;
    CdtMatrix3<m_type, m_row, m_col> operator +(const CdtMatrix<m_row, m_col, m_type>& m) const;
    CdtMatrix3<m_type, m_row, m_col> operator -(const CdtMatrix<m_row, m_col, m_type>& m) const;
    CdtMatrix3<m_type, m_row, m_col> operator *(const m_type s) const;
    CdtMatrix3<m_type, m_row, m_col> operator /(const m_type s) const;

    template <uint16_t col>
    CdtMatrix<m_row, col, m_type> operator *(const CdtMatrix<m_col, col, m_type>& m) const;         // RotMat * matrix
    CdtMatrix3<m_type, m_row, m_col> operator *(const CdtMatrix3<m_type, m_row, m_col>& m) const;   // RotMat * matrix
    CdtRotation operator *(const CdtRotation& m) const;                                             // RotMat * RotMat
    CdtVector<m_row, m_type> operator *(const CdtVector<m_col, m_type>& v) const;                   // RotMat * vector
    CdtVector3<m_type, m_row> operator *(const CdtVector3<m_type, m_col>& v) const;                 // RotMat * vector
    CdtRotation operator &(const CdtVector<m_col, m_type>& v) const;                                // RotMat * [v]x, []x is skew-symmetric matrix
    CdtRotation operator &(const CdtVector3<m_type, m_col>& v) const;                               // RotMat * [v]x, []x is skew-symmetric matrix

    /* Comparison operators */
    bool operator ==(const CdtRotation& m) const;                       // (true or false) matrix == matrix
    bool operator !=(const CdtRotation& m) const;                       // (true or false) matrix != matrix
    bool operator ==(const CdtMatrix3<m_type, m_row, m_col>& m) const;  // (true or false) matrix == matrix
    bool operator !=(const CdtMatrix3<m_type, m_row, m_col>& m) const;  // (true or false) matrix != matrix
    bool operator ==(const CdtMatrix<m_row, m_col, m_type>& m) const;   // (true or false) matrix == matrix
    bool operator !=(const CdtMatrix<m_row, m_col, m_type>& m) const;   // (true or false) matrix != matrix

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row, uint16_t col> friend class CdtRotation;
    template <typename type, uint16_t row, uint16_t col> friend class CdtMatrix3;
    template <typename type, uint16_t row, uint16_t col> friend class CdtTransform;
    template <uint16_t row, uint16_t col, typename type> friend class CdtMatrix;

    template <typename type, uint16_t row> friend class CdtQuaternion;

    /* Friend template function */
    template <typename type, uint16_t row, uint16_t col>
    friend CdtMatrix3<type, row, col> operator *(const type s, const CdtRotation<type, row, col>& m); // scalar * RotMat
};

#include "dtRotation.ipp"
