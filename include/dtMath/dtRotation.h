/*!
\file       dtRotation.h
\brief      dtMath, Rotation matrix class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTROTATION_H_
#define DTMATH_DTROTATION_H_

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
template <typename t_type, uint16_t t_row> class Quaternion;
template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Matrix3;

template <typename t_type = float, uint16_t t_row = 3, uint16_t t_col = 3>
class Rotation
{
private:
    t_type m_tolerance = std::numeric_limits<t_type>::epsilon();
    t_type m_elem[t_row * t_col];
    Rotation(const t_type *element);
    inline void Euler2RotMat(const uint16_t order, const t_type *e);
    inline void Quat2RotMat(const t_type *q);

public:
    Rotation();
    Rotation(const t_type *element, const size_t n_byte);
    Rotation(const char c, const t_type *element, const size_t n_byte);
    Rotation(
        const t_type m00, const t_type m01, const t_type m02,
        const t_type m10, const t_type m11, const t_type m12,
        const t_type m20, const t_type m21, const t_type m22);
    Rotation(const uint16_t order, const t_type angle);
    Rotation(const uint16_t order, const t_type angle1, const t_type angle2);
    Rotation(const uint16_t order, const t_type angle1, const t_type angle2, const t_type angle3);
    Rotation(const Rotation &m);
    Rotation(const Matrix3<t_type, t_row, t_col> &m);
    Rotation(const Matrix<t_row, t_col, t_type> &m);
    Rotation(const Matrix<0, 0, t_type> &m);
    Rotation(const uint16_t order, const Vector3<t_type, 3> &e);
    Rotation(const uint16_t order, const Vector<3, t_type> &e);
    Rotation(const uint16_t order, const Vector<0, t_type> &e);
    Rotation(const Quaternion<t_type, 4> &q);
    ~Rotation() {}

    void SetZero();
    void SetIdentity();
    void SetDiagonal(const t_type d1, const t_type d2, const t_type d3);
    void SetDiagonal(const t_type *element, const size_t n_byte);
    void SetDiagonal(const Vector<t_row, t_type> &v);
    void SetDiagonal(const Vector<0, t_type> &v);
    void SetDiagonal(const Vector3<t_type, t_row> &v);
    void SetFill(const t_type value);
    void SetElement(const t_type *element, const size_t n_byte);
    void SetElement(
        const t_type m00, const t_type m01, const t_type m02,
        const t_type m10, const t_type m11, const t_type m12,
        const t_type m20, const t_type m21, const t_type m22);

    void SetElement(const uint16_t order, const t_type angle);
    void SetElement(const uint16_t order, const t_type angle1, const t_type angle2);
    void SetElement(const uint16_t order, const t_type angle1, const t_type angle2, const t_type angle3);

    void SetElement(const Rotation &m);
    void SetElement(const Matrix3<t_type, t_row, t_col> &m);
    void SetElement(const Matrix<t_row, t_col, t_type> &m);
    void SetElement(const Matrix<0, 0, t_type> &m);

    void SetElement(const uint16_t order, const Vector3<t_type, 3> &e);
    void SetElement(const uint16_t order, const Vector<3, t_type> &e);
    void SetElement(const uint16_t order, const Vector<0, t_type> &e);
    void SetElement(const uint16_t order, const t_type *e);

    void SetElement(const Quaternion<t_type, 4> &q);
    void SetElement(const t_type *q);
    void SetElement(const t_type w, const t_type x, const t_type y, const t_type z);

    void SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2);
    void SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2);

    const t_type *const GetElementsAddr() const;
    Vector3<t_type, 3> GetRowVec(const uint16_t idxRow) const;
    Vector3<t_type, 3> GetColVec(const uint16_t idxCol) const;
    int8_t GetRowVec(const uint16_t idxRow, Vector3<t_type, 3> &v) const;
    int8_t GetColVec(const uint16_t idxCol, Vector3<t_type, 3> &v) const;
    Vector3<t_type, 3> GetEulerAngles(uint16_t order) const;
    Rotation Transpose() const;
    Rotation log() const;                                // log(R) = ln(R) : SO(3) -> so(3), SO(3) is Special Orthogonal Group, so(3) is the set of skew-symmetric 3x3 matrices
    Vector3<t_type, 3> Log() const;                      // Log(R) = u*phi : SO(3) -> R3
    Rotation ode(t_type wx, t_type wy, t_type wz) const; // dR/dt = R*[w]x
    Rotation ode(t_type *w) const;                       // dR/dt = R*[w]x
    Rotation ode(Vector3<t_type, 3> w) const;            // dR/dt = R*[w]x
    Rotation ode(Vector<3, t_type> w) const;             // dR/dt = R*[w]x
    Rotation ode(Vector<0, t_type> w) const;             // dR/dt = R*[w]x
    Rotation Inv() const;

    /* Member access operators */
    t_type &operator()(uint16_t irow, uint16_t icol);             // returns a row of modifiable elements
    const t_type &operator()(uint16_t irow, uint16_t icol) const; // returns a row of non-modifiable elements

    /* Assignment operators */
    Rotation &operator=(const Rotation &m);                      // matrix = matrix
    CommaInit<t_row * t_col, t_type> operator<<(const t_type s); // Init first matrix elements

    /* Arithmetic operators */
    Rotation operator-() const; // minus sign
    Matrix3<t_type, t_row, t_col> operator+(const Rotation &m) const;
    Matrix3<t_type, t_row, t_col> operator-(const Rotation &m) const;
    Matrix3<t_type, t_row, t_col> operator+(const Matrix3<t_type, t_row, t_col> &m) const;
    Matrix3<t_type, t_row, t_col> operator-(const Matrix3<t_type, t_row, t_col> &m) const;
    Matrix3<t_type, t_row, t_col> operator+(const Matrix<t_row, t_col, t_type> &m) const;
    Matrix3<t_type, t_row, t_col> operator-(const Matrix<t_row, t_col, t_type> &m) const;
    Matrix3<t_type, t_row, t_col> operator+(const Matrix<0, 0, t_type> &m) const;
    Matrix3<t_type, t_row, t_col> operator-(const Matrix<0, 0, t_type> &m) const;
    Matrix3<t_type, t_row, t_col> operator*(const t_type s) const;
    Matrix3<t_type, t_row, t_col> operator/(const t_type s) const;

    template <uint16_t col>
    Matrix<t_row, col, t_type> operator*(const Matrix<t_col, col, t_type> &m) const;       // RotMat * matrix
    Matrix<0, 0, t_type> operator*(const Matrix<0, 0, t_type> &m) const;                   // RotMat * matrix
    Matrix3<t_type, t_row, t_col> operator*(const Matrix3<t_type, t_row, t_col> &m) const; // RotMat * matrix
    Rotation operator*(const Rotation &m) const;                                           // RotMat * RotMat
    Vector<t_row, t_type> operator*(const Vector<t_col, t_type> &v) const;                 // RotMat * vector
    Vector<t_row, t_type> operator*(const Vector<0, t_type> &v) const;                     // RotMat * vector
    Vector3<t_type, t_row> operator*(const Vector3<t_type, t_col> &v) const;               // RotMat * vector
    Rotation operator&(const Vector<t_col, t_type> &v) const;                              // RotMat * [v]x, []x is skew-symmetric matrix
    Rotation operator&(const Vector<0, t_type> &v) const;                                  // RotMat * [v]x, []x is skew-symmetric matrix
    Rotation operator&(const Vector3<t_type, t_col> &v) const;                             // RotMat * [v]x, []x is skew-symmetric matrix

    /* Comparison operators */
    bool operator==(const Rotation &m) const;                      // (true or false) matrix == matrix
    bool operator!=(const Rotation &m) const;                      // (true or false) matrix != matrix
    bool operator==(const Matrix3<t_type, t_row, t_col> &m) const; // (true or false) matrix == matrix
    bool operator!=(const Matrix3<t_type, t_row, t_col> &m) const; // (true or false) matrix != matrix
    bool operator==(const Matrix<t_row, t_col, t_type> &m) const;  // (true or false) matrix == matrix
    bool operator!=(const Matrix<t_row, t_col, t_type> &m) const;  // (true or false) matrix != matrix
    bool operator==(const Matrix<0, 0, t_type> &m) const;          // (true or false) matrix == matrix
    bool operator!=(const Matrix<0, 0, t_type> &m) const;          // (true or false) matrix != matrix

    void Print(const char endChar = 0);

    /* Friend classes */
    template <typename type, uint16_t row, uint16_t col> friend class Rotation;
    template <typename type, uint16_t row, uint16_t col> friend class Matrix3;
    template <typename type, uint16_t row, uint16_t col> friend class Transform;
    template <uint16_t row, uint16_t col, typename type> friend class Matrix;

    template <typename type, uint16_t row> friend class Quaternion;

    /* Friend template function */
    template <typename type, uint16_t row, uint16_t col>
    friend Matrix3<type, row, col> operator*(const type s, const Rotation<type, row, col> &m); // scalar * RotMat
};

} // namespace Math
} // namespace dt

#include "dtRotation.tpp"

#endif // DTMATH_DTROTATION_H_
