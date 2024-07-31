/*!
\file       dtFullPivLU.h
\brief      dtMath, LU Decomposition with partial pivoting(Doolittle form) class for Dynamic Allocation
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 05. 14
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTFULL_PIV_LU0_H_
#define DTMATH_DTFULL_PIV_LU0_H_

#include "dtDefine.h"
#include "dtFullPivLU.h"

#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
#include <stdint.h>
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
class FullPivLU<0, 0, t_type>
{
private:
    uint16_t m_row, m_col, m_size;
    t_type *m_elem;
    t_type *m_inv;
    int *m_pivot, *m_pivotCol, *m_permCol;
    int8_t m_isOk;
    int8_t Compute(); // Compute Lower/Upper Triangular Matrix, Doolittle form

public:
    FullPivLU();
    FullPivLU(const uint16_t row, const uint16_t col, const t_type *element);
    template <uint16_t row, uint16_t col>
    FullPivLU(const Matrix<row, col, t_type> &m);
    FullPivLU(const Matrix3<t_type, 3, 3> &m);
    FullPivLU(const Matrix<0, 0, t_type> &m);
    ~FullPivLU();

    int8_t Compute(const uint16_t row, const uint16_t col, const t_type *element); // Compute Lower/Upper Triangular Matrix, Doolittle form
    template <uint16_t row, uint16_t col>
    int8_t Compute(const Matrix<row, col, t_type> &m); // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const Matrix3<t_type, 3, 3> &m);    // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const Matrix<0, 0, t_type> &m);     // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t IsOk() { return m_isOk; }

    Matrix<0, 0, t_type> GetMatrix() const;  // return matrix A including L/U matrix
    Matrix<0, 0, t_type> GetMatrixL() const; // return Lower Triangular matrix 'L'
    Matrix<0, 0, t_type> GetMatrixU() const; // return Upper Triangular matrix 'U'
    Matrix<0, 0, t_type> GetMatrixP() const; // return Permutation matrix 'P'
    Matrix<0, 0, t_type> GetMatrixQ() const; // return Permutation matrix 'Q'

    template <uint16_t row, uint16_t col>
    int8_t Solve(const Vector<row, t_type> &b, Vector<col, t_type> &x); // Solve x = (LU)^-1 * b
    int8_t Solve(const Vector<0, t_type> &b, Vector<0, t_type> &x);     // Solve x = (LU)^-1 * b
    template <uint16_t row>
    Vector<0, t_type> Solve(const Vector<row, t_type> &b, int8_t *isOk = nullptr); // Solve x = (LU)^-1 * b
    Vector<0, t_type> Solve(const Vector<0, t_type> &b, int8_t *isOk = nullptr);   // Solve x = (LU)^-1 * b

    template <uint16_t row, uint16_t col>
    int8_t Inverse(Matrix<row, col, t_type> &inv);        // Inverse matrix of LU matrix
    int8_t Inverse(Matrix<0, 0, t_type> &inv);            // Inverse matrix of LU matrix
    int8_t Inverse(Matrix3<t_type, 3, 3> &inv);           // Inverse matrix of LU matrix
    Matrix<0, 0, t_type> Inverse(int8_t *isOk = nullptr); // Inverse matrix of LU matrix
};

} // namespace Math
} // namespace dt

#include "dtFullPivLU0.tpp"

#endif // DTMATH_DTFULL_PIV_LU0_H_
