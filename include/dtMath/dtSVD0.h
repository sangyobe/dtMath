/*!
\file       dtSVD0.h
\brief      dtMath, Singular Value Decomposition solver class for dynamic memory allocation
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 05. 22
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTSVD0_H_
#define DTMATH_DTSVD0_H_

#include "dtDefine.h"
#include "dtSVD.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <string.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

/*
This procedure programmed here is based on the method of Golub and Reinsch as given on
pages 134 - 151 of the "Handbook for Automatic Computation vol II - Linear Algebra"
edited by Wilkinson and Reinsch and published by Springer - Verlag, 1971.
*/

/*
This routine decomposes an mxn matrix A, with m >= n,
into a product of the three matrices U, S, and VT,

i.e. A = U S VT,

where U is an mxn matrix whose columns are mutually orthogonal,
      S is an nxn diagonal matrix,
      V is an nxn orthogonal matrix.
      VT denotes the transpose of V.

If m < n, then the procedure may be used for the matrix AT.
The singular values of A are the diagonal elements of the diagonal matrix S
and correspond to the positive square roots of the eigenvalues of the matrix ATA.
*/

namespace dt
{
namespace Math
{

template <typename t_type>
class SVD<0, 0, t_type>
{

private:
    uint16_t m_row, m_col, m_size;

    const int MAX_ITERATION_COUNT = 30;
    t_type *m_A;
    t_type *m_U;
    t_type *m_S;
    t_type *m_V;
    // t_type m_superDiagonal[(t_row >= t_col) ? t_col : t_row]; // This array is used to store the super-diagonal elements resulting from the Householder reduction of the matrix A to bidiagonal form.
    t_type *m_superDiagonal; // This array is used to store the super-diagonal elements resulting from the Householder reduction of the matrix A to bidiagonal form.
    t_type *m_inv;
    int8_t m_isOk;

public:
    SVD();
    SVD(const uint16_t row, const uint16_t col, const t_type *element);
    template <uint16_t row, uint16_t col>
    SVD(const Matrix<row, col, t_type> &m);
    SVD(const Matrix3<t_type, 3, 3> &m);
    SVD(const Matrix<0, 0, t_type> &m);
    ~SVD();

    int8_t Compute(const uint16_t row, const uint16_t col, const t_type *element); // Compute Lower/Upper Triangular Matrix, Doolittle form
    template <uint16_t row, uint16_t col>
    int8_t Compute(const Matrix<row, col, t_type> &m); // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const Matrix3<t_type, 3, 3> &m);    // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t Compute(const Matrix<0, 0, t_type> &m);     // Compute Lower/Upper Triangular Matrix, Doolittle form
    int8_t IsOk() { return m_isOk; }

    Matrix<0, 0, t_type> GetMatrixU() const; // return matrix A including L/U matrix
    Matrix<0, 0, t_type> GetMatrixS() const; // return Lower Triangular matrix
    Matrix<0, 0, t_type> GetMatrixV() const; // return Upper Triangular matrix

    int8_t Solve(const Matrix<0, 0, t_type> &b, Matrix<0, 0, t_type> &x, t_type tolerance = 0); // Solve x = (USVT)^-1 * b
    int8_t Solve(const Vector<0, t_type> &b, Vector<0, t_type> &x, t_type tolerance = 0);       // Solve x = (USVT)^-1 * b

    Matrix<0, 0, t_type> Solve(const Matrix<0, 0, t_type> &b, int8_t *isOk = nullptr, t_type tolerance = 0); // Solve x = (USVT)^-1 * b
    template <uint16_t row>
    Vector<0, t_type> Solve(const Vector<row, t_type> &b, int8_t *isOk = nullptr, t_type tolerance = 0); // Solve x = (LU)^-1 * b
    Vector<0, t_type> Solve(const Vector<0, t_type> &b, int8_t *isOk = nullptr, t_type tolerance = 0);   // Solve x = (USVT)^-1 * b

    int8_t Inverse(Matrix<0, 0, t_type> &inv, t_type tolerance = 0); // Inverse matrix of USVT matrix
    template <uint16_t row, uint16_t col>
    int8_t Inverse(Matrix<row, col, t_type> &inv, t_type tolerance = 0);        // Inverse matrix of USVT matrix
    Matrix<0, 0, t_type> Inverse(int8_t *isOk = nullptr, t_type tolerance = 0); // Inverse matrix of USVT matrix

private:
    void HouseholdersReductionToBidiagonalForm_Mxn0();
    int8_t GivensReductionToDiagonalForm_Mxn0();
    void SortByDecreasingSingularValues_Mxn0();

    void HouseholdersReductionToBidiagonalForm_mxN0();
    int8_t GivensReductionToDiagonalForm_mxN0();
    void SortByDecreasingSingularValues_mxN0();
};

} // namespace Math
} // namespace dt

#include "dtSVD0.tpp"

#endif // DTMATH_DTSVD0_H_