/*!
\file       dtSVD.h
\brief      dtMath, Singular Value Decomposition solver class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTSVD_H_
#define DTMATH_DTSVD_H_

#include "dtDefine.h"

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

template <uint16_t t_row, typename t_type> class Vector;
template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;
template <typename t_type, uint16_t t_row, uint16_t t_col> class Matrix3;

template <uint16_t t_row, uint16_t t_col, typename t_type = float>
class SVD
{
private:
    const int MAX_ITERATION_COUNT = 30;
    t_type m_A[t_row * t_col];
    t_type m_U[(t_row >= t_col) ? t_row * t_col : t_col * t_col];
    t_type m_S[(t_row >= t_col) ? t_col : t_row];
    t_type m_V[(t_row >= t_col) ? t_col * t_col : t_row * t_col];
    // t_type m_superDiagonal[(t_row >= t_col) ? t_col : t_row]; // This array is used to store the super-diagonal elements resulting from the Householder reduction of the matrix A to bidiagonal form.
    t_type m_superDiagonal[(t_row >= t_col) ? t_row : t_col]; // This array is used to store the super-diagonal elements resulting from the Householder reduction of the matrix A to bidiagonal form.
    t_type m_inv[t_row * t_col];
    int8_t m_isOk;  // SVD solved
    int8_t m_isInv; // inverse matrix generated

public:
    SVD();
    SVD(const t_type *element, const size_t n_byte);
    SVD(const Matrix<t_row, t_col, t_type> &m);
    SVD(const Matrix3<t_type, t_row, t_col> &m);

    int8_t Compute(const t_type *element, const size_t n_byte); // Compute Singular Value Decomposition
    int8_t Compute(const Matrix<t_row, t_col, t_type> &m);      // Compute Singular Value Decomposition
    int8_t Compute(const Matrix3<t_type, t_row, t_col> &m);     // Compute Singular Value Decomposition
    int8_t IsOk() { return m_isOk; }

    Matrix<t_row, t_row, t_type> GetMatrixU() const; // return U matrix, m x m, left singular vectors of A
    Matrix<t_row, t_col, t_type> GetMatrixS() const; // return S matrix, m x n, singular values of A
    Matrix<t_col, t_col, t_type> GetMatrixV() const; // return V matrix, n x n, right singular vectors of A

    template <uint16_t col>
    int8_t Solve(const Matrix<t_row, col, t_type> &b, Matrix<t_col, col, t_type> &x, t_type tolerance = 0); // Solve x = (USVT)^-1 * b
    int8_t Solve(const Vector<t_row, t_type> &b, Vector<t_col, t_type> &x, t_type tolerance = 0);           // Solve x = (USVT)^-1 * b
    template <uint16_t col>
    Matrix<t_col, col, t_type> Solve(const Matrix<t_row, col, t_type> &b, int8_t *isOk = nullptr, t_type tolerance = 0); // Solve x = (USVT)^-1 * b
    Vector<t_col, t_type> Solve(const Vector<t_row, t_type> &b, int8_t *isOk = nullptr, t_type tolerance = 0);           // Solve x = (USVT)^-1 * b

    int8_t Inverse(Matrix<t_col, t_row, t_type> &inv, t_type tolerance = 0);            // Inverse matrix of USVT matrix
    int8_t Inverse(Matrix3<t_type, t_col, t_row> &inv, t_type tolerance = 0);           // Inverse matrix of USVT matrix
    Matrix<t_col, t_row, t_type> Inverse(int8_t *isOk = nullptr, t_type tolerance = 0); // Inverse matrix of USVT matrix

    int8_t InverseArray(t_type *inv, t_type tolerance = 0);             // Inverse array of USVT matrix
    t_type *InverseArray(int8_t *isOk = nullptr, t_type tolerance = 0); // Inverse array of USVT matrix

private:
    void HouseholdersReductionToBidiagonalForm_Mxn();
    int8_t GivensReductionToDiagonalForm_Mxn();
    void SortByDecreasingSingularValues_Mxn();

    void HouseholdersReductionToBidiagonalForm_mxN();
    int8_t GivensReductionToDiagonalForm_mxN();
    void SortByDecreasingSingularValues_mxN();
};

} // namespace Math
} // namespace dt

#include "dtSVD.tpp"
#include "dtSVD0.h"

#endif // DTMATH_DTSVD_H_
