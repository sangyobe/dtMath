/*!
\file      dtCscAMD.h
\brief      dtMath, Cholesky decomposition(L*L^T form) for Sparse Matrix Class
\author     Muhammad Zahak Jamal, zahakj@gmail.com
\author     Who is next author?
\date      2023. 06. 10
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTCSC_AMD_H_
#define DTMATH_DTCSC_AMD_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
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

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class CscAMD
{
private:
    int m_elemNum = 0;            // no. of non zero elements
    m_type m_elem[m_row * m_col]; // non zero elements of the sparse matrix
    int m_rowIdx[m_row * m_col];  // row indices
    int m_colPtr[m_col + 1];      // column index pointer

    int flag = 0;
    int tdfs(int j, int k, int *head, int *next, int *post, int *stck);
    int P[m_col + 1] = {
        0,
    }; // Permutation vector
    int Pinv[m_col] = {
        0,
    }; // Inverse permutation vector
    int wclear(int mark, int lemax, int *w);
    int Compute(); // Compute AMD for P vector
    int flip(const int &i) { return -i - 2; }
    void returnPinv(); // Compute the Pinv vector

public:
    CscAMD();                                         // Default constructor
    CscAMD(const CscMatrix<m_row, m_col, m_type> &m); // Constructor with sparse matrix as argument

    int AMDOrder(const CscMatrix<m_row, m_col, m_type> &m); // Find P and Pinv Vector

    Vector<m_col, int> GetP();                                              // Get P in  Vector
    Vector<m_col, int> GetPinv();                                           // Get Pinv in  Vector
    Vector<m_col, m_type> PermuteVector(const Vector<m_col, m_type> &b);    // Permute vector (b_perm = P * b)
    Vector<m_col, m_type> InvPermuteVector(const Vector<m_col, m_type> &b); // Inverse Permute the vector (b_invperm = Pinv * b)
    CscMatrix<m_row, m_col, m_type> GetPermutedMatrix();                    // Apply permutation to square matrix (B = PAPT)

    template <uint16_t row, uint16_t col, typename type> friend class CscLLT; // Target Friend Classes
    // template <uint16_t row, uint16_t col, typename type> friend class  CscLDLT;
};

} // namespace Math
} // namespace dt

#include "dtCscAMD.tpp"

#endif // DTMATH_DTCSC_AMD_H_