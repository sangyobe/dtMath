/*!
\file       dtUpperTriangular.hpp
\brief      dtMath, Upper triangular matrix solver class
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
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

template <uint16_t m_row, uint16_t m_col, typename m_type = float>
class CdtUpperTriangular
{
public:
    CdtUpperTriangular() {}

    // for General Upper Triangular Matrix
    static int8_t Solve(CdtMatrix<m_row, m_col, m_type> &U, CdtVector<m_row, m_type> &b, CdtVector<m_col, m_type> &x);
    static CdtVector<m_col, m_type> Solve(CdtMatrix<m_row, m_col, m_type> &U, CdtVector<m_row, m_type> &b, int8_t *isOk = nullptr);
    static int8_t Inverse(CdtMatrix<m_row, m_col, m_type> U, CdtMatrix<m_row, m_col, m_type>& invU);
    static CdtMatrix<m_row, m_col, m_type> Inverse(CdtMatrix<m_row, m_col, m_type> U, int8_t *isOk = nullptr);

    // for Unit Upper Triangular Matrix
    static int8_t SolveUnit(CdtMatrix<m_row, m_col, m_type> &U, CdtVector<m_row, m_type> &b, CdtVector<m_col, m_type> &x);
    static CdtVector<m_col, m_type> SolveUnit(CdtMatrix<m_row, m_col, m_type> &U, CdtVector<m_row, m_type> &b, int8_t *isOk = nullptr);
    static int8_t InverseUnit(CdtMatrix<m_row, m_col, m_type> U, CdtMatrix<m_row, m_col, m_type>& invU);
    static CdtMatrix<m_row, m_col, m_type> InverseUnit(CdtMatrix<m_row, m_col, m_type> U, int8_t *isOk = nullptr);
};

#include "dtUpperTriangular.ipp"
