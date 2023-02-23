/*!
\file       djBezier.h
\brief      Bezier class
\author     Dong-jin Hyun
\editor     Dong-hyun Lee, phenom8305@gmail.com
\date       2021. 04. 19
\version    1.0.0
*/

#ifndef DTMATH_DJBEZIER_H_
#define DTMATH_DJBEZIER_H_

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <limits>

#include <dtMath/dtMath.h>

// ncol: space dimension
// mrow: the number of control points

//template <int ncol, typename m_type = float>
//struct BezierOut
//{
//    CdtVector<ncol, m_type> pos;
//    CdtVector<ncol, m_type> vel;
//};

template <uint16_t mrow, uint16_t ncol, uint16_t order = 8, typename m_type = float>
class CdjBezier
{
public:
    CdtVector<ncol, m_type> pos;
    CdtVector<ncol, m_type> vel;

    m_type BinomialCoeff(int n, int k);
    m_type BernsteinPoly(int n, int k, m_type ctrlParam);
    void BezierPoly(CdtMatrix<order, ncol, m_type> &ctrlPoints, m_type ctrlParam, m_type period);
    void BezierInterp(CdtMatrix<mrow, ncol, m_type> &ctrlPoints, m_type ctrlParam, m_type period);
};

#include "djBezier.tpp"

#endif // DTMATH_DJBEZIER_H_
