/*!
\file       dtSwingBezier.h
\brief      Bezier trajactory used during the swing period of the legged robot.
\author     Dongjin Hyun, mecjin@gmail.com
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2021. 07. 01
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTSWING_BEZIER_H_
#define DTMATH_DTSWING_BEZIER_H_

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <dtMath/dtMath.h>

#include <cmath>
#include <limits>

template <typename m_type = float>
class CdtSwingBezier
{
public:
    CdtVector3<m_type> pos_x;   // unit:[x]
    CdtVector3<m_type> vel_xps; // unit:[x/sec]

private:
    m_type m_smpTime_ms;
    m_type m_period_ms;
    m_type m_elapsed_ms;            // sum of sampling time
    m_type m_ctrlParam;             // m_elapsed_ms / m_period_ms, range: 0 ~ 1
    m_type m_binomialCoeff15[16];   // pre-calculated binomial coefficient about n = 15, k is index of the array
    m_type m_binomialCoeff14[15];   // pre-calculated binomial coefficient about n = 14, k is index of the array
    CdtVector3<m_type> m_ctrlPts[16];

public:
    CdtSwingBezier();
    int8_t Init(const uint32_t sampleTime_us, const m_type period_ms);
    int8_t SetPeriod(const m_type period_ms);
    int8_t SetTarget(
        CdtVector3<m_type> &startPos, CdtVector3<m_type> &endPos,
        CdtVector3<m_type> &startVel, CdtVector3<m_type> &endVel,
        const m_type height, const m_type period_ms);
    int8_t SetTarget(
        CdtVector3<m_type> &startPos, CdtVector3<m_type> &endPos,
        CdtVector3<m_type> &startVel, CdtVector3<m_type> &endVel,
        const m_type height);
    CdtVector3<m_type>& CtrlPts(const uint8_t idx) { return m_ctrlPts[idx]; }

    void Compute(); // Compute Bezier Poly with internel elapsed time and period
    void Compute(m_type ctrlParam); // Compute Bezier Poly with ctrlParam and period

private:
    m_type BinomialCoeff(int n, int k);
    m_type BernsteinPoly(int n, int k, m_type ctrlParam);
};

#include "dtSwingBezier.tpp"

#endif // DTMATH_DTSWING_BEZIER_H_
