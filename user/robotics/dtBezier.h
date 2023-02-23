
#ifndef DTMATH_DTBEZIER_H_
#define DTMATH_DTBEZIER_H_

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

template <uint16_t m_ctrlNum, typename m_type = float>
class CdtBezier
{
public:
    CdtVector3<m_type> pos_x;   // unit:[x]
    CdtVector3<m_type> vel_xps; // unit:[x/sec]

private:
    m_type m_smpTime_ms;
    m_type m_period_ms;
    m_type m_elapsed_ms;                    // sum of sampling time
    m_type m_ctrlParam;                     // m_elapsed_ms / m_period_ms, range: 0 ~ 1
    m_type m_binomialCoeffN1[m_ctrlNum];    // pre-calculated binomial coefficient about n = 15, k is index of the array
    m_type m_binomialCoeffN2[m_ctrlNum-1];  // pre-calculated binomial coefficient about n = 14, k is index of the array
    CdtVector3<m_type> m_ctrlPts[m_ctrlNum];

public:
    CdtBezier();
    int8_t Init(const uint32_t sampleTime_us, const m_type period_ms);
    int8_t SetPeriod(const m_type period_ms);
    int8_t SetCtrlPoints(const CdtVector3<m_type> *ctrlPts);
    void Compute(); // Compute Bezier Poly with internel elapsed time and period
    void Compute(m_type ctrlParam); // Compute Bezier Poly with ctrlParam and period

private:
    m_type BinomialCoeff(int n, int k);
    m_type BernsteinPoly(int n, int k, m_type ctrlParam);
};

#include "dtBezier.tpp"

#endif // DTMATH_DTBEZIER_H_
