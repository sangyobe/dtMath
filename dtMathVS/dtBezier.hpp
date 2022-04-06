#pragma once

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include "./dtMath/dtMath.h"
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

template<uint16_t m_ctrlNum, typename m_type>
inline CdtBezier<m_ctrlNum, m_type>::CdtBezier()
{
    m_smpTime_ms = 1;
    m_period_ms = -1;
    m_elapsed_ms = 2;
    m_ctrlParam = 2;

    for (int i = 0; i < (m_ctrlNum-1); i++)
    {
        m_binomialCoeffN1[i] = BinomialCoeff(m_ctrlNum-1, i);
        m_binomialCoeffN2[i] = BinomialCoeff(m_ctrlNum-2, i);
    }
    m_binomialCoeffN1[m_ctrlNum-1] = BinomialCoeff(m_ctrlNum - 1, m_ctrlNum - 1);
}

template<uint16_t m_ctrlNum, typename m_type>
inline int8_t CdtBezier<m_ctrlNum, m_type>::Init(const uint32_t sampleTime_us, const m_type period_ms)
{
    if (sampleTime_us <= 0)
    {
        m_smpTime_ms = 1;
        return -1;
    }

    if (period_ms <= 0)
    {
        m_period_ms = -1;
        return -1;
    }

    m_smpTime_ms = sampleTime_us / 1000.0f;
    m_period_ms = period_ms;
    m_elapsed_ms = 2;
    m_ctrlParam = 2;

    return 0;
}

template<uint16_t m_ctrlNum, typename m_type>
inline int8_t CdtBezier<m_ctrlNum, m_type>::SetPeriod(const m_type period_ms)
{
    if (period_ms <= 0)
    {
        m_period_ms = -1;
        return -1;
    }

    m_period_ms = period_ms;

    return 0;
}

template<uint16_t m_ctrlNum, typename m_type>
inline int8_t CdtBezier<m_ctrlNum, m_type>::SetCtrlPoints(const CdtVector3<m_type>* ctrlPts)
{
    for (uint8_t i = 0; i < m_ctrlNum; i++)
    {
        m_ctrlPts[i] = ctrlPts[i];
    }
    
    return 0;
}

template<uint16_t m_ctrlNum, typename m_type>
inline void CdtBezier<m_ctrlNum, m_type>::Compute()
{ // Compute BezierPoly with elapsed time and period
    m_ctrlParam = m_elapsed_ms / m_period_ms;

    if (m_ctrlParam > 1 || m_ctrlParam < 0) return;

    pos_x.SetZero();
    vel_xps.SetZero();

    for (int i = 0; i < (m_ctrlNum - 1); i++)
    {
        pos_x += BernsteinPoly(m_ctrlNum - 1, i, m_ctrlParam) * m_ctrlPts[i];
        vel_xps += BernsteinPoly(m_ctrlNum - 2, i, m_ctrlParam) * (m_ctrlPts[i + 1] - m_ctrlPts[i]);
    }

    pos_x += BernsteinPoly(m_ctrlNum - 1, m_ctrlNum - 1, m_ctrlParam) * m_ctrlPts[m_ctrlNum - 1];
    vel_xps *= ((m_ctrlNum - 1) * 1000 / m_period_ms);

    m_elapsed_ms += m_smpTime_ms;
}

template<uint16_t m_ctrlNum, typename m_type>
inline void CdtBezier<m_ctrlNum, m_type>::Compute(m_type ctrlParam)
{ // Compute BezierPoly with control parameter
    if (ctrlParam > 1 || ctrlParam < 0) return;

    pos_x.SetZero();
    vel_xps.SetZero();

    for (int i = 0; i < (m_ctrlNum - 1); i++)
    {
        pos_x += BernsteinPoly(m_ctrlNum - 1, i, ctrlParam) * m_ctrlPts[i];
        vel_xps += BernsteinPoly(m_ctrlNum - 2, i, ctrlParam) * (m_ctrlPts[i + 1] - m_ctrlPts[i]);
    }

    pos_x += BernsteinPoly(m_ctrlNum - 1, m_ctrlNum - 1, ctrlParam) * m_ctrlPts[m_ctrlNum - 1];
    vel_xps *= ((m_ctrlNum - 1) * 1000 / m_period_ms);
}

template<uint16_t m_ctrlNum, typename m_type>
inline m_type CdtBezier<m_ctrlNum, m_type>::BinomialCoeff(int n, int k)
{
    // Computing the value of binomial coefficients
    // (n, k)b = n^{\underset{-}{k}} / k!        if 2k <= n
    //         = n^{\underset{-}{n-k}} / (n-k)!  if 2k > n

    if (n < k) return 0;

    int min, max;
    m_type bc = 1;

    if ((n - k) >= k)
    {
        max = n - k;
        min = k;
    }
    else
    {
        max = k;
        min = n - k;
    }

    for (int i = 1; i <= min; i++)
    {
        bc *= ((m_type)(max + i) / (m_type)i);
    }

    return bc;
}

template<uint16_t m_ctrlNum, typename m_type>
inline m_type CdtBezier<m_ctrlNum, m_type>::BernsteinPoly(int n, int k, m_type x)
{
    // Bernstein basis polynomials
    // b_{k,n}(x) = (n, k)b * x^{k} * (1-x)^{n-k}
    // (,)b is binomial coefficient

    //return BinomialCoeff(n, k) * std::pow(x, k) * std::pow(1 - x, n - k);

    switch (n)
    {
    case (m_ctrlNum - 2):
        return m_binomialCoeffN2[k] * std::pow(x, k) * std::pow(1 - x, n - k);
    case (m_ctrlNum - 1):
    default:
        return m_binomialCoeffN1[k] * std::pow(x, k) * std::pow(1 - x, n - k);
    }
}
