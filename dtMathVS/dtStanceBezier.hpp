/*!
\file       dtStanceBezier.hpp
\brief      Bezier trajactory used during the stance period of the legged robot.
\author     Dongjin Hyun, mecjin@gmail.com
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2021. 07. 01
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

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

template <typename m_type = float>
class CdtStanceBezier
{
public:
    CdtVector3<m_type> pos_x;   // unit:[x]
    CdtVector3<m_type> vel_xps; // unit:[x/sec]

private:
    m_type m_smpTime_ms;
    m_type m_period_ms;
    m_type m_elapsed_ms;            // sum of sampling time
    m_type m_ctrlParam;             // m_elapsed_ms / m_period_ms, range: 0 ~ 1
    m_type m_binomialCoeff7[8];   // pre-calculated binomial coefficient about n = 7, k is index of the array
    m_type m_binomialCoeff6[7];   // pre-calculated binomial coefficient about n = 6, k is index of the array
    CdtVector3<m_type> m_ctrlPts[8];

public:
    CdtStanceBezier();
    int8_t Init(const uint32_t sampleTime_us, const m_type period_ms);
    int8_t SetPeriod(const m_type period_ms);
    int8_t SetTarget(
        CdtVector3<m_type> &startPos, CdtVector3<m_type> &endPos,
        CdtVector3<m_type> &startVel, CdtVector3<m_type> &endVel,
        const m_type period_ms);
    int8_t SetTarget(
        CdtVector3<m_type> &startPos, CdtVector3<m_type> &endPos,
        CdtVector3<m_type> &startVel, CdtVector3<m_type> &endVel);

    CdtVector3<m_type>& CtrlPts(const uint8_t idx) { return m_ctrlPts[idx]; }

    void Compute(); // Compute Bezier Poly with internel elapsed time and period
    void Compute(m_type ctrlParam); // Compute Bezier Poly with ctrlParam and period

private:
    m_type BinomialCoeff(int n, int k);
    m_type BernsteinPoly(int n, int k, m_type ctrlParam);
};

template<typename m_type>
inline CdtStanceBezier<m_type>::CdtStanceBezier()
{
    m_smpTime_ms = 1;
    m_period_ms = -1;
    m_elapsed_ms = 2;
    m_ctrlParam = 2;

    for (int i = 0; i < 7; i++)
    {
        m_binomialCoeff7[i] = BinomialCoeff(7, i);
        m_binomialCoeff6[i] = BinomialCoeff(6, i);
    }
    m_binomialCoeff7[7] = BinomialCoeff(7, 7);

    CdtVector3<m_type> zeroValue;
    SetTarget(zeroValue, zeroValue, zeroValue, zeroValue);
}

template<typename m_type>
inline int8_t CdtStanceBezier<m_type>::Init(const uint32_t sampleTime_us, const m_type period_ms)
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

    CdtVector3<m_type> zeroValue;
    SetTarget(zeroValue, zeroValue, zeroValue, zeroValue);

    return 0;
}

template<typename m_type>
inline int8_t CdtStanceBezier<m_type>::SetPeriod(const m_type period_ms)
{
    if (period_ms <= 0)
    {
        m_period_ms = -1;
        return -1;
    }

    m_period_ms = period_ms;

    return 0;
}

template<typename m_type>
inline int8_t CdtStanceBezier<m_type>::SetTarget(
    CdtVector3<m_type>& startPos, CdtVector3<m_type>& endPos,
    CdtVector3<m_type>& startVel, CdtVector3<m_type>& endVel,
    const m_type period_ms)
{ // Make Control Points
    if (period_ms <= std::numeric_limits<m_type>::epsilon())
    {
        m_period_ms = -1;
        m_elapsed_ms = 2;
        return -1;
    }

    m_period_ms = period_ms;

    CdtVector3<m_type> startDel;
    CdtVector3<m_type> center;
    CdtVector3<m_type> endDel;

    m_ctrlPts[0] = startPos;

    startDel = startVel * (m_period_ms / 7000.0f); // unit change from ms to sec
    m_ctrlPts[1] = startPos + startDel;
    m_ctrlPts[2] = startPos + (startDel * 2.0f);

    center = (endPos + startPos) * 0.5f;
    m_ctrlPts[3] = center;
    m_ctrlPts[4] = center;

    endDel = endVel * (m_period_ms / 7000.0f); // unit change from ms to sec
    m_ctrlPts[5] = endPos - (endDel * 2.0f);
    m_ctrlPts[6] = endPos - endDel;

    m_ctrlPts[7] = endPos;

    m_elapsed_ms = 0;

    return 0;
}

template<typename m_type>
inline int8_t CdtStanceBezier<m_type>::SetTarget(
    CdtVector3<m_type>& startPos, CdtVector3<m_type>& endPos,
    CdtVector3<m_type>& startVel, CdtVector3<m_type>& endVel)
{ // Make Control Points
    if (m_period_ms < 0)
    {
        m_elapsed_ms = 2;
        return -1;
    }

    CdtVector3<m_type> startDel;
    CdtVector3<m_type> center;
    CdtVector3<m_type> endDel;

    m_ctrlPts[0] = startPos;

    startDel = startVel * (m_period_ms / 7000.0f); // unit change from ms to sec
    m_ctrlPts[1] = startPos + startDel;
    m_ctrlPts[2] = startPos + (startDel * 2.0f);

    center = (endPos + startPos) * 0.5f;

    m_ctrlPts[3] = center;
    m_ctrlPts[4] = center;

    endDel = endVel * (m_period_ms / 7000.0f); // unit change from ms to sec
    m_ctrlPts[5] = endPos - (endDel * 2.0f);
    m_ctrlPts[6] = endPos - endDel;

    m_ctrlPts[7] = endPos;

    m_elapsed_ms = 0;

    return 0;
}

template<typename m_type>
inline void CdtStanceBezier<m_type>::Compute()
{ // Compute BezierPoly with elapsed time and period
    m_ctrlParam = m_elapsed_ms / m_period_ms;

    if (m_ctrlParam > 1 || m_ctrlParam < 0) return;

    pos_x.SetZero();
    vel_xps.SetZero();

    for (int i = 0; i < 7; i++)
    {
        pos_x += BernsteinPoly(7, i, m_ctrlParam) * m_ctrlPts[i];
        vel_xps += BernsteinPoly(6, i, m_ctrlParam) * (m_ctrlPts[i + 1] - m_ctrlPts[i]);
    }

    pos_x += BernsteinPoly(7, 7, m_ctrlParam) * m_ctrlPts[7];
    vel_xps *= (7000 / m_period_ms);

    m_elapsed_ms += m_smpTime_ms;
}

template<typename m_type>
inline void CdtStanceBezier<m_type>::Compute(m_type ctrlParam)
{ // Compute BezierPoly with control parameter
    if (ctrlParam > 1 || ctrlParam < 0) return;

    pos_x.SetZero();
    vel_xps.SetZero();

    for (int i = 0; i < 7; i++)
    {
        pos_x += BernsteinPoly(7, i, ctrlParam) * m_ctrlPts[i];
        vel_xps += BernsteinPoly(6, i, ctrlParam) * (m_ctrlPts[i + 1] - m_ctrlPts[i]);
    }

    pos_x += BernsteinPoly(7, 7, ctrlParam) * m_ctrlPts[7];
    vel_xps *= (7000 / m_period_ms);
}

template<typename m_type>
inline m_type CdtStanceBezier<m_type>::BinomialCoeff(int n, int k)
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

template<typename m_type>
inline m_type CdtStanceBezier<m_type>::BernsteinPoly(int n, int k, m_type x)
{
    // Bernstein basis polynomials
    // b_{k,n}(x) = (n, k)b * x^{k} * (1-x)^{n-k}
    // (,)b is binomial coefficient

    //return BinomialCoeff(n, k) * std::pow(x, k) * std::pow(1 - x, n - k);

    switch (n)
    {
    case 6:
        return m_binomialCoeff6[k] * std::pow(x, k) * std::pow(1 - x, n - k);
    case 7:
    default:
        return m_binomialCoeff7[k] * std::pow(x, k) * std::pow(1 - x, n - k);
    }
}
