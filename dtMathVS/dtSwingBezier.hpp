/*!
\file       dtSwingBezier.hpp
\brief      Bezier trajactory used during the swing period of the legged robot.
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

template<typename m_type>
inline CdtSwingBezier<m_type>::CdtSwingBezier()
{
    m_smpTime_ms = 1;
    m_period_ms = -1;
    m_elapsed_ms = 2;
    m_ctrlParam = 2;

    for (int i = 0; i < 15; i++)
    {
        m_binomialCoeff15[i] = BinomialCoeff(15, i);
        m_binomialCoeff14[i] = BinomialCoeff(14, i);
    }
    m_binomialCoeff15[15] = BinomialCoeff(15, 15);

    CdtVector3<m_type> zeroValue;
    SetTarget(zeroValue, zeroValue, zeroValue, zeroValue, 0);
}

template<typename m_type>
inline int8_t CdtSwingBezier<m_type>::Init(const uint32_t sampleTime_us, const m_type period_ms)
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
    SetTarget(zeroValue, zeroValue, zeroValue, zeroValue, 0);

    return 0;
}

template<typename m_type>
inline int8_t CdtSwingBezier<m_type>::SetPeriod(const m_type period_ms)
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
inline int8_t CdtSwingBezier<m_type>::SetTarget(
    CdtVector3<m_type>& startPos, CdtVector3<m_type>& endPos,
    CdtVector3<m_type>& startVel, CdtVector3<m_type>& endVel,
    const m_type height, const m_type period_ms)
{ // Make Control Points
    if (period_ms <= std::numeric_limits<m_type>::epsilon())
    {
        m_period_ms = -1;
        m_elapsed_ms = 2;
        return -1;
    }

    m_period_ms = period_ms;

    m_type heightBase = (startPos(2) > endPos(2)) ? startPos(2) : endPos(2);
    CdtVector3<m_type> startDel;
    CdtVector3<m_type> delta;
    CdtVector3<m_type> trail;
    CdtVector3<m_type> center;
    CdtVector3<m_type> retraction;
    CdtVector3<m_type> endDel;

    m_ctrlPts[0] = startPos;

    startDel = startVel * (m_period_ms / 15000.0f); // unit change from ms to sec
    m_ctrlPts[1] = startPos + startDel;
    m_ctrlPts[2] = startPos + (startDel * 2.0f);

    delta = (endPos - startPos);
    trail = startPos - (0.1f * delta); // drag
    trail(2) = startPos(2) + (0.6f * height);
    m_ctrlPts[3] = trail;
    m_ctrlPts[4] = trail;

    trail(2) = startPos(2) + height;
    m_ctrlPts[5] = trail; // trail2

    center = (endPos + startPos) * 0.5f;
    center(2) = heightBase + height;
    m_ctrlPts[6] = center;
    m_ctrlPts[9] = center;

    center(2) = heightBase + 1.2f * height;
    m_ctrlPts[7] = center; // center peak
    m_ctrlPts[8] = center; // center peak

    retraction = endPos + (0.3f * delta);
    retraction(2) = endPos(2) + (0.5f * height);
    m_ctrlPts[11] = retraction;
    m_ctrlPts[12] = retraction;

    retraction(2) = endPos(2) + (1.1f * height);
    m_ctrlPts[10] = retraction; // retraction2

    endDel = endVel * (m_period_ms / 15000.0f); // unit change from ms to sec
    m_ctrlPts[13] = endPos - (endDel * 2.0f);
    m_ctrlPts[14] = endPos - endDel;

    m_ctrlPts[15] = endPos;

    m_elapsed_ms = 0;

    return 0;
}

template<typename m_type>
inline int8_t CdtSwingBezier<m_type>::SetTarget(CdtVector3<m_type>& startPos, CdtVector3<m_type>& endPos, CdtVector3<m_type>& startVel, CdtVector3<m_type>& endVel, const m_type height)
{ // Make Control Points
    if (m_period_ms < 0)
    {
        m_elapsed_ms = 2;
        return -1;
    }

    m_type heightBase = (startPos(2) > endPos(2)) ? startPos(2) : endPos(2);
    CdtVector3<m_type> startDel;
    CdtVector3<m_type> delta;
    CdtVector3<m_type> trail;
    CdtVector3<m_type> center;
    CdtVector3<m_type> retraction;
    CdtVector3<m_type> endDel;

    m_ctrlPts[0] = startPos;

    startDel = startVel * (m_period_ms / 15000.0f); // unit change from ms to sec
    m_ctrlPts[1] = startPos + startDel;
    m_ctrlPts[2] = startPos + (startDel * 2.0f);

    delta = (endPos - startPos);
    trail = startPos - (0.1f * delta); // drag
    trail(2) = startPos(2) + (0.6f * height);
    m_ctrlPts[3] = trail;
    m_ctrlPts[4] = trail;

    trail(2) = startPos(2) + height;
    m_ctrlPts[5] = trail; // trail2

    center = (endPos + startPos) * 0.5f;
    center(2) = heightBase + height;
    m_ctrlPts[6] = center;
    m_ctrlPts[9] = center;

    center(2) = heightBase + 1.2f * height;
    m_ctrlPts[7] = center; // center peak
    m_ctrlPts[8] = center; // center peak

    retraction = endPos + (0.3f * delta);
    retraction(2) = endPos(2) + (0.5f * height);
    m_ctrlPts[11] = retraction;
    m_ctrlPts[12] = retraction;

    retraction(2) = endPos(2) + (1.1f * height);
    m_ctrlPts[10] = retraction; // retraction2

    endDel = endVel * (m_period_ms / 15000.0f); // unit change from ms to sec
    m_ctrlPts[13] = endPos - (endDel * 2.0f);
    m_ctrlPts[14] = endPos - endDel;

    m_ctrlPts[15] = endPos;

    m_elapsed_ms = 0;

    return 0;
}

template<typename m_type>
inline void CdtSwingBezier<m_type>::Compute()
{ // Compute BezierPoly with elapsed time and period
    m_ctrlParam = m_elapsed_ms / m_period_ms;

    if (m_ctrlParam > 1 || m_ctrlParam < 0) return;

    pos_x.SetZero();
    vel_xps.SetZero();

    for (int i = 0; i < 15; i++)
    {
        pos_x += BernsteinPoly(15, i, m_ctrlParam) * m_ctrlPts[i];
        vel_xps += BernsteinPoly(14, i, m_ctrlParam) * (m_ctrlPts[i + 1] - m_ctrlPts[i]);
    }

    pos_x += BernsteinPoly(15, 15, m_ctrlParam) * m_ctrlPts[15];
    vel_xps *= (15000 / m_period_ms);

    m_elapsed_ms += m_smpTime_ms;
}

template<typename m_type>
inline void CdtSwingBezier<m_type>::Compute(m_type ctrlParam)
{ // Compute BezierPoly with control parameter
    if (ctrlParam > 1 || ctrlParam < 0) return;

    pos_x.SetZero();
    vel_xps.SetZero();

    for (int i = 0; i < 15; i++)
    {
        pos_x += BernsteinPoly(15, i, ctrlParam) * m_ctrlPts[i];
        vel_xps += BernsteinPoly(14, i, ctrlParam) * (m_ctrlPts[i + 1] - m_ctrlPts[i]);
    }

    pos_x += BernsteinPoly(15, 15, ctrlParam) * m_ctrlPts[15];
    vel_xps *= (15000 / m_period_ms);
}

template<typename m_type>
inline m_type CdtSwingBezier<m_type>::BinomialCoeff(int n, int k)
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
inline m_type CdtSwingBezier<m_type>::BernsteinPoly(int n, int k, m_type x)
{
    // Bernstein basis polynomials
    // b_{k,n}(x) = (n, k)b * x^{k} * (1-x)^{n-k}
    // (,)b is binomial coefficient

    //return BinomialCoeff(n, k) * std::pow(x, k) * std::pow(1 - x, n - k);

    switch (n)
    {
    case 14:
        return m_binomialCoeff14[k] * std::pow(x, k) * std::pow(1 - x, n - k);
    case 15:
    default:
        return m_binomialCoeff15[k] * std::pow(x, k) * std::pow(1 - x, n - k);
    }
}
