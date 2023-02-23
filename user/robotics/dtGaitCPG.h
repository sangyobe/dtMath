
#ifndef DTMATH_DTGAIT_CPG_H_
#define DTMATH_DTGAIT_CPG_H_

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

#define STANCE  0
#define SWING   1
#define TOTAL   2

template <typename m_type = float>
class CdtGaitCPG
{
public:
    uint8_t refState, flwState;         // state 0 is stance and 1 is swing
    m_type refCtrlParam, flwCtrlParam;  // ctrlparam range [0 1]

private:
    m_type m_smpTime_ms;
    m_type m_elapsed_ms;
    m_type m_period_ms[3];      // [0]:stance time, [1]:swing time, [2]:stance + swing time
    m_type m_refPeriod_ms[3];   // [0]:stance time, [1]:swing time, [2]:stance + swing time
    m_type m_flwPeriod_ms[3];   // [0]:stance time, [1]:swing time, [2]:stance + swing time
    m_type m_phaseDelay_pu;     // unit:[percent unit], 0 ~ 1

public:
    CdtGaitCPG();
    int8_t Init(const uint32_t sampleTime_us, const m_type stanceTime_ms, const m_type swingTime_ms);
    int8_t SetSampleTime(const uint32_t sampleTime_us);
    int8_t SetPeriod(const m_type stanceTime_ms, const m_type swingTime_ms);
    void Compute();

private:
    void ReferenceCPG();
    void FollowCPG();
};

#include "dtGaitCPG.tpp"

#endif // DTMATH_DTGAIT_CPG_H_
