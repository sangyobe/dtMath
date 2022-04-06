/*!
\file       dtGaitPlaner.hpp
\brief      Gait Planner using Swing/Stance Bezier trajactory for the legged robot.
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
#include "dtSwingBezier.hpp"
#include "dtStanceBezier.hpp"

#include <cmath>
#include <limits>

#define SWING   0
#define STANCE  1
#define TOTAL   2

template <uint8_t m_num, typename m_type = float>
class CdtGaitPlaner
{
public:
    uint8_t state[m_num];       // state 0 is stance and 1 is swing
    m_type ctrlParam[m_num];    // ctrlparam range [0 1]

private:
    m_type m_smpTime_ms;
    m_type m_newPeriod_ms[3];       // [0]:swing time, [1]:stance time, [2]:stance + swing time

    // Reference CPG variables
    m_type m_refElapseTime_ms;
    uint8_t m_refState;             // state 0 is stance and 1 is swing
    m_type m_refCtrlParam;          // ctrlparam range [0 1]
    m_type m_refPeriod_ms[3];       // [0]:swing time, [1]:stance time, [2]:stance + swing time
    uint8_t m_refUpdatePeriod;      // new period will be updated when the swing period starts.

    // Follow CPG variables
    m_type m_newPhaseDelay_ms[m_num];   // not yet updated phase delay, unit[ms]
    m_type m_phaseDelay_ms[m_num];      // followCPG phase delay, unit[ms]
    m_type m_elapseTime_ms[m_num];      // followCPG elapsed time
    m_type m_period_ms[m_num][3];       // [0]:swing time, [1]:stance time, [2]:stance + swing time
    m_type m_stancePeriod_ms[m_num];    // temporarily save the stance period.
    uint8_t m_updatePeriod[m_num];      // activated after ref CPG period updated, new period will be updated when the swing period starts.
    uint8_t m_updatePhase[m_num];       // followCPG phase update, 0: all updated, 1: not yet updated, 2: modifying stance period for phase delay

    // Trajectory variables
    CdtSwingBezier<m_type>  m_swingBezier[m_num];
    CdtStanceBezier<m_type> m_stanceBezier[m_num];

    CdtVector3<m_type> m_zeroVec3;
    CdtVector3<m_type> *m_pStartPos[m_num];   // address of current foot position
    CdtVector3<m_type> *m_pStartVel[m_num];   // address of current foot velocity

    m_type m_swingHeight[m_num];                // the highest height of the foot during the swing period
    CdtVector3<m_type> m_swingEndPos[m_num];    // swing target position
    CdtVector3<m_type> m_swingEndVel[m_num];    // swing target velocity, typically zero
    CdtVector3<m_type> m_stanceEndPos[m_num];   // stance target position
    CdtVector3<m_type> m_stanceEndVel[m_num];   // stance target velocity, typically zero

    // Gait stop variables
    uint8_t m_stop;                 // m_stop >= m_num, all leg completely stop
    uint8_t m_updateStop[m_num];    // 0 is completely stop, 1: moving to zero pos, 2: it has moved to the zero pos
                                    
    
public:
    CdtGaitPlaner();

    int8_t Init(const uint32_t sampleTime_us, const m_type swingTime_ms, const m_type stanceTime_ms, const m_type *phaseDelay_pu);
    int8_t BindData(const uint8_t idx, CdtVector3<m_type> *pStartPos, CdtVector3<m_type> *pStartVel);
    int8_t SetParam(const m_type stanceTime_ms, const m_type swingTime_ms, const m_type *phaseDelay_pu);
    int8_t SetPeriod(const m_type stanceTime_ms, const m_type swingTime_ms);
    int8_t SetPhase(const uint8_t idx, const m_type phase_pu);

    int8_t SetStep(CdtVector3<m_type> &step, const m_type height = -1);
    int8_t SetStep(CdtVector3<m_type> &step, CdtVector3<m_type> &endVel, const m_type height = -1);
    int8_t SetStep(const uint8_t idx, CdtVector3<m_type> &step, const m_type height = -1);
    int8_t SetStep(const uint8_t idx, CdtVector3<m_type> &step, CdtVector3<m_type> &endVel, const m_type height = -1);
    int8_t SetTarget(const uint8_t idx,
        CdtVector3<m_type> &swingEndPos, CdtVector3<m_type> &swingEndVel,
        CdtVector3<m_type> &stanceEndPos, CdtVector3<m_type> &stanceEndVel,
        const m_type height = -1);

    CdtVector3<m_type>& Pos(uint8_t idx);   // @ local coordinate
    CdtVector3<m_type>& Vel(uint8_t idx);   // @ local coordiante

    void Start();
    void Stop();
    void Compute();

private:
    void ReferenceCPG();
    void FollowCPG(int idx);
};

template<uint8_t m_num, typename m_type>
inline CdtGaitPlaner<m_num, m_type>::CdtGaitPlaner()
{
    m_smpTime_ms = 1;
    m_newPeriod_ms[SWING] = 300;
    m_newPeriod_ms[STANCE] = 1000;
    m_newPeriod_ms[TOTAL] = 1300;

    m_refElapseTime_ms = 0;
    m_refState = STANCE;
    m_refCtrlParam = 2;
    m_refPeriod_ms[SWING] = m_newPeriod_ms[SWING];
    m_refPeriod_ms[STANCE] = m_newPeriod_ms[STANCE];
    m_refPeriod_ms[TOTAL] = m_newPeriod_ms[TOTAL];
    m_refUpdatePeriod = 0;

    for (int i = 0; i < m_num; i++)
    {
        m_newPhaseDelay_ms[i] = 0;
        m_phaseDelay_ms[i] = 0;
        m_elapseTime_ms[i] = 0;
        state[i] = STANCE;
        ctrlParam[i] = 2;
        m_period_ms[i][SWING] = m_newPeriod_ms[SWING];
        m_period_ms[i][STANCE] = m_newPeriod_ms[STANCE];
        m_period_ms[i][TOTAL] = m_newPeriod_ms[TOTAL];
        m_stancePeriod_ms[i] = m_newPeriod_ms[STANCE];
        m_updatePeriod[i] = 0;
        m_updatePhase[i] = 0;

        m_pStartPos[i] = nullptr;
        m_pStartVel[i] = nullptr;
        m_swingHeight[i] = 0.3f; // unit[m]
        m_updateStop[i] = 0;
    }

    m_stop = m_num;
}

template<uint8_t m_num, typename m_type>
inline int8_t CdtGaitPlaner<m_num, m_type>::Init(const uint32_t sampleTime_us, const m_type swingTime_ms, const m_type stanceTime_ms, const m_type *phaseDelay_pu)
{
    if (sampleTime_us <= 0)
    {
        m_smpTime_ms = 1;
        return -1;
    }

    m_smpTime_ms = sampleTime_us / 1000.0f;
    m_newPeriod_ms[SWING] = swingTime_ms;
    m_newPeriod_ms[STANCE] = stanceTime_ms;
    m_newPeriod_ms[TOTAL] = swingTime_ms + stanceTime_ms;

    m_refElapseTime_ms = 0;
    m_refState = STANCE;
    m_refCtrlParam = 2;
    m_refPeriod_ms[SWING] = m_newPeriod_ms[SWING];
    m_refPeriod_ms[STANCE] = m_newPeriod_ms[STANCE];
    m_refPeriod_ms[TOTAL] = m_newPeriod_ms[TOTAL];
    m_refUpdatePeriod = 0;

    for (int i = 0; i < m_num; i++)
    {
        m_newPhaseDelay_ms[i] = m_newPeriod_ms[TOTAL] * phaseDelay_pu[i];
        m_phaseDelay_ms[i] = m_newPhaseDelay_ms[i];
        m_elapseTime_ms[i] = m_newPeriod_ms[TOTAL] - m_phaseDelay_ms[i];
        state[i] = STANCE;
        ctrlParam[i] = 2;
        m_period_ms[i][SWING] = m_newPeriod_ms[SWING];
        m_period_ms[i][STANCE] = m_newPeriod_ms[STANCE];
        m_period_ms[i][TOTAL] = m_newPeriod_ms[TOTAL];
        m_stancePeriod_ms[i] = m_newPeriod_ms[STANCE];
        m_updatePeriod[i] = 0;
        m_updatePhase[i] = 0;

        m_pStartPos[i] = nullptr;
        m_pStartVel[i] = nullptr;
        m_swingHeight[i] = 0.3f; // unit[m]
        m_updateStop[i] = 0;
    }

    m_stop = m_num;

    return 0;
}

template<uint8_t m_num, typename m_type>
inline int8_t CdtGaitPlaner<m_num, m_type>::BindData(const uint8_t idx, CdtVector3<m_type>* pStartPos, CdtVector3<m_type>* pStartVel)
{
    if (pStartPos == nullptr) return -1;
    if (pStartVel == nullptr) return -1;

    m_pStartPos[idx] = pStartPos;
    m_pStartVel[idx] = pStartVel;

    return 0;
}

template<uint8_t m_num, typename m_type>
inline int8_t CdtGaitPlaner<m_num, m_type>::SetParam(const m_type stanceTime_ms, const m_type swingTime_ms, const m_type * phaseDelay_pu)
{
    if (m_refUpdatePeriod) return -1;

    for (uint8_t i = 0; i < m_num; i++)
    {
        if (m_updatePeriod[i] != 0) return -1;
        if (m_updatePhase[i] != 0) return -1;
    }

    // new period
    m_newPeriod_ms[STANCE] = stanceTime_ms;
    m_newPeriod_ms[SWING] = swingTime_ms;
    m_newPeriod_ms[TOTAL] = stanceTime_ms + swingTime_ms;

    // new phase
    for (uint8_t i = 0; i < m_num; i++)
    {
        m_newPhaseDelay_ms[i] = m_newPeriod_ms[TOTAL] * phaseDelay_pu[i];
        m_updatePhase[i] = 1;
    }

    m_refUpdatePeriod = 1;

    return 0;
}

template<uint8_t m_num, typename m_type>
inline int8_t CdtGaitPlaner<m_num, m_type>::SetPeriod(const m_type swingTime_ms, const m_type stanceTime_ms)
{
    if (m_refUpdatePeriod) return -1;

    for (uint8_t i = 0; i < m_num; i++)
    {
        if (m_updatePeriod[i] != 0) return -1;
        if (m_updatePhase[i] != 0) return -1;
    }

    // new period
    m_newPeriod_ms[STANCE] = stanceTime_ms;
    m_newPeriod_ms[SWING] = swingTime_ms;
    m_newPeriod_ms[TOTAL] = stanceTime_ms + swingTime_ms;

    // new phase
    for (uint8_t i = 0; i < m_num; i++)
    {
        m_newPhaseDelay_ms[i] = m_newPeriod_ms[TOTAL] * m_phaseDelay_ms[i];
        m_updatePhase[i] = 1;
    }

    m_refUpdatePeriod = 1;

    return 0;
}

template<uint8_t m_num, typename m_type>
inline int8_t CdtGaitPlaner<m_num, m_type>::SetPhase(const uint8_t idx, const m_type phaseDelay_pu)
{
    if (idx >= m_num) return -1;
    if (m_updatePhase[idx] != 0) return -1;

    m_newPhaseDelay_ms[idx] = m_newPeriod_ms[TOTAL] * phaseDelay_pu;
    m_updatePhase[idx] = 1;

    return 0;
}

template<uint8_t m_num, typename m_type>
inline int8_t CdtGaitPlaner<m_num, m_type>::SetStep(CdtVector3<m_type>& step, const m_type height)
{
    for (uint8_t i = 0; i < m_num; i++)
    {
        m_swingEndPos[i] = step;
        m_swingEndVel[i] = m_zeroVec3;

        m_stanceEndPos[i] = -step;
        m_stanceEndVel[i] = m_zeroVec3;

        if (height >= 0)
            m_swingHeight[i] = height;
    }

    return 0;
}

template<uint8_t m_num, typename m_type>
inline int8_t CdtGaitPlaner<m_num, m_type>::SetStep(CdtVector3<m_type>& step, CdtVector3<m_type>& endVel, const m_type height)
{
    for (uint8_t i = 0; i < m_num; i++)
    {
        m_swingEndPos[i] = step;
        m_swingEndVel[i] = m_zeroVec3;

        m_stanceEndPos[i] = -step;
        m_stanceEndVel[i] = endVel;

        if (height >= 0)
            m_swingHeight[i] = height;
    }

    return 0;
}

template<uint8_t m_num, typename m_type>
inline int8_t CdtGaitPlaner<m_num, m_type>::SetStep(const uint8_t idx, CdtVector3<m_type>& step, const m_type height)
{
    if (idx >= m_num) return -1;

    m_swingEndPos[idx] = step;
    m_swingEndVel[idx] = m_zeroVec3;

    m_stanceEndPos[idx] = -step;
    m_stanceEndVel[idx] = m_zeroVec3;

    if (height >= 0)
        m_swingHeight[idx] = height;

    return 0;
}

template<uint8_t m_num, typename m_type>
inline int8_t CdtGaitPlaner<m_num, m_type>::SetStep(const uint8_t idx, CdtVector3<m_type>& step, CdtVector3<m_type>& endVel, const m_type height)
{
    if (idx >= m_num) return -1;

    m_swingEndPos[idx] = step;
    m_swingEndVel[idx] = m_zeroVec3;

    m_stanceEndPos[idx] = -step;
    m_stanceEndVel[idx] = endVel;

    if (height >= 0)
        m_swingHeight[idx] = height;

    return 0;
}

template<uint8_t m_num, typename m_type>
inline int8_t CdtGaitPlaner<m_num, m_type>::SetTarget(const uint8_t idx,
    CdtVector3<m_type>& swingEndPos, CdtVector3<m_type>& swingEndVel,
    CdtVector3<m_type>& stanceEndPos, CdtVector3<m_type>& stanceEndVel, const m_type height)
{
    if (idx >= m_num) return -1;
    
    m_swingEndPos[idx] = swingEndPos;
    m_swingEndVel[idx] = swingEndVel;

    m_stanceEndPos[idx] = stanceEndPos;
    m_stanceEndVel[idx] = stanceEndVel;

    if (height >= 0)
        m_swingHeight[idx] = height;

    return 0;
}

template<uint8_t m_num, typename m_type>
inline CdtVector3<m_type>& CdtGaitPlaner<m_num, m_type>::Pos(uint8_t idx)
{
    if (idx >= m_num) return m_stanceBezier[idx].pos_x;

    if (state[idx] == SWING) return m_swingBezier[idx].pos_x;
    else return m_stanceBezier[idx].pos_x;
}

template<uint8_t m_num, typename m_type>
inline CdtVector3<m_type>& CdtGaitPlaner<m_num, m_type>::Vel(uint8_t idx)
{
    if (idx >= m_num) return m_stanceBezier[idx].vel_xps;

    if (state[idx] == SWING) return m_swingBezier[idx].vel_xps;
    else return m_stanceBezier[idx].vel_xps;
}

template<uint8_t m_num, typename m_type>
inline void CdtGaitPlaner<m_num, m_type>::Start()
{
    m_stop = 0;
    m_refElapseTime_ms = 0;

    for (uint8_t i = 0; i < m_num; i++)
    {
        m_elapseTime_ms[i] = 0;
        m_updateStop[i] = 0;
    }        
}

template<uint8_t m_num, typename m_type>
inline void CdtGaitPlaner<m_num, m_type>::Stop()
{
    for (uint8_t i=0; i<m_num; i++)
        m_updateStop[i] = 1;
}

template<uint8_t m_num, typename m_type>
inline void CdtGaitPlaner<m_num, m_type>::Compute()
{
    // receive(or update)
    // start pos
    // start vel

    // Generate 
    ReferenceCPG();

    for (uint8_t i = 0; i < m_num; i++)
    {
        FollowCPG(i);
        switch (state[i])
        {
        case SWING:
            m_swingBezier[i].Compute(ctrlParam[i]);
            break;
        case STANCE:
            m_stanceBezier[i].Compute(ctrlParam[i]);
            break;
        }
    }
}

template<uint8_t m_num, typename m_type>
inline void CdtGaitPlaner<m_num, m_type>::ReferenceCPG()
{
    if (m_stop >= m_num) return;

    // Start of swing period
    if (m_refElapseTime_ms == 0 && m_refUpdatePeriod)
    { // update period
        m_refPeriod_ms[SWING] = m_newPeriod_ms[SWING];
        m_refPeriod_ms[STANCE] = m_newPeriod_ms[STANCE];
        m_refPeriod_ms[TOTAL] = m_newPeriod_ms[TOTAL];
        m_refUpdatePeriod = 0;

        // reference period updated and than we will update FollowCPG period
        for (uint8_t i = 0; i< m_num; i++)
            m_updatePeriod[i] = 1;
    }

    // Make Control Parameter
    if (m_refElapseTime_ms < m_refPeriod_ms[SWING])
    {
        m_refState = SWING;
        m_refCtrlParam = m_refElapseTime_ms / m_refPeriod_ms[SWING];
    }
    else
    {
        m_refState = STANCE;
        m_refCtrlParam = (m_refElapseTime_ms - m_refPeriod_ms[SWING]) / m_refPeriod_ms[STANCE];
    }

    m_refElapseTime_ms += m_smpTime_ms;
    if (m_refElapseTime_ms >= m_refPeriod_ms[TOTAL])
        m_refElapseTime_ms = 0;
}

template<uint8_t m_num, typename m_type>
inline void CdtGaitPlaner<m_num, m_type>::FollowCPG(int idx)
{
    if (m_stop >= m_num) return;

    if (m_elapseTime_ms[idx] < m_smpTime_ms)
    {// Start of swing period
        switch (m_updateStop[idx])
        {// update stop
        case 1:// step1. move to zero
            m_swingBezier[idx].SetTarget(*m_pStartPos[idx], m_zeroVec3, *m_pStartVel[idx], m_zeroVec3, m_swingHeight[idx]);
            m_updateStop[idx] = 2;
            break;
        case 2:// step2. stop completely
            m_stanceBezier[idx].SetTarget(*m_pStartPos[idx], *m_pStartPos[idx], m_zeroVec3, m_zeroVec3);
            state[idx] = STANCE;
            ctrlParam[idx] = 0;
            m_updateStop[idx] = 0;
            m_stop++;
            return;
        }

        if (m_updateStop[idx] == 0)
        {// update step
            m_swingBezier[idx].SetTarget(*m_pStartPos[idx], m_swingEndPos[idx], *m_pStartVel[idx], m_swingEndVel[idx], m_swingHeight[idx]);
        }

        if (m_updatePeriod[idx] && m_updateStop[idx] == 0)
        {// update period
            m_period_ms[idx][SWING] = m_newPeriod_ms[SWING];
            m_period_ms[idx][STANCE] = m_newPeriod_ms[STANCE];
            m_period_ms[idx][TOTAL] = m_newPeriod_ms[TOTAL];
            m_stancePeriod_ms[idx] = m_newPeriod_ms[STANCE];
            m_updatePeriod[idx] = 0;

            // update bezier period
            m_swingBezier[idx].SetPeriod(m_period_ms[idx][SWING]);
            m_stanceBezier[idx].SetPeriod(m_period_ms[idx][STANCE]);
        }
    }
    else if (std::abs(m_elapseTime_ms[idx] - m_period_ms[idx][SWING]) < m_smpTime_ms)
    {// Start of stance period
        switch (m_updateStop[idx])
        {// update stop
        case 1:// step1. move to zero
            m_stanceBezier[idx].SetTarget(*m_pStartPos[idx], m_zeroVec3, *m_pStartVel[idx], m_zeroVec3);
            m_updateStop[idx] = 2;
            break;
        case 2:// step2. stop completely
            m_stanceBezier[idx].SetTarget(*m_pStartPos[idx], *m_pStartPos[idx], m_zeroVec3, m_zeroVec3);
            state[idx] = STANCE;
            ctrlParam[idx] = 0;
            m_updateStop[idx] = 0;
            m_stop++;
            return;
        }

        if (m_updateStop[idx] == 0)
        {// update step
            m_stanceBezier[idx].SetTarget(*m_pStartPos[idx], m_stanceEndPos[idx], *m_pStartVel[idx], m_stanceEndVel[idx]);
        }

        if (m_updatePhase[idx] == 1 && m_updateStop[idx] == 0)
        {// update stance period for controlling phase
            m_updatePhase[idx] = 2;
            m_stancePeriod_ms[idx] = m_period_ms[idx][STANCE];
            m_period_ms[idx][STANCE] += m_newPhaseDelay_ms[idx] - m_phaseDelay_ms[idx];

            if (m_period_ms[idx][STANCE] <= std::numeric_limits<m_type>::epsilon())
                m_period_ms[idx][STANCE] = m_newPeriod_ms[TOTAL];
            
            m_phaseDelay_ms[idx] = m_newPhaseDelay_ms[idx];

            // update stance bezier period
            m_stanceBezier[idx].SetPeriod(m_period_ms[idx][STANCE]);
        }
    }

    // Make Control Parameter
    if (m_elapseTime_ms[idx] < m_period_ms[idx][SWING])
    {
        state[idx] = SWING;
        ctrlParam[idx] = m_elapseTime_ms[idx] / m_period_ms[idx][SWING];

        if (m_updatePhase[idx] == 2)
        {// recovery stance period
            m_updatePhase[idx] = 0;
            m_period_ms[idx][STANCE] = m_stancePeriod_ms[idx];

            // update stance bezier period
            m_stanceBezier[idx].SetPeriod(m_period_ms[idx][STANCE]);
        }
    }
    else
    {
        state[idx] = STANCE;
        ctrlParam[idx] = (m_elapseTime_ms[idx] - m_period_ms[idx][SWING]) / m_period_ms[idx][STANCE];
    }

    m_elapseTime_ms[idx] += m_smpTime_ms;
    if (m_elapseTime_ms[idx] >= m_period_ms[idx][TOTAL])
        m_elapseTime_ms[idx] = 0;
}
