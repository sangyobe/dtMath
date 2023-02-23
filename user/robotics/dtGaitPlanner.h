/*!
\file       dtGaitPlaner.h
\brief      Gait Planner using Swing/Stance Bezier trajactory for the legged robot.
\author     Dongjin Hyun, mecjin@gmail.com
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2021. 07. 01
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTGAIT_PLANNER_H_
#define DTMATH_DTGAIT_PLANNER_H_

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <dtMath/dtMath.h>

#include "dtSwingBezier.h"
#include "dtStanceBezier.h"

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

#include "dtGaitPlanner.tpp"

#endif // DTMATH_DTGAIT_PLANNER_H_
