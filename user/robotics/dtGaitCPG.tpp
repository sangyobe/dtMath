
template<typename m_type>
inline CdtGaitCPG<m_type>::CdtGaitCPG()
{
    m_smpTime_ms = 1;
    m_elapsed_ms = 0;
    
    m_period_ms[STANCE] = 1000;
    m_period_ms[SWING]  = 300;
    m_period_ms[TOTAL]  = 1300;

    m_refPeriod_ms[STANCE] = m_period_ms[0];
    m_refPeriod_ms[SWING]  = m_period_ms[1];
    m_refPeriod_ms[TOTAL]  = m_period_ms[2];

    m_flwPeriod_ms[STANCE] = m_period_ms[0];
    m_flwPeriod_ms[SWING]  = m_period_ms[1];
    m_flwPeriod_ms[TOTAL]  = m_period_ms[2];
    
}

template<typename m_type>
inline int8_t CdtGaitCPG<m_type>::Init(const uint32_t sampleTime_us, const m_type stanceTime_ms, const m_type swingTime_ms)
{
    if (sampleTime_us <= 0)
    {
        m_smpTime_ms = 1;
        return -1;
    }

    m_smpTime_ms = sampleTime_us / 1000.0f;

    m_elapsed_ms = 0;

    m_period_ms[STANCE] = stanceTime_ms;
    m_period_ms[SWING]  = swingTime_ms;
    m_period_ms[TOTAL]  = stanceTime_ms + swingTime_ms;

    m_refPeriod_ms[STANCE] = m_period_ms[0];
    m_refPeriod_ms[SWING]  = m_period_ms[1];
    m_refPeriod_ms[TOTAL]  = m_period_ms[2];

    m_flwPeriod_ms[STANCE] = m_period_ms[0];
    m_flwPeriod_ms[SWING]  = m_period_ms[1];
    m_flwPeriod_ms[TOTAL]  = m_period_ms[2];

    return 0;
}

template<typename m_type>
inline int8_t CdtGaitCPG<m_type>::SetSampleTime(const uint32_t sampleTime_us)
{
    if (sampleTime_us <= 0)
    {
        m_smpTime_ms = 1;
        return -1;
    }

    m_smpTime_ms = sampleTime_us / 1000.0f;

    return 0;
}

template<typename m_type>
inline int8_t CdtGaitCPG<m_type>::SetPeriod(const m_type stanceTime_ms, const m_type swingTime_ms)
{
    m_period_ms[STANCE] = stanceTime_ms;
    m_period_ms[SWING]  = swingTime_ms;
    m_period_ms[TOTAL]  = stanceTime_ms + swingTime_ms;

    return 0;
}

template<typename m_type>
inline void CdtGaitCPG<m_type>::Compute()
{
    ReferenceCPG();
    FollowCPG();

    m_elapsed_ms += m_smpTime_ms;

    if (m_elapsed_ms > m_refPeriod_ms[TOTAL])
    {
        m_elapsed_ms = 0;
    }
}

template<typename m_type>
inline void CdtGaitCPG<m_type>::ReferenceCPG()
{
    m_type mod = m_elapsed_ms % (m_refPeriod_ms[TOTAL]);

    if (mod > m_refPeriod_ms[STANCE])
    {
        refState = SWING;   // swing
        refCtrlParam = (mod - m_refPeriod_ms[STANCE]) / m_refPeriod_ms[SWING];
    }
    else
    {
        refState = STANCE;  // stance
        refCtrlParam = mod / m_refPeriod_ms[STANCE];

        // update period
        m_refPeriod_ms[STANCE] = m_period_ms[STANCE];  // stance
        m_refPeriod_ms[SWING]  = m_period_ms[SWING];   // swing
        m_refPeriod_ms[TOTAL]  = m_period_ms[TOTAL];   // stance + swing
    }
}

template<typename m_type>
inline void CdtGaitCPG<m_type>::FollowCPG()
{
    // making the dependency of follow cpg on ref cpg
    m_type period_pu, stanceRatio_pu;

    if (refState == SWING)
        period_pu = (refCtrlParam * m_flwPeriod_ms[SWING] + m_flwPeriod_ms[STANCE]) / m_flwPeriod_ms[TOTAL];
    else
        period_pu = (refCtrlParam * m_flwPeriod_ms[STANCE]) / m_flwPeriod_ms[TOTAL];

    period_pu -= m_phaseDelay_pu;
    if (period_pu < 0) period_pu += 1;

    stanceRatio_pu = m_flwPeriod_ms[STANCE] / m_flwPeriod_ms[TOTAL]; // ratio of stance time

    if (period_pu <= stanceRatio_pu)
    {
        flwState = STANCE; // Stance
        flwCtrlParam = period_pu / stanceRatio_pu;

        // update period
        m_flwPeriod_ms[STANCE] = m_period_ms[STANCE];  // stance
        m_flwPeriod_ms[SWING]  = m_period_ms[SWING];   // swing
        m_flwPeriod_ms[TOTAL]  = m_period_ms[TOTAL];   // stance + swing
    }
    else
    {
        flwState = SWING; // Swing
        flwCtrlParam = (period_pu - stanceRatio_pu) / (1 - stanceRatio_pu);
    }
}
