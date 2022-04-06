#pragma once

template<uint16_t m_dimN, uint16_t m_dimM, typename m_type>
inline CdtCARE<m_dimN, m_dimM, m_type>::CdtCARE()
{
    m_mEye.SetIdentity();
    m_maxIteration = 100;
    m_tolerance = std::numeric_limits<m_type>::epsilon();
    //m_tolerance = 1e-9;
}

template<uint16_t m_dimN, uint16_t m_dimM, typename m_type>
inline void CdtCARE<m_dimN, m_dimM, m_type>::SetSolveCriteria(m_type tolerance, uint16_t maxIteration)
{
    m_tolerance = tolerance;
    m_maxIteration = maxIteration;
}

template<uint16_t m_dimN, uint16_t m_dimM, typename m_type>
inline CdtMatrix<m_dimN, m_dimN, m_type> CdtCARE<m_dimN, m_dimM, m_type>::Solve(
    const CdtMatrix<m_dimN, m_dimN, m_type>& A,
    const CdtMatrix<m_dimN, m_dimM, m_type>& B,
    const CdtMatrix<m_dimN, m_dimN, m_type>& Q,
    const CdtMatrix<m_dimM, m_dimM, m_type>& R)
{
    // Is Q Symmetric matrix?
    if (Q != Q.Transpose()) return CdtMatrix<m_dimN, m_dimN, m_type>();

    m_type relative_norm;
    uint16_t iteration = 0;

    //m_mH.SetBlock(0, 0, A);
    //m_mH.SetBlock(0, m_dimN, B * R.LLT().Solve(B.Transpose()));
    //m_mH.SetBlock(m_dimN, 0, Q);
    //m_mH.SetBlock(m_dimN, m_dimN, -A.Transpose());
    //
    //m_mZ = m_mH;

    m_mZ.SetBlock(0, 0, A);
    m_mZ.SetBlock(0, m_dimN, B * R.LLT().Solve(B.Transpose()));
    m_mZ.SetBlock(m_dimN, 0, Q);
    m_mZ.SetBlock(m_dimN, m_dimN, -A.Transpose());

    do {
        m_mZold = m_mZ;
        // R. Byers. Solving the algebraic Riccati equation with the matrix sign
        // function. Linear Algebra Appl., 85:267?279, 1987
        // Added determinant scaling to improve convergence (converges in rough half
        // the iterations with this)
        m_type ck = std::pow(std::abs(m_mZ.Determinant()), m_power);
        m_mZ *= ck;
        m_mZ = m_mZ - (m_mZ - m_mZ.Inv()) * 0.5;
        relative_norm = (m_mZ - m_mZold).GetNorm();
        iteration++;
    } while ((iteration < m_maxIteration) && (relative_norm > m_tolerance));

    m_mZ.GetBlock(0, 0, m_mW11);
    m_mZ.GetBlock(0, m_dimN, m_mW12);
    m_mZ.GetBlock(m_dimN, 0, m_mW21);
    m_mZ.GetBlock(m_dimN, m_dimN, m_mW22);

    m_mLhs.SetBlock(0, 0, m_mW12);
    m_mLhs.SetBlock(m_dimN, 0, m_mW22 + m_mEye);

    m_mRhs.SetBlock(0, 0, m_mW11 + m_mEye);
    m_mRhs.SetBlock(m_dimN, 0, m_mW21);

    return m_mLhs.SVD().Solve(m_mRhs);
}

template<uint16_t m_dimN, uint16_t m_dimM, typename m_type>
inline int8_t CdtCARE<m_dimN, m_dimM, m_type>::Solve(
    const CdtMatrix<m_dimN, m_dimN, m_type>& A,
    const CdtMatrix<m_dimN, m_dimM, m_type>& B,
    const CdtMatrix<m_dimN, m_dimN, m_type>& Q,
    const CdtMatrix<m_dimM, m_dimM, m_type>& R,
    CdtMatrix<m_dimN, m_dimN, m_type>& X)
{
    // Is Q Symmetric matrix?
    if (Q != Q.Transpose()) return -1;

    m_type relative_norm;
    uint16_t iteration = 0;

    //m_mH.SetBlock(0, 0, A);
    //m_mH.SetBlock(0, m_dimN, B * R.LLT().Solve(B.Transpose()));
    //m_mH.SetBlock(m_dimN, 0, Q);
    //m_mH.SetBlock(m_dimN, m_dimN, -A.Transpose());
    //
    //m_mZ = m_mH;

    m_mZ.SetBlock(0, 0, A);
    m_mZ.SetBlock(0, m_dimN, B * R.LLT().Solve(B.Transpose()));
    m_mZ.SetBlock(m_dimN, 0, Q);
    m_mZ.SetBlock(m_dimN, m_dimN, -A.Transpose());

    do {
        m_mZold = m_mZ;
        // R. Byers. Solving the algebraic Riccati equation with the matrix sign
        // function. Linear Algebra Appl., 85:267?279, 1987
        // Added determinant scaling to improve convergence (converges in rough half
        // the iterations with this)
        m_type ck = std::pow(std::abs(m_mZ.Determinant()), m_power);
        m_mZ *= ck;
        m_mZ = m_mZ - (m_mZ - m_mZ.Inv()) * 0.5;
        relative_norm = (m_mZ - m_mZold).GetNorm();
        iteration++;
    } while ((iteration < m_maxIteration) && (relative_norm > m_tolerance));

    m_mZ.GetBlock(0, 0, m_mW11);
    m_mZ.GetBlock(0, m_dimN, m_mW12);
    m_mZ.GetBlock(m_dimN, 0, m_mW21);
    m_mZ.GetBlock(m_dimN, m_dimN, m_mW22);

    m_mLhs.SetBlock(0, 0, m_mW12);
    m_mLhs.SetBlock(m_dimN, 0, m_mW22 + m_mEye);

    m_mRhs.SetBlock(0, 0, m_mW11 + m_mEye);
    m_mRhs.SetBlock(m_dimN, 0, m_mW21);

    return m_mLhs.SVD().Solve(m_mRhs, X);
}
