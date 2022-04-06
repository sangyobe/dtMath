/*!
\file	    dtQuadProg.hpp
\brief	    Quadratic Programming Solver
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date	    2019. 04. 03
\see        file QuadProg++.h and QuadProg++.cc
\see        https://github.com/liuq/QuadProgpp
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

/*
 File $Id: QuadProg++.hh 232 2007-06-21 12:29:00Z digasper $

 The quadprog_solve() function implements the algorithm of Goldfarb and Idnani
 for the solution of a (convex) Quadratic Programming problem
 by means of an active-set dual method.

The problem is in the form:

min 0.5 * x^t * G * x + g0^t * x
s.t.
    CE^t * x + ce0 = 0
    CI^t * x + ci0 >= 0

 The matrix and vectors dimensions are as follows:
     G: n * n
        g0: n

    CE: n * p
       ce0: p

    CI: n * m
       ci0: m

         x: n

 The function will return the cost of the solution written in the x vector or
 std::numeric_limits::infinity() if the problem is infeasible. In the latter case
 the value of the x vector is not correct.

 References: D. Goldfarb, A. Idnani. A numerically stable dual method for solving
             strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.

 Notes:
  1. pay attention in setting up the vectors ce0 and ci0.
       If the constraints of your problem are specified in the form
       A^T x = b and C^T x >= d, then you should set ce0 = -b and ci0 = -d.
  2. The matrices have column dimension equal to MATRIX_DIM,
     a constant set to 20 in this file (by means of a #define macro).
     If the matrices are bigger than 20 x 20 the limit could be
         increased by means of a -DMATRIX_DIM=n on the compiler command line.
  3. The matrix G is modified within the function since it is used to compute
     the G = L^T L cholesky factorization for further computations inside the function.
     If you need the original matrix G you should make a copy of it and pass the copy
     to the function.

 Author: Luca Di Gaspero
             DIEGM - University of Udine, Italy
                 luca.digaspero@uniud.it
                 http://www.diegm.uniud.it/digaspero/

 The author will be grateful if the researchers using this software will
 acknowledge the contribution of this function in their research papers.

 Copyright (c) 2007-2016 Luca Di Gaspero

 This software may be modified and distributed under the terms
 of the MIT license.  See the LICENSE file for details.
*/

#pragma once

#if defined(_WIN32) || defined(__linux__)
#include <iostream>
#include <cstdlib>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <cmath>
#include <algorithm>
#include <limits>

template <uint16_t m_row, uint16_t m_col, typename m_type> class CdtMatrix;
template <uint16_t m_row, typename m_type> class CdtVector;

//#define DT_QP_DEBUG

#ifdef DT_QP_DEBUG
using namespace std;
#endif

template <int m_dimN, int m_dimM, int m_dimP, typename m_type = float>
class CdtQuadProg
{
private:
    m_type m_inf;
    //CdtMatrix<m_dimN, m_dimN, m_type> m_mG;	// Positive definite symmetric Matrix of the cost function
    CdtVector<m_dimN, m_type> m_vG;             // Vector of the cost function
    CdtMatrix<m_dimN, m_dimN, m_type> m_mL;	    // Cholesky decomposed matrix of the matrix mG
    m_type m_mCi[m_dimN * m_dimM];              // Matrix of the inequality constraints
    m_type m_vCi[m_dimM];                       // Vector of the inequality constraints
    CdtMatrix<m_dimN, m_dimN, m_type> m_mJ0, m_mJ, m_mR;
    CdtVector<m_dimN, m_type> m_vX0;

    CdtVector<m_dimN, m_type> m_vZ, m_vD, m_vNp, m_vOldX;
    CdtVector<m_dimM + m_dimP, m_type> m_vS, m_vR, m_vU, m_vOldU;
    CdtVector<m_dimM + m_dimP, int> m_vA, m_vOldA; // Active Set A
    CdtVector<m_dimM + m_dimP, int> m_vIaIncl;
    CdtVector<m_dimM, int> m_vConstIaIncl;
    CdtVector<m_dimM + m_dimP, bool> m_vIaExcl;

    m_type m_fx, m_fx0;
    m_type m_c1; // c1 is trace of G matrix, c1 * c2 is an estimate for cond(G)
    m_type m_c2; // c2 is trace of J matrix, c1 * c2 is an estimate for cond(G)
    m_type m_normR;

    m_type m_t;			// t is the step lenght, which is the minimum of
    m_type m_t1, m_t2;	// the partial step length t1 and the full step length t2

    int m_iter;

public:
    CdtQuadProg();
    ~CdtQuadProg() {}
    int Iteration() { return m_iter; }
    m_type GetCost() { return m_fx; }

    int8_t SetFixedParam(
        const CdtMatrix<m_dimN, m_dimN, m_type> &mG, const CdtVector<m_dimN, m_type>& vG,
        const CdtMatrix<m_dimN, m_dimM, m_type>& mCi, const CdtVector<m_dimM, m_type>& vCi);

    int8_t Solve(
        const CdtMatrix<m_dimN, m_dimP, m_type>& mCe, const CdtVector<m_dimP, m_type>& vCe,
        CdtVector<m_dimN, m_type>& vX);

    int8_t SetFixedParam(
        const CdtMatrix<m_dimN, m_dimN, m_type> &mG,
        const CdtMatrix<m_dimN, m_dimM, m_type>& mCi, const CdtVector<m_dimM, m_type>& vCi);

    int8_t Solve(
        const CdtVector<m_dimN, m_type>& vG,
        const CdtMatrix<m_dimN, m_dimP, m_type>& mCe, const CdtVector<m_dimP, m_type>& vCe,
        CdtVector<m_dimN, m_type>& vX);

    int8_t SetFixedParam(
        const CdtMatrix<m_dimN, m_dimN, m_type> &mG, const CdtVector<m_dimN, m_type>& vG);

    int8_t Solve(
        const CdtMatrix<m_dimN, m_dimP, m_type>& mCe, const CdtVector<m_dimP, m_type>& vCe,
        const CdtMatrix<m_dimN, m_dimM, m_type>& mCi, const CdtVector<m_dimM, m_type>& vCi,
        CdtVector<m_dimN, m_type>& vX);

    int8_t Solve(
        const CdtMatrix<m_dimN, m_dimN, m_type>& mG, const CdtVector<m_dimN, m_type>& vG,
        const CdtMatrix<m_dimN, m_dimP, m_type>& mCe, const CdtVector<m_dimP, m_type>& vCe,
        const CdtMatrix<m_dimN, m_dimM, m_type>& mCi, const CdtVector<m_dimM, m_type>& vCi,
        CdtVector<m_dimN, m_type>& vX);

private:
    //m_type Distance(m_type a, m_type b); //the euclidean distance between two numbers

    //m_type DotProduct(const CdtVector<m_dimN, m_type>& x, const CdtVector<m_dimN, m_type>& y);
    int8_t CholeskyDecomposition(CdtMatrix<m_dimN, m_dimN, m_type>& A);
    void ForwardElimination(const CdtMatrix<m_dimN, m_dimN, m_type>& L, CdtVector<m_dimN, m_type>& y, const CdtVector<m_dimN, m_type>& b);
    void BackwardElimination(const CdtMatrix<m_dimN, m_dimN, m_type>& U, CdtVector<m_dimN, m_type>& x, const CdtVector<m_dimN, m_type>& y);
    void CholeskySolve(const CdtMatrix<m_dimN, m_dimN, m_type>& L, CdtVector<m_dimN, m_type>& x, const CdtVector<m_dimN, m_type>& b);

    void ComputeVecD(CdtVector<m_dimN, m_type>& vD, const CdtMatrix<m_dimN, m_dimN, m_type>& mJ, const CdtVector<m_dimN, m_type>& vNp); // d = J^T * np
    void UpdateVecZ(CdtVector<m_dimN, m_type>& vZ, const CdtMatrix<m_dimN, m_dimN, m_type>& mJ, const CdtVector<m_dimN, m_type>& vD, const int iq); // z = J2 * d2
    void UpdateVecR(const CdtMatrix<m_dimN, m_dimN, m_type>& mR, CdtVector<m_dimP + m_dimM, m_type>& vR, const CdtVector<m_dimN, m_type>& vD, const int iq); // r = R^-1 d

    m_type Distance(m_type a, m_type b);
    bool AddConstraint(CdtMatrix<m_dimN, m_dimN, m_type>& mR, CdtMatrix<m_dimN, m_dimN, m_type>& mJ, CdtVector<m_dimN, m_type>& vD, int &iq, m_type& R_norm);
    int8_t DelConstraint(CdtMatrix<m_dimN, m_dimN, m_type>& mR, CdtMatrix<m_dimN, m_dimN, m_type>& mJ, CdtVector<m_dimP + m_dimM, int>& vA, CdtVector<m_dimP + m_dimM, m_type>& vU, int &iq, const int l);
};

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::CdtQuadProg() :
    m_fx(0), m_c1(0), m_c2(0), m_normR(0),
    m_t(0), m_t1(0), m_t2(0)
{
    if (std::numeric_limits<m_type>::has_infinity)
        m_inf = std::numeric_limits<m_type>::infinity();
    else
        m_inf = (m_type)3.40282e+038;

    for (int i = 0; i < m_dimM; i++)
        m_vConstIaIncl(i) = i;
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline int8_t CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::SetFixedParam(
    const CdtMatrix<m_dimN, m_dimN, m_type>& mG, const CdtVector<m_dimN, m_type>& vG,
    const CdtMatrix<m_dimN, m_dimM, m_type>& mCi, const CdtVector<m_dimM, m_type>& vCi)
{
    register int cnt, i, j;

    m_mL = mG;
    m_vG = vG;
    memcpy(m_mCi, mCi.GetElementsAddr(), sizeof(m_mCi));
    memcpy(m_vCi, vCi.GetElementsAddr(), sizeof(m_vCi));

    /* Preprocessing phase -------------------------------------------------------------------------*/
    /* compute the trace of the original matrix G */
    m_c1 = mG.Trace();

    /* decompose the matrix G in the form LLT */
    if (CholeskyDecomposition(m_mL)) return -1;

    /*compute the inverse of the factorized matrix inv(G), this is the initial value for H(or J) */
    m_c2 = 0;
    for (i = 0; i < m_dimN; i++)
    {
        m_vD(i) = 1.0;

        ForwardElimination(m_mL, m_vZ, m_vD);

        for (cnt = m_dimN >> 2, j = 0; cnt > 0; cnt--, j += 4)
        {
            m_mJ0(i, j) = m_vZ(j);
            m_mJ0(i, j + 1) = m_vZ(j + 1);
            m_mJ0(i, j + 2) = m_vZ(j + 2);
            m_mJ0(i, j + 3) = m_vZ(j + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, j++)
            m_mJ0(i, j) = m_vZ(j);

        m_c2 += m_vZ(i); // compute the trace of the matrix J
        m_vD(i) = 0;
    }
    /* c1 * c2 is an estimate for cond(G) */

    /* Step 0: Find the unconstrained minimum ------------------------------------------------------*/
    /* Find the unconstrained minimizer of the quadratic form 1/2 * xt * G * x + xt * g0
     * this is a feasible point in the dual space
     * G * x = -g0
     * x = - (inv(G) * g0)
     * and compute the current solution value
     * f = 1/2 * xt * G * x + xt * g0
     * f = 1/2 * xt * g0
    */
    CholeskySolve(m_mL, m_vX0, vG);
    m_vX0 *= -1;
    m_fx0 = 0.5f * (m_vX0.dot(vG));

    return 0;
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline int8_t CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::Solve(
    const CdtMatrix<m_dimN, m_dimP, m_type>& mCe, const CdtVector<m_dimP, m_type>& vCe,
    CdtVector<m_dimN, m_type>& vX)
{
    /* m is the number of inequality constraints */
    /* p is the number of equality constraints */

    int cnt, i, j, k, l; // indices
    int ip; // this is the index of the constraint to be added to the active set
    int iq;
    m_iter = 0;
    m_type ss;
    m_type psi; // this value will contain the sum of all infeasibilities
    m_type sum;

    /* initialize the matrix R */
    m_mR.SetZero();
    m_normR = 1.0; /* this variable will hold the norm of the matrix R */

    /*compute the inverse of the factorized matrix inv(G), this is the initial value for H */
    m_mJ = m_mJ0;

    /* Step 0: Find the unconstrained minimum ------------------------------------------------------*/
    /* Find the unconstrained minimizer of the quadratic form 1/2 * xt * G * x + xt * g0
     * this is a feasible point in the dual space
     * G * x = -g0
     * x = - (inv(G) * g0)
     * and compute the current solution value
     * f = 1/2 * xt * G * x + xt * g0
     * f = 1/2 * xt * g0
    */
    vX = m_vX0;
    m_fx = m_fx0;
    
    /* Add equality constraints to the working set vector A */
    iq = 0;
    for (i = 0; i < m_dimP; i++)
    {
        // Get colum vector from CE Matrix
        mCe.GetColVec(i, m_vNp);

        // compute d, d = J^t * np;
        ComputeVecD(m_vD, m_mJ, m_vNp);

        // z = J2 * d2;
        UpdateVecZ(m_vZ, m_mJ, m_vD, iq);

        // r = R^-1 d, R is Upper Triangle Matrix of QR decompostion
        UpdateVecR(m_mR, m_vR, m_vD, iq);

        // compute full step length t2: i.e., the minimum step in primal space s.t. the contraint becomes feasible
        m_t2 = 0;
        if (std::abs(m_vZ.dot(m_vZ)) > std::numeric_limits<m_type>::epsilon()) // i.e. z != 0
            m_t2 = (-(m_vNp.dot(vX)) - vCe(i)) / (m_vZ.dot(m_vNp));

        /* set x = x + t2 * z */
        for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
        {
            vX(k) += m_t2 * m_vZ(k);
            vX(k + 1) += m_t2 * m_vZ(k + 1);
            vX(k + 2) += m_t2 * m_vZ(k + 2);
            vX(k + 3) += m_t2 * m_vZ(k + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
        {
            vX(k) += m_t2 * m_vZ(k);
        }

        /* set u = u+ */
        m_vU(iq) = m_t2;
        for (cnt = iq >> 2, k = 0; cnt > 0; cnt--, k += 4)
        {
            m_vU(k) -= m_t2 * m_vR(k);
            m_vU(k + 1) -= m_t2 * m_vR(k + 1);
            m_vU(k + 2) -= m_t2 * m_vR(k + 2);
            m_vU(k + 3) -= m_t2 * m_vR(k + 3);
        }
        for (cnt = iq % 4; cnt > 0; cnt--, k++)
        {
            m_vU(k) -= m_t2 * m_vR(k);
        }

        /* compute the new solution value */
        m_fx += 0.5f * (m_t2 * m_t2) * (m_vZ.dot(m_vNp));
        m_vA(i) = -i - 1;

        if (!AddConstraint(m_mR, m_mJ, m_vD, iq, m_normR))
        {
            // Equality constraints are linearly dependent
            //throw std::runtime_error("Constraints are linearly dependent");
            //return m_fx;
            return -1;
        }
    }
    /* End: Add equality constraints to the working set A */

    /* set vIaIncl = K \ A */
    m_vIaIncl.SetBlock(0, m_vConstIaIncl);

    /* Step 1: choose a violated constraint --------------------------------------------------------*/
Label_1:
    m_iter++;

    for (i = m_dimP; i < iq; i++)
    {
        ip = m_vA(i);
        m_vIaIncl(ip) = -1; // vector vIaIncl(m_dimP + m_dimM)
    }

    /* compute s[x] = Ci^T * x + ci0 for all elements of K \ A */
    psi = 0;	// this value will contain the sum of all infeasibilities
    ip = 0;		// ip will be the index of the chosen violated constraint
    ss = 0;

    for (i = 0; i < m_dimM; i++)
    {
        m_vIaExcl(i) = true;
        m_vS(i) = 0;

        for (cnt = m_dimN >> 2, j = 0; cnt > 0; cnt--, j += 4)
        {
            m_vS(i) += m_mCi[j * m_dimM + i] * vX(j);
            m_vS(i) += m_mCi[(j + 1) * m_dimM + i] * vX(j + 1);
            m_vS(i) += m_mCi[(j + 2) * m_dimM + i] * vX(j + 2);
            m_vS(i) += m_mCi[(j + 3) * m_dimM + i] * vX(j + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, j++)
        {
            m_vS(i) += m_mCi[j * m_dimM + i] * vX(j);
        }

        m_vS(i) += m_vCi[i];
        psi += std::min(0, sum);
        //psi += (m_vS(i) < 0) ? m_vS(i) : 0;
    }

    /* numerically there are not infeasibilities anymore */
    if (std::abs(psi) <= m_dimM * std::numeric_limits<m_type>::epsilon() * m_c1 * m_c2 * 100.0)
    {
        return 0;
    }

    /* save old values for u, A and x*/
    for (cnt = iq >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        m_vOldU(i) = m_vU(i); m_vOldA(i) = m_vA(i);
        m_vOldU(i + 1) = m_vU(i + 1); m_vOldA(i + 1) = m_vA(i + 1);
        m_vOldU(i + 2) = m_vU(i + 2); m_vOldA(i + 2) = m_vA(i + 2);
        m_vOldU(i + 3) = m_vU(i + 3); m_vOldA(i + 3) = m_vA(i + 3);
    }
    for (cnt = iq % 4; cnt > 0; cnt--, i++)
    {
        m_vOldU(i) = m_vU(i); m_vOldA(i) = m_vA(i);
    }

    m_vOldX = vX;

    /* Step 2: check for feasibility and determine a new S-pair ------------------------------------*/
Label_2:
    for (cnt = m_dimM >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        if (m_vS(i) < ss && m_vIaIncl(i) != -1 && m_vIaExcl(i))
        {
            ss = m_vS(i); ip = i;
        }
        if (m_vS(i + 1) < ss && m_vIaIncl(i + 1) != -1 && m_vIaExcl(i + 1))
        {
            ss = m_vS(i + 1); ip = i + 1;
        }
        if (m_vS(i + 2) < ss && m_vIaIncl(i + 2) != -1 && m_vIaExcl(i + 2))
        {
            ss = m_vS(i + 2); ip = i + 2;
        }
        if (m_vS(i + 3) < ss && m_vIaIncl(i + 3) != -1 && m_vIaExcl(i + 3))
        {
            ss = m_vS(i + 3); ip = i + 3;
        }
    }
    for (cnt = m_dimM % 4; cnt > 0; cnt--, i++)
    {
        if (m_vS(i) < ss && m_vIaIncl(i) != -1 && m_vIaExcl(i))
        {
            ss = m_vS(i); ip = i;
        }
    }

    if (ss >= 0)
    {
        return 0;
    }

    /* set np = n[ip] */
    for (cnt = m_dimN >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        m_vNp(i) = m_mCi[i * m_dimM + ip];
        m_vNp(i + 1) = m_mCi[(i + 1) * m_dimM + ip];
        m_vNp(i + 2) = m_mCi[(i + 2) * m_dimM + ip];
        m_vNp(i + 3) = m_mCi[(i + 3) * m_dimM + ip];
    }
    for (cnt = m_dimN % 4; cnt > 0; cnt--, i++)
    {
        m_vNp(i) = m_mCi[i * m_dimM + ip];
    }

    /* set u = [u 0]^T */
    m_vU(iq) = 0;

    /* add ip to the active set A */
    m_vA(iq) = ip;

    /* Step 2a: determine step direction ------------------------------------------------------------*/
Label_2a:
    /* compute z = H np: the step direction in the primal space (through J, see the paper) */
    ComputeVecD(m_vD, m_mJ, m_vNp);
    UpdateVecZ(m_vZ, m_mJ, m_vD, iq);

    /* compute r = N* np (if q > 0): the negative of the step direction in the dual space */
    UpdateVecR(m_mR, m_vR, m_vD, iq);

    /* Step 2b: compute step length ----------------------------------------------------------------*/
    l = 0;

    /* (i) Compute t1: partial step length (maximum step in dual space without violating dual feasibility) */
    m_t1 = m_inf; // +inf
    // find the index l s.t. it reaches the minimum of u+[x] / r
    for (k = m_dimP; k < iq; k++)
    {
        if (m_vR(k) > 0)
        {
            if (m_vU(k) / m_vR(k) < m_t1)
            {
                m_t1 = m_vU(k) / m_vR(k);
                l = m_vA(k);
            }
        }
    }

    /* (ii) Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible) */
    if (std::abs(m_vZ.dot(m_vZ)) > std::numeric_limits<m_type>::epsilon()) // i.e. z != 0
    {
        m_t2 = -m_vS(ip) / (m_vZ.dot(m_vNp));
        if (m_t2 < 0) // patch suggested by Takano Akio for handling numerical inconsistencies
            m_t2 = m_inf;
    }
    else
        m_t2 = m_inf; // +inf

    /* (iii) Step length, t: The step is chosen as the minimum of t1 and t2 */
    m_t = std::min(m_t1, m_t2);
    //m_t = (m_t1 < m_t2) ? m_t1 : m_t2;

    /* Step 2c: Determine new S-pair and take step: ------------------------------------------------*/
    /* case (i): No step in primal or dual space */
    if (m_t >= m_inf)
    {
        /* QPP are infeasible */
        // FIXME: unbounded to raise
        m_fx = m_inf;
        return -1;
    }

    /* case (ii): Step in dual space */
    if (m_t2 >= m_inf)
    {
        /* set u = u +  t * [-r 1] and drop constraint l from the active set A */
        for (cnt = iq >> 2, k = 0; cnt > 0; cnt--, k += 4)
        {
            m_vU(k) -= m_t * m_vR(k);
            m_vU(k + 1) -= m_t * m_vR(k + 1);
            m_vU(k + 2) -= m_t * m_vR(k + 2);
            m_vU(k + 3) -= m_t * m_vR(k + 3);
        }
        for (cnt = iq % 4; cnt > 0; cnt--, k++)
        {
            m_vU(k) -= m_t * m_vR(k);
        }

        m_vU(iq) += m_t;
        m_vIaIncl(l) = l;
        DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, l);

        goto Label_2a;
    }

    /* case (iii): Step in primal and dual space */
    /* set x = x + t * z */
    for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        vX(k) += m_t * m_vZ(k);
        vX(k + 1) += m_t * m_vZ(k + 1);
        vX(k + 2) += m_t * m_vZ(k + 2);
        vX(k + 3) += m_t * m_vZ(k + 3);
    }
    for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
    {
        vX(k) += m_t * m_vZ(k);
    }

    // update the solution value F(x)
    m_fx += m_t * (m_vZ.dot(m_vNp)) * (0.5f * m_t + m_vU(iq));

    // u = u + t * [-r 1]
    for (cnt = iq >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        m_vU(k) -= m_t * m_vR(k);
        m_vU(k + 1) -= m_t * m_vR(k + 1);
        m_vU(k + 2) -= m_t * m_vR(k + 2);
        m_vU(k + 3) -= m_t * m_vR(k + 3);
    }
    for (cnt = iq % 4; cnt > 0; cnt--, k++)
    {
        m_vU(k) -= m_t * m_vR(k);
    }

    m_vU(iq) += m_t;

    /* If t == t2(full step), u = (u^+), add constraint update H and N* and go to Step1 */
    if (std::abs(m_t - m_t2) < std::numeric_limits<m_type>::epsilon())
    {
        // full step has taken
        // add constraint ip to the active set
        if (!AddConstraint(m_mR, m_mJ, m_vD, iq, m_normR))
        {
            m_vIaExcl(ip) = false;
            DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, ip);

            m_vIaIncl.SetBlock(0, m_vConstIaIncl);

            for (i = m_dimP; i < iq; i++)
            {
                m_vA(i) = m_vOldA(i);
                m_vU(i) = m_vOldU(i);
                m_vIaIncl(m_vA(i)) = -1;
            }

            vX = m_vOldX;

            goto Label_2; /* go to step 2 */
        }
        else
            m_vIaIncl(ip) = -1;

        goto Label_1;
    }

    /* a patial step has taken */
    /* If t == t1(partial step), drop constraint l, update H and N*, and go to Step_2a */

    // drop constraint l
    m_vIaIncl(l) = l;
    DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, l);

    // update s[ip] = CI * x + ci0 */
    m_vS(ip) = 0;
    for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        m_vS(ip) += m_mCi[k * m_dimM + ip] * vX(k);
        m_vS(ip) += m_mCi[(k + 1) * m_dimM + ip] * vX(k + 1);
        m_vS(ip) += m_mCi[(k + 2) * m_dimM + ip] * vX(k + 2);
        m_vS(ip) += m_mCi[(k + 3) * m_dimM + ip] * vX(k + 3);
    }
    for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
        m_vS(ip) += m_mCi[k * m_dimM + ip] * vX(k);

    m_vS(ip) += m_vCi[ip];

    goto Label_2a;
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline int8_t CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::SetFixedParam(
    const CdtMatrix<m_dimN, m_dimN, m_type> &mG,
    const CdtMatrix<m_dimN, m_dimM, m_type>& mCi, const CdtVector<m_dimM, m_type>& vCi)
{
    register int cnt, i, j;

    m_mL = mG;
    memcpy(m_mCi, mCi.GetElementsAddr(), sizeof(m_mCi));
    memcpy(m_vCi, vCi.GetElementsAddr(), sizeof(m_vCi));

    /* Preprocessing phase -------------------------------------------------------------------------*/
    /* compute the trace of the original matrix G */
    m_c1 = mG.Trace();

    /* decompose the matrix G in the form LLT */
    if (CholeskyDecomposition(m_mL)) return -1;

    /*compute the inverse of the factorized matrix U, this is the initial value for H(or J) */
    for (i = 0; i < m_dimN; i++)
    {
        m_vD(i) = 1.0;

        ForwardElimination(m_mL, m_vZ, m_vD);

        for (cnt = m_dimN >> 2, j = 0; cnt > 0; cnt--, j += 4)
        {
            m_mJ0(i, j) = m_vZ(j);
            m_mJ0(i, j + 1) = m_vZ(j + 1);
            m_mJ0(i, j + 2) = m_vZ(j + 2);
            m_mJ0(i, j + 3) = m_vZ(j + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, j++)
            m_mJ0(i, j) = m_vZ(j);

        m_c2 += m_vZ(i); // compute the trace of the matrix J
        m_vD(i) = 0;
    }
    /* c1 * c2 is an estimate for cond(G) */

    return 0;
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline int8_t CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::Solve(
    const CdtVector<m_dimN, m_type>& vG,
    const CdtMatrix<m_dimN, m_dimP, m_type>& mCe, const CdtVector<m_dimP, m_type>& vCe,
    CdtVector<m_dimN, m_type>& vX)
{
    /* m is the number of inequality constraints */
    /* p is the number of equality constraints */

    register int cnt, i, j, k, l; // indices
    int ip; // this is the index of the constraint to be added to the active set
    int iq;
    m_iter = 0;
    m_type ss;
    m_type psi; // this value will contain the sum of all infeasibilities

    /* initialize the matrix R */
    m_mR.SetZero();
    m_normR = 1.0; /* this variable will hold the norm of the matrix R */

    /*compute the inverse of the factorized matrix inv(G), this is the initial value for H */
    m_mJ = m_mJ0;

    /* Step 0: Find the unconstrained minimum ------------------------------------------------------*/
    /* Find the unconstrained minimizer of the quadratic form 1/2 * xt * G * x + xt * g0
     * this is a feasible point in the dual space
     * G * x = -g0
     * x = - (inv(G) * g0)
     * and compute the current solution value
     * f = 1/2 * xt * G * x + xt * g0
     * f = 1/2 * xt * g0
    */
    CholeskySolve(m_mL, vX, vG);
    vX *= -1;
    m_fx = 0.5f * (vX.dot(vG));

    /* Add equality constraints to the working set vector A */
    iq = 0;
    for (i = 0; i < m_dimP; i++)
    {
        // Get colum vector from CE Matrix
        mCe.GetColVec(i, m_vNp);

        // compute d, d = J^t * np;
        ComputeVecD(m_vD, m_mJ, m_vNp);

        // z = J2 * d2;
        UpdateVecZ(m_vZ, m_mJ, m_vD, iq);

        // r = R^-1 d, R is Upper Triangle Matrix of QR decompostion
        UpdateVecR(m_mR, m_vR, m_vD, iq);

        // compute full step length t2: i.e., the minimum step in primal space s.t. the contraint becomes feasible
        m_t2 = 0;
        if (std::abs(m_vZ.dot(m_vZ)) > std::numeric_limits<m_type>::epsilon()) // i.e. z != 0
            m_t2 = (-(m_vNp.dot(vX)) - vCe(i)) / (m_vZ.dot(m_vNp));

        /* set x = x + t2 * z */
        for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
        {
            vX(k) += m_t2 * m_vZ(k);
            vX(k + 1) += m_t2 * m_vZ(k + 1);
            vX(k + 2) += m_t2 * m_vZ(k + 2);
            vX(k + 3) += m_t2 * m_vZ(k + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
        {
            vX(k) += m_t2 * m_vZ(k);
        }

        /* set u = u+ */
        m_vU(iq) = m_t2;
        for (k = 0; k < iq; k++)
            m_vU(k) -= m_t2 * m_vR(k);

        /* compute the new solution value */
        m_fx += 0.5f * (m_t2 * m_t2) * (m_vZ.dot(m_vNp));
        m_vA(i) = -i - 1;

        if (!AddConstraint(m_mR, m_mJ, m_vD, iq, m_normR))
        {
            // Equality constraints are linearly dependent
            //throw std::runtime_error("Constraints are linearly dependent");
            //return m_fx;
            if (cost != nullptr) *cost = m_fx;
            return -1;
        }
    }
    /* End: Add equality constraints to the working set A */

    /* set vIaIncl = K \ A */
    m_vIaIncl.SetBlock(0, m_vConstIaIncl);

    /* Step 1: choose a violated constraint --------------------------------------------------------*/
Label_1:
    m_iter++;

    for (i = m_dimP; i < iq; i++)
    {
        ip = m_vA(i);
        m_vIaIncl(ip) = -1; // vector vIaIncl(m_dimP + m_dimM)
    }

    /* compute s[x] = Ci^T * x + ci0 for all elements of K \ A */
    psi = 0;	// this value will contain the sum of all infeasibilities
    ip = 0;		// ip will be the index of the chosen violated constraint
    ss = 0;

    for (i = 0; i < m_dimM; i++)
    {
        m_vIaExcl(i) = true;
        m_vS(i) = 0;

        for (cnt = m_dimN >> 2, j = 0; cnt > 0; cnt--, j += 4)
        {
            m_vS(i) += m_mCi[j * m_dimM + i] * vX(j);
            m_vS(i) += m_mCi[(j + 1) * m_dimM + i] * vX(j + 1);
            m_vS(i) += m_mCi[(j + 2) * m_dimM + i] * vX(j + 2);
            m_vS(i) += m_mCi[(j + 3) * m_dimM + i] * vX(j + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, j++)
        {
            m_vS(i) += m_mCi[j * m_dimM + i] * vX(j);
        }

        m_vS(i) += m_vCi[i];
        psi += std::min(0, sum);
        //psi += (m_vS(i) < 0) ? m_vS(i) : 0;
    }

    /* numerically there are not infeasibilities anymore */
    if (std::abs(psi) <= m_dimM * std::numeric_limits<m_type>::epsilon() * m_c1 * m_c2 * 100.0)
    {
        return 0;
    }

    /* save old values for u, A and x*/
    for (cnt = iq >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        m_vOldU(i) = m_vU(i); m_vOldA(i) = m_vA(i);
        m_vOldU(i + 1) = m_vU(i + 1); m_vOldA(i + 1) = m_vA(i + 1);
        m_vOldU(i + 2) = m_vU(i + 2); m_vOldA(i + 2) = m_vA(i + 2);
        m_vOldU(i + 3) = m_vU(i + 3); m_vOldA(i + 3) = m_vA(i + 3);
    }
    for (cnt = iq % 4; cnt > 0; cnt--, i++)
    {
        m_vOldU(i) = m_vU(i); m_vOldA(i) = m_vA(i);
    }

    m_vOldX = vX;

    /* Step 2: check for feasibility and determine a new S-pair ------------------------------------*/
Label_2:
    for (cnt = m_dimM >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        if (m_vS(i) < ss && m_vIaIncl(i) != -1 && m_vIaExcl(i))
        {
            ss = m_vS(i); ip = i;
        }
        if (m_vS(i + 1) < ss && m_vIaIncl(i + 1) != -1 && m_vIaExcl(i + 1))
        {
            ss = m_vS(i + 1); ip = i + 1;
        }
        if (m_vS(i + 2) < ss && m_vIaIncl(i + 2) != -1 && m_vIaExcl(i + 2))
        {
            ss = m_vS(i + 2); ip = i + 2;
        }
        if (m_vS(i + 3) < ss && m_vIaIncl(i + 3) != -1 && m_vIaExcl(i + 3))
        {
            ss = m_vS(i + 3); ip = i + 3;
        }
    }
    for (cnt = m_dimM % 4; cnt > 0; cnt--, i++)
    {
        if (m_vS(i) < ss && m_vIaIncl(i) != -1 && m_vIaExcl(i))
        {
            ss = m_vS(i); ip = i;
        }
    }

    if (ss >= 0)
    {
        if (cost != nullptr) *cost = m_fx;
        return 0;
    }

    /* set np = n[ip] */
    for (cnt = m_dimN >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        m_vNp(i) = m_mCi[i * m_dimM + ip];
        m_vNp(i + 1) = m_mCi[(i + 1) * m_dimM + ip];
        m_vNp(i + 2) = m_mCi[(i + 2) * m_dimM + ip];
        m_vNp(i + 3) = m_mCi[(i + 3) * m_dimM + ip];
    }
    for (cnt = m_dimN % 4; cnt > 0; cnt--, i++)
    {
        m_vNp(i) = m_mCi[i * m_dimM + ip];
    }

    /* set u = [u 0]^T */
    m_vU(iq) = 0;

    /* add ip to the active set A */
    m_vA(iq) = ip;

    /* Step 2a: determine step direction ------------------------------------------------------------*/
Label_2a:
    /* compute z = H np: the step direction in the primal space (through J, see the paper) */
    ComputeVecD(m_vD, m_mJ, m_vNp);
    UpdateVecZ(m_vZ, m_mJ, m_vD, iq);

    /* compute r = N* np (if q > 0): the negative of the step direction in the dual space */
    UpdateVecR(m_mR, m_vR, m_vD, iq);

    /* Step 2b: compute step length ----------------------------------------------------------------*/
    l = 0;

    /* (i) Compute t1: partial step length (maximum step in dual space without violating dual feasibility) */
    m_t1 = m_inf; // +inf
    // find the index l s.t. it reaches the minimum of u+[x] / r
    for (k = m_dimP; k < iq; k++)
    {
        if (m_vR(k) > 0)
        {
            if (m_vU(k) / m_vR(k) < m_t1)
            {
                m_t1 = m_vU(k) / m_vR(k);
                l = m_vA(k);
            }
        }
    }

    /* (ii) Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible) */
    if (std::abs(m_vZ.dot(m_vZ)) > std::numeric_limits<m_type>::epsilon()) // i.e. z != 0
    {
        m_t2 = -m_vS(ip) / (m_vZ.dot(m_vNp));
        if (m_t2 < 0) // patch suggested by Takano Akio for handling numerical inconsistencies
            m_t2 = m_inf;
    }
    else
        m_t2 = m_inf; // +inf

    /* (iii) Step length, t: The step is chosen as the minimum of t1 and t2 */
    m_t = std::min(m_t1, m_t2);
    //m_t = (m_t1 < m_t2) ? m_t1 : m_t2;

    /* Step 2c: Determine new S-pair and take step: ------------------------------------------------*/
    /* case (i): No step in primal or dual space */
    if (m_t >= m_inf)
    {
        /* QPP are infeasible */
        // FIXME: unbounded to raise
        m_fx = m_inf;
        return -1;
    }

    /* case (ii): Step in dual space */
    if (m_t2 >= m_inf)
    {
        /* set u = u +  t * [-r 1] and drop constraint l from the active set A */
        for (cnt = iq >> 2, k = 0; cnt > 0; cnt--, k += 4)
        {
            m_vU(k) -= m_t * m_vR(k);
            m_vU(k + 1) -= m_t * m_vR(k + 1);
            m_vU(k + 2) -= m_t * m_vR(k + 2);
            m_vU(k + 3) -= m_t * m_vR(k + 3);
        }
        for (cnt = iq % 4; cnt > 0; cnt--, k++)
        {
            m_vU(k) -= m_t * m_vR(k);
        }

        m_vU(iq) += m_t;
        m_vIaIncl(l) = l;
        DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, l);

        goto Label_2a;
    }

    /* case (iii): Step in primal and dual space */
    /* set x = x + t * z */
    for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        vX(k) += m_t * m_vZ(k);
        vX(k + 1) += m_t * m_vZ(k + 1);
        vX(k + 2) += m_t * m_vZ(k + 2);
        vX(k + 3) += m_t * m_vZ(k + 3);
    }
    for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
    {
        vX(k) += m_t * m_vZ(k);
    }

    // update the solution value F(x)
    m_fx += m_t * (m_vZ.dot(m_vNp)) * (0.5f * m_t + m_vU(iq));

    // u = u + t * [-r 1]
    for (cnt = iq >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        m_vU(k) -= m_t * m_vR(k);
        m_vU(k + 1) -= m_t * m_vR(k + 1);
        m_vU(k + 2) -= m_t * m_vR(k + 2);
        m_vU(k + 3) -= m_t * m_vR(k + 3);
    }
    for (cnt = iq % 4; cnt > 0; cnt--, k++)
    {
        m_vU(k) -= m_t * m_vR(k);
    }

    m_vU(iq) += m_t;

    /* If t == t2(full step), u = (u^+), add constraint update H and N* and go to Step1 */
    if (std::abs(m_t - m_t2) < std::numeric_limits<m_type>::epsilon())
    {
        // full step has taken
        // add constraint ip to the active set
        if (!AddConstraint(m_mR, m_mJ, m_vD, iq, m_normR))
        {
            m_vIaExcl(ip) = false;
            DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, ip);

            m_vIaIncl.SetBlock(0, m_vConstIaIncl);

            for (i = m_dimP; i < iq; i++)
            {
                m_vA(i) = m_vOldA(i);
                m_vU(i) = m_vOldU(i);
                m_vIaIncl(m_vA(i)) = -1;
            }

            vX = m_vOldX;

            goto Label_2; /* go to step 2 */
        }
        else
            m_vIaIncl(ip) = -1;

        goto Label_1;
    }

    /* a patial step has taken */
    /* If t == t1(partial step), drop constraint l, update H and N*, and go to Step_2a */

    // drop constraint l
    m_vIaIncl(l) = l;
    DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, l);

    // update s[ip] = CI * x + ci0 */
    m_vS(ip) = 0;
    for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        m_vS(ip) += m_mCi[k * m_dimM + ip] * vX(k);
        m_vS(ip) += m_mCi[(k + 1) * m_dimM + ip] * vX(k + 1);
        m_vS(ip) += m_mCi[(k + 2) * m_dimM + ip] * vX(k + 2);
        m_vS(ip) += m_mCi[(k + 3) * m_dimM + ip] * vX(k + 3);
    }
    for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
        m_vS(ip) += m_mCi[k * m_dimM + ip] * vX(k);

    m_vS(ip) += m_vCi[ip];

    goto Label_2a;
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline int8_t CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::SetFixedParam(
    const CdtMatrix<m_dimN, m_dimN, m_type>& mG, const CdtVector<m_dimN, m_type>& vG)
{
    register int cnt, i, j;

    m_mL = mG;
    m_vG = vG;

    /* Preprocessing phase -------------------------------------------------------------------------*/
    /* compute the trace of the original matrix G */
    m_c1 = mG.Trace();

    /* decompose the matrix G in the form LLT */
    if (CholeskyDecomposition(m_mL)) return -1;

    /*compute the inverse of the factorized matrix inv(G), this is the initial value for H(or J) */
    m_c2 = 0;
    for (i = 0; i < m_dimN; i++)
    {
        m_vD(i) = 1.0;

        ForwardElimination(m_mL, m_vZ, m_vD);

        for (cnt = m_dimN >> 2, j = 0; cnt > 0; cnt--, j += 4)
        {
            m_mJ0(i, j) = m_vZ(j);
            m_mJ0(i, j + 1) = m_vZ(j + 1);
            m_mJ0(i, j + 2) = m_vZ(j + 2);
            m_mJ0(i, j + 3) = m_vZ(j + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, j++)
            m_mJ0(i, j) = m_vZ(j);

        m_c2 += m_vZ(i); // compute the trace of the matrix J
        m_vD(i) = 0;
    }
    /* c1 * c2 is an estimate for cond(G) */

    /* Step 0: Find the unconstrained minimum ------------------------------------------------------*/
    /* Find the unconstrained minimizer of the quadratic form 1/2 * xt * G * x + xt * g0
     * this is a feasible point in the dual space
     * G * x = -g0
     * x = - (inv(G) * g0)
     * and compute the current solution value
     * f = 1/2 * xt * G * x + xt * g0
     * f = 1/2 * xt * g0
    */
    CholeskySolve(m_mL, m_vX0, vG);
    m_vX0 *= -1;
    m_fx0 = 0.5f * (m_vX0.dot(vG));

    return 0;
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline int8_t CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::Solve(
    const CdtMatrix<m_dimN, m_dimP, m_type>& mCe, const CdtVector<m_dimP, m_type>& vCe,
    const CdtMatrix<m_dimN, m_dimM, m_type>& mCi, const CdtVector<m_dimM, m_type>& vCi,
    CdtVector<m_dimN, m_type>& vX)
{
    /* m is the number of inequality constraints */
    /* p is the number of equality constraints */

    int cnt, i, j, k, l; // indices
    int ip; // this is the index of the constraint to be added to the active set
    int iq;
    m_iter = 0;
    m_type ss;
    m_type psi; // this value will contain the sum of all infeasibilities

    m_mJ = m_mJ0;
    m_mR.SetZero();
    m_normR = 1.0; /* this variable will hold the norm of the matrix R */

    /* Step 0: Find the unconstrained minimum ------------------------------------------------------*/
    /* Find the unconstrained minimizer of the quadratic form 1/2 * xt * G * x + xt * g0
     * this is a feasible point in the dual space
     * G * x = -g0
     * x = - (inv(G) * g0)
     * and compute the current solution value
     * f = 1/2 * xt * G * x + xt * g0
     * f = 1/2 * xt * g0
    */
    vX = m_vX0;
    m_fx = m_fx0;

    /* Add equality constraints to the working set vector A */
    iq = 0;
    for (i = 0; i < m_dimP; i++)
    {
        // Get colum vector from CE Matrix
        mCe.GetColVec(i, m_vNp);

        // compute d, d = J^t * np;
        ComputeVecD(m_vD, m_mJ, m_vNp);

        // z = J2 * d2;
        UpdateVecZ(m_vZ, m_mJ, m_vD, iq);

        // r = R^-1 d, R is Upper Triangle Matrix of QR decompostion
        UpdateVecR(m_mR, m_vR, m_vD, iq);

        // compute full step length t2: i.e., the minimum step in primal space s.t. the contraint becomes feasible
        m_t2 = 0;
        if (std::abs(m_vZ.dot(m_vZ)) > std::numeric_limits<m_type>::epsilon()) // i.e. z != 0
            m_t2 = (-(m_vNp.dot(vX)) - vCe(i)) / (m_vZ.dot(m_vNp));

        /* set x = x + t2 * z */
        for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
        {
            vX(k) += m_t2 * m_vZ(k);
            vX(k + 1) += m_t2 * m_vZ(k + 1);
            vX(k + 2) += m_t2 * m_vZ(k + 2);
            vX(k + 3) += m_t2 * m_vZ(k + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
        {
            vX(k) += m_t2 * m_vZ(k);
        }

        /* set u = u+ */
        m_vU(iq) = m_t2;
        for (k = 0; k < iq; k++)
            m_vU(k) -= m_t2 * m_vR(k);

        /* compute the new solution value */
        m_fx += 0.5f * (m_t2 * m_t2) * (m_vZ.dot(m_vNp));
        m_vA(i) = -i - 1;

        if (!AddConstraint(m_mR, m_mJ, m_vD, iq, m_normR))
        {
            // Equality constraints are linearly dependent
            //throw std::runtime_error("Constraints are linearly dependent");
            //return m_fx;
            return -1;
        }
    }
    /* End: Add equality constraints to the working set A */

    /* set vIaIncl = K \ A */
    m_vIaIncl.SetBlock(0, m_vConstIaIncl);

    /* Step 1: choose a violated constraint --------------------------------------------------------*/
Label_1:
    m_iter++;

    for (i = m_dimP; i < iq; i++)
    {
        ip = m_vA(i);
        m_vIaIncl(ip) = -1; // vector vIaIncl(m_dimP + m_dimM)
    }

    /* compute s[x] = Ci^T * x + ci0 for all elements of K \ A */
    psi = 0;	// this value will contain the sum of all infeasibilities
    ip = 0;		// ip will be the index of the chosen violated constraint
    ss = 0;

    for (i = 0; i < m_dimM; i++)
    {
        m_vIaExcl(i) = true;
        m_vS(i) = 0;

        for (cnt = m_dimN >> 2, j = 0; cnt > 0; cnt--, j += 4)
        {
            m_vS(i) += mCi(j, i) * vX(j);
            m_vS(i) += mCi(j + 1, i) * vX(j + 1);
            m_vS(i) += mCi(j + 2, i) * vX(j + 2);
            m_vS(i) += mCi(j + 3, i) * vX(j + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, j++)
        {
            m_vS(i) += mCi(j, i) * vX(j);
        }

        m_vS(i) += vCi(i);
        psi += std::min(0, m_vS(i));
        //psi += (m_vS(i) < 0) ? m_vS(i) : 0;
    }

    /* numerically there are not infeasibilities anymore */
    if (std::abs(psi) <= m_dimM * std::numeric_limits<m_type>::epsilon() * m_c1 * m_c2 * 100.0)
    {
        return 0;
    }

    /* save old values for u, A and x*/
    for (cnt = iq >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        m_vOldU(i) = m_vU(i); m_vOldA(i) = m_vA(i);
        m_vOldU(i + 1) = m_vU(i + 1); m_vOldA(i + 1) = m_vA(i + 1);
        m_vOldU(i + 2) = m_vU(i + 2); m_vOldA(i + 2) = m_vA(i + 2);
        m_vOldU(i + 3) = m_vU(i + 3); m_vOldA(i + 3) = m_vA(i + 3);
    }
    for (cnt = iq % 4; cnt > 0; cnt--, i++)
    {
        m_vOldU(i) = m_vU(i); m_vOldA(i) = m_vA(i);
    }

    m_vOldX = vX;

    /* Step 2: check for feasibility and determine a new S-pair ------------------------------------*/
Label_2:
    for (cnt = m_dimM >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        if (m_vS(i) < ss && m_vIaIncl(i) != -1 && m_vIaExcl(i))
        {
            ss = m_vS(i); ip = i;
        }
        if (m_vS(i + 1) < ss && m_vIaIncl(i + 1) != -1 && m_vIaExcl(i + 1))
        {
            ss = m_vS(i + 1); ip = i + 1;
        }
        if (m_vS(i + 2) < ss && m_vIaIncl(i + 2) != -1 && m_vIaExcl(i + 2))
        {
            ss = m_vS(i + 2); ip = i + 2;
        }
        if (m_vS(i + 3) < ss && m_vIaIncl(i + 3) != -1 && m_vIaExcl(i + 3))
        {
            ss = m_vS(i + 3); ip = i + 3;
        }
    }
    for (cnt = m_dimM % 4; cnt > 0; cnt--, i++)
    {
        if (m_vS(i) < ss && m_vIaIncl(i) != -1 && m_vIaExcl(i))
        {
            ss = m_vS(i); ip = i;
        }
    }

    if (ss >= 0)
    {
        return 0;
    }

    /* set np = n[ip] */
    mCi.GetColVec(ip, m_vNp);

    /* set u = [u 0]^T */
    m_vU(iq) = 0;

    /* add ip to the active set A */
    m_vA(iq) = ip;

    /* Step 2a: determine step direction ------------------------------------------------------------*/
Label_2a:
    /* compute z = H np: the step direction in the primal space (through J, see the paper) */
    ComputeVecD(m_vD, m_mJ, m_vNp);
    UpdateVecZ(m_vZ, m_mJ, m_vD, iq);

    /* compute r = N* np (if q > 0): the negative of the step direction in the dual space */
    UpdateVecR(m_mR, m_vR, m_vD, iq);

    /* Step 2b: compute step length ----------------------------------------------------------------*/
    l = 0;

    /* (i) Compute t1: partial step length (maximum step in dual space without violating dual feasibility) */
    m_t1 = m_inf; // +inf
    // find the index l s.t. it reaches the minimum of u+[x] / r
    for (k = m_dimP; k < iq; k++)
    {
        if (m_vR(k) > 0)
        {
            if (m_vU(k) / m_vR(k) < m_t1)
            {
                m_t1 = m_vU(k) / m_vR(k);
                l = m_vA(k);
            }
        }
    }

    /* (ii) Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible) */
    if (std::abs(m_vZ.dot(m_vZ)) > std::numeric_limits<m_type>::epsilon()) // i.e. z != 0
    {
        m_t2 = -m_vS(ip) / (m_vZ.dot(m_vNp));
        if (m_t2 < 0) // patch suggested by Takano Akio for handling numerical inconsistencies
            m_t2 = m_inf;
    }
    else
        m_t2 = m_inf; // +inf

    /* (iii) Step length, t: The step is chosen as the minimum of t1 and t2 */
    m_t = std::min(m_t1, m_t2);
    //m_t = (m_t1 < m_t2) ? m_t1 : m_t2;

    /* Step 2c: Determine new S-pair and take step: ------------------------------------------------*/
    /* case (i): No step in primal or dual space */
    if (m_t >= m_inf)
    {
        /* QPP are infeasible */
        // FIXME: unbounded to raise
        m_fx = m_inf;
        return -1;
    }

    /* case (ii): Step in dual space */
    if (m_t2 >= m_inf)
    {
        /* set u = u +  t * [-r 1] and drop constraint l from the active set A */
        for (cnt = iq >> 2, k = 0; cnt > 0; cnt--, k += 4)
        {
            m_vU(k) -= m_t * m_vR(k);
            m_vU(k + 1) -= m_t * m_vR(k + 1);
            m_vU(k + 2) -= m_t * m_vR(k + 2);
            m_vU(k + 3) -= m_t * m_vR(k + 3);
        }
        for (cnt = iq % 4; cnt > 0; cnt--, k++)
        {
            m_vU(k) -= m_t * m_vR(k);
        }

        m_vU(iq) += m_t;
        m_vIaIncl(l) = l;
        DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, l);

        goto Label_2a;
    }

    /* case (iii): Step in primal and dual space */
    /* set x = x + t * z */
    for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        vX(k) += m_t * m_vZ(k);
        vX(k + 1) += m_t * m_vZ(k + 1);
        vX(k + 2) += m_t * m_vZ(k + 2);
        vX(k + 3) += m_t * m_vZ(k + 3);
    }
    for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
    {
        vX(k) += m_t * m_vZ(k);
    }

    // update the solution value F(x)
    m_fx += m_t * (m_vZ.dot(m_vNp)) * (0.5f * m_t + m_vU(iq));

    // u = u + t * [-r 1]
    for (cnt = iq >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        m_vU(k) -= m_t * m_vR(k);
        m_vU(k + 1) -= m_t * m_vR(k + 1);
        m_vU(k + 2) -= m_t * m_vR(k + 2);
        m_vU(k + 3) -= m_t * m_vR(k + 3);
    }
    for (cnt = iq % 4; cnt > 0; cnt--, k++)
    {
        m_vU(k) -= m_t * m_vR(k);
    }

    m_vU(iq) += m_t;

    /* If t == t2(full step), u = (u^+), add constraint update H and N* and go to Step1 */
    if (std::abs(m_t - m_t2) < std::numeric_limits<m_type>::epsilon())
    {
        // full step has taken
        // add constraint ip to the active set
        if (!AddConstraint(m_mR, m_mJ, m_vD, iq, m_normR))
        {
            m_vIaExcl(ip) = false;
            DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, ip);

            m_vIaIncl.SetBlock(0, m_vConstIaIncl);

            for (i = m_dimP; i < iq; i++)
            {
                m_vA(i) = m_vOldA(i);
                m_vU(i) = m_vOldU(i);
                m_vIaIncl(m_vA(i)) = -1;
            }

            vX = m_vOldX;

            goto Label_2; /* go to step 2 */
        }
        else
            m_vIaIncl(ip) = -1;

        goto Label_1;
    }

    /* a patial step has taken */
    /* If t == t1(partial step), drop constraint l, update H and N*, and go to Step_2a */

    // drop constraint l
    m_vIaIncl(l) = l;
    DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, l);

    // update s[ip] = CI * x + ci0 */
    m_vS(ip) = 0;
    for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        m_vS(ip) += mCi(k, ip) * vX(k);
        m_vS(ip) += mCi(k + 1, ip) * vX(k + 1);
        m_vS(ip) += mCi(k + 2, ip) * vX(k + 2);
        m_vS(ip) += mCi(k + 3, ip) * vX(k + 3);
    }
    for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
        m_vS(ip) += mCi(k, ip) * vX(k);

    m_vS(ip) += vCi(ip);

    goto Label_2a;
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline int8_t CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::Solve(
    const CdtMatrix<m_dimN, m_dimN, m_type>& mG,  const CdtVector<m_dimN, m_type>& vG,
    const CdtMatrix<m_dimN, m_dimP, m_type>& mCe, const CdtVector<m_dimP, m_type>& vCe,
    const CdtMatrix<m_dimN, m_dimM, m_type>& mCi, const CdtVector<m_dimM, m_type>& vCi,
    CdtVector<m_dimN, m_type>& vX)
{
    /* m is the number of inequality constraints */
    /* p is the number of equality constraints */

    int cnt, i, j, k, l; // indices
    int ip; // this is the index of the constraint to be added to the active set
    int iq;
    m_iter = 0;
    m_type ss;
    m_type psi; // this value will contain the sum of all infeasibilities
    m_type sum;

    /* Preprocessing phase -------------------------------------------------------------------------*/
    /* compute the trace of the original matrix G */
    m_c1 = mG.Trace();

    /* decompose the matrix G in the form U^T U */
    m_mL = mG;
    if (CholeskyDecomposition(m_mL)) return -1;
#ifdef DT_QP_DEBUG
    cout << "/* Preprocessing phase */" << endl;
    cout << "// Cholesky Decomposition(G = L*L^t)" << endl;
    cout << "mat L =" << endl; m_mL.Print('\n');
    cout << endl;
#endif

    /* initialize the matrix R */
    m_mR.SetZero();
    m_normR = 1.0; /* this variable will hold the norm of the matrix R */

    /*compute the inverse of the factorized matrix inv(G), this is the initial value for H */
    m_c2 = 0;
    for (i = 0; i < m_dimN; i++)
    {
        m_vD(i) = 1.0;

        ForwardElimination(m_mL, m_vZ, m_vD);

        for (cnt = m_dimN >> 2, j = 0; cnt > 0; cnt--, j += 4)
        {
            m_mJ(i, j) = m_vZ(j);
            m_mJ(i, j + 1) = m_vZ(j + 1);
            m_mJ(i, j + 2) = m_vZ(j + 2);
            m_mJ(i, j + 3) = m_vZ(j + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, j++)
            m_mJ(i, j) = m_vZ(j);

        m_c2 += m_vZ(i); // compute the trace of the matrix J
        m_vD(i) = 0;
    }
    /* c1 * c2 is an estimate for cond(G) */

#ifdef DT_QP_DEBUG
    cout << "// Inverse of L is the initial value for H(=J)" << endl;
    cout << "mat H(=J) =" << endl; m_mJ.Print('\n');
    cout << "// Trace of matrix G and H(=J)" << endl;
    cout << "trace of mat G = " << m_c1 << endl;
    cout << "trace of mat J = " << m_c2 << endl;
    cout << endl << endl;
#endif

    /* Step 0: Find the unconstrained minimum ------------------------------------------------------*/
    /* Find the unconstrained minimizer of the quadratic form 1/2 * xt * G * x + xt * g0
     * this is a feasible point in the dual space
     * G * x = -g0
     * x = - (inv(G) * g0)
     * and compute the current solution value
     * f = 1/2 * xt * G * x + xt * g0
     * f = 1/2 * xt * g0
    */
    CholeskySolve(m_mL, vX, vG);
    vX *= -1;
    m_fx = 0.5f * (vX.dot(vG));
#ifdef DT_QP_DEBUG
    cout << "/* Step 0: Find the unconstrained minimum */" << endl;
    cout << "solution f(x) = " << m_fx << endl;
    cout << "vec x =" << endl; vX.Print('\n');
    cout << endl;
#endif

    /* Add equality constraints to the working set vector A */
#ifdef DT_QP_DEBUG
    cout << "// Add equality constraints to the working set vector A" << endl;
#endif
    iq = 0;
    for (i = 0; i < m_dimP; i++)
    {
        // Get colum vector from CE Matrix
        mCe.GetColVec(i, m_vNp);

        // compute d, d = J^t * np;
        ComputeVecD(m_vD, m_mJ, m_vNp);

        // z = J2 * d2;
        UpdateVecZ(m_vZ, m_mJ, m_vD, iq);

        // r = R^-1 d, R is Upper Triangle Matrix of QR decompostion
        UpdateVecR(m_mR, m_vR, m_vD, iq);

        // compute full step length t2: i.e., the minimum step in primal space s.t. the contraint becomes feasible
        m_t2 = 0;
        if (std::abs(m_vZ.dot(m_vZ)) > std::numeric_limits<m_type>::epsilon()) // i.e. z != 0
            m_t2 = (-(m_vNp.dot(vX)) - vCe(i)) / (m_vZ.dot(m_vNp));

        /* set x = x + t2 * z */
        for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
        {
            vX(k) += m_t2 * m_vZ(k);
            vX(k + 1) += m_t2 * m_vZ(k + 1);
            vX(k + 2) += m_t2 * m_vZ(k + 2);
            vX(k + 3) += m_t2 * m_vZ(k + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
        {
            vX(k) += m_t2 * m_vZ(k);
        }

        /* set u = u+ */
        m_vU(iq) = m_t2;
        for (cnt = iq >> 2, k = 0; cnt > 0; cnt--, k += 4)
        {
            m_vU(k) -= m_t2 * m_vR(k);
            m_vU(k + 1) -= m_t2 * m_vR(k + 1);
            m_vU(k + 2) -= m_t2 * m_vR(k + 2);
            m_vU(k + 3) -= m_t2 * m_vR(k + 3);
        }
        for (cnt = iq % 4; cnt > 0; cnt--, k++)
        {
            m_vU(k) -= m_t2 * m_vR(k);
        }

        /* compute the new solution value */
        m_fx += 0.5f * (m_t2 * m_t2) * (m_vZ.dot(m_vNp));
        m_vA(i) = -i - 1;

        if (!AddConstraint(m_mR, m_mJ, m_vD, iq, m_normR))
        {
            // Equality constraints are linearly dependent
            //throw std::runtime_error("Constraints are linearly dependent");
            //return m_fx;
            return -1;
        }
    }
    /* End: Add equality constraints to the working set A */

    /* set vIaIncl = K \ A */
    m_vIaIncl.SetBlock(0, m_vConstIaIncl);

    /* Step 1: choose a violated constraint --------------------------------------------------------*/
Label_1:
    m_iter++;
#ifdef DT_QP_DEBUG
    cout << "/* Step 1: Choose a violated Constraint */" << endl;
    cout << "iter = " << m_iter << endl;
    cout << "vec x =" << endl; vX.Print('\n');
    cout << endl;
#endif

    for (i = m_dimP; i < iq; i++)
    {
        ip = m_vA(i);
        m_vIaIncl(ip) = -1; // vector vIaIncl(m_dimP + m_dimM)
    }

    /* compute s[x] = Ci^T * x + ci0 for all elements of K \ A */
    psi = 0;	// this value will contain the sum of all infeasibilities
    ip = 0;		// ip will be the index of the chosen violated constraint
    ss = 0;

    for (i = 0; i < m_dimM; i++)
    {
        m_vIaExcl(i) = true;
        sum = 0;

        for (cnt = m_dimN >> 2, j = 0; cnt > 0; cnt--, j += 4)
        {
            sum += mCi(j, i) * vX(j);
            sum += mCi(j + 1, i) * vX(j + 1);
            sum += mCi(j + 2, i) * vX(j + 2);
            sum += mCi(j + 3, i) * vX(j + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, j++)
        {
            sum += mCi(j, i) * vX(j);
        }

        m_vS(i) = sum + vCi(i);
        psi += std::min(0, sum);
        //psi += (m_vS(i) < 0) ? m_vS(i) : 0;
    }
#ifdef DT_QP_DEBUG
    cout << "// compute s[x] = Ci^T * x + ci0 for all elements of K \\ A" << endl;
    cout << "vec s =" << endl; m_vS.Print('\n');
    cout << "sum of all infeasibilites, psi = " << psi << endl;
    cout << endl;
#endif

    /* numerically there are not infeasibilities anymore */
    if (std::abs(psi) <= m_dimM * std::numeric_limits<m_type>::epsilon() * m_c1 * m_c2 * 100.0)
    {
#ifdef DT_QP_DEBUG
        cout << "// numerically there are not infeasibilities anymore" << endl;
        cout << "f(x) = " << m_fx << endl;
        cout << "!! dhQP Bye !!" << endl;
        cout << endl << endl;
#endif
        //if (cost != nullptr) *cost = m_fx;
        return 0;
    }

    /* save old values for u, A and x*/
    for (cnt = iq >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        m_vOldU(i) = m_vU(i); m_vOldA(i) = m_vA(i);
        m_vOldU(i + 1) = m_vU(i + 1); m_vOldA(i + 1) = m_vA(i + 1);
        m_vOldU(i + 2) = m_vU(i + 2); m_vOldA(i + 2) = m_vA(i + 2);
        m_vOldU(i + 3) = m_vU(i + 3); m_vOldA(i + 3) = m_vA(i + 3);
    }
    for (cnt = iq % 4; cnt > 0; cnt--, i++)
    {
        m_vOldU(i) = m_vU(i); m_vOldA(i) = m_vA(i);
    }

    m_vOldX = vX;

    /* Step 2: check for feasibility and determine a new S-pair ------------------------------------*/
Label_2:
    for (cnt = m_dimM >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        if (m_vS(i) < ss && m_vIaIncl(i) != -1 && m_vIaExcl(i))
        {
            ss = m_vS(i); ip = i;
        }
        if (m_vS(i + 1) < ss && m_vIaIncl(i + 1) != -1 && m_vIaExcl(i + 1))
        {
            ss = m_vS(i + 1); ip = i + 1;
        }
        if (m_vS(i + 2) < ss && m_vIaIncl(i + 2) != -1 && m_vIaExcl(i + 2))
        {
            ss = m_vS(i + 2); ip = i + 2;
        }
        if (m_vS(i + 3) < ss && m_vIaIncl(i + 3) != -1 && m_vIaExcl(i + 3))
        {
            ss = m_vS(i + 3); ip = i + 3;
        }
    }
    for (cnt = m_dimM % 4; cnt > 0; cnt--, i++)
    {
        if (m_vS(i) < ss && m_vIaIncl(i) != -1 && m_vIaExcl(i))
        {
            ss = m_vS(i); ip = i;
        }
    }

    if (ss >= 0)
    {
        //if (cost != nullptr) *cost = m_fx;
        return 0;
    }

    /* set np = n[ip] */
    mCi.GetColVec(ip, m_vNp);
    
    /* set u = [u 0]^T */
    m_vU(iq) = 0;

    /* add ip to the active set A */
    m_vA(iq) = ip;

#ifdef DT_QP_DEBUG
    cout << "/* Step 2: check for feasibility and determine a new S-pair */" << endl;
    cout << "Trying with constraint, ip = " << ip << endl;
    cout << "vec np =" << endl; m_vNp.Print('\n');
    cout << endl;
#endif

    /* Step 2a: determine step direction ------------------------------------------------------------*/
Label_2a:
    /* compute z = H np: the step direction in the primal space (through J, see the paper) */
    ComputeVecD(m_vD, m_mJ, m_vNp);
    UpdateVecZ(m_vZ, m_mJ, m_vD, iq);

    /* compute r = N* np (if q > 0): the negative of the step direction in the dual space */
    UpdateVecR(m_mR, m_vR, m_vD, iq);

#ifdef DT_QP_DEBUG
    cout << "/* Step 2a: determine step direction */" << endl;
    cout << "Step direction z" << endl;
    cout << "vec z =" << endl; m_vZ.Print('\n');
    cout << "vec r =" << endl; m_vR.Print('\n');
    cout << "vec u =" << endl; m_vU.Print('\n');
    cout << "vec d =" << endl; m_vD.Print('\n');
    cout << "vec A =" << endl; m_vA.Print('\n');
    cout << endl << endl;
#endif

    /* Step 2b: compute step length ----------------------------------------------------------------*/
    l = 0;

    /* (i) Compute t1: partial step length (maximum step in dual space without violating dual feasibility) */
    m_t1 = m_inf; // +inf
    // find the index l s.t. it reaches the minimum of u+[x] / r
    for (k = m_dimP; k < iq; k++)
    {
        if (m_vR(k) > 0)
        {
            if (m_vU(k) / m_vR(k) < m_t1)
            {
                m_t1 = m_vU(k) / m_vR(k);
                l = m_vA(k);
            }
        }
    }

    /* (ii) Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible) */
    if (std::abs(m_vZ.dot(m_vZ)) > std::numeric_limits<m_type>::epsilon()) // i.e. z != 0
    {
        m_t2 = -m_vS(ip) / (m_vZ.dot(m_vNp));
        if (m_t2 < 0) // patch suggested by Takano Akio for handling numerical inconsistencies
            m_t2 = m_inf;
    }
    else
        m_t2 = m_inf; // +inf
#ifdef DT_QP_DEBUG
    cout << "// (ii) Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible)" << endl;
    cout << "t2 = " << m_t2 << endl;
    cout << "*ref: inf = " << m_inf << endl;
    cout << endl;
#endif

    /* (iii) Step length, t: The step is chosen as the minimum of t1 and t2 */
    m_t = std::min(m_t1, m_t2);
    //m_t = (m_t1 < m_t2) ? m_t1 : m_t2;
#ifdef DT_QP_DEBUG
    cout << "// (iii) Step length, t: The step is chosen as the minimum of t1 and t2" << endl;
    cout << "t = " << m_t << " (t1 = " << m_t1 << ", t2 = " << m_t2 << ") " << endl;
    cout << endl << endl;
#endif

    /* Step 2c: Determine new S-pair and take step: ------------------------------------------------*/
    /* case (i): No step in primal or dual space */
    if (m_t >= m_inf)
    {
        /* QPP are infeasible */
        // FIXME: unbounded to raise
#ifdef DT_QP_DEBUG
        cout << "/* Step 2c: Determine new S-pair and take step */" << endl;
        cout << "// Case (i): No step in primal or dual space" << endl;
        cout << "f(x) = " << m_fx << endl;
        cout << "t = inf" << endl;
        cout << "!! dtQP Bye !!" << endl;
        cout << endl << endl;
#endif
        //if (cost != nullptr) *cost = m_inf;
        //return -1;
        m_fx = m_inf;
        return 0;
    }
#ifdef DT_QP_DEBUG
    cout << "/* Step 2c: Determine new S-pair and take step */" << endl;
    cout << "// (t < inf) -> PASS:(i) No step in primal or dual space" << endl;
    cout << endl;
#endif

    /* case (ii): Step in dual space */
    if (m_t2 >= m_inf)
    {
        /* set u = u +  t * [-r 1] and drop constraint l from the active set A */
        for (cnt = iq >> 2, k = 0; cnt > 0; cnt--, k += 4)
        {
            m_vU(k) -= m_t * m_vR(k);
            m_vU(k + 1) -= m_t * m_vR(k + 1);
            m_vU(k + 2) -= m_t * m_vR(k + 2);
            m_vU(k + 3) -= m_t * m_vR(k + 3);
        }
        for (cnt = iq % 4; cnt > 0; cnt--, k++)
        {
            m_vU(k) -= m_t * m_vR(k);
        }

        m_vU(iq) += m_t;
        m_vIaIncl(l) = l;
        DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, l);
#ifdef DT_QP_DEBUG
        cout << "// Case (ii): Step in dual space" << endl;
        cout << "f(x) = " << m_fx << endl;
        cout << "vec x =" endl; vX.Print('\n');
        cout << "vec A =" endl; m_vA.Print('\n');
        cout << "vec u =" endl; m_vU.Print('\n');
        cout << "mat R =" endl; m_mR.Print('\n');
        cout << "mat J =" endl; m_mJ.Print('\n');
        cout << endl;
#endif
        goto Label_2a; // go to Step 2a: determine step direction
    }
#ifdef DT_QP_DEBUG
    cout << "// (t2 < inf) -> PASS:(ii) Step in dual space" << endl;
    cout << endl << endl;
#endif

    /* case (iii): Step in primal and dual space */
    /* set x = x + t * z */
    for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        vX(k) += m_t * m_vZ(k);
        vX(k + 1) += m_t * m_vZ(k + 1);
        vX(k + 2) += m_t * m_vZ(k + 2);
        vX(k + 3) += m_t * m_vZ(k + 3);
    }
    for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
    {
        vX(k) += m_t * m_vZ(k);
    }

    // update the solution value F(x)
    m_fx += m_t * (m_vZ.dot(m_vNp)) * (0.5f * m_t + m_vU(iq));

    // u = u + t * [-r 1]
    for (cnt = iq >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        m_vU(k) -= m_t * m_vR(k);
        m_vU(k + 1) -= m_t * m_vR(k + 1);
        m_vU(k + 2) -= m_t * m_vR(k + 2);
        m_vU(k + 3) -= m_t * m_vR(k + 3);
    }
    for (cnt = iq % 4; cnt > 0; cnt--, k++)
    {
        m_vU(k) -= m_t * m_vR(k);
    }

    m_vU(iq) += m_t;
#ifdef DT_QP_DEBUG
    cout << "// Case (iii): Step in primal and dual space (both)" << endl;
    cout << "f(x) = " << m_fx << endl;
    cout << "vec x =" endl; vX.Print('\n');
    cout << "vec A =" endl; m_vA.Print('\n');
    cout << "vec u =" endl; m_vU.Print('\n');
    cout << "mat R =" endl; m_mR.Print('\n');
    cout << "mat J =" endl; m_mJ.Print('\n');
    cout << endl;
#endif

    /* If t == t2(full step), u = (u^+), add constraint update H and N* and go to Step1 */
    if (std::abs(m_t - m_t2) < std::numeric_limits<m_type>::epsilon())
    {
#ifdef DT_QP_DEBUG
        cout << "(t == t2) -> Full step : t = " << m_t << endl;
        cout << "vec x =" << endl; vX.Print('\n');
#endif
        // full step has taken
        // add constraint ip to the active set
        if (!AddConstraint(m_mR, m_mJ, m_vD, iq, m_normR))
        {
            m_vIaExcl(ip) = false;
            DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, ip);

            m_vIaIncl.SetBlock(0, m_vConstIaIncl);

            for (i = m_dimP; i < iq; i++)
            {
                m_vA(i) = m_vOldA(i);
                m_vU(i) = m_vOldU(i);
                m_vIaIncl(m_vA(i)) = -1;
            }

            vX = m_vOldX;

            goto Label_2; /* go to step 2 */
        }
        else
            m_vIaIncl(ip) = -1;

        goto Label_1;
    }

    /* a patial step has taken */
    /* If t == t1(partial step), drop constraint l, update H and N*, and go to Step_2a */
#ifdef DT_QP_DEBUG
    cout << "(t == t1) -> Partial step : t = " << m_t << endl;
    cout << "vec x =" << endl; vX.Print('\n');
    cout << endl;
#endif

    // drop constraint l
    m_vIaIncl(l) = l;
    DelConstraint(m_mR, m_mJ, m_vA, m_vU, iq, l);

    // update s[ip] = CI * x + ci0 */
    sum = 0;
    for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k += 4)
    {
        sum += mCi(k, ip) * vX(k);
        sum += mCi(k + 1, ip) * vX(k + 1);
        sum += mCi(k + 2, ip) * vX(k + 2);
        sum += mCi(k + 3, ip) * vX(k + 3);
    }
    for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
        sum += mCi(k, ip) * vX(k);

    m_vS(ip) = sum + vCi(ip);

#ifdef DT_QP_DEBUG
    cout << "updated vec s =" << endl; m_vS.Print('\n');
    cout << "GO TO STEP 2a!!" << endl;
    cout << endl << endl;
#endif
    goto Label_2a; // go to Step 2a: determine step direction
}

/*
template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline m_type CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::DotProduct(const CdtVector<m_dimN, m_type>& x, const CdtVector<m_dimN, m_type>& y)
{
    register int i;
    register m_type sum;

    sum = 0;
    for (i = 0; i < m_dimN; i++)
        sum += x(i) * y(i);
    return sum;
}
*/

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline int8_t CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::CholeskyDecomposition(CdtMatrix<m_dimN, m_dimN, m_type>& A)
{
    int i, j, k;
    m_type sum;

    for (i = 0; i < m_dimN; i++)
    {
        for (j = i; j < m_dimN; j++)
        {
            sum = A(i, j);
            for (k = i - 1; k >= 0; k--)
                sum -= A(i, k) * A(j, k);
            if (i == j)
            {
                if (sum <= 0)
                { // error
                    //std::ostringstream os;
                    //os << "Error in cholesky decomposition, sum: " << sum;
                    //throw std::logic_error(os.str());
                    //exit(-1);
                    return -1; //
                }
                A(i, i) = std::sqrt(sum);
            }
            else
                A(j, i) = sum / A(i, i);
        }
        for (k = i + 1; k < m_dimN; k++)
            A(i, k) = A(k, i);
    }

    return 0;
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline void CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::ForwardElimination(const CdtMatrix<m_dimN, m_dimN, m_type>& L, CdtVector<m_dimN, m_type>& y, const CdtVector<m_dimN, m_type>& b)
{
    int i, j;

    y(0) = b(0) / L(0, 0);
    for (i = 1; i < m_dimN; i++)
    {
        y(i) = b(i);
        for (j = 0; j < i; j++)
            y(i) -= L(i, j) * y(j);
        y(i) = y(i) / L(i, i);
    }
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline void CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::BackwardElimination(const CdtMatrix<m_dimN, m_dimN, m_type>& U, CdtVector<m_dimN, m_type>& x, const CdtVector<m_dimN, m_type>& y)
{
    int i, j;

    x(m_dimN - 1) = y(m_dimN - 1) / U(m_dimN - 1, m_dimN - 1);
    for (i = m_dimN - 2; i >= 0; i--)
    {
        x(i) = y(i);
        for (j = i + 1; j < m_dimN; j++)
            x(i) -= U(i, j) * x(j);
        x(i) = x(i) / U(i, i);
    }
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline void CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::CholeskySolve(const CdtMatrix<m_dimN, m_dimN, m_type>& L, CdtVector<m_dimN, m_type>& x, const CdtVector<m_dimN, m_type>& b)
{
    CdtVector<m_dimN, m_type> y;
    /* Solve L * y = b */
    ForwardElimination(L, y, b);
    /* Solve L^T * x = y */
    BackwardElimination(L, x, y);
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline void CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::ComputeVecD(CdtVector<m_dimN, m_type>& vD, const CdtMatrix<m_dimN, m_dimN, m_type>& mJ, const CdtVector<m_dimN, m_type>& vNp)
{
    int i, j, cnt;
    m_type sum;

    /* compute d = H^T * np */
    for (i = 0; i < m_dimN; i++)
    {
        sum = 0;
        for (cnt = m_dimN >> 2, j = 0; cnt > 0; cnt--, j += 4)
        {
            sum += mJ(j, i) * vNp(j);
            sum += mJ(j + 1, i) * vNp(j + 1);
            sum += mJ(j + 2, i) * vNp(j + 2);
            sum += mJ(j + 3, i) * vNp(j + 3);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, j++)
            sum += mJ(j, i) * vNp(j);

        vD(i) = sum;
    }
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline void CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::UpdateVecZ(CdtVector<m_dimN, m_type>& vZ, const CdtMatrix<m_dimN, m_dimN, m_type>& mJ, const CdtVector<m_dimN, m_type>& vD, const int iq)
{
    int cnt, i, j;

    // setting of z = J2 * d2 = H * d
    for (i = 0; i < m_dimN; i++)
    {
        vZ(i) = 0;
        // iq: 0 ~ (m_dimP-1);
        for (cnt = (m_dimN - iq) >> 2, j = iq; cnt > 0; cnt--, j += 4)
        {
            vZ(i) += mJ(i, j) * vD(j);
            vZ(i) += mJ(i, j + 1) * vD(j + 1);
            vZ(i) += mJ(i, j + 2) * vD(j + 2);
            vZ(i) += mJ(i, j + 3) * vD(j + 3);
        }
        for (cnt = (m_dimN - iq) % 4; cnt > 0; cnt--, j++)
        {
            vZ(i) += mJ(i, j) * vD(j);
        }
    }
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline void CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::UpdateVecR(const CdtMatrix<m_dimN, m_dimN, m_type>& mR, CdtVector<m_dimP + m_dimM, m_type>& vR, const CdtVector<m_dimN, m_type>& vD, const int iq)
{
    // R is Upper Triangle Matrix of QR decompostion
    int cnt, i, j;
    m_type sum;

    /* setting of r = R^-1 d */
    for (i = iq - 1; i >= 0; i--)
    {
        sum = 0;
        //for (j = i + 1; j < iq; j++)
        //    sum += mR(i, j) * vR(j);
        for (cnt = (iq - i - 1) >> 2, j = i + 1; cnt > 0; cnt--, j += 4)
        {
            sum += mR(i, j) * vR(j);
            sum += mR(i, j + 1) * vR(j + 1);
            sum += mR(i, j + 2) * vR(j + 2);
            sum += mR(i, j + 3) * vR(j + 3);
        }
        for (cnt = (iq - i + 1) % 4; cnt > 0; cnt--, j++)
            sum += mR(i, j) * vR(j);

        vR(i) = (vD(i) - sum) / mR(i, i);
    }
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline m_type CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::Distance(m_type a, m_type b)
{
    // Increases numerical accuracy
    m_type absA, absB, t;

    absA = std::abs(a);
    absB = std::abs(b);

    if (absA > absB)
    {
        t = (absB / absA);
        return absA * sqrt(1.0 + t * t);
    }
    else if (absB > absA)
    {
        t = (absA / absB);
        return absB * sqrt(1.0 + t * t);
    }

    return absA * std::sqrt(2.0);
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline bool CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::AddConstraint(CdtMatrix<m_dimN, m_dimN, m_type>& mR, CdtMatrix<m_dimN, m_dimN, m_type>& mJ, CdtVector<m_dimN, m_type>& vD, int &iq, m_type& R_norm)
{
#ifdef DT_QP_DEBUG
    cout << "->  Come in Func. AddConstraint(R, J, d, iq, norm R), with iq(" << iq << ")" << endl;
#endif
    int i, j, k, cnt;
    m_type cc, ss, h, xny;
    m_type t1[4], t2[4];

    /* we have to find the Givens rotation which will reduce the element d(j) to zero.
      if it is already zero we don't have to do anything, except of decreasing j
     */
    for (j = m_dimN - 1; j >= iq + 1; j--)
    {
        /* The Givens rotation is done with the matrix (cc cs, cs -cc).
        If cc is one, then element (j) of d is zero compared with element
        (j - 1). Hence we don't have to do anything.
        If cc is zero, then we just have to switch column (j) and column (j - 1)
        of J. Since we only switch columns in J, we have to be careful how we
        update d depending on the sign of gs.
        Otherwise we have to apply the Givens rotation to these columns.
        The i - 1 element of d has to be updated to h. */
        cc = vD(j - 1);
        ss = vD(j);
        h = Distance(cc, ss);

        if (std::abs(h) < std::numeric_limits<m_type>::epsilon()) // h == 0
            continue;

        vD(j) = 0;
        ss = ss / h;
        cc = cc / h;

        if (cc < 0)
        {
            cc = -cc;
            ss = -ss;
            vD(j - 1) = -h;
        }
        else
            vD(j - 1) = h;

        xny = ss / (1.0 + cc);

        //for (k = 0; k < m_dimN; k++)
        //{
        //    t1 = mJ(k, j - 1);
        //    t2 = mJ(k, j);
        //    mJ(k, j - 1) = t1 * cc + t2 * ss;
        //    mJ(k, j) = xny * (t1 + mJ(k, j - 1)) - t2;
        //}
        for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k+=4)
        {
            t1[0] = mJ(k, j - 1);
            t2[0] = mJ(k, j);
            mJ(k, j - 1) = t1[0] * cc + t2[0] * ss;
            mJ(k, j) = xny * (t1[0] + mJ(k, j - 1)) - t2[0];

            t1[1] = mJ(k + 1, j - 1);
            t2[1] = mJ(k + 1, j);
            mJ(k + 1, j - 1) = t1[1] * cc + t2[1] * ss;
            mJ(k + 1, j) = xny * (t1[1] + mJ(k + 1, j - 1)) - t2[1];

            t1[2] = mJ(k + 2, j - 1);
            t2[2] = mJ(k + 2, j);
            mJ(k + 2, j - 1) = t1[2] * cc + t2[2] * ss;
            mJ(k + 2, j) = xny * (t1[2] + mJ(k + 2, j - 1)) - t2[2];

            t1[3] = mJ(k + 3, j - 1);
            t2[3] = mJ(k + 3, j);
            mJ(k + 3, j - 1) = t1[3] * cc + t2[3] * ss;
            mJ(k + 3, j) = xny * (t1[3] + mJ(k + 3, j - 1)) - t2[3];
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
        {
            t1[0] = mJ(k, j - 1);
            t2[0] = mJ(k, j);
            mJ(k, j - 1) = t1[0] * cc + t2[0] * ss;
            mJ(k, j) = xny * (t1[0] + mJ(k, j - 1)) - t2[0];
        }
    }

    /* update the number of constraints added*/
    iq++;

    /* To update R we have to put the iq components of the d vector
       into column iq - 1 of R
     */
    for (i = 0; i < iq; i++)
        mR(i, iq - 1) = vD(i);

#ifdef DT_QP_DEBUG
    cout << "iq = " << iq << endl;
    cout << "mat J =" << endl; mJ.Print('\n');
    cout << "vec d =" << endl; vD.Print('\n');
    cout << "mat R =" << endl; mR.Print('\n');
    cout << "norm of mat R = " << R_norm << endl;
    cout << endl;
#endif
    if (std::abs(vD(iq - 1)) <= std::numeric_limits<m_type>::epsilon()* R_norm)
    {
        // problem degenerate
#ifdef DT_QP_DEBUG
        cout << "!! Fail Add Constraint !!" << endl;
#endif
        return false;
    }

    R_norm = std::max<m_type>(R_norm, std::abs(vD(iq - 1)));

    return true;
}

template<int m_dimN, int m_dimM, int m_dimP, typename m_type>
inline int8_t CdtQuadProg<m_dimN, m_dimM, m_dimP, m_type>::DelConstraint(CdtMatrix<m_dimN, m_dimN, m_type> &mR, CdtMatrix<m_dimN, m_dimN, m_type> &mJ, CdtVector<m_dimP + m_dimM, int> &vA, CdtVector<m_dimP + m_dimM, m_type> &vU, int &iq, const int l)
{
#ifdef DT_QP_DEBUG
    cout << "->  Come in Func. DelConstraint(R, J, A, u, iq, l), with iq(" << iq << "), l(" << l << ")" << endl;
#endif
    int i, j, k, cnt, qq = 0; // just to prevent warnings from smart compilers
    m_type cc, ss, h, xny;
    m_type t1[4], t2[4];
    bool found = false;

    /* Find the index qq for active constraint l to be removed */
    for (i = m_dimP; i < iq; i++)
    {
        if (vA(i) == l)
        {
            qq = i;
            found = true;
            break;
        }
    }

    if (!found)
    {
        //std::ostringstream os;
        //os << "Attempt to delete non existing constraint, constraint: " << l;
        //throw std::invalid_argument(os.str());
        return -1;
    }

    /* remove the constraint from the active set and the duals */
    for (i = qq; i < iq - 1; i++)
    {
        vA(i) = vA(i + 1);
        vU(i) = vU(i + 1);
        //for (j = 0; j < m_dimN; j++)
        //    mR(j, i) = mR(j, i + 1);
        for (cnt = m_dimN >> 2, j = 0; cnt > 0; cnt--, j += 4)
        {
            mR(j, i) = mR(j, i + 1);
            mR(j + 1, i) = mR(j + 1, i + 1);
            mR(j + 2, i) = mR(j + 2, i + 1);
            mR(j + 3, i) = mR(j + 3, i + 1);
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, j++)
        {
            mR(j, i) = mR(j, i + 1);
        }
    }

    vA(iq - 1) = vA(iq);
    vU(iq - 1) = vU(iq);
    vA(iq) = 0;
    vU(iq) = 0;
    for (j = 0; j < iq; j++)
        mR(j, iq - 1) = 0;

    /* constraint has been fully removed */
    iq--;

    if (iq == 0)
    {
#ifdef DT_QP_DEBUG
        cout << "iq == 0, just return" << endl;
        cout << "vec A =" << endl; vA.Print('\n');
        cout << "vec u =" << endl; vU.Print('\n');
        cout << "mat R =" << endl; mR.Print('\n');
        cout << "<- Come out Func. DelConstraint(R, J, A, u, iq, l), with iq(" << iq << "), l(" << l << ")" << endl;
        cout << endl;
#endif
        return 0;
    }

    for (j = qq; j < iq; j++)
    {
        cc = mR(j, j);
        ss = mR(j + 1, j);
        h = distance(cc, ss);

        if (std::abs(h) < std::numeric_limits<m_type>::epsilon()) // h == 0
            continue;

        cc = cc / h;
        ss = ss / h;
        mR(j + 1, j) = 0;

        if (cc < 0)
        {
            mR(j, j) = -h;
            cc = -cc;
            ss = -ss;
        }
        else
            mR(j, j) = h;

        xny = ss / (1.0f + cc);
        //for (k = j + 1; k < iq; k++)
        //{
        //    t1[0] = mR(j, k);
        //    t2[0] = mR(j + 1, k);
        //    mR(j, k) = t1[0] * cc + t2[0] * ss;
        //    mR(j + 1, k) = xny * (t1[0] + mR(j, k)) - t2[0];
        //}
        for (cnt = (iq - j - 1) >> 2, k = j + 1; cnt > 0; cnt--, k += 4)
        {
            t1[0] = mR(j, k);
            t2[0] = mR(j + 1, k);
            mR(j, k) = t1[0] * cc + t2[0] * ss;
            mR(j + 1, k) = xny * (t1[0] + mR(j, k)) - t2[0];

            t1[1] = mR(j, k + 1);
            t2[1] = mR(j + 1, k + 1);
            mR(j, k + 1) = t1[1] * cc + t2[1] * ss;
            mR(j + 1, k + 1) = xny * (t1[1] + mR(j, k + 1)) - t2[1];

            t1[2] = mR(j, k + 2);
            t2[2] = mR(j + 1, k + 2);
            mR(j, k + 2) = t1[2] * cc + t2[2] * ss;
            mR(j + 1, k + 2) = xny * (t1[2] + mR(j, k + 2)) - t2[2];

            t1[3] = mR(j, k + 3);
            t2[3] = mR(j + 1, k + 3);
            mR(j, k + 3) = t1[3] * cc + t2[3] * ss;
            mR(j + 1, k + 3) = xny * (t1[3] + mR(j, k + 3)) - t2[3];
        }
        for (cnt = (iq - j - 1) % 4; cnt > 0; cnt--, k++)
        {
            t1[0] = mR(j, k);
            t2[0] = mR(j + 1, k);
            mR(j, k) = t1[0] * cc + t2[0] * ss;
            mR(j + 1, k) = xny * (t1[0] + mR(j, k)) - t2[0];
        }

        //for (k = 0; k < m_dimN; k++)
        //{
        //    t1 = mJ(k, j);
        //    t2 = mJ(k, j + 1);
        //    mJ(k, j) = t1 * cc + t2 * ss;
        //    mJ(k, j + 1) = xny * (mJ(k, j) + t1) - t2;
        //}
        for (cnt = m_dimN >> 2, k = 0; cnt > 0; cnt--, k+=4)
        {
            t1[0] = mJ(k, j);
            t2[0] = mJ(k, j + 1);
            mJ(k, j) = t1[0] * cc + t2[0] * ss;
            mJ(k, j + 1) = xny * (mJ(k, j) + t1[0]) - t2[0];

            t1[1] = mJ(k + 1, j);
            t2[1] = mJ(k + 1, j + 1);
            mJ(k + 1, j) = t1[1] * cc + t2[1] * ss;
            mJ(k + 1, j + 1) = xny * (mJ(k + 1, j) + t1[1]) - t2[1];

            t1[2] = mJ(k + 2, j);
            t2[2] = mJ(k + 2, j + 1);
            mJ(k + 2, j) = t1[2] * cc + t2[2] * ss;
            mJ(k + 2, j + 1) = xny * (mJ(k + 2, j) + t1[2]) - t2[2];

            t1[3] = mJ(k + 3, j);
            t2[3] = mJ(k + 3, j + 1);
            mJ(k + 3, j) = t1[3] * cc + t2[3] * ss;
            mJ(k + 3, j + 1) = xny * (mJ(k + 3, j) + t1[3]) - t2[3];
        }
        for (cnt = m_dimN % 4; cnt > 0; cnt--, k++)
        {
            t1[0] = mJ(k, j);
            t2[0] = mJ(k, j + 1);
            mJ(k, j) = t1[0] * cc + t2[0] * ss;
            mJ(k, j + 1) = xny * (mJ(k, j) + t1[0]) - t2[0];
        }
    }

    return 0;
}
