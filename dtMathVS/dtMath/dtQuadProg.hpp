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
public:
    CdtQuadProg();
    ~CdtQuadProg() {}

    int8_t SetObjectFunc(const CdtMatrix<m_dimN, m_dimN, m_type> &mG, const CdtVector<m_dimN, m_type>& vG);
    int8_t UpdateMatrixG(const CdtMatrix<m_dimN, m_dimN, m_type> &mG);
    int8_t UpdateVectorG(const CdtVector<m_dimN, m_type>& vG);
    int8_t UpdateObjectFunc(const CdtMatrix<m_dimN, m_dimN, m_type> &mG, const CdtVector<m_dimN, m_type>& vG);

    // General QuadProg: Consider both equality and inequality constraints.
    int8_t Solve(
        const CdtMatrix<m_dimN, m_dimP, m_type>& mCe, const CdtVector<m_dimP, m_type>& vCe,
        const CdtMatrix<m_dimN, m_dimM, m_type>& mCi, const CdtVector<m_dimM, m_type>& vCi,
        CdtVector<m_dimN, m_type>& vX);

    // Experiment! Modified QuadProg: Consider only inequality constraint.
    int8_t Solve(
        const CdtMatrix<m_dimN, m_dimM, m_type>& mCi, const CdtVector<m_dimM, m_type>& vCi,
        CdtVector<m_dimN, m_type>& vX);

    int GetIteration() { return m_iter; }
    m_type GetObjectValue() { return m_fx; }

private:
    m_type m_inf;
    //CdtMatrix<m_dimN, m_dimN, m_type> m_mG;	// Positive definite symmetric Matrix of the cost function
    CdtVector<m_dimN, m_type> m_vG;             // Vector of the cost function
    CdtMatrix<m_dimN, m_dimN, m_type> m_mL;	    // Cholesky decomposed matrix of the matrix mG
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
    m_type m_psiThreshold;

    int m_iter;

private:
    //m_type DotProduct(const CdtVector<m_dimN, m_type>& x, const CdtVector<m_dimN, m_type>& y);
    int8_t CholeskyDecomposition(CdtMatrix<m_dimN, m_dimN, m_type>& A);
    void ForwardElimination(const CdtMatrix<m_dimN, m_dimN, m_type>& L, CdtVector<m_dimN, m_type>& y, const CdtVector<m_dimN, m_type>& b);
    void BackwardElimination(const CdtMatrix<m_dimN, m_dimN, m_type>& U, CdtVector<m_dimN, m_type>& x, const CdtVector<m_dimN, m_type>& y);
    void CholeskySolve(const CdtMatrix<m_dimN, m_dimN, m_type>& L, CdtVector<m_dimN, m_type>& x, const CdtVector<m_dimN, m_type>& b);

    void ComputeVecD(CdtVector<m_dimN, m_type>& vD, const CdtMatrix<m_dimN, m_dimN, m_type>& mJ, const CdtVector<m_dimN, m_type>& vNp); // d = J^T * np
    void UpdateVecZ(CdtVector<m_dimN, m_type>& vZ, const CdtMatrix<m_dimN, m_dimN, m_type>& mJ, const CdtVector<m_dimN, m_type>& vD, const int iq); // z = J2 * d2
    void UpdateVecR(const CdtMatrix<m_dimN, m_dimN, m_type>& mR, CdtVector<m_dimP + m_dimM, m_type>& vR, const CdtVector<m_dimN, m_type>& vD, const int iq); // r = R^-1 d

    m_type Distance(m_type a, m_type b); //the euclidean distance between two numbers
    bool AddConstraint(CdtMatrix<m_dimN, m_dimN, m_type>& mR, CdtMatrix<m_dimN, m_dimN, m_type>& mJ, CdtVector<m_dimN, m_type>& vD, int &iq, m_type& R_norm);
    int8_t DelConstraint(CdtMatrix<m_dimN, m_dimN, m_type>& mR, CdtMatrix<m_dimN, m_dimN, m_type>& mJ, CdtVector<m_dimP + m_dimM, int>& vA, CdtVector<m_dimP + m_dimM, m_type>& vU, int &iq, const int l);
};

#include "dtQuadProg.ipp"
