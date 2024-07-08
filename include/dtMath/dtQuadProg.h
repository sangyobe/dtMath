/*!
\file	    dtQuadProg.h
\brief	    Quadratic Programming Solver
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
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

#ifndef DTMATH_DTQUAD_PROG_H_
#define DTMATH_DTQUAD_PROG_H_

#if defined(_WIN32) || defined(__linux__)
#include <cstdlib>
#include <iostream>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type> class Matrix;
template <uint16_t t_row, typename t_type> class Vector;

// #define DT_QP_DEBUG

#ifdef DT_QP_DEBUG
using namespace std;
#endif

template <int m_dimN, int m_dimM, int m_dimP, typename t_type = float>
class dtQuadProg
{
public:
    dtQuadProg();
    ~dtQuadProg() {}

    int8_t SetObjectFunc(const Matrix<m_dimN, m_dimN, t_type> &mG, const Vector<m_dimN, t_type> &vG);
    int8_t UpdateMatrixG(const Matrix<m_dimN, m_dimN, t_type> &mG);
    int8_t UpdateVectorG(const Vector<m_dimN, t_type> &vG);
    int8_t UpdateObjectFunc(const Matrix<m_dimN, m_dimN, t_type> &mG, const Vector<m_dimN, t_type> &vG);

    // General QuadProg: Consider both equality and inequality constraints.
    int8_t Solve(
        const Matrix<m_dimN, m_dimP, t_type> &mCe, const Vector<m_dimP, t_type> &vCe,
        const Matrix<m_dimN, m_dimM, t_type> &mCi, const Vector<m_dimM, t_type> &vCi,
        Vector<m_dimN, t_type> &vX);

    // Experiment! Modified QuadProg: Consider only inequality constraint.
    int8_t Solve(
        const Matrix<m_dimN, m_dimM, t_type> &mCi, const Vector<m_dimM, t_type> &vCi,
        Vector<m_dimN, t_type> &vX);

    int GetIteration() { return m_iter; }
    t_type GetObjectValue() { return m_fx; }

private:
    t_type m_inf;
    // Matrix<m_dimN, m_dimN, t_type> m_mG;	// Positive definite symmetric Matrix of the cost function
    Vector<m_dimN, t_type> m_vG;         // Vector of the cost function
    Matrix<m_dimN, m_dimN, t_type> m_mL; // Cholesky decomposed matrix of the matrix mG
    Matrix<m_dimN, m_dimN, t_type> m_mJ0, m_mJ, m_mR;
    Vector<m_dimN, t_type> m_vX0;

    Vector<m_dimN, t_type> m_vZ, m_vD, m_vNp, m_vOldX;
    Vector<m_dimM + m_dimP, t_type> m_vS, m_vR, m_vU, m_vOldU;
    Vector<m_dimM + m_dimP, int> m_vA, m_vOldA; // Active Set A
    Vector<m_dimM + m_dimP, int> m_vIaIncl;
    Vector<m_dimM, int> m_vConstIaIncl;
    Vector<m_dimM + m_dimP, bool> m_vIaExcl;

    t_type m_fx, m_fx0;
    t_type m_c1; // c1 is trace of G matrix, c1 * c2 is an estimate for cond(G)
    t_type m_c2; // c2 is trace of J matrix, c1 * c2 is an estimate for cond(G)
    t_type m_normR;

    t_type m_t;        // t is the step lenght, which is the minimum of
    t_type m_t1, m_t2; // the partial step length t1 and the full step length t2
    t_type m_psiThreshold;

    int m_iter;

private:
    // t_type DotProduct(const Vector<m_dimN, t_type>& x, const Vector<m_dimN, t_type>& y);
    int8_t CholeskyDecomposition(Matrix<m_dimN, m_dimN, t_type> &A);
    void ForwardElimination(const Matrix<m_dimN, m_dimN, t_type> &L, Vector<m_dimN, t_type> &y, const Vector<m_dimN, t_type> &b);
    void BackwardElimination(const Matrix<m_dimN, m_dimN, t_type> &U, Vector<m_dimN, t_type> &x, const Vector<m_dimN, t_type> &y);
    void CholeskySolve(const Matrix<m_dimN, m_dimN, t_type> &L, Vector<m_dimN, t_type> &x, const Vector<m_dimN, t_type> &b);

    void ComputeVecD(Vector<m_dimN, t_type> &vD, const Matrix<m_dimN, m_dimN, t_type> &mJ, const Vector<m_dimN, t_type> &vNp);                      // d = J^T * np
    void UpdateVecZ(Vector<m_dimN, t_type> &vZ, const Matrix<m_dimN, m_dimN, t_type> &mJ, const Vector<m_dimN, t_type> &vD, const int iq);          // z = J2 * d2
    void UpdateVecR(const Matrix<m_dimN, m_dimN, t_type> &mR, Vector<m_dimP + m_dimM, t_type> &vR, const Vector<m_dimN, t_type> &vD, const int iq); // r = R^-1 d

    t_type Distance(t_type a, t_type b); // the euclidean distance between two numbers
    bool AddConstraint(Matrix<m_dimN, m_dimN, t_type> &mR, Matrix<m_dimN, m_dimN, t_type> &mJ, Vector<m_dimN, t_type> &vD, int &iq, t_type &R_norm);
    int8_t DelConstraint(Matrix<m_dimN, m_dimN, t_type> &mR, Matrix<m_dimN, m_dimN, t_type> &mJ, Vector<m_dimP + m_dimM, int> &vA, Vector<m_dimP + m_dimM, t_type> &vU, int &iq, const int l);
};

} // namespace Math
} // namespace dt

#include "dtQuadProg.tpp"

#endif // DTMATH_DTQUAD_PROG_H_
