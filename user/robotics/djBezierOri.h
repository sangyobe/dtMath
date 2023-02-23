
#ifndef DTMATH_DJBEZIER_ORI_H_
#define DTMATH_DJBEZIER_ORI_H_

#ifndef __Bezier_H__
#define __Bezier_H__

#include <dtMath/dtMath.h>

#include <float.h>
#include <iostream>
using namespace std;

// ncol: space dimension
// mrow: the number of control points

template <uint16_t ncol>
struct output
{
    CdtVector<ncol> pos;
    CdtVector<ncol> vel;
};

template <uint16_t mrow, uint16_t ncol>
class Bezier
{
public:

    float BinomialCoeff(const uint16_t n, const uint16_t k);
    float BernsteinPoly(const uint16_t n, const uint16_t k, float ctrlParam);

    output<ncol> BezierPoly(CdtMatrix<mrow, ncol> ctrlPoints, float ctrlParam, float period);
    output<ncol> BezierInterp(CdtMatrix<mrow, ncol> ctrlPoints, float ctrlParam, float period);


};

#include "djBezierOri.tpp"


#endif

#endif // DTMATH_DJBEZIER_ORI_H_

