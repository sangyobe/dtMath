#if defined(_WIN32) || defined(__linux__)
#include <stdio.h>
#include <iostream>
#include "dhTimeCheck.hpp"
#include "dhTerm.h"
using namespace std;
#endif /* _WIN32 || __linux__ */

#if defined(_WIN32)
#include "./dtMath/dtGnuPlot.hpp"
#endif

#include "testPrint.h"
#include "testQuadProg.h"
#include "./dtMath/dtMath.h"

void Test_QuadProg()
{
    PrintTitle("Test Quadratic Programming Problem");
    dtMathQuadProg();
}

void dtMathQuadProg()
{
    PrintHeading("dtMath Quadratic Programming Solve ");

    double cost = 0;

    CdtQuadProg<2, 3, 1, double> qp;
    CdtMatrix<2, 2, double> mG;
    CdtVector<2, double> vG0;
    CdtMatrix<2, 1, double> mCe;
    CdtVector<1, double> vCe0;
    CdtMatrix<2, 3, double> mCi;
    CdtVector<3, double> vCi0;
    CdtVector<2, double> vX;

    mG << 4.0, -2.0, -2.0, 4.0;
    vG0 << 6.0, 0.0;
    mCe << 1.0, 1.0;
    vCe0 << -3.0;
    mCi << 1.0, 0.0, 1.0, 0.0, 1.0, 1.0;
    vCi0 << 0.0, 0.0, -2.0;    

    // Solution is 
    // F(x): 12
    // X: [1, 2]^T
    qp.SetObjectFunc(mG, vG0);
    if (qp.Solve(mCe, vCe0, mCi, vCi0, vX))
    {
        Printf("dtMath QuadProg Solve Error!\n");
    }
    else
    {
        Printf("Result:\n");
        Printf("f(x) = %f\n", cost);
        Printf("solution X =\n");
        vX.Print();
        Println;
    }
}

#if defined(_WIN32) || defined(__linux__)
void Test_QuadrupedRobotQP()
{
    PrintTitle("Test Quadruped Robot QP");
    //dtMathQuadrupedRobotQP();
    //dtMathQuadrupedRobotQP_DJ();
    //dtMathQuadrupedRobotQP_DH();
    dtMathMoBedQP();
}

void dtMathQuadrupedRobotQP()
{
    PrintHeading("Quadruped Robot Example using dtMath ");
    int curX, curY;
    bool fRun = true;
    char ch = 0;

    CdhTimeCheck time;
    int8_t rtn = 0;
    double elapseTime = 0;
    double maxTime = 0;
    double minTime = 1000;

    double mu = 0.5;
    double maxF = 100;

    // n = 18, m = 24, p = 6
    CdtVector<18, double> vX;       // Vector x = [ fFL, fBL, fFR, fBR, A*f-b ]^T, x:R(12+6)

    CdtMatrix<18, 18, double> mG;   // mG = [ W1, 0 ; 0 W2 ]
    CdtVector<18, double> vG;       // vG0 is 0

    CdtMatrix<3, 6, double> mSubCI;
    CdtMatrix<18, 24, double> mCI;  // diag = [CIFL, CIBL, CIFR, CIBR, 0]
    CdtVector<24, double> vCi;

    CdtVector3<double> rFR, rBL, rFL, rBR;
    CdtMatrix<6, 12, double> mA;    // mA = [  I3   I3   I3   I3  ]
                                    //      [ rFLx rBLx rFRx rBRx ]    
    CdtMatrix<6, 6, double> mI6;
    CdtVector<6, double> vB;        // Vector b = [ fd, td ]^T, b:R(6)
    CdtMatrix<18, 6, double> mCE;   // = [A -I]^T
    CdtVector<6, double> vCe;       // = -vB

    CdtQuadProg<18, 24, 6, double> qp;
    CdtVector<12, double> fc;

    //double w[18] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 10, 10, 100, 1, 1, 1 };
    double w[18] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 100, 100, 100, 50, 50, 50 };
    mG.SetDiagonal(w, sizeof(w));
    //cout << "mG:\n" << mG << endl;

    vG.SetZero();
    //cout << "vG:\n" << vG << endl << endl;

    mSubCI <<
        1, -1, 0, 0, 0, 0,
        0, 0, 1, -1, 0, 0,
        mu, mu, mu, mu, 1, -1;
    //cout << "mSubCI:\n" << mSubCI << endl << endl;

    mCI.SetBlock(0, 0, mSubCI);
    mCI.SetBlock(3, 6, mSubCI);
    mCI.SetBlock(6, 12, mSubCI);
    mCI.SetBlock(9, 18, mSubCI);
    //cout << "mCI:\n" << mCI << endl << endl;

    vCi << 0, 0, 0, 0, 0, maxF, 0, 0, 0, 0, 0, maxF, 0, 0, 0, 0, 0, maxF, 0, 0, 0, 0, 0, maxF;
    //cout << "vCi:\n" << vCi << endl << endl;

    rFR << 2, -1, -1;
    rBL << -2, 1, -1;
    rFL << 2, 1, -1;
    rBR << -2, -1, -1;

    mA <<
        1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
        0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
        0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
        0, -rFL(2), rFL(1), 0, -rBL(2), rBL(1), 0, -rFR(2), rFR(1), 0, -rBR(2), rBR(1),
        rFL(2), 0, -rFL(0), rBL(2), 0, -rBL(0), rFR(2), 0, -rFR(0), rBR(2), 0, -rBR(0),
        -rFL(1), rFL(0), 0, -rBL(1), rBL(0), 0, -rFR(1), rFR(0), 0, -rBR(1), rBR(0), 0;
    //cout << "mA:\n" << mA << endl << endl;

    mI6.SetIdentity();
    mCE.SetBlock(0, 0, mA.Transpose());
    mCE.SetBlock(12, 0, -mI6);
    //cout << "mCE:\n" << mCI << endl << endl;

    vB << 0, 0, 100, 0, 0, 0;
    //cout << "vB:\n" << vB << endl << endl;

    vCe = -vB;
    //cout << "vCe:\n" << vCe << endl << endl;

    qp.SetObjectFunc(mG, vG);

    dhTerm::GetCurPos(&curY, &curX);
    dhTerm::print(curY + 1, 0, "rFL = [%6.3f %6.3f %6.3f]", rFL(0), rFL(1), rFL(2));
    dhTerm::print(curY + 2, 0, "rBL = [%6.3f %6.3f %6.3f]", rBL(0), rBL(1), rBL(2));
    dhTerm::print(curY + 3, 0, "rFR = [%6.3f %6.3f %6.3f]", rFR(0), rFR(1), rFR(2));
    dhTerm::print(curY + 4, 0, "rBR = [%6.3f %6.3f %6.3f]", rBR(0), rBR(1), rBR(2));
    dhTerm::print(curY + 5, 0, "  b = [%6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f]", vB(0), vB(1), vB(2), vB(3), vB(4), vB(5));

    while (fRun)
    {
        if (dhTerm::kbhit())
        {
            ch = getchar();
            switch (ch)
            {
            case 'q': rFL(0) += 0.1; break;  case 'w': rFL(0) -= 0.1; break;
            case 'a': rFL(1) += 0.1; break;  case 's': rFL(1) -= 0.1; break;
            case 'z': rFL(2) += 0.1; break;  case 'x': rFL(2) -= 0.1; break;

            case 'e': rFR(0) += 0.1; break;  case 'r': rFR(0) -= 0.1; break;
            case 'd': rFR(1) += 0.1; break;  case 'f': rFR(1) -= 0.1; break;
            case 'c': rFR(2) += 0.1; break;  case 'v': rFR(2) -= 0.1; break;

            case 't': rBL(0) += 0.1; break;  case 'y': rBL(0) -= 0.1; break;
            case 'g': rBL(1) += 0.1; break;  case 'h': rBL(1) -= 0.1; break;
            case 'b': rBL(2) += 0.1; break;  case 'n': rBL(2) -= 0.1; break;

            case 'u': rBR(0) += 0.1; break;  case 'i': rBR(0) -= 0.1; break;
            case 'j': rBR(1) += 0.1; break;  case 'k': rBR(1) -= 0.1; break;
            case 'm': rBR(2) += 0.1; break;  case ',': rBR(2) -= 0.1; break;

            case 'Q': vB(0) += 5; break;  case 'W': vB(0) -= 5; break;
            case 'A': vB(1) += 5; break;  case 'S': vB(1) -= 5; break;
            case 'Z': vB(2) += 5; break;  case 'X': vB(2) -= 5; break;

            case 'E': vB(3) += 5; break;  case 'R': vB(3) -= 5; break;
            case 'D': vB(4) += 5; break;  case 'F': vB(4) -= 5; break;
            case 'C': vB(5) += 5; break;  case 'V': vB(5) -= 5; break;

            case '!': fRun = false; break;
            }
            dhTerm::print(curY + 1, 0, "rFL = [%6.3f %6.3f %6.3f]", rFL(0), rFL(1), rFL(2));
            dhTerm::print(curY + 2, 0, "rBL = [%6.3f %6.3f %6.3f]", rBL(0), rBL(1), rBL(2));
            dhTerm::print(curY + 3, 0, "rFR = [%6.3f %6.3f %6.3f]", rFR(0), rFR(1), rFR(2));
            dhTerm::print(curY + 4, 0, "rBR = [%6.3f %6.3f %6.3f]", rBR(0), rBR(1), rBR(2));
            dhTerm::print(curY + 5, 0, "  b = [%6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f]", vB(0), vB(1), vB(2), vB(3), vB(4), vB(5));
        }

        /// Update Matrix
        mCE.SetBlock(0, 3, rFL.GetSkew().Transpose());
        mCE.SetBlock(3, 3, rBL.GetSkew().Transpose());
        mCE.SetBlock(6, 3, rFR.GetSkew().Transpose());
        mCE.SetBlock(9, 3, rBR.GetSkew().Transpose());
        vCe = -vB;

        /// Solve QP
        time.Start();
        rtn = qp.Solve(mCE, vCe, mCI, vCi, vX);
        time.Stop();

        /// Print Data
        dhTerm::print(curY + 6, 0, "Solver Result: %s", (rtn) ? "Fail   " : "Success");
        dhTerm::print(curY + 7, 0, "Cost Func Result: %f", qp.GetObjectValue());
        dhTerm::print(curY + 8, 0, "Xf_FL = [%6.3f %6.3f %6.3f]", vX(0), vX(1), vX(2));
        dhTerm::print(curY + 9, 0, "Xf_BL = [%6.3f %6.3f %6.3f]", vX(3), vX(4), vX(5));
        dhTerm::print(curY + 10, 0, "Xf_FR = [%6.3f %6.3f %6.3f]", vX(6), vX(7), vX(8));
        dhTerm::print(curY + 11, 0, "Xf_BR = [%6.3f %6.3f %6.3f]", vX(9), vX(10), vX(11));
        dhTerm::print(curY + 12, 0, "Xd_TR = [%6.3f %6.3f %6.3f]", vX(12), vX(13), vX(14));
        dhTerm::print(curY + 13, 0, "Xd_OR = [%6.3f %6.3f %6.3f]", vX(15), vX(16), vX(17));
        dhTerm::print(curY + 14, 0, "elapse time: %6.3f [ms]", time.GetElapsedTime_msec());

        vX.GetBlock(0, fc);
        (mA * fc - vB).Print();
    }
    Println;
    Println;
}

void dtMathQuadrupedRobotQP_DJ()
{// ver (f-f_prev) of DJ
    PrintHeading("Quadruped Robot Example using dtMath, ver (f-f_prev) of DJ");
    int curX, curY;
    bool fRun = true;
    char ch = 0;

    CdhTimeCheck time;
    int8_t rtn = 0;
    double elapseTime = 0;
    double maxTime = 0;
    double minTime = 100;

    double mu = 0.5;
    double maxF = 1000;

    // n = 18, m = 24, p = 6
    CdtVector3<double> vXprev[4];
    CdtVector<18, double> vX;       // Vector x = [ fFL, fBL, fFR, fBR, A*f-b ]^T, x:R(12+6)

    CdtMatrix<18, 18, double> mG;   // mG = [ W1, 0 ; 0 W2 ]
    CdtVector<18, double> vG;       // vG0 is 0

    CdtMatrix<3, 6, double> mSubCI;
    CdtMatrix<18, 24, double> mCI;  // diag = [CIFL, CIBL, CIFR, CIBR, 0]
    CdtVector<24, double> vCi;

    CdtVector3<double> rFR, rBL, rFL, rBR;
    CdtMatrix<6, 12, double> mA;    // mA = [  I3   I3   I3   I3  ]
                                    //      [ rFLx rBLx rFRx rBRx ]    
    CdtMatrix<6, 6, double> mI6;
    CdtVector<6, double> vB;        // Vector b = [ fd, td ]^T, b:R(6)
    CdtMatrix<18, 6, double> mCE;   // = [A -I]^T
    CdtVector<6, double> vCe;       // = -vB

    CdtQuadProg<18, 24, 6, double> qp;

    double Gc = 0.1;
    double Gd = 0.1;
    //double w[18] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 10, 10, 100, 1, 1, 1 };
    double w[18] = { Gc + Gd, Gc + Gd, Gc + Gd, Gc + Gd, Gc + Gd, Gc + Gd, Gc + Gd, Gc + Gd, Gc + Gd, Gc + Gd, Gc + Gd, Gc + Gd, 100, 100, 100, 50, 50, 50 };
    mG.SetDiagonal(w, sizeof(w));
    //cout << "mG:\n" << mG << endl;

    vG.SetZero();
    //cout << "vG:\n" << vG << endl << endl;

    mSubCI <<
        1, -1, 0, 0, 0, 0,
        0, 0, 1, -1, 0, 0,
        mu, mu, mu, mu, 1, -1;
    //cout << "mSubCI:\n" << mSubCI << endl << endl;

    mCI.SetBlock(0, 0, mSubCI);
    mCI.SetBlock(3, 6, mSubCI);
    mCI.SetBlock(6, 12, mSubCI);
    mCI.SetBlock(9, 18, mSubCI);
    //cout << "mCI:\n" << mCI << endl << endl;

    vCi << 0, 0, 0, 0, 0, maxF, 0, 0, 0, 0, 0, maxF, 0, 0, 0, 0, 0, maxF, 0, 0, 0, 0, 0, maxF;
    //cout << "vCi:\n" << vCi << endl << endl;

    rFL << 2, 1, -1;
    rBL << -2, 1, -1;
    rFR << 2, -1, -1;
    rBR << -2, -1, -1;

    mA <<
        1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
        0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
        0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
        0, -rFL(2), rFL(1), 0, -rBL(2), rBL(1), 0, -rFR(2), rFR(1), 0, -rBR(2), rBR(1),
        rFL(2), 0, -rFL(0), rBL(2), 0, -rBL(0), rFR(2), 0, -rFR(0), rBR(2), 0, -rBR(0),
        -rFL(1), rFL(0), 0, -rBL(1), rBL(0), 0, -rFR(1), rFR(0), 0, -rBR(1), rBR(0), 0;
    //cout << "mA:\n" << mA << endl << endl;

    mI6.SetIdentity();
    mCE.SetBlock(0, 0, mA.Transpose());
    mCE.SetBlock(12, 0, -mI6);
    //cout << "mCE:\n" << mCI << endl << endl;

    vB << 0, 0, 100, 0, 0, 0;
    //cout << "vB:\n" << vB << endl << endl;

    vCe = -vB;
    //cout << "vCe:\n" << vCe << endl << endl;

    qp.SetObjectFunc(mG, vG);

    dhTerm::GetCurPos(&curY, &curX);
    dhTerm::print(curY + 1, 0, "rFL = [%6.3f %6.3f %6.3f]", rFL(0), rFL(1), rFL(2));
    dhTerm::print(curY + 2, 0, "rBL = [%6.3f %6.3f %6.3f]", rBL(0), rBL(1), rBL(2));
    dhTerm::print(curY + 3, 0, "rFR = [%6.3f %6.3f %6.3f]", rFR(0), rFR(1), rFR(2));
    dhTerm::print(curY + 4, 0, "rBR = [%6.3f %6.3f %6.3f]", rBR(0), rBR(1), rBR(2));
    dhTerm::print(curY + 5, 0, "  b = [%6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f]", vB(0), vB(1), vB(2), vB(3), vB(4), vB(5));

    while (fRun)
    {
        if (dhTerm::kbhit())
        {
            ch = getchar();
            switch (ch)
            {
            case 'q': rFL(0) += 0.1; break;  case 'w': rFL(0) -= 0.1; break;
            case 'a': rFL(1) += 0.1; break;  case 's': rFL(1) -= 0.1; break;
            case 'z': rFL(2) += 0.1; break;  case 'x': rFL(2) -= 0.1; break;

            case 'e': rFR(0) += 0.1; break;  case 'r': rFR(0) -= 0.1; break;
            case 'd': rFR(1) += 0.1; break;  case 'f': rFR(1) -= 0.1; break;
            case 'c': rFR(2) += 0.1; break;  case 'v': rFR(2) -= 0.1; break;

            case 't': rBL(0) += 0.1; break;  case 'y': rBL(0) -= 0.1; break;
            case 'g': rBL(1) += 0.1; break;  case 'h': rBL(1) -= 0.1; break;
            case 'b': rBL(2) += 0.1; break;  case 'n': rBL(2) -= 0.1; break;

            case 'u': rBR(0) += 0.1; break;  case 'i': rBR(0) -= 0.1; break;
            case 'j': rBR(1) += 0.1; break;  case 'k': rBR(1) -= 0.1; break;
            case 'm': rBR(2) += 0.1; break;  case ',': rBR(2) -= 0.1; break;

            case 'Q': vB(0) += 5; break;  case 'W': vB(0) -= 5; break;
            case 'A': vB(1) += 5; break;  case 'S': vB(1) -= 5; break;
            case 'Z': vB(2) += 5; break;  case 'X': vB(2) -= 5; break;

            case 'E': vB(3) += 5; break;  case 'R': vB(3) -= 5; break;
            case 'D': vB(4) += 5; break;  case 'F': vB(4) -= 5; break;
            case 'C': vB(5) += 5; break;  case 'V': vB(5) -= 5; break;

            case '!': fRun = false; break;
            }
            dhTerm::print(curY + 1, 0, "rFL = [%6.3f %6.3f %6.3f]", rFL(0), rFL(1), rFL(2));
            dhTerm::print(curY + 2, 0, "rBL = [%6.3f %6.3f %6.3f]", rBL(0), rBL(1), rBL(2));
            dhTerm::print(curY + 3, 0, "rFR = [%6.3f %6.3f %6.3f]", rFR(0), rFR(1), rFR(2));
            dhTerm::print(curY + 4, 0, "rBR = [%6.3f %6.3f %6.3f]", rBR(0), rBR(1), rBR(2));
            dhTerm::print(curY + 5, 0, "  b = [%6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f]", vB(0), vB(1), vB(2), vB(3), vB(4), vB(5));
        }

        /* Update Matrix */
        mCE.SetBlock(0, 3, rFL.GetSkew().Transpose());
        mCE.SetBlock(3, 3, rBL.GetSkew().Transpose());
        mCE.SetBlock(6, 3, rFR.GetSkew().Transpose());
        mCE.SetBlock(9, 3, rBR.GetSkew().Transpose());
        vCe = -vB;

        vX.GetBlockVec3(0, vXprev[0]);
        vX.GetBlockVec3(3, vXprev[1]);
        vX.GetBlockVec3(6, vXprev[2]);
        vX.GetBlockVec3(9, vXprev[3]);

        vG.SetBlock(0, -2 * Gd*vXprev[0]);
        vG.SetBlock(3, -2 * Gd*vXprev[1]);
        vG.SetBlock(6, -2 * Gd*vXprev[2]);
        vG.SetBlock(9, -2 * Gd*vXprev[3]);

        /* Solve QP */
        time.Start();
        //rtn = qp.Solve(mG, vG, mCE, vCe, mCI, vCi, vX, &cost);
        qp.UpdateVectorG(vG);
        rtn = qp.Solve(mCE, vCe, mCI, vCi, vX);
        time.Stop();

        /* Print Data */
        dhTerm::print(curY + 6, 0, "Solver Result: %s", (rtn) ? "Fail   " : "Success");
        dhTerm::print(curY + 7, 0, "Cost Func Result: %f", qp.GetObjectValue());
        dhTerm::print(curY + 8, 0, "Xf_FL = [%6.3f %6.3f %6.3f]", vX(0), vX(1), vX(2));
        dhTerm::print(curY + 9, 0, "Xf_BL = [%6.3f %6.3f %6.3f]", vX(3), vX(4), vX(5));
        dhTerm::print(curY + 10, 0, "Xf_FR = [%6.3f %6.3f %6.3f]", vX(6), vX(7), vX(8));
        dhTerm::print(curY + 11, 0, "Xf_BR = [%6.3f %6.3f %6.3f]", vX(9), vX(10), vX(11));
        dhTerm::print(curY + 12, 0, "Xd_TR = [%6.3f %6.3f %6.3f]", vX(12), vX(13), vX(14));
        dhTerm::print(curY + 13, 0, "Xd_OR = [%6.3f %6.3f %6.3f]", vX(15), vX(16), vX(17));
        dhTerm::print(curY + 14, 0, "elapse time: %6.3f [ms]", time.GetElapsedTime_msec());
        if (rtn != 0) break;
    }
    Println;
    Println;
}

void dtMathQuadrupedRobotQP_DH()
{// ver (f-f_prev) of DJ
    PrintHeading("Quadruped Robot Example using dtMath ver (f-f_prev) of DH");
    int curX, curY;
    bool fRun = true;
    char ch = 0;

    CdhTimeCheck time;
    int8_t rtn = 0;
    double elapseTime = 0;
    double maxTime = 0;
    double minTime = 100;

    double mu = 0.5;
    double maxF = 1;

    // n = 30, m = 24, p = 18
    CdtVector<30, double> vX;       // Vector x = [ fFR, fBL, fFL, fBR, F-Fprev, A*f-b]^T, x:R(30=12+12+6)
    CdtVector<12, double> vFprev;

    CdtMatrix<30, 30, double> mG;   // mG = [ Gc, 0 0; 0 Gd 0 ; 0 0 Ge ], mG:R(30x30)
    CdtVector<30, double> vG;       // vG0 is 0

    CdtMatrix<3, 6, double> mSubCI;
    CdtMatrix<30, 24, double> mCI;  // diag = [CIFR, CIBL, CIFL, CIBR, 0]
    CdtVector<24, double> vCi;

    CdtVector3<double> rFR, rBL, rFL, rBR;
    CdtMatrix<6, 12, double> mA;    // mA = [  I3   I3   I3   I3  ]
                                    //      [ rFRx rBLx rFLx rBLx ]
    CdtMatrix<6, 6, double> mI6;
    CdtMatrix<12, 12, double> mI12;
    CdtVector<6, double> vB;        // Vector b = [ fd, td ]^T, b:R(6)
    CdtMatrix<30, 18, double> mCE;  // = [A  0 -I]^T
                                    //   [I -I  0]
    CdtVector<18, double> vCe;      // = -[vB; vFprev]

    CdtQuadProg<30, 24, 18, double> qp;

    double Gc = 0.1;
    double Gd = 10.0;
    double w[30] = {
        Gc, Gc, Gc, Gc, Gc, Gc, Gc, Gc, Gc, Gc, Gc, Gc,
        Gd, Gd, Gd, Gd, Gd, Gd, Gd, Gd, Gd, Gd, Gd, Gd,
        100, 100, 100, 50, 50, 50 };
    mG.SetDiagonal(w, sizeof(w));
    //cout << "mG:\n" << mG << endl;

    vG.SetZero();
    //cout << "vG:\n" << vG << endl << endl;

    mSubCI <<
         1, -1,  0,  0, 0, 0,
         0,  0,  1, -1, 0, 0,
        mu, mu, mu, mu, 1,-1;
    //cout << "mSubCI:\n" << mSubCI << endl << endl;

    mCI.SetBlock(0, 0, mSubCI);
    mCI.SetBlock(3, 6, mSubCI);
    mCI.SetBlock(6, 12, mSubCI);
    mCI.SetBlock(9, 18, mSubCI);
    //cout << "mCI:\n" << mCI << endl << endl;

    vCi << 0, 0, 0, 0, 0, maxF, 0, 0, 0, 0, 0, maxF, 0, 0, 0, 0, 0, maxF, 0, 0, 0, 0, 0, maxF;
    //cout << "vCi:\n" << vCi << endl << endl;

    rFR << 2, -1, -1;
    rBL << -2, 1, -1;
    rFL << 2, 1, -1;
    rBR << -2, -1, -1;

    mA <<
        1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
        0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
        0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
        0, -rFR(2), rFR(1), 0, -rBL(2), rBL(1), 0, -rFL(2), rFL(1), 0, -rBR(2), rBR(1),
        rFR(2), 0, -rFR(0), rBL(2), 0, -rBL(0), rFL(2), 0, -rFL(0), rBR(2), 0, -rBR(0),
        -rFR(1), rFR(0), 0, -rBL(1), rBL(0), 0, -rFL(1), rFL(0), 0, -rBR(1), rBR(0), 0;
    //cout << "mA:\n" << mA << endl << endl;

    mI6.SetIdentity();
    mI12.SetIdentity();
    mCE.SetBlock(0, 0, mA.Transpose());
    mCE.SetBlock(0, 6, mI12);
    mCE.SetBlock(12, 6, -mI12);
    mCE.SetBlock(24, 0, -mI6);
    //cout << "mCE:\n" << mCI << endl << endl;

    vB << 0, 0, 100, 0, 0, 0;
    //cout << "vB:\n" << vB << endl << endl;

    vCe.SetBlock(0, -vB);       // = -[    vB]
    vCe.SetBlock(6, -vFprev);   //    [vFprev]
    //cout << "vCe:\n" << vCe << endl << endl;

    qp.SetObjectFunc(mG, vG);

    dhTerm::GetCurPos(&curY, &curX);
    dhTerm::print(curY + 1, 0, "rFL = [%6.3f %6.3f %6.3f]", rFL(0), rFL(1), rFL(2));
    dhTerm::print(curY + 2, 0, "rBL = [%6.3f %6.3f %6.3f]", rBL(0), rBL(1), rBL(2));
    dhTerm::print(curY + 3, 0, "rFR = [%6.3f %6.3f %6.3f]", rFR(0), rFR(1), rFR(2));
    dhTerm::print(curY + 4, 0, "rBR = [%6.3f %6.3f %6.3f]", rBR(0), rBR(1), rBR(2));
    dhTerm::print(curY + 5, 0, "  b = [%6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f]", vB(0), vB(1), vB(2), vB(3), vB(4), vB(5));

    while (fRun)
    {
        if (dhTerm::kbhit())
        {
            ch = getchar();
            switch (ch)
            {
            case 'q': rFL(0) += 0.1; break;  case 'w': rFL(0) -= 0.1; break;
            case 'a': rFL(1) += 0.1; break;  case 's': rFL(1) -= 0.1; break;
            case 'z': rFL(2) += 0.1; break;  case 'x': rFL(2) -= 0.1; break;

            case 'e': rFR(0) += 0.1; break;  case 'r': rFR(0) -= 0.1; break;
            case 'd': rFR(1) += 0.1; break;  case 'f': rFR(1) -= 0.1; break;
            case 'c': rFR(2) += 0.1; break;  case 'v': rFR(2) -= 0.1; break;

            case 't': rBL(0) += 0.1; break;  case 'y': rBL(0) -= 0.1; break;
            case 'g': rBL(1) += 0.1; break;  case 'h': rBL(1) -= 0.1; break;
            case 'b': rBL(2) += 0.1; break;  case 'n': rBL(2) -= 0.1; break;

            case 'u': rBR(0) += 0.1; break;  case 'i': rBR(0) -= 0.1; break;
            case 'j': rBR(1) += 0.1; break;  case 'k': rBR(1) -= 0.1; break;
            case 'm': rBR(2) += 0.1; break;  case ',': rBR(2) -= 0.1; break;

            case 'Q': vB(0) += 5; break;  case 'W': vB(0) -= 5; break;
            case 'A': vB(1) += 5; break;  case 'S': vB(1) -= 5; break;
            case 'Z': vB(2) += 5; break;  case 'X': vB(2) -= 5; break;

            case 'E': vB(3) += 5; break;  case 'R': vB(3) -= 5; break;
            case 'D': vB(4) += 5; break;  case 'F': vB(4) -= 5; break;
            case 'C': vB(5) += 5; break;  case 'V': vB(5) -= 5; break;

            case '!': fRun = false; break;
            }
            dhTerm::print(curY + 1, 0, "rFL = [%6.3f %6.3f %6.3f]", rFL(0), rFL(1), rFL(2));
            dhTerm::print(curY + 2, 0, "rBL = [%6.3f %6.3f %6.3f]", rBL(0), rBL(1), rBL(2));
            dhTerm::print(curY + 3, 0, "rFR = [%6.3f %6.3f %6.3f]", rFR(0), rFR(1), rFR(2));
            dhTerm::print(curY + 4, 0, "rBR = [%6.3f %6.3f %6.3f]", rBR(0), rBR(1), rBR(2));
            dhTerm::print(curY + 5, 0, "  b = [%6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f]", vB(0), vB(1), vB(2), vB(3), vB(4), vB(5));
        }
        
        /// Update Matrix
        vX.GetBlock(0, vFprev);
        mCE.SetBlock(0, 3, rFL.GetSkew().Transpose());
        mCE.SetBlock(3, 3, rBL.GetSkew().Transpose());
        mCE.SetBlock(6, 3, rFR.GetSkew().Transpose());
        mCE.SetBlock(9, 3, rBR.GetSkew().Transpose());
        vCe.SetBlock(0, -vB);       // = -[    vB]
        vCe.SetBlock(6, -vFprev);   //    [vFprev]

        /// Solve QP
        time.Start();
        //rtn = qp.Solve(mG, vG, mCE, vCe, mCI, vCi, vX);
        //rtn = qp.Solve(mCE, vCe, mCI, vCi, vX);
        //rtn = qp.Solve(-(mG*vX), mCE, vCe, vX);
        rtn = qp.Solve(mCE, vCe, mCI, vCi, vX);
        time.Stop();

        /// Print Data
        dhTerm::print(curY + 6, 0, "Solver Result: %s", (rtn) ? "Fail   " : "Success");
        dhTerm::print(curY + 7, 0, "Cost Func Result: %f", qp.GetObjectValue());
        dhTerm::print(curY + 8, 0, "Xf_FL = [%6.3f %6.3f %6.3f]", vX(0), vX(1), vX(2));
        dhTerm::print(curY + 9, 0, "Xf_BL = [%6.3f %6.3f %6.3f]", vX(3), vX(4), vX(5));
        dhTerm::print(curY + 10, 0, "Xf_FR = [%6.3f %6.3f %6.3f]", vX(6), vX(7), vX(8));
        dhTerm::print(curY + 11, 0, "Xf_BR = [%6.3f %6.3f %6.3f]", vX(9), vX(10), vX(11));
        dhTerm::print(curY + 12, 0, "Xd_TR = [%6.3f %6.3f %6.3f]", vX(24), vX(25), vX(26));
        dhTerm::print(curY + 13, 0, "Xd_OR = [%6.3f %6.3f %6.3f]", vX(27), vX(28), vX(29));
        dhTerm::print(curY + 14, 0, "elapse time: %6.3f [ms]", time.GetElapsedTime_msec());
    }
    Println;
    Println;
}

void dtMathMoBedQP()
{
    PrintHeading("Quadruped Robot Example using dtMath ver (f-f_prev) of DH");
    int curX, curY;
    bool fRun = true;
    char ch = 0;

    CdhTimeCheck time;
    int8_t rtn = 0;
    double elapseTime = 0;
    double maxTime = 0;
    double minTime = 100;

    double mu = 0.5;
    double maxF = 500;

    // n = 30, m = 24, p = 18
    CdtVector<18, double> vX;       // Vector x = [ fFR, fBL, fFL, fBR, F-Fprev, A*f-b]^T, x:R(30=12+12+6)
    CdtVector<12, double> vX2;

    CdtMatrix<18, 18, double> mG;   // mG = [ Gc, 0 0; 0 Gd 0 ; 0 0 Ge ], mG:R(30x30)
    CdtVector<18, double> vG;       // vG0 is 0

    CdtMatrix<3, 6, double> mSubCI;
    CdtMatrix<18, 8, double> mCI;  // diag = [CIFR, CIBL, CIFL, CIBR, 0]
    CdtVector<8, double> vCi;

    CdtVector3<double> rFR, rBL, rFL, rBR;
    CdtMatrix<6, 12, double> mA;    // mA = [  I3   I3   I3   I3  ]
                                    //      [ rFRx rBLx rFLx rBLx ]
    CdtMatrix<6, 6, double> mI6;
    CdtVector<6, double> vB;        // Vector b = [ fd, td ]^T, b:R(6)
    CdtMatrix<18, 6, double> mCE;  // = [A  0 -I]^T
                                    //   [I -I  0]
    CdtVector<6, double> vCe;      // = -[vB; vFprev]

    CdtQuadProg<18, 8, 6, double> qp;

    double Gc = 0.1;
    double Gd = 10.0;
    double w[18] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 100, 100, 100, 100, 100, 100 };
    mG.SetDiagonal(w, sizeof(w));
    //cout << "mG:\n" << mG << endl;

    vG.SetZero();
    //cout << "vG:\n" << vG << endl << endl;

    mSubCI << 0, 0, 0, 0, 1, -1;
    //cout << "mSubCI:\n" << mSubCI << endl << endl;

    mCI.SetBlock(0, 0, mSubCI);
    mCI.SetBlock(3, 2, mSubCI);
    mCI.SetBlock(6, 4, mSubCI);
    mCI.SetBlock(9, 6, mSubCI);
    //cout << "mCI:\n" << mCI << endl << endl;

    vCi << 10, maxF, 10, maxF, 10, maxF, 10, maxF;
    //cout << "vCi:\n" << vCi << endl << endl;

    rFL << 0.38, 0.32, -0.18;
    rBL << -0.38, 0.32, -0.18;
    rFR << 0.38, -0.32, -0.18;
    rBR << -0.38, -0.32, -0.18;

    mA <<
        1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
        0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
        0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
        0, -rFL(2), rFL(1), 0, -rBL(2), rBL(1), 0, -rFR(2), rFR(1), 0, -rBR(2), rBR(1),
        rFL(2), 0, -rFL(0), rBL(2), 0, -rBL(0), rFR(2), 0, -rFR(0), rBR(2), 0, -rBR(0),
        -rFL(1), rFL(0), 0, -rBL(1), rBL(0), 0, -rFR(1), rFR(0), 0, -rBR(1), rBR(0), 0;
    //cout << "mA:\n" << mA << endl << endl;

    mI6.SetIdentity();
    mCE.SetBlock(0, 0, mA.Transpose());
    mCE.SetBlock(12, 0, -mI6);
    //cout << "mCE:\n" << mCI << endl << endl;

    vB << 0, 0, 100, 0, 0, 0;
    //cout << "vB:\n" << vB << endl << endl;

    vCe.SetBlock(0, -vB);       // = -[    vB]
    //vCe.SetBlock(6, -vFprev);   //    [vFprev]
    //cout << "vCe:\n" << vCe << endl << endl;

    qp.SetObjectFunc(mG, vG);

    dhTerm::GetCurPos(&curY, &curX);
    dhTerm::print(curY + 1, 0, "rFL = [%6.3f %6.3f %6.3f]", rFL(0), rFL(1), rFL(2));
    dhTerm::print(curY + 2, 0, "rBL = [%6.3f %6.3f %6.3f]", rBL(0), rBL(1), rBL(2));
    dhTerm::print(curY + 3, 0, "rFR = [%6.3f %6.3f %6.3f]", rFR(0), rFR(1), rFR(2));
    dhTerm::print(curY + 4, 0, "rBR = [%6.3f %6.3f %6.3f]", rBR(0), rBR(1), rBR(2));
    dhTerm::print(curY + 5, 0, "  b = [%6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f]", vB(0), vB(1), vB(2), vB(3), vB(4), vB(5));

    while (fRun)
    {
        if (dhTerm::kbhit())
        {
            ch = getchar();
            switch (ch)
            {
            case 'q': rFL(0) += 0.01; break;  case 'w': rFL(0) -= 0.01; break;
            case 'a': rFL(1) += 0.01; break;  case 's': rFL(1) -= 0.01; break;
            case 'z': rFL(2) += 0.01; break;  case 'x': rFL(2) -= 0.01; break;

            case 'e': rFR(0) += 0.01; break;  case 'r': rFR(0) -= 0.01; break;
            case 'd': rFR(1) += 0.01; break;  case 'f': rFR(1) -= 0.01; break;
            case 'c': rFR(2) += 0.01; break;  case 'v': rFR(2) -= 0.01; break;

            case 't': rBL(0) += 0.01; break;  case 'y': rBL(0) -= 0.01; break;
            case 'g': rBL(1) += 0.01; break;  case 'h': rBL(1) -= 0.01; break;
            case 'b': rBL(2) += 0.01; break;  case 'n': rBL(2) -= 0.01; break;

            case 'u': rBR(0) += 0.01; break;  case 'i': rBR(0) -= 0.01; break;
            case 'j': rBR(1) += 0.01; break;  case 'k': rBR(1) -= 0.01; break;
            case 'm': rBR(2) += 0.01; break;  case ',': rBR(2) -= 0.01; break;

            case 'Q': vB(0) += 5; break;  case 'W': vB(0) -= 5; break;
            case 'A': vB(1) += 5; break;  case 'S': vB(1) -= 5; break;
            case 'Z': vB(2) += 5; break;  case 'X': vB(2) -= 5; break;

            case 'E': vB(3) += 5; break;  case 'R': vB(3) -= 5; break;
            case 'D': vB(4) += 5; break;  case 'F': vB(4) -= 5; break;
            case 'C': vB(5) += 5; break;  case 'V': vB(5) -= 5; break;

            case '!': fRun = false; break;
            }
            dhTerm::print(curY + 1, 0, "rFL = [%6.3f %6.3f %6.3f]", rFL(0), rFL(1), rFL(2));
            dhTerm::print(curY + 2, 0, "rBL = [%6.3f %6.3f %6.3f]", rBL(0), rBL(1), rBL(2));
            dhTerm::print(curY + 3, 0, "rFR = [%6.3f %6.3f %6.3f]", rFR(0), rFR(1), rFR(2));
            dhTerm::print(curY + 4, 0, "rBR = [%6.3f %6.3f %6.3f]", rBR(0), rBR(1), rBR(2));
            dhTerm::print(curY + 5, 0, "  b = [%6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f]", vB(0), vB(1), vB(2), vB(3), vB(4), vB(5));
        }

        /// Update Matrix
        mCE.SetBlock(0, 3, rFL.GetSkew().Transpose());
        mCE.SetBlock(3, 3, rBL.GetSkew().Transpose());
        mCE.SetBlock(6, 3, rFR.GetSkew().Transpose());
        mCE.SetBlock(9, 3, rBR.GetSkew().Transpose());
        vCe.SetBlock(0, -vB);       // = -[    vB]
        mA.SetBlock(3, 0, rFL.GetSkew());
        mA.SetBlock(3, 3, rBL.GetSkew());
        mA.SetBlock(3, 6, rFR.GetSkew());
        mA.SetBlock(3, 9, rBR.GetSkew());

        /// Solve QP
        time.Start();
        //rtn = qp.Solve(mG, vG, mCE, vCe, mCI, vCi, vX);
        //rtn = qp.Solve(mCE, vCe, mCI, vCi, vX);
        //rtn = qp.Solve(-(mG*vX), mCE, vCe, vX);
        rtn = qp.Solve(mCE, vCe, mCI, vCi, vX);
        vX2 = (mA.Transpose() * (mA * mA.Transpose()).Inv()) * vB;
        time.Stop();

        /// Print Data
        dhTerm::print(curY + 6, 0, "Solver Result: %s", (rtn) ? "Fail   " : "Success");
        dhTerm::print(curY + 7, 0, "Cost Func Result: %f", qp.GetObjectValue());
        dhTerm::print(curY + 8, 0, "Xf_FL = [%6.3f %6.3f %6.3f]", vX(0), vX(1), vX(2));
        dhTerm::print(curY + 9, 0, "Xf_BL = [%6.3f %6.3f %6.3f]", vX(3), vX(4), vX(5));
        dhTerm::print(curY + 10, 0, "Xf_FR = [%6.3f %6.3f %6.3f]", vX(6), vX(7), vX(8));
        dhTerm::print(curY + 11, 0, "Xf_BR = [%6.3f %6.3f %6.3f]", vX(9), vX(10), vX(11));
        dhTerm::print(curY + 12, 0, "Xd_TR = [%6.3f %6.3f %6.3f]", vX(12), vX(13), vX(14));
        dhTerm::print(curY + 13, 0, "Xd_OR = [%6.3f %6.3f %6.3f]", vX(15), vX(16), vX(17));
        dhTerm::print(curY + 14, 0, "elapse time: %6.3f [ms]", time.GetElapsedTime_msec());

        dhTerm::print(curY + 15, 0, "Xf_FL = [%6.3f %6.3f %6.3f]", vX2(0), vX2(1), vX2(2));
        dhTerm::print(curY + 16, 0, "Xf_BL = [%6.3f %6.3f %6.3f]", vX2(3), vX2(4), vX2(5));
        dhTerm::print(curY + 17, 0, "Xf_FR = [%6.3f %6.3f %6.3f]", vX2(6), vX2(7), vX2(8));
        dhTerm::print(curY + 18, 0, "Xf_BR = [%6.3f %6.3f %6.3f]", vX2(9), vX2(10), vX2(11));
    }
    Println;
    Println;
}
#endif /* defined(_WIN32) || defined(__linux__) */