#if defined(_WIN32) || defined(__linux__)
#include <stdio.h>
#include <iostream>
using namespace std;
#endif /* _WIN32 || __linux__ */

#include "testPrint.h"
#include "testContiAlgebraicRiccatiEq.h"
#include "./dtMath/dtMath.h"

void Test_ContiAlgebraicRiccatiEq()
{
    PrintTitle("Test Continuous Algebraic Riccati Equastion");
    dtContiAlgebraicRiccatiEq();
}

void dtContiAlgebraicRiccatiEq()
{
    PrintHeading("dtMath Continuous Algebraic Riccati Eq Solve ");

    CdtMatrix<2, 2> A1, A2, Q1, Q2;
    CdtMatrix<2, 1> B1, B2;
    CdtMatrix<1, 1> R1, R2;

    A1 << -3, 2, 1, 1;
    B1 << 0, 1;
    Q1 << 3, 0, 0, 3;
    //Q1 << 1,-1,-1, 1;
    R1 << 3;

    //A2 << 0, 1, 10, 0;
    //B2 << 0, 1;
    //Q2 << 1, 0, 0, 1;
    //R2 << 1;

    Printf("matrix A =\n");
    A1.Print();
    Println;

    Printf("matrix B =\n");
    B1.Print();
    Println;

    Printf("matrix Q =\n");
    Q1.Print();
    Println;

    Printf("matrix R =\n");
    R1.Print();
    Println;

    CdtCARE<2, 1> care1, care2;

    Printf("solution X =\n");
    care1.Solve(A1, B1, Q1, R1).Print();
    Println;
    //care2.Solve(A2, B2, Q2, R2);
}