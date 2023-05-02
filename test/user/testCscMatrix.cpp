#include "testCscMatrix.h"
#include "testPrint.h"

#include <dtMath/dtMath.h>

using dtMath::dtCscMatrix;
using dtMath::dtMatrix;

void Test_CscMatrix()
{
    PrintTitle("Test CSC Matrix");

    PrintHeading("Matrix Member Functions ");
    Printf("/* Class Create: with dtMatrix */\n");
    dtMatrix<6, 5> mat65;
    mat65 << 0, 0, 0, 0, 2,
        0, 1, 0, 0, 0,
        0, 0, 3, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 4,
        0, 0, 5, 0, 0;
    dtCscMatrix<6, 5> cscMat65(mat65);

    cscMat65.Print();
    Println;

    Printf("/* Function: GetVec() */\n");
    Printf("cscMat65.GetRowVec(1) = \n");
    cscMat65.GetRowVec(1).Print('\n');
    cscMat65.GetRowVec(5).Print('\n');
    Println;

    Printf("cscMat65.GetColVec(2) = \n");
    cscMat65.GetColVec(2).Print();
    Println;

    cscMat65.GetDenseMat().Print('\n');

    dtMatrix<4, 3> mat43;
    mat43 << 1, 0, 0,
        0, 2, 3,
        0, 0, 0,
        0, 4, 0;
    dtCscMatrix<4, 3> cscMat43(mat43);
    cscMat43.Transpose().GetDenseMat().Print('\n');
}