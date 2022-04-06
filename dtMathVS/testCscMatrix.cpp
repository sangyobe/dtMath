#include "testPrint.h"
#include "testCscMatrix.h"
#include "./dtMath/dtMath.h"

void Test_CscMatrix()
{
    PrintTitle("Test CSC Matrix");

    PrintHeading("Matrix Member Functions ");
    Printf("/* Class Create: with CdtMatrix */\n");
    CdtMatrix<6, 5> mat65;    
    mat65 << 
        0,0,0,0,2,
        0,1,0,0,0,
        0,0,3,0,0,
        0,0,0,0,0,
        0,0,0,0,4,
        0,0,5,0,0;
    CdtCscMatrix<6, 5> cscMat65(mat65);

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

    CdtMatrix<4, 3> mat43;
    mat43 <<
        1, 0, 0,
        0, 2, 3,
        0, 0, 0,
        0, 4, 0;
    CdtCscMatrix<4, 3> cscMat43(mat43);
    cscMat43.Transpose().GetDenseMat().Print('\n');


}