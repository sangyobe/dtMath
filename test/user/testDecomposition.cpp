#include "testPrint.h"
#include "testDecomposition.h"

#include <dtMath/dtMath.h>

#include <utility/dhTimeCheck.h>

void Test_Decomposition()
{
    PrintTitle("Test Decomposition");
    LowerTriangular();
    UnitLowerTriangular();
    UpperTriangular();
    UnitUpperTriangular();
    NoPivotLU();
    NoPivotLUMat3();
    PartialPivotLU();
    CholeskyLLT();
    CholeskyLDLT();
    HouseholderQR();
    SingularValue();
    SingularValueMat3();
}

void LowerTriangular()
{
    PrintHeading("Lower Triangular Matrix ");
    double b[4] = { 1, 2, 7, 3 };
    CdtVector<4, double> vecB(b, sizeof(b));
    CdtVector<4, double> vecX;

    /* Lower Triangular Matrix Test */
    double L[16] = { 9, 0, 0, 0, 1, 8, 0, 0, 2, 1, 7, 0, 1, 1, 4, 5 };
    CdtMatrix<4, 4, double> matL(L, sizeof(L));
    CdtMatrix<4, 4, double> invMatL;

    Printf("Target Matrix A:\n");
    matL.Print();

    Printf("\nSolve1\n");
    CdtLowerTriangular<4, 4, double>::Solve(matL, vecB, vecX);
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    vecX.Print();

    Printf("\nSolve2\n");
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    CdtLowerTriangular<4, 4, double>::Solve(matL, vecB).Print();

    Printf("\nInverse1:\n");
    CdtLowerTriangular<4, 4, double>::Inverse(matL, invMatL);
    invMatL.Print();

    Printf("\nInverse2:\n");
    invMatL = CdtLowerTriangular<4, 4, double>::Inverse(matL);
    invMatL.Print();

    Println;
}

void UnitLowerTriangular()
{
    PrintHeading("Unit Lower Triangular Matrix ");
    double b[4] = { 1, 2, 7, 3 };
    CdtVector<4, double> vecB(b, sizeof(b));
    CdtVector<4, double> vecX;

    /* Unit Lower Triangular Matrix Test */
    double uintL[16] = { 1, 0, 0, 0, 3, 1, 0, 0, 2, 7, 1, 0, 9, 5, 4, 1 };
    CdtMatrix<4, 4, double> matUintL(uintL, sizeof(uintL));
    CdtMatrix<4, 4, double> invMatUnitL;

    Printf("Target Matrix A:\n");
    matUintL.Print();

    Printf("\nSolve1\n");
    CdtLowerTriangular<4, 4, double>::Solve(matUintL, vecB, vecX);
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    vecX.Print();

    Printf("\nSolve2\n");
    vecX = CdtLowerTriangular<4, 4, double>::Solve(matUintL, vecB);
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    vecX.Print();

    Printf("\nInverse1:\n");
    CdtLowerTriangular<4, 4, double>::InverseUnit(matUintL, invMatUnitL);
    invMatUnitL.Print();

    Printf("\nInverse2:\n");
    invMatUnitL = CdtLowerTriangular<4, 4, double>::InverseUnit(matUintL);
    invMatUnitL.Print();

    Println;
}

void UpperTriangular()
{
    PrintHeading("Upper Triangular Matrix ");
    double b[4] = { 1, 2, 7, 3 };
    CdtVector<4, double> vecB(b, sizeof(b));
    CdtVector<4, double> vecX;

    /* Upper Triangular Matrix Test */
    double U[16] = { 9, 1, 4, 5, 0, 8, 1, 3, 0, 0, 7, 2, 0, 0, 0, 5 };
    CdtMatrix<4, 4, double> matU(U, sizeof(U));
    CdtMatrix<4, 4, double> invMatU;
    vecX.SetFill(0);

    Printf("Target Matrix A:\n");
    matU.Print();

    Printf("\nSolve1\n");
    CdtUpperTriangular<4, 4, double>::Solve(matU, vecB, vecX);
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    vecX.Print();

    Printf("\nSolve2\n");
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    CdtUpperTriangular<4, 4, double>::Solve(matU, vecB).Print();

    Printf("\nInverse1:\n");
    CdtUpperTriangular<4, 4, double>::Inverse(matU, invMatU);
    invMatU.Print();

    Printf("\nInverse2:\n");
    invMatU = CdtUpperTriangular<4, 4, double>::Inverse(matU);
    invMatU.Print();

    Println;
}

void UnitUpperTriangular()
{
    PrintHeading("Unit Upper Triangular Matrix ");

    double b[4] = { 1, 2, 7, 3 };
    CdtVector<4, double> vecB(b, sizeof(b));
    CdtVector<4, double> vecX;

    /* Unit Upper Triangular Matrix Test */
    double unitU[16] = { 1, 3, 4, 5, 0, 1, 7, 3, 0, 0, 1, 9, 0, 0, 0, 1 };
    CdtMatrix<4, 4, double> matUnitU(unitU, sizeof(unitU));
    CdtMatrix<4, 4, double> invMatUnitU;

    Printf("Target Matrix A:\n");
    matUnitU.Print();

    Printf("\nSolve1\n");
    CdtUpperTriangular<4, 4, double>::Solve(matUnitU, vecB, vecX);
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    vecX.Print();

    Printf("\nSolve2\n");
    vecX = CdtUpperTriangular<4, 4, double>::Solve(matUnitU, vecB);
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    vecX.Print();

    Printf("\nInverse1:\n");
    CdtUpperTriangular<4, 4, double>::InverseUnit(matUnitU, invMatUnitU);
    invMatUnitU.Print();

    Printf("\nInverse2:\n");
    invMatUnitU = CdtUpperTriangular<4, 4, double>::InverseUnit(matUnitU);
    invMatUnitU.Print();

    Println;
}

void NoPivotLU()
{
    PrintHeading("No Pivoting LU Decomposition ");
    /* No Pivoting LU Decomposition Test */
    double A[16] = { 9,1,2,4,1,8,2,3,2,1,7,3,1,3,5,19 };
    double b[4] = { 1, 2, 7, 3 };
    
    CdtMatrix<4, 4, double> matA(A, sizeof(A));
    CdtMatrix<4, 4, double> matLU;
    CdtMatrix<4, 4, double> invLU;
    CdtVector<4, double> vecB(b, sizeof(b));
    CdtVector<4, double> vecX;

    Printf("Target Matrix A:\n");
    matA.Print();

    Printf("\nLU:\n");
    matA.NoPivLU().GetMatrix().Print();

    Printf("\nL:\n");
    matA.NoPivLU().GetMatrixL().Print();

    Printf("\nU:\n");
    matA.NoPivLU().GetMatrixU().Print();

    Printf("\nSolve1\n");
    matA.NoPivLU().Solve(vecB, vecX);
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    vecX.Print();

    Printf("\nSolve2\n");
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    matA.NoPivLU().Solve(vecB).Print();

    Printf("\nInverse1:\n");
    matA.NoPivLU().Inverse(invLU);
    invLU.Print();

    Printf("\nInverse2:\n");
    matA.NoPivLU().Inverse().Print();

    Println;
}

void NoPivotLUMat3()
{
    PrintHeading("No Pivoting LU Decomposition ");
    /* No Pivoting LU Decomposition Test */
    double A[9] = { 9,1,2, 1,8,2, 2,1,7 };
    double b[4] = { 1, 2, 7 };

    CdtMatrix3<double> matA(A, sizeof(A));
    CdtMatrix3<double> matLU;
    CdtMatrix3<double> invLU;
    CdtVector<3, double> vecB(b, sizeof(b));
    CdtVector<3, double> vecX;

    Printf("Target Matrix A:\n");
    matA.Print();

    Printf("\nLU:\n");
    matA.NoPivLU().GetMatrix().Print();

    Printf("\nL:\n");
    matA.NoPivLU().GetMatrixL().Print();

    Printf("\nU:\n");
    matA.NoPivLU().GetMatrixU().Print();

    Printf("\nSolve1\n");
    matA.NoPivLU().Solve(vecB, vecX);
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    vecX.Print();

    Printf("\nSolve2\n");
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    matA.NoPivLU().Solve(vecB).Print();

    Printf("\nInverse1:\n");
    matA.NoPivLU().Inverse(invLU);
    invLU.Print();

    Printf("\nInverse2:\n");
    matA.NoPivLU().Inverse().Print();

    matA.Inv().Print();
    Println;
}

void PartialPivotLU()
{
    PrintHeading("Partial Pivoting LU Decomposition ");
    /* Partial Pivoting LU Decomposition Test */
    double A[16] = {
        19, 1,21, 4,
         1, 9, 2, 3,
        22, 1, 7, 3,
         1,31, 5, 3 };
    double b[4] = { 1, 2, 7, 3 };
    
    CdtMatrix<4, 4, double> matA(A, sizeof(A));
    CdtMatrix<4, 4, double> invA;
    CdtVector<4, double> vecB(b, sizeof(b));
    CdtVector<4, double> vecX;

    Printf("Target Matrix A:\n");
    matA.Print();

    Printf("\nLU:\n");
    matA.PartialPivLU().GetMatrix().Print();

    Printf("\nL:\n");
    matA.PartialPivLU().GetMatrixL().Print();

    Printf("\nU:\n");
    matA.PartialPivLU().GetMatrixU().Print();

    Printf("\nP:\n");
    matA.PartialPivLU().GetMatrixP().Print();

    Printf("\nSolve1\n");
    matA.PartialPivLU().Solve(vecB, vecX);
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    vecX.Print();

    Printf("\nSolve2:\n");
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    matA.PartialPivLU().Solve(vecB).Print();

    Printf("\nInverse1:\n");
    matA.PartialPivLU().Inverse(invA);
    invA.Print();

    Printf("\nInverse2:\n");
    matA.PartialPivLU().Inverse().Print();

    Println;
}

void CholeskyLLT()
{
    PrintHeading("Cholesky Decomposition LLT ");
    /* Cholesky Decomposition A */
    double A[16] = {
        81, 9,18, 9,
         9,65,10, 9,
        18,10,54,31,
         9, 9,31,43 };
    double b[4] = { 1, 2, 7, 3 };

    CdtMatrix<4, 4, double> matA(A, sizeof(A));
    CdtMatrix<4, 4, double> invA;
    CdtVector<4, double> vecB(b, sizeof(b));
    CdtVector<4, double> vecX;

    Printf("Target Matrix A:\n");
    matA.Print();

    Printf("\nCholesky LLT:\n");
    matA.LLT().GetMatrix().Print();

    Printf("\nCholesky L:\n");
    matA.LLT().GetMatrixL().Print();

    Printf("\nCholesky LT:\n");
    matA.LLT().GetMatrixU().Print();

    Printf("\nSolve1\n");
    matA.LLT().Solve(vecB, vecX);
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    vecX.Print();

    Printf("\nSolve2\n");
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    matA.LLT().Solve(vecB).Print();

    Printf("\nInverse1:\n");
    matA.LLT().Inverse(invA);
    invA.Print();

    Printf("\nInverse2:\n");
    matA.LLT().Inverse().Print();

    Println;
}

void CholeskyLDLT()
{
    PrintHeading("Cholesky Decomposition LDLT ");

    double A[16] = {
        81, 9,18, 9,
         9,65,10, 9,
        18,10,54,31,
         9, 9,31,43 };
    double b[4] = { 1, 2, 7, 3 };

    CdtMatrix<4, 4, double> matA(A, sizeof(A));
    CdtMatrix<4, 4, double> invA;
    CdtVector<4, double> vecB(b, sizeof(b));
    CdtVector<4, double> vecX;

    Printf("Target Matrix A:\n");
    matA.Print();

    Printf("\nCholesky LDLT:\n");
    matA.LDLT().GetMatrix().Print();

    Printf("\nCholesky L:\n");
    matA.LDLT().GetMatrixL().Print();

    Printf("\nCholesky D:\n");
    matA.LDLT().GetMatrixD().Print();

    Printf("\nCholesky U(LT):\n");
    matA.LDLT().GetMatrixU().Print();

    Printf("\nSolve1\n");
    matA.LDLT().Solve(vecB, vecX);
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    vecX.Print();

    Printf("\nSolve2\n");
    Printf("b:\n");
    vecB.Print();
    Printf("x:\n");
    matA.LDLT().Solve(vecB).Print();

    Printf("\nInverse1:\n");
    matA.LDLT().Inverse(invA);
    invA.Print();

    Printf("\nInverse2:\n");
    matA.LDLT().Inverse().Print();

    Println;
}

void HouseholderQR()
{
    PrintHeading("Householder QR Decomposition of Matrix ");
    double A43[4 * 3] = { -1, -1, 1, 4, 3, 2, -1, -1, 5, 1, 3, 7 };
    double A34[3 * 4] = { -1, -1, 1, 4, 3, 2, -1, -1, 5, 1, 3, 7 };

    CdtMatrix<4, 3, double> matA43(A43, sizeof(A43));
    CdtMatrix<3, 4, double> matA34(A34, sizeof(A34));
    CdtQR<4, 3, double> QR43(matA43);
    CdtQR<3, 4, double> QR34(matA34);

    Printf("/* Skinny Matrix A (4 x 3) - Rank 3 */\n");
    Printf("Target Matrix:\n");
    matA43.Print();

    Printf("\nQR Q:\n");
    //matA43.QR().GetMatrixQ().Print();
    QR43.GetMatrixQ().Print();
    Printf("\nQR R:\n");
    //matA43.QR().GetMatrixR().Print();
    QR43.GetMatrixR().Print();
    Printf("\nQ*R:\n");
    (QR43.GetMatrixQ() * QR43.GetMatrixR()).Print();
    Printf("\n\n");

    Printf("/* Fat Matrix A (3 x 4) - Rank 3 */\n");
    Printf("Target Matrix:\n");
    matA34.Print();

    Printf("\nQR Q:\n");
    //matA34.QR().GetMatrixQ().Print();
    QR34.GetMatrixQ().Print();
    Printf("\nQR R:\n");
    //matA34.QR().GetMatrixR().Print();
    QR34.GetMatrixR().Print();
    Printf("\nQ*R:\n");
    (QR34.GetMatrixQ() * QR34.GetMatrixR()).Print();

    Println;
}

void SingularValue()
{
    PrintHeading("Singular Value Decomposition of Matrix ");

    double A43[4 * 3] = { 8,1,3,8,2,1,9,4,7,7,5,2 };  // rank 3
    double B43[4 * 3] = { 8,1,2,8,2,4,9,4,8,7,5,10 }; // rank 2
    double x4[4] = { 1, 2, 7, 3 };
    double x3[3] = { 1, 2, 7 };

    CdtMatrix<4, 3, double> matA43(A43, sizeof(A43));
    CdtMatrix<4, 3, double> matB43(B43, sizeof(B43));
    CdtMatrix<3, 4, double> invA43;
    CdtMatrix<3, 4, double> invB43;
    CdtVector<4, double> vecB4;
    CdtVector<3, double> vecX3(x3, sizeof(x3));

    CdtMatrix<3, 4, double> matA34(matA43.Transpose());
    CdtMatrix<3, 4, double> matB34(matB43.Transpose());
    CdtMatrix<4, 3, double> invA34;
    CdtMatrix<4, 3, double> invB34;
    CdtVector<3, double> vecB3;
    CdtVector<4, double> vecX4(x4, sizeof(x4));
    
    Printf("/* Skinny Matrix A (4 x 3) - Rank 3 */\n");
    Printf("Target Matrix:\n");
    matA43.Print();

    Printf("\nSVD U:\n");
    matA43.SVD().GetMatrixU().Print();
    Printf("\nSVD S:\n");
    matA43.SVD().GetMatrixS().Print();
    Printf("\nSVD V:\n");
    matA43.SVD().GetMatrixV().Print();
    Printf("\nU*S*VT:\n");
    (matA43.SVD().GetMatrixU() * matA43.SVD().GetMatrixS() * matA43.SVD().GetMatrixV().Transpose()).Print();

    Printf("\nVector b:\n");
    vecB4 = matA43 * vecX3;
    vecB4.Print();

    Printf("\nSolution Vector X:\n");
    vecX3.Print();

    Printf("\nSolve1 Vector X:\n");
    matA43.SVD().Solve(vecB4, vecX3);
    vecX3.Print();

    Printf("\nSolve2 Vector X:\n");
    matA43.SVD().Solve(vecB4).Print();

    Printf("\nInverse1:\n");
    matA43.SVD().Inverse(invA43);
    invA43.Print();

    Printf("\nInverse2:\n");
    matA43.SVD().Inverse().Print();

    Printf("\n\n");

    Printf("/* Skinny Matrix A (4 x 3) - Rank 2 */\n");
    Printf("Target Matrix:\n");
    matB43.Print();

    Printf("\nSVD U:\n");
    matB43.SVD().GetMatrixU().Print();
    Printf("\nSVD S:\n");
    matB43.SVD().GetMatrixS().Print();
    Printf("\nSVD V:\n");
    matB43.SVD().GetMatrixV().Print();
    Printf("\nU*S*VT:\n");
    (matB43.SVD().GetMatrixU() * matB43.SVD().GetMatrixS() * matB43.SVD().GetMatrixV().Transpose()).Print();

    Printf("\nVector b:\n");
    vecB4 = matB43 * vecX3;
    vecB4.Print();

    Printf("\nSolution Vector X:\n");
    vecX3.Print();

    Printf("\nSolve1 Vector X:\n");
    matB43.SVD().Solve(vecB4, vecX3);
    vecX3.Print();

    Printf("\nSolve2 Vector X:\n");
    matB43.SVD().Solve(vecB4).Print();

    Printf("\nInverse1:\n");
    matB43.SVD().Inverse(invA43);
    invA43.Print();

    Printf("\nInverse2:\n");
    matB43.SVD().Inverse().Print();

    Printf("\n\n");

    Printf("/* Fat Matrix A (3 x 4) - Rank 3 */\n");
    Printf("Target Matrix:\n");
    matA34.Print();

    Printf("\nSVD U:\n");
    matA34.SVD().GetMatrixU().Print();
    Printf("\nSVD S:\n");
    matA34.SVD().GetMatrixS().Print();
    Printf("\nSVD V:\n");
    matA34.SVD().GetMatrixV().Print();
    Printf("\nU*S*VT:\n");
    (matA34.SVD().GetMatrixU() * matA34.SVD().GetMatrixS() * matA34.SVD().GetMatrixV().Transpose()).Print();

    vecX4.SetElement(x4, sizeof(x4));
    Printf("\nVector b:\n");
    vecB3 = matA34 * vecX4;
    vecB3.Print();

    Printf("\nSolution Vector X:\n");
    vecX4.Print();

    Printf("\nSolve1 Vector X:\n");
    matA34.SVD().Solve(vecB3, vecX4);
    vecX4.Print();

    Printf("\nSolve2 Vector X:\n");
    matA34.SVD().Solve(vecB3).Print();

    Printf("\nInverse1:\n");
    matA34.SVD().Inverse(invA34);
    invA34.Print();

    Printf("\nInverse2:\n");
    matA34.SVD().Inverse().Print();

    Printf("\n\n");

    Printf("/* Fat Matrix A (3 x 4) - Rank 2 */\n");
    Printf("Target Matrix:\n");
    matB34.Print();

    Printf("\nSVD U:\n");
    matB34.SVD().GetMatrixU().Print();
    Printf("\nSVD S:\n");
    matB34.SVD().GetMatrixS().Print();
    Printf("\nSVD V:\n");
    matB34.SVD().GetMatrixV().Print();
    Printf("\nU*S*VT:\n");
    (matB34.SVD().GetMatrixU() * matB34.SVD().GetMatrixS() * matB34.SVD().GetMatrixV().Transpose()).Print();

    vecX4.SetElement(x4, sizeof(x4));
    Printf("\nVector b:\n");
    vecB3 = matB34 * vecX4;
    vecB3.Print();

    Printf("\nSolution Vector X:\n");
    vecX4.Print();

    Printf("\nSolve1 Vector X:\n");
    matB34.SVD().Solve(vecB3, vecX4);
    vecX4.Print();

    Printf("\nSolve2 Vector X:\n");
    matB34.SVD().Solve(vecB3).Print();

    Printf("\nInverse1:\n");
    matB34.SVD().Inverse(invA34);
    invA34.Print();

    Printf("\nInverse2:\n");
    matB34.SVD().Inverse().Print();

    Println;
}

void SingularValueMat3()
{
    PrintHeading("Singular Value Decomposition of Matrix3 ");

    double A[3 * 3] = { 1, 12, 7, 3, 8, 9, 12, 3, 19 }; // rank 3
    double B[3 * 3] = { 1, 12, 7, 3, 8, 9, 6, 16, 18 }; // rank 2
    double C[3 * 3] = { 1, 12, 2, 3, 8, 6, 12, 3, 24 }; // rank 2
    double x[3] = { 1, 5, 2 };

    CdtMatrix3<double> matA(A, sizeof(A));
    CdtMatrix3<double> matB(B, sizeof(B));
    CdtMatrix3<double> matC(C, sizeof(C));
    CdtMatrix3<double> invA;
    CdtMatrix3<double> invB;
    CdtMatrix3<double> invC;
    CdtVector<3, double> vecB;
    CdtVector<3, double> vecX(x, sizeof(x));

    Printf("/* Matrix(3 x 3) - Rank 3*/\n");
    Printf("Target Matrix A:\n");
    matA.Print();

    Printf("\nSVD U:\n");
    matA.SVD().GetMatrixU().Print();
    Printf("\nSVD S:\n");
    matA.SVD().GetMatrixS().Print();
    Printf("\nSVD V:\n");
    matA.SVD().GetMatrixV().Print();
    Printf("\nU*S*VT:\n");
    (matA.SVD().GetMatrixU() * matA.SVD().GetMatrixS() * matA.SVD().GetMatrixV().Transpose()).Print();

    Printf("\nVector b:\n");
    vecB = matA * vecX;
    vecB.Print();

    Printf("\nSolution Vector X:\n");
    vecX.Print();

    Printf("\nSolve1 Vector X:\n");
    matA.SVD().Solve(vecB, vecX);
    vecX.Print();

    Printf("\nSolve2 Vector X:\n");
    matA.SVD().Solve(vecB).Print();

    Printf("\nInverse1:\n");
    matA.SVD().Inverse(invA);
    invA.Print();

    Printf("\nInverse2:\n");
    matA.SVD().Inverse().Print();

    Printf("\n\n");

    Printf("/* Matrix(3 x 3) - Row Rank 2*/\n");
    Printf("Target Matrix A:\n");
    matB.Print();

    Printf("\nSVD U:\n");
    matB.SVD().GetMatrixU().Print();
    Printf("\nSVD S:\n");
    matB.SVD().GetMatrixS().Print();
    Printf("\nSVD V:\n");
    matB.SVD().GetMatrixV().Print();
    Printf("\nU*S*VT:\n");
    (matB.SVD().GetMatrixU() * matB.SVD().GetMatrixS() * matB.SVD().GetMatrixV().Transpose()).Print();

    vecX.SetElement(x, sizeof(x));
    Printf("\nVector b:\n");
    vecB = matB * vecX;
    vecB.Print();

    Printf("\nSolution Vector X:\n");
    vecX.Print();

    Printf("\nSolve1 Vector X:\n");
    matB.SVD().Solve(vecB, vecX);
    vecX.Print();

    Printf("\nSolve2 Vector X:\n");
    matB.SVD().Solve(vecB).Print();

    Printf("\nInverse1:\n");
    matB.SVD().Inverse(invA);
    invA.Print();

    Printf("\nInverse2:\n");
    matB.SVD().Inverse().Print();

    Printf("\n\n");

    Printf("/* Matrix(3 x 3) - Col Rank 2*/\n");
    Printf("Target Matrix A:\n");
    matC.Print();

    Printf("\nSVD U:\n");
    matC.SVD().GetMatrixU().Print();
    Printf("\nSVD S:\n");
    matC.SVD().GetMatrixS().Print();
    Printf("\nSVD V:\n");
    matC.SVD().GetMatrixV().Print();
    Printf("\nU*S*VT:\n");
    (matC.SVD().GetMatrixU() * matC.SVD().GetMatrixS() * matC.SVD().GetMatrixV().Transpose()).Print();

    vecX.SetElement(x, sizeof(x));
    Printf("\nVector b:\n");
    vecB = matC * vecX;
    vecB.Print();

    Printf("\nSolution Vector X:\n");
    vecX.Print();

    Printf("\nSolve1 Vector X:\n");
    matC.SVD().Solve(vecB, vecX);
    vecX.Print();

    Printf("\nSolve2 Vector X:\n");
    matC.SVD().Solve(vecB).Print();

    Printf("\nInverse1:\n");
    matC.SVD().Inverse(invA);
    invA.Print();

    Printf("\nInverse2:\n");
    matC.SVD().Inverse().Print();

    Println;
}
