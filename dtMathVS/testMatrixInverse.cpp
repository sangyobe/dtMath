#include "testPrint.h"
#include "testMatrixInverse.h"
#include "./dtMath/dtMath.h"

void Test_MatrixInverse()
{
    PrintTitle("Test Inverse of a matrix");
    MatrixInverse();
    Matrix3Inverse();
    RotationInverse();
    TransformInverse();
}

void MatrixInverse()
{
    PrintHeading("Inverse of tMatrix ");
    
    int8_t isOk = 0;
    double M5[5*5] = { 13,11,7,16,5, 3,37,9,7,8, 12,14,31,17,3, 6,19,12,21,15, 7,10,11,13,33 }; // rank 5 (full rank)
    double M4[5*5] = { 13,11,7,16,5, 3,37,9,7,8, 26,22,14,32,10, 6,19,12,21,15, 7,10,11,13,33 }; // rank 4

    double VF[4*3] = { 8,1,3, 8,2,1, 9,4,7, 7,5,2 };    // full rank
    double VC[4*3] = { 8,1,2, 8,2,4, 9,4,8, 7,5,10 };   // col rank 2
    double VR[4*3] = { 8,1,3, 8,2,1,16,2,6, 16,4,2 };   // row rank 2

    double HF[3*4] = { 8,8,9,7, 1,2,4,5, 3,1,7,2 };     // full rank
    double HR[3*4] = { 8,8,9,7, 1,2,4,5, 2,4,8,10 };    // row rank 2
    double HC[3*4] = { 8,8,16,16, 1,2,2,4, 3,1,6,2 };    // col rank 2

    CdtMatrix<5, 5, double> matM5(M5, sizeof(M5));
    CdtMatrix<5, 5, double> matM4(M4, sizeof(M4));

    CdtMatrix<4, 3, double> matVF(VF, sizeof(VF));
    CdtMatrix<4, 3, double> matVC(VC, sizeof(VC));
    CdtMatrix<4, 3, double> matVR(VR, sizeof(VR));

    CdtMatrix<3, 4, double> matHF(HF, sizeof(HF));
    CdtMatrix<3, 4, double> matHC(HC, sizeof(HC));
    CdtMatrix<3, 4, double> matHR(HR, sizeof(HR));

    Printf("/* Inverse of a full-rank square matrix(5x5) */\n");
    Printf("Target Matrix:\n");
    matM5.Print();
    Printf("Inverse:\n");
    matM5.Inv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");
    Printf("invM * M:\n");
    (matM5.Inv() * matM5).Print();
    Println;
    
    Printf("/* Inverse of a rank4 square matrix(5x5) */\n");
    Printf("Target Matrix:\n");
    matM4.Print();
    Printf("Inverse:\n");
    matM4.Inv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");

    Printf("/* Pseudo Inverse of a full-rank square matrix(5x5) */\n");
    Printf("Target Matrix:\n");
    matM5.Print();
    Printf("Pseudo Inverse:\n");
    matM5.PInv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");
    Printf("PInvM * M:\n");
    (matM5.PInv() * matM5).Print();
    Println;
    
    Printf("/* Pseudo Inverse of a rank4 square matrix(5x5) */\n");
    Printf("Target Matrix:\n");
    matM4.Print();
    Printf("Pseudo Inverse:\n");
    matM4.PInv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");
    Printf("PInvM * M:\n");
    (matM4.PInv() * matM4).Print();
    Println;

    Printf("/* Pseudo Inverse of a full-rank rectangle matrix(4x3) */\n");
    Printf("Target Matrix:\n");
    matVF.Print();
    Printf("Pseudo Inverse:\n");
    matVF.PInv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");
    Printf("PInvM * M:\n");
    (matVF.PInv() * matVF).Print();
    Println;
    
    Printf("/* Pseudo Inverse of a row rank 2 rectangle matrix(4x3) */\n");
    Printf("Target Matrix:\n");
    matVR.Print();
    Printf("Pseudo Inverse:\n");
    matVR.PInv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");
    Printf("PInvM * M:\n");
    (matVR.PInv() * matVR).Print();
    Println;

    Printf("/* Pseudo Inverse of a col rank 2 rectangle matrix(4x3) */\n");
    Printf("Target Matrix:\n");
    matVC.Print();
    Printf("Pseudo Inverse:\n");
    matVC.PInv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");
    Printf("PInvM * M:\n");
    (matVC.PInv() * matVC).Print();
    Println;

    Printf("/* Pseudo Inverse of a full-rank rectangle matrix(3x4) */\n");
    Printf("Target Matrix:\n");
    matHF.Print();
    Printf("Pseudo Inverse:\n");
    matHF.PInv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");
    Printf("PInvM * M:\n");
    (matHF.PInv() * matHF).Print();
    Println;

    Printf("/* Pseudo Inverse of a row rank 2 rectangle matrix(3x4) */\n");
    Printf("Target Matrix:\n");
    matHR.Print();
    Printf("Pseudo Inverse:\n");
    matHR.PInv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");
    Printf("PInvM * M:\n");
    (matHR.PInv() * matHR).Print();
    Println;

    Printf("/* Pseudo Inverse of a col rank 2 rectangle matrix(3x4) */\n");
    Printf("Target Matrix:\n");
    matHC.Print();
    Printf("Pseudo Inverse:\n");
    matHC.PInv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");
    Printf("PInvM * M:\n");
    (matHC.PInv() * matHC).Print();
    Println;
}

void Matrix3Inverse()
{
    PrintHeading("Inverse of tMatrix3 ");

    int8_t isOk = 0;
    double A3[3*3] = { 13,11,7, 3,37,9, 12,14,31 }; // rank 3
    double A2[3*3] = { 13,11,7, 3,37,9, 26,22,14 }; // rank 2

    CdtMatrix3<double> matA3(A3, sizeof(A3));
    CdtMatrix3<double> matA2(A2, sizeof(A2));

    Printf("/* Inverse of a full-rank square matrix(3x3) */\n");
    Printf("Target Matrix:\n");
    matA3.Print();
    Printf("Inverse:\n");
    matA3.Inv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");

    Printf("/* Inverse of a row rank2 square matrix(3x3) */\n");
    Printf("Target Matrix:\n");
    matA2.Print();
    Printf("Inverse:\n");
    matA2.Inv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");

    Printf("/* Pseudo Inverse of a full-rank square matrix(3x3) */\n");
    Printf("Target Matrix\n");
    matA3.Print();
    Printf("Pseudo Inverse:\n");
    matA3.PInv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");

    Printf("/* Pseudo Inverse of a row rank2 square matrix(3x3) */\n");
    Printf("Target Matrix:\n");
    matA2.Print();
    Printf("Pseudo Inverse:\n");
    matA2.PInv(&isOk).Print();
    Printf("isOk:%s\n\n", isOk ? "Yes" : "No");
}

void RotationInverse()
{
    PrintHeading("Inverse of a rotation matrix ");
    Printf("/* rot - euler zyx(10, 20, 30) */\n");
    CdtRotation<double> rot(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 10 * DEG2RADd, 20 * DEG2RADd, 30 * DEG2RADd);
    rot.Print();
    Printf("rot.Inv()\n");
    rot.Inv().Print();
    Println;
}

void TransformInverse()
{
    PrintHeading("Inverse of Transformation matrix ");
    Printf("/* Transform matrix(rot - euler zyx(10, 20, 30), pos - 1, 2, 3) */\n");
    CdtRotation<double> rot(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 10 * DEG2RADd, 20 * DEG2RADd, 30 * DEG2RADd);
    CdtVector3<double> p(1, 2, 3);
    CdtTransform<double> tf(rot, p);
    tf.Print();
    Printf("tf.Inv()\n");
    tf.Inv().Print();
    Println;
}