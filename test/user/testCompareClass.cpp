#include "testCompareClass.h"
#include "testPrint.h"

#include <dtMath/dtMath.h>

#include <utility/dhTimeCheck.h>

using dtMath::dtMat3;
using dtMath::dtMatrix;
using dtMath::dtVector;

void MatAdd();
void MatScalar();
void MatMat();
void MatVec();
void MatTranspose();
void MatInv();
void MatPInv();

void Test_CompareClass()
{
    PrintTitle("Test Compare Class");

    PrintHeading("Matrix + Matrix ");
    MatAdd();

    PrintHeading("Matrix * Scalar ");
    MatScalar();

    PrintHeading("Matrix x Matrix ");
    MatMat();

    PrintHeading("Matrix x Vector ");
    MatVec();

    PrintHeading("Matrix Transpose ");
    MatTranspose();

    PrintHeading("Matrix Inverse ");
    MatInv();

    PrintHeading("Matrix Pseudo Inverse ");
    MatPInv();
}

void MatAdd()
{
    float a[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    float b[9] = {9, 8, 7, 6, 5, 4, 3, 2, 1};

    dtMatrix<3, 3, float> matA(a, sizeof(a));
    dtMatrix<3, 3, float> matB(b, sizeof(b));
    dtMatrix<3, 3, float> matC;

    dtMat3 mat3A(a, sizeof(a));
    dtMat3 mat3B(b, sizeof(b));
    dtMat3 mat3C;

    CdhTimeCheck elapseTime;

    elapseTime.Start();
    for (int i = 0; i < 1000000; i++)
        matC = matA + matB;
    elapseTime.Stop();
    Printf("// dtMath dtMatrix Mat + Mat\n");
    Printf("matA =\n");
    matA.Print();
    Printf("matB =\n");
    matB.Print();
    Printf("Result matC is\n");
    matC.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 1000000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 1000000);
#endif

    elapseTime.Start();
    for (int i = 0; i < 1000000; i++)
        mat3C = mat3A + mat3B;
    elapseTime.Stop();
    Printf("// dtMath dtMatrix3 Mat + Mat\n");
    Printf("mat3A =\n");
    mat3A.Print();
    Printf("mat3B =\n");
    mat3B.Print();
    Printf("Result mat3C is\n");
    mat3C.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 1000000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 1000000);
#endif
}

void MatScalar()
{
    float a[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    dtMatrix<3, 3, float> matA(a, sizeof(a));
    dtMatrix<3, 3, float> matC;

    dtMat3 mat3A(a, sizeof(a));
    dtMat3 mat3C;

    CdhTimeCheck elapseTime;

    elapseTime.Start();
    for (int i = 0; i < 1000000; i++)
        matC = matA * 100.0f;
    elapseTime.Stop();
    Printf("// dtMath dtMatrix Mat x Scalar\n");
    Printf("matA =\n");
    matA.Print();
    Printf("Result matA * 100 is\n");
    matC.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 1000000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 1000000);
#endif

    elapseTime.Start();
    for (int i = 0; i < 1000000; i++)
        mat3C = mat3A * 100.0f;
    elapseTime.Stop();
    Printf("// dtMath dtMatrix3 Mat x Scalar\n");
    Printf("mat3A =\n");
    mat3A.Print();
    Printf("Result mat3A * 100 is\n");
    mat3C.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 1000000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 1000000);
#endif
}

void MatMat()
{
    float a[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    float b[9] = {9, 8, 7, 6, 5, 4, 3, 2, 1};

    dtMatrix<3, 3, float> matA(a, sizeof(a));
    dtMatrix<3, 3, float> matB(b, sizeof(b));
    dtMatrix<3, 3, float> matC;

    dtMat3 mat3A(a, sizeof(a));
    dtMat3 mat3B(b, sizeof(b));
    dtMat3 mat3C;

    CdhTimeCheck elapseTime;

    elapseTime.Start();
    for (int i = 0; i < 100000; i++)
        matC = matA * matB;
    elapseTime.Stop();
    Printf("// dtMath dtMatrix Mat x Mat\n");
    Printf("matA =\n");
    matA.Print();
    Printf("matB =\n");
    matB.Print();
    Printf("Result matA * matB is \n");
    matC.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#endif

    elapseTime.Start();
    for (int i = 0; i < 100000; i++)
        mat3C = mat3A * mat3B;
    elapseTime.Stop();
    Printf("// dtMath dtMatrix3 Mat x Mat\n");
    Printf("mat3A =\n");
    mat3A.Print();
    Printf("mat3B =\n");
    mat3B.Print();
    Printf("Result mat3A * mat3B is \n");
    mat3C.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#endif
}

void MatVec()
{
    float a[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    dtMatrix<3, 3, float> matA(a, sizeof(a));
    dtMat3 mat3A(a, sizeof(a));
    dtVector<3, float> vecA;
    dtVector<3, float> vecB;

    vecA << 1, 2, 3;

    CdhTimeCheck elapseTime;

    elapseTime.Start();
    for (int i = 0; i < 100000; i++)
        vecB = matA * vecA;
    elapseTime.Stop();
    Printf("// dtMath dtMatrix Mat x Vec\n");
    Printf("matA =\n");
    matA.Print();
    Printf("vecA =\n");
    vecA.Print();
    Printf("Result matA * vecA is \n");
    vecB.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#endif

    elapseTime.Start();
    for (int i = 0; i < 100000; i++)
        vecB = mat3A * vecA;
    elapseTime.Stop();
    Printf("// dtMath dtMatrix3 Mat x Vec\n");
    Printf("mat3A =\n");
    mat3A.Print();
    Printf("vecA =\n");
    vecA.Print();
    Printf("Result mat3A * vecA is \n");
    vecB.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#endif
}

void MatTranspose()
{
    float a[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    dtMatrix<3, 3, float> matA(a, sizeof(a));
    dtMatrix<3, 3, float> matC;

    dtMat3 mat3A(a, sizeof(a));
    dtMat3 mat3C;

    CdhTimeCheck elapseTime;

    elapseTime.Start();
    for (int i = 0; i < 100000; i++)
        matC = matA.Transpose();
    elapseTime.Stop();
    Printf("// dtMath dtMatrix Mat Transpose\n");
    Printf("Result matA.Transpose is \n");
    matC.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#endif

    elapseTime.Start();
    for (int i = 0; i < 100000; i++)
        mat3C = mat3A.Transpose();
    elapseTime.Stop();
    Printf("// dtMath dtMatrix3 Mat Transpose\n");
    Printf("Result mat3A.Transpose is \n");
    mat3C.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#endif
}

void MatInv()
{
    float a[9] = {3, 2, 1, 1, 6, 5, 2, 3, 9}; // full rank

    dtMatrix<3, 3, float> matA(a, sizeof(a));
    dtMatrix<3, 3, float> matC;

    dtMat3 mat3A(a, sizeof(a));
    dtMat3 mat3C;

    CdhTimeCheck elapseTime;

    elapseTime.Start();
    for (int i = 0; i < 100000; i++)
        matC = matA.Inv();
    elapseTime.Stop();
    Printf("// dtMath dtMatrix Mat Inverse\n");
    Printf("matA =\n");
    matA.Print();
    Printf("Result is \n");
    matC.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#endif

    elapseTime.Start();
    for (int i = 0; i < 100000; i++)
        mat3C = mat3A.Inv();
    elapseTime.Stop();
    Printf("// dtMath dtMatrix3 Mat Inverse\n");
    Printf("mat3A =\n");
    mat3A.Print();
    Printf("Result is\n");
    mat3C.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#endif
}

void MatPInv()
{
    float a[9] = {3, 2, 1, 1, 6, 5, 2, 12, 10}; // rank2

    dtMatrix<3, 3, float> matA(a, sizeof(a));
    dtMatrix<3, 3, float> matC;

    dtMat3 mat3A(a, sizeof(a));
    dtMat3 mat3C;

    CdhTimeCheck elapseTime;

    elapseTime.Start();
    for (int i = 0; i < 100000; i++)
        matC = matA.PInv();
    elapseTime.Stop();
    Printf("// dtMath dtMatrix Mat Pseudo Inverse\n");
    Printf("matA =\n");
    matA.Print();
    Printf("Result is \n");
    matC.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#endif

    elapseTime.Start();
    for (int i = 0; i < 100000; i++)
        mat3C = mat3A.PInv();
    elapseTime.Stop();
    Printf("// dtMath dtMatrix3 Mat Pseudo Inverse\n");
    Printf("mat3A =\n");
    mat3A.Print();
    Printf("Result is \n");
    mat3C.Print();
#if defined(_WIN32) || defined(__linux__)
    Printf("\033[33mTotal elapsed time is %fms\033[0m\n", elapseTime.GetElapsedTime_msec());
    Printf("\033[33mAvg elapsed time is %fms\033[0m\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#else
    Printf("Total elapsed time is %fms\n", elapseTime.GetElapsedTime_msec());
    Printf("Avg elapsed time is %fms\n\n", elapseTime.GetElapsedTime_msec() / 100000);
#endif
}