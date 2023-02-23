#include "testPrint.h"
#include "testMatrix.h"

#include <dtMath/dtMath.h>


CdtMatrix<3, 3> mat33;
CdtMatrix<3, 4> mat34;
CdtMatrix<4, 3> mat43;
CdtMatrix<4, 4> mat44;
CdtMatrix<6, 3> mat63;
CdtMatrix<6, 4> mat64;
CdtMatrix<6, 6> mat66;

CdtMatrix3<> mat3;
CdtRotation<> rot;
CdtTransform<> trf;

CdtVector<3> vec31;
CdtVector<4> vec41;
CdtVector<6> vec61;
CdtVector3<> vec3;
CdtVector4<> vec4;
CdtVector6<> vec6;

void Test_Matrix()
{
    PrintTitle("Test Matrix");
    MatInit();
    MatMemberFunc();
    MatMemberAccessOperator();
    MatArithmeticSum();
    MatArithmeticScalarProduct();
    MatArithmeticMatProduct();
    MatArithmeticVecProduct();
    MatComparisonOperator();
}

void MatInit()
{
    PrintHeading("Initialize Matrix ");
    Printf("CommaInit, mat << 1, 2, 3, 4, 5, 6, 7, 8, 9;\n");
    Printf("mat = \n");

    CdtMatrix<3, 3> mat;
    mat << 1, 2, 3,
        4, 5, 6,
        7, 8, 9;

    mat.Print();
    Println;
}

void MatMemberFunc()
{
    PrintHeading("Matrix Member Functions ");
    Printf("/* Class Create: with array */\n");
    float a[2 * 2] = { 1,2,3,4 };
    CdtMatrix<2, 2> a22(a, sizeof(a));
    Printf("a[2*2] = {1,2,3,4}\n");
    Printf("a22 = \n");
    a22.Print();
    Println;

    Printf("/* Class Create: with diagonal */\n");
    float d[3] = { 1,2,3 };
    CdtMatrix<3, 3> d33('d', d, sizeof(d));
    Printf("d3 = {1,2,3}\n");
    Printf("d33 = \n");
    d33.Print();
    Println;

    Printf("/* Class Create: with matrix */\n");
    CdtMatrix<3, 3> m33(d33);
    Printf("arg mat is d33\n");
    Printf("m33 = \n");
    m33.Print();
    Println;

    Printf("/* Function: SetZero() */\n");
    Printf("mat34 = \n");
    mat34.SetZero();
    mat34.Print();
    Println;

    Printf("/* Function: SetIdentity() */\n");
    Printf("mat33 = \n");
    mat33.SetIdentity();
    mat33.Print();
    Println;

    Printf("mat34 = \n");
    mat34.SetIdentity();
    mat34.Print();
    Println;

    Printf("mat43 = \n");
    mat43.SetIdentity();
    mat43.Print();
    Println;

    Printf("/* Function: SetDiagonal() */\n");
    Printf("diag = [1,2,3,4,5,6]T\n");
    float dia[6] = { 1,2,3,4,5,6 };
    Printf("mat33 = \n");
    mat33.SetDiagonal(dia, sizeof(dia));
    mat33.Print();
    Println;

    Printf("mat34 = \n");
    mat34.SetDiagonal(dia, sizeof(dia));
    mat34.Print();
    Println;

    Printf("mat43 = \n");
    mat43.SetDiagonal(dia, sizeof(dia));
    mat43.Print();
    Println;

    Printf("/* Function: SetFill() */\n");
    Printf("value = 100\n");
    float value = 100;
    Printf("mat33 = \n");
    mat33.SetFill(value);
    mat33.Print();
    Println;

    Printf("/* Function: SetElement() */\n");
    Printf("elem[9] = {9,8,7,6,5,4,3,2,1};\n");
    float elem33[3 * 3] = { 9,8,7,6,5,4,3,2,1 };
    Printf("mat33 = \n");
    mat33.SetElement(elem33, sizeof(elem33));
    mat33.Print();
    Println;

    Printf("mat34 = \n");
    Printf("elem[12] = {0,1,2,3,4,5,6,7,8,9,10,11};\n");
    float elem34[3 * 4] = { 0,1,2,3,4,5,6,7,8,9,10,11 };
    mat34.SetElement(elem34, sizeof(elem34));
    mat34.Print();
    Println;

    Printf("mat43 = \n");
    Printf("elem[12] = {0,1,2,3,4,5,6,7,8,9,10,11};\n");
    float elem43[4 * 3] = { 0,1,2,3,4,5,6,7,8,9,10,11 };
    mat43.SetElement(elem43, sizeof(elem43));
    mat43.Print();
    Println;

    Printf("/* Function: SetBlock() */\n");
    Printf("mat34.SetBlock(0,0,mat33) = \n");
    mat34.SetBlock(0, 0, mat33);
    mat34.Print();
    Println;

    Printf("mat43.SetBlock(1,0,mat33) = \n");
    mat43.SetBlock(1, 0, mat33);
    mat43.Print();
    Println;

    Printf("mat63.SetBlock(0, 0, mat3); \n");
    Printf("mat63.SetBlock(3, 0, Rz(30)); = \n");
    mat3.SetFill(10);
    rot.SetElement(Z_AXIS, 30 * DEG2RADf);
    mat63.SetBlock(0, 0, mat3);
    mat63.SetBlock(3, 0, rot);
    mat63.Print();
    Println;

    Printf("/* Function: SetSwapVec() */\n");
    Printf("mat63.SetSwapRowVec(0, 5) = \n");
    mat63.SetSwapRowVec(0, 5);
    mat63.Print();
    Println;

    Printf("mat63.SetSwapColVec(1, 2) = \n");
    mat63.SetSwapColVec(1, 2);
    mat63.Print();
    Println;

    Printf("/* Function: GetSize() */\n");
    Printf("mat63.GetRowSize() = ");
    Printf("%d\n", mat63.GetRowSize());
    Printf("mat63.GetColSize() = ");
    Printf("%d\n", mat63.GetColSize());
    Println;

    Printf("/* Function: GetBlock() */\n");
    Printf("mat63.GetBlock<3,3>(3,0) = \n");
    mat33 = mat63.GetBlock<3, 3>(3, 0);
    mat33.Print();
    Println;

    Printf("mat63.GetBlock<2,2>(0,0) = \n");
    mat63.GetBlock<2, 2>(0, 0).Print();
    Println;

    Printf("/* Function: GetVec() */\n");
    Printf("mat63.GetRowVec(1) = \n");
    mat63.GetRowVec(1).Print();
    Println;

    Printf("mat63.GetColVec(2) = \n");
    mat63.GetColVec(2).Print();
    Println;

    Printf("/* Function: GetCscMat() = \n");
    mat63.GetCscMat().Print();
    Println;

    Printf("/* Function: Transpose() */\n");
    mat63.Transpose().Print();
    Println;
}

void MatMemberAccessOperator()
{
    PrintHeading("Matrix Access operator ");
    Printf("mat63.Print() = \n");
    mat63.Print();
    Println;

    Printf("/* Operator: () */\n");
    Printf("mat63(i,j) = 100; & Printf(""%%f"", mat63(i,j));\n");
    for (uint16_t i = 0; i < 6; i++)
    {
        for (uint16_t j = 0; j < 3; j++)
        {
            mat63(i, j) = 100;
            Printf("%7.3f ", mat63(i, j));
        }
        Printf("\n");
    }
    Println;
}

void MatArithmeticSum()
{
    PrintHeading("Matrix Arithmetic operators - Sum ");
    Printf("mat63.SetFill(1)\n");
    mat63.SetFill(1);
    Printf("/* Operator: -() */\n");
    Printf("(-mat63)\n");
    (-mat63).Print();
    Println;

    Printf("/* Operator: +(CdtMatrix &m) */\n");
    Printf("mat63 + (-mat63)\n");
    (mat63 + (-mat63)).Print();
    Println;

    Printf("/* Operator: -(CdtMatrix &m) */\n");
    Printf("mat63 - mat63\n");
    (mat63 - mat63).Print();
    Println;

    Printf("/* Operator: +(CdtMatrix3 &m) */\n");
    Printf("mat33.SetFill(1)\n");
    Printf("mat3.SetFill(-1)\n");
    Printf("mat33 + mat3 =\n");
    mat33.SetFill(1);
    mat3.SetFill(-1);
    (mat33 + mat3).Print();
    Println;

    Printf("/* Operator: -(CdtMatrix3 &m) */\n");
    Printf("mat33.SetFill(1)\n");
    Printf("mat3.SetFill(-1)\n");
    Printf("mat33 - (-mat3) =\n");
    mat33.SetFill(1);
    mat3.SetFill(-1);
    (mat33 - (-mat3)).Print();
    Println;

    Printf("/* Operator: +(CdtRotation &m) */\n");
    Printf("mat33.SetFill(1)\n");
    Printf("rot.SetFill(-1)\n");
    Printf("mat33 + rot = \n");
    mat33.SetFill(1);
    rot.SetFill(-1);
    (mat33 + rot).Print();
    Println;

    Printf("/* Operator: -(CdtRotation &m) */\n");
    Printf("mat33.SetFill(1)\n");
    Printf("rot.SetFill(-1)\n");
    Printf("mat33 - (-rot) = \n");
    mat33.SetFill(1);
    rot.SetFill(-1);
    (mat33 - (-rot)).Print();
    Println;

    Printf("/* Operator: +(CdtTransform &m) */\n");
    Printf("mat44.SetIdentity()\n");
    Printf("trf.SetIdentity()\n");
    Printf("mat44 + trf = \n");
    mat44.SetIdentity();
    trf.SetIdentity();
    (mat44 + trf).Print();
    Println;

    Printf("/* Operator: -(CdtTransform &m) */\n");
    Printf("mat44.SetIdentity()\n");
    Printf("trf.SetIdentity()\n");
    Printf("mat44 - trf = \n");
    mat44.SetIdentity();
    trf.SetIdentity();
    (mat44 - trf).Print();
    Println;
}

void MatArithmeticScalarProduct()
{
    PrintHeading("Matrix Arithmetic operators - Scalar product ");

    Printf("/* Operator: *(m_type s) and /(m_type &s) */\n");
    Printf("mat44.SetFill(1)\n");
    Printf("(mat44 * 0.5) - (mat44 / 2.0) =\n");
    mat44.SetFill(1);
    ((mat44 * 0.5) - (mat44 / 2.0)).Print();
    Println;

    Printf("/* Operator: *(m_type s, CdtMatrix &m) */\n");
    Printf("mat44.SetFill(1)\n");
    Printf("(0.5f * mat44) - (mat44 / 2.0) =\n");
    mat44.SetFill(1);
    ((0.5f * mat44) - (mat44 / 2.0)).Print();
    Println;
}

void MatArithmeticMatProduct()
{
    PrintHeading("Matrix Arithmetic operators - Mat x Mat ");
    Printf("/* Matrix Product: operator *(const CdtMatrix<> &m) */\n");
    float m63[18];
    float m34[12];
    for (int i = 0; i < 18; i++) m63[i] = i + 1.0f;
    for (int i = 0; i < 12; i++) m34[i] = i + 1.0f;
    mat63.SetElement(m63, sizeof(m63));
    mat34.SetElement(m34, sizeof(m34));
    Printf("mat63 = \n");
    mat63.Print();
    Printf("mat34 = \n");
    mat34.Print();
    Printf("mat64 = mat63 * mat34\n");
    mat64 = mat63 * mat34;
    mat64.Print();
    Println;

    Printf("/* Matrix Product: operator *(const CdtMatrix3<> &m) */\n");
    float m43[12];
    float m33[9];
    for (int i = 0; i < 12; i++) m43[i] = i + 1.0f;
    for (int i = 0; i < 9; i++) m33[i] = i + 1.0f;
    mat43.SetElement(m43, sizeof(m43));
    mat3.SetElement(m33, sizeof(m33));
    Printf("mat43 =\n");
    mat43.Print();
    Printf("mat3 =\n");
    mat3.Print();
    Printf("mat43 = mat43 * mat3\n");
    mat43 = mat43 * mat3;
    mat43.Print();
    Println;

    Printf("/* Matrix Product: operator *(const CdtRotation<> &m) */\n");
    for (int i = 0; i < 12; i++) m43[i] = i + 1.0f;
    mat43.SetElement(m43, sizeof(m43));
    rot.SetElement(AXIS3(2, 1, 0), 10.0f * DEG2RADf, 20.0f * DEG2RADf, 30.0f * DEG2RADf);
    Printf("mat43 =\n");
    mat43.Print();
    Printf("rot = Rz(10) * Ry(20) * Rx(30) =\n");
    rot.Print();
    Printf("mat43 = mat43 * rot\n");
    mat43 = mat43 * rot;
    mat43.Print();
    Println;

    Printf("/* Matrix Product: operator *(const CdtTransform<> &m) */\n");
    float m64[24];
    for (int i = 0; i < 24; i++) m64[i] = i + 1.0f;
    mat64.SetElement(m64, sizeof(m64));
    trf.SetIdentity();
    Printf("mat64 =\n");
    mat64.Print();
    Printf("trf = Identity =\n");
    Printf("mat64 = mat64 * trf\n");
    mat64 = mat64 * trf;
    mat64.Print();
    Println;
}

void MatArithmeticVecProduct()
{
    PrintHeading("Matrix Arithmetic operators - Mat x Vec ");
    float m63[18];
    float m64[24];
    float m66[36];

    Printf("/* Vector Product: operator *(const CdtVector<> &v) */\n");
    for (int i = 0; i < 18; i++) m63[i] = i + 1.0f;
    mat63.SetElement(m63, sizeof(m63));
    vec31 << 1, 2, 3;
    Printf("mat63 =\n");
    mat63.Print();
    Printf("vec31 =\n");
    vec31.Print();
    Printf("vec61 = mat63 * vec31\n");
    vec61 = mat63 * vec31;
    vec61.Print();
    Println;

    Printf("/* Vector Product: operator *(const CdtVector3<> &v) */\n");
    for (int i = 0; i < 18; i++) m63[i] = i + 1.0f;
    mat63.SetElement(m63, sizeof(m63));
    vec3.SetElement(1.0, 2.0, 3.0);
    Printf("mat63 =\n");
    mat63.Print();
    Printf("vec3 =\n");
    vec3.Print();
    Printf("vec61 = mat63 * vec3\n");
    vec61 = mat63 * vec3;
    vec61.Print();
    Println;

    Printf("/* Vector Product: operator *(const CdtVector4<> &v) */\n");
    for (int i = 0; i < 24; i++) m64[i] = i + 1.0f;
    mat64.SetElement(m64, sizeof(m64));
    vec4.SetElement(1.0, 2.0, 3.0, 4.0);
    Printf("mat64 =\n");
    mat64.Print();
    Printf("vec4 =\n");
    vec4.Print();
    Printf("vec61 = mat64 * vec4\n");
    vec61 = mat64 * vec4;
    vec61.Print();
    Println;

    Printf("/* Vector Product: operator *(const CdtVector6<> &v) */\n");
    for (int i = 0; i < 36; i++) m66[i] = i + 1.0f;
    mat66.SetElement(m66, sizeof(m66));
    vec6.SetElement(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    Printf("mat66 =\n");
    mat66.Print();
    Printf("vec6 =\n");
    vec6.Print();
    Printf("vec61 = mat66 * vec6\n");
    vec61 = mat66 * vec6;
    vec61.Print();
    Println;
}

void MatComparisonOperator()
{
    PrintHeading("Matrix Comparison operators ");
    Printf("/* Operator: ==(CdtMatrix &m) and !=(CdtMatrix &m) */\n");

    Printf("mat63 == mat63? -> ");
    if (mat63 == mat63) Printf("true");
    else Printf("false");
    Println;

    Printf("mat63 != mat63? -> ");
    if (mat63 != mat63) Printf("true");
    else Printf("false");
    Println;

    Printf("/* Operator: ==(CdtMatrix3 &m) and !=(CdtMatrix3 &m) */\n");
    mat3.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    mat33 = mat3;

    Printf("mat33 = mat3\n");
    Printf("mat33 == mat3? -> ");
    if (mat33 == mat3) Printf("true");
    else Printf("false");
    Println;

    Printf("mat33 != mat3? -> ");
    if (mat33 != mat3) Printf("true");
    else Printf("false");
    Println;

    Printf("/* Operator: ==(CdtRotation &m) and !=(CdtRotation &m) */\n");
    rot.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    mat33 = rot;

    Printf("mat33 = rot\n");
    Printf("mat33 == rot? -> ");
    if (mat33 == rot) Printf("true");
    else Printf("false");
    Println;

    Printf("mat33 != rot? -> ");
    if (mat33 != rot) Printf("true");
    else Printf("false");
    Println;

    Printf("/* Operator: ==(CdtTransform &m) and !=(CdtTransform &m) */\n");
    trf.SetElement(rot, vec3);
    mat44 = trf;

    Printf("mat33 = trf\n");
    Printf("mat33 == trf? -> ");
    if (mat44 == trf) Printf("true");
    else Printf("false");
    Println;

    Printf("mat33 == trf? -> ");
    if (mat44 != trf) Printf("true");
    else Printf("false");
    Println;
    Println;
}
