#include "testPrint.h"
#include "testVector.h"
#include "./dtMath/dtMath.h"

void Test_Vector()
{
    PrintTitle("Test Vector");
    VecInit();
    VecMemberFunc();
    VecMemberAccessOperator();
    VecArithmeticSum();
    VecArithmeticProduct();
    VecComparisonOperator();
}

void VecInit()
{
    PrintHeading("Initialize Vector ");
    Printf("CommaInit, vec << 1, 2, 3, 4, 5, 6;\n");
    Printf("vec = \n");

    CdtVector<6> vec;
    vec << 1, 2, 3, 4, 5, 6;

    vec.Print();
    Println;
}

void VecMemberFunc()
{
    PrintHeading("Vector Member Functions ");

    CdtVector<6> vec61;
    CdtVector3<> vec3;
    CdtVector4<> vec4;
    CdtVector6<> vec6;

    Printf("/* Class Create: with array */\n");
    float v[3] = { 1,2,3 };
    CdtVector<3> v31(v, sizeof(v));
    Printf("array[3] = {1,2,3}\n");
    Printf("vec =\n");
    v31.Print();
    Println;

    Printf("/* Class Create: with CdtVector */\n");
    CdtVector<3> v31_2(v31);
    Printf("vec = [1,2,3]T\n");
    Printf("vec =\n");
    v31_2.Print();
    Println;

    Printf("/* Class Create: with CdtMatrix */\n");
    CdtMatrix<3, 1> m(v, sizeof(v));
    CdtVector<3> v31_3(m);
    Printf("matrix31 = [1,2,3]T\n");
    Printf("vec =\n");
    v31_3.Print();
    Println;

    Printf("/* Function: SetZero() */\n");
    v31.SetZero();
    Printf("vec =\n");
    v31.Print();
    Println;

    Printf("/* Function: SetFill() */\n");
    Printf("value is 100\n");
    Printf("vec =\n");
    v31.SetFill(100);
    v31.Print();
    Println;

    Printf("/* Function: SetElement(const m_type *element) */\n");
    Printf("element[3] = {1,2,3}\n");
    Printf("vec =\n");
    v31.SetElement(v, sizeof(v));
    v31.Print();
    Println;

    Printf("/* Function: SetElement(const CdtVector &v) */\n");
    Printf("vec = [1,2,3]T\n");
    v31.SetElement(v31_2);
    v31.Print();
    Println;

    Printf("/* Function: SetElement(const CdtMatrix<m_row, 1, m_type> &v) */\n");
    Printf("matrix = [1,2,3]T\n");
    v31.SetElement(m);
    v31.Print();
    Println;

    Printf("/* Function: SetBlock(const uint16_t idxRow, const CdtVector<row, m_type> &v) */\n");
    vec61.SetBlock(0, v31);
    vec61.SetBlock(3, v31_2);
    Printf("v1 =\n");
    v31.Print();
    Printf("v2 =\n");
    v31_2.Print();
    Printf("SetBlock(0,v1); SetBlock(3,v2);\n");
    Printf("vec =\n");
    vec61.Print();
    Println;

    Printf("/* Function: SetBlock(const uint16_t idxRow, const CdtVector3<m_type, 3> &v) */\n");
    vec3.SetElement(1, 2, 3);
    vec61.SetBlock(0, vec3);
    vec61.SetBlock(3, vec3);
    Printf("v3 =\n");
    vec3.Print();
    Printf("SetBlock(0,v3); SetBlock(3,v3);\n");
    Printf("vec =\n");
    vec61.Print();
    Println;

    Printf("/* Function: SetBlock(const uint16_t idxRow, const CdtVector4<m_type, 4> &v) */\n");
    vec4.SetElement(1, 2, 3, 4);
    vec61.SetBlock(0, vec4);
    vec61.SetBlock(3, vec4);
    Printf("v4 =\n");
    vec4.Print();
    Printf("SetBlock(0,v4); SetBlock(3,v4);\n");
    Printf("vec =\n");
    vec61.Print();
    Println;

    Printf("/* Function: SetBlock(const uint16_t idxRow, const CdtVector4<m_type, 4> &v) */\n");
    vec6.SetElement(1, 2, 3, 4, 5, 6);
    vec61.SetBlock(0, vec6);
    Printf("v6 =\n");
    vec6.Print();
    Printf("SetBlock(0,v6);\n");
    Printf("vec =\n");
    vec61.Print();
    Println;

    Printf("/* Function: SetBlock(const uint16_t idxRow, const CdtMatrix<row, 1, m_type> &v */\n");
    float m41[4] = { 4, 3, 2, 1 };
    CdtMatrix<4, 1> mat41(m41, sizeof(m41));
    vec61.SetBlock(0, mat41);
    Printf("mat41 =\n");
    mat41.Print();
    Printf("SetBlock(0,mat41);\n");
    Printf("vec =\n");
    vec61.Print();
    Println;

    Printf("/* Function: SetSwap(const uint16_t i, const uint16_t j) */\n");
    Printf("vec =\n");
    vec61.Print();
    vec61.SetSwap(3, 0);
    vec61.SetSwap(1, 2);
    Printf("SetSwap(3,0); SetSwap(1,2);\n");
    Printf("vec =\n");
    vec61.Print();
    Println;

    Printf("/* Function: SetNormalize() */\n");
    Printf("vec =\n");
    vec61.Print();
    vec61.SetNormalize();
    Printf("Normalized vec =\n");
    vec61.Print();
    Println;

    Printf("/* Function: GetDim() */\n");
    Printf("vec =\n");
    vec61.Print();
    Printf("Dimension is %d\n", vec61.GetDim());
    Println;

    Printf("/* Function: GetNorm() */\n");
    vec61.SetBlock(0, vec6);
    Printf("vec =\n");
    vec61.Print();
    Printf("Norm is %.6f\n", vec61.GetNorm());
    Println;

    Printf("/* Function: GetSqNorm() */\n");
    Printf("vec =\n");
    vec61.Print();
    Printf("Squared Norm is %.6f\n", vec61.GetSqNorm());
    Println;

    Printf("/* Function: GetSum() */\n");
    Printf("vec =\n");
    vec61.Print();
    Printf("Sum of elements is %.6f\n", vec61.GetSum());
    Println;

    Printf("/* Function: GetNormalized() */\n");
    Printf("vec =\n");
    vec61.Print();
    Printf("Normalized vec =\n");
    vec61.GetNormalized().Print();
    Println;

    Printf("/* Function: Transpose() */\n");
    Printf("vec =\n");
    vec61.Print();
    Printf("Transposed vec =\n");
    vec61.Transpose().Print();
    Println;
}

void VecMemberAccessOperator()
{
    CdtVector<6> vec61;

    PrintHeading("Vector Member Access operator ");
    Printf("/* Operator: () */\n");
    Printf("vec(i) = i; Printf(""%%f"", vec(i));\n");

    for (uint16_t i = 0; i < 6; i++)
    {
        vec61(i) = i;
        Printf("%7.3f ", vec61(i));
    }

    Println;
    Println;
}

void VecArithmeticSum()
{
    PrintHeading("Vector Arithmetic operators - Sum ");
    CdtVector<6> vec61;
    CdtVector<3> vec31;
    CdtVector<4> vec41;
    CdtVector3<> vec3;
    CdtVector4<> vec4;
    CdtVector6<> vec6;
    CdtMatrix<6, 1> mat61;
    float mat[6] = { 1,2,3,4,5,6 };

    vec61 << 1, 2, 3, 4, 5, 6;
    vec31 << 1, 2, 3;
    vec41 << 1, 2, 3, 4;
    vec3.SetElement(1, 2, 3);
    vec4.SetElement(1, 2, 3, 4);
    vec6.SetElement(1, 2, 3, 4, 5, 6);
    mat61.SetElement(mat, sizeof(mat));

    Printf("vec =\n");
    vec61.Print();
    Println;

    Printf("/* Operator: -() */\n");
    Printf("-vec =\n");
    (-vec61).Print();
    Println;

    Printf("/* Operator: +(CdtVector &v) */\n");
    Printf("vec + (-vec) =\n");
    (vec61 + (-vec61)).Print();
    Println;

    Printf("/* Operator: -(CdtVector &v) */\n");
    Printf("vec - vec =\n");
    (vec61 - vec61).Print();
    Println;

    Printf("/* Operator: +(CdtVector3 &v) */\n");
    Printf("vec = [1,2,3]T\n");
    Printf("vec3 = [1,2,3]T\n");
    Printf("vec + (-vec3) =\n");
    (vec31 + (-vec3)).Print();
    Println;

    Printf("/* Operator: -(CdtVector3 &v) */\n");
    Printf("vec - vec3 =\n");
    (vec31 - vec3).Print();
    Println;

    Printf("/* Operator: +(CdtVector4 &v) */\n");
    Printf("vec = [1,2,3,4]T\n");
    Printf("vec4 = [1,2,3,4]T\n");
    Printf("vec + (-vec4) =\n");
    (vec41 + (-vec4)).Print();
    Println;

    Printf("/* Operator: -(CdtVector4 &v) */\n");
    Printf("vec = [1,2,3,4]T\n");
    Printf("vec4 = [1,2,3,4]T\n");
    Printf("vec - vec4 =\n");
    (vec41 - vec4).Print();
    Println;

    Printf("/* Operator: +(CdtVector6 &v) */\n");
    Printf("vec = [1,2,3,4,5,6]T\n");
    Printf("vec6 = [1,2,3,4,5,6]T\n");
    Printf("vec + (-vec6) =\n");
    (vec61 + (-vec6)).Print();
    Println;

    Printf("/* Operator: -(CdtVector6 &v) */\n");
    Printf("vec = [1,2,3,4,5,6]T\n");
    Printf("vec6 = [1,2,3,4,5,6]T\n");
    Printf("vec - vec6 =\n");
    (vec61 - vec6).Print();
    Println;

    Printf("/* Operator: +(CdtMatrix<m_row, 1, m_type> &v) */\n");
    Printf("vec = [1,2,3,4,5,6]T\n");
    Printf("mat61 = [1,2,3,4,5,6]T\n");
    Printf("vec + (-mat61) =\n");
    (vec61 + (-mat61)).Print();
    Println;

    Printf("/* Operator: -(CdtMatrix<m_row, 1, m_type> &v) */\n");
    Printf("vec = [1,2,3,4,5,6]T\n");
    Printf("mat61 = [1,2,3,4,5,6]T\n");
    Printf("vec - mat61 =\n");
    (vec61 - mat61).Print();
    Println;
}

void VecArithmeticProduct()
{
    CdtVector<6> vec61;
    CdtVector<3> vec31;
    CdtVector<4> vec41;
    CdtVector3<> vec3;
    CdtVector4<> vec4;
    CdtVector6<> vec6;
    CdtMatrix<6, 1> mat61;
    float mat[6] = { 1,1,1,1,1,1 };

    vec61.SetFill(1);
    vec31.SetFill(1);
    vec41.SetFill(1);
    vec3.SetFill(1);
    vec4.SetFill(1);
    vec6.SetFill(1);
    mat61.SetElement(mat, sizeof(mat));

    PrintHeading("Vector Arithmetic operators - Product ");

    Printf("/* Operator: *(const m_type s) */\n");
    Printf("vec = [1,1,1,1,1,1]T\n");
    Printf("vec * 0.1 =\n");
    (vec61 * 0.1f).Print();
    Println;

    Printf("/* Operator: *(const m_type s) */\n");
    Printf("vec = [1,1,1,1,1,1,]T\n");
    Printf("vec * 0.1 =\n");
    (vec61 / 10.0f).Print();
    Println;

    Printf("/* Operator: *(const m_type s, CdtVector& v) */\n");
    Printf("vec = [1,1,1,1,1,1]T\n");
    Printf("10 * vec =\n");
    (10.0f * vec61).Print();
    Println;

    Printf("/* Operator: *(const CdtMatrix<1, col, m_type>& m) */\n");
    Printf("vec = [1,1,1,1,1,1]T, mat16 = [1,1,1,1,1,1]\n");
    Printf("vec * mat16 =\n");
    CdtMatrix<1, 6> vec16;
    vec16.SetFill(1);
    (vec61 * vec16).Print();
    Println;

    Printf("/* Operator: dot(const CdtVector& v) */\n");
    Printf("vec = [1,1,1,1,1,1]T\n");
    Printf("vec dot vec = %f\n", vec61.dot(vec61));
    Println;

    Printf("/* Operator: dot(const CdtVector3& v) */\n");
    Printf("vec = [1,1,1]T, vec3 = [1,1,1]T\n");
    Printf("vec dot vec3 = %f\n", vec31.dot(vec3));
    Println;

    Printf("/* Operator: dot(const CdtVector4& v) */\n");
    Printf("vec = [1,1,1,1]T, vec4 = [1,1,1,1]T\n");
    Printf("vec dot vec4 = %f\n", vec41.dot(vec4));
    Println;

    Printf("/* Operator: dot(const CdtVector6& v) */\n");
    Printf("vec = [1,1,1,1,1,1]T, vec6 = [1,1,1,1,1,1]T\n");
    Printf("vec dot vec6 = %f\n", vec61.dot(vec6));
    Println;

    Printf("/* Operator: dot(const CdtMatrix<m_row, 1, m_type>& v) */\n");
    Printf("vec = [1,1,1,1,1,1]T, mat61 = [1,1,1,1,1,1]T\n");
    Printf("vec dot mat61 = %f\n", vec61.dot(mat61));
    Println;
}

void VecComparisonOperator()
{
    CdtVector<6> vec61;
    CdtVector<3> vec31;
    CdtVector<4> vec41;
    CdtVector3<> vec3;
    CdtVector4<> vec4;
    CdtVector6<> vec6;
    CdtMatrix<6, 1> mat61;
    float mat[6] = { 1,2,3,4,5,6 };

    vec61 << 1., 2., 3., 4., 5., 6.;
    vec31 << 1., 2., 3.;
    vec41 << 1., 2., 3., 4.;
    vec3.SetElement(1, 2, 3);
    vec4.SetElement(1, 2, 3, 4);
    vec6.SetElement(1, 2, 3, 4, 5, 6);
    mat61.SetElement(mat, sizeof(mat));

    PrintHeading("Vector Comparison operators ");
    Printf("/* Operator: ==(CdtVector &v) and !=(CdtVector &v) */\n");

    Printf("vec61 == vec61? -> ");
    if (vec61 == vec61) Printf("true");
    else Printf("false");
    Println;

    Printf("vec61 != vec61? -> ");
    if (vec61 != vec61) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(CdtVector3 &v) and !=(CdtVector3 &v) */\n");

    Printf("vec31 = vec3\n");
    Printf("vec31 == vec3? -> ");
    if (vec31 == vec3) Printf("true");
    else Printf("false");
    Println;

    Printf("vec31 != vec3? -> ");
    if (vec31 != vec3) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(CdtVector4 &v) and !=(CdtVector4 &v) */\n");
    
    Printf("vec41 = vec4\n");
    Printf("vec41 == vec4? -> ");
    if (vec41 == vec4) Printf("true");
    else Printf("false");
    Println;

    Printf("vec41 != vec4? -> ");
    if (vec41 != vec4) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(CdtVector6 &v) and !=(CdtVector6 &v) */\n");
    
    Printf("vec61 = vec6\n");
    Printf("vec61 == vec6? -> ");
    if (vec61 == vec6) Printf("true");
    else Printf("false");
    Println;

    Printf("vec61 != vec6? -> ");
    if (vec61 != vec6) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(CdtMatrix &v) and !=(CdtMatrix &v) */\n");
    
    Printf("vec61 = mat61\n");
    Printf("vec61 == mat61? -> ");
    if (vec61 == mat61) Printf("true");
    else Printf("false");
    Println;

    Printf("vec61 != mat61? -> ");
    if (vec61 != mat61) Printf("false");
    else Printf("true");
    Println;
    Println;
}