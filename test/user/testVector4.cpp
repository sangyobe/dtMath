#include "testPrint.h"
#include "testVector4.h"

#include <dtMath/dtMath.h>


void Test_Vector4()
{
    PrintTitle("Test Vector4");
    Vec4Init();
    Vec4MemberFunc();
    Vec4MemberAccessOperator();
    Vec4ArithmeticSum();
    Vec4ArithmeticProduct();
    Vec4ComparisonOperator();
}

void Vec4Init()
{
    PrintHeading("Initialize Vector4 ");
    Printf("CommaInit, vec << 1, 2, 3, 4;\n");
    Printf("vec = \n");

    CdtVector4<> vec;
    vec << 1, 2, 3, 4;

    vec.Print();
    Println;
}

void Vec4MemberFunc()
{
    PrintHeading("Vector4 Member Functions ");

    float v[4] = { 1,2,3,4 };

    Printf("/* Class Create: () */\n");
    CdtVector4<> v1;
    v1.Print();
    Println;

    Printf("/* Class Create: with array */\n");
    CdtVector4<> v2(v, sizeof(v));
    Printf("array[4] = {1,2,3,4}\n");
    Printf("vec =\n");
    v2.Print();
    Println;

    Printf("/* Class Create: with element */\n");
    CdtVector4<> v3(1, 2, 3, 4);
    Printf("elements are 1, 2, 3, 4\n");
    Printf("vec =\n");
    v3.Print();
    Println;

    Printf("/* Class Create: with CdtVector4 */\n");
    CdtVector4<> v4(v3);
    Printf("vec = [1,2,3,4]T\n");
    Printf("vec =\n");
    v4.Print();
    Println;

    Printf("/* Class Create: with CdtVector */\n");
    CdtVector<4> v41(v, sizeof(v));
    CdtVector4<> v5(v41);
    Printf("vec = [1,2,3,4]T\n");
    Printf("vec =\n");
    v5.Print();
    Println;

    Printf("/* Class Create: with CdtMatrix */\n");
    CdtMatrix<4, 1> m41(v, sizeof(v));
    CdtVector4<> v6(m41);
    Printf("mat41 = [1,2,3,4]T\n");
    Printf("vec =\n");
    v6.Print();
    Println;

    Printf("/* Function: SetZero() */\n");
    Printf("vec =\n");
    v2.SetZero();
    v2.Print();
    Println;

    Printf("/* Function: SetFill(value) */\n");
    Printf("value is 1\n");
    Printf("vec =\n");
    v2.SetFill(1);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(* element) */\n");
    Printf("element[4] = {1,2,3,4}\n");
    Printf("vec =\n");
    v2.SetElement(v, sizeof(v));
    v2.Print();
    Println;

    Printf("/* Function: SetElement(v0, v1, v2, v3) */\n");
    Printf("v0=1, v1=2, v2=3, v3=4\n");
    Printf("vec =\n");
    v2.SetElement(1, 2, 3, 4);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(CdtVector4) */\n");
    Printf("vec = [1,2,3,4]T\n");
    Printf("vec =\n");
    v2.SetElement(v3);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(CdtVector) */\n");
    Printf("vec = [1,2,3,4]T\n");
    Printf("vec =\n");
    v2.SetElement(v41);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(CdtMatrix) */\n");
    Printf("mat41 = [1,2,3,4]T\n");
    Printf("vec =\n");
    v2.SetElement(m41);
    v2.Print();
    Println;

    Printf("/* Function: SetSwap(i, j) */\n");
    Printf("vec =\n");
    v2.Print();
    Printf("SetSwap(0,2);\n");
    Printf("vec =\n");
    v2.SetSwap(0, 2);
    v2.Print();
    Println;

    Printf("/* Function: SetNormalize() */\n");
    v2.SetFill(1);
    Printf("vec =\n");
    v2.Print();
    v2.SetNormalize();
    Printf("Normalized vec =\n");
    v2.Print();
    Println;

    Printf("/* Function: GetNorm() */\n");
    v2.SetFill(1);
    Printf("vec =\n");
    v2.Print();
    Printf("Norm is %f\n", v2.GetNorm());
    Println;

    Printf("/* Function: GetSqNorm() */\n");
    v2.SetFill(1);
    Printf("vec =\n");
    v2.Print();
    Printf("Squared Norm is %f\n", v2.GetSqNorm());
    Println;

    Printf("/* Function: GetSum() */\n");
    v2.SetFill(1);
    Printf("vec =\n");
    v2.Print();
    Printf("Sum of elements is %f\n", v2.GetSum());
    Println;

    Printf("/* Function: GetNormalized() */\n");
    v2.SetFill(1);
    Printf("vec =\n");
    v2.Print();
    Printf("Normalized vec =\n");
    v2.GetNormalized().Print();
    Println;

    Printf("/* Function: Transpose() */\n");
    v2.SetElement(1, 2, 3, 4);
    Printf("vec =\n");
    v2.Print();
    Printf("Transposed vec =\n");
    v2.Transpose().Print();
    Println;
}

void Vec4MemberAccessOperator()
{
    CdtVector4<> v1;

    PrintHeading("Vector4 Access operator ");
    Printf("/* Operator: () */\n");
    Printf("vec(i) = i; Printf(""%%f"", vec(i));\n");

    for (uint16_t i = 0; i < 4; i++)
    {
        v1(i) = i;
        Printf("%7.3f ", v1(i));
    }

    Println;
    Println;
}

void Vec4ArithmeticSum()
{
    PrintHeading("Vector4 Arithmetic operators - Sum ");

    float v[4] = { 1,2,3,4 };
    CdtVector4<> v1(1, 2, 3, 4);
    CdtVector<4> v41(v, sizeof(v));
    CdtMatrix<4, 1> m41(v, sizeof(v));

    Printf("vec =\n");
    v1.Print();
    Println;

    Printf("/* Operator: -() */\n");
    Printf("-vec =\n");
    (-v1).Print();
    Println;

    Printf("/* Operator: +(CdtVector4) */\n");
    Printf("vec + (-vec) =\n");
    (v1 + (-v1)).Print();
    Println;

    Printf("/* Operator: -(CdtVector4) */\n");
    Printf("vec - vec =\n");
    (v1 - v1).Print();
    Println;

    Printf("/* Operator: +(CdtVector) */\n");
    Printf("vec = [1,2,3,4]T, vec41 = [1,2,3,4]T\n");
    Printf("vec + (-vec41) =\n");
    (v1 + (-v41)).Print();
    Println;

    Printf("/* Operator: -(CdtVector) */\n");
    Printf("vec = [1,2,3,4]T, vec41 = [1,2,3,4]T\n");
    Printf("vec - vec41) =\n");
    (v1 - v41).Print();
    Println;

    Printf("/* Operator: +(CdtMatrix) */\n");
    Printf("vec = [1,2,3,4]T, mat41 = [1,2,3,4]T\n");
    Printf("vec + (-mat41) =\n");
    (v1 + (-m41)).Print();
    Println;

    Printf("/* Operator: -(CdtMatrix) */\n");
    Printf("vec = [1,2,3,4]T, mat41 = [1,2,3,4]T\n");
    Printf("vec - mat41 =\n");
    (v1 - m41).Print();
    Println;
}

void Vec4ArithmeticProduct()
{
    PrintHeading("Vector4 Arithmetic operators - Product ");

    float v[4] = { 1,2,3,4 };
    CdtVector4<> v1(1, 2, 3, 4);
    CdtVector4<> v2(v, sizeof(v));
    CdtVector<4> v41(v, sizeof(v));
    CdtMatrix<4, 1> m41(v, sizeof(v));
    CdtMatrix<1, 4> m14(v, sizeof(v));

    Printf("/* Operator: *(scalar) */\n");
    Printf("vec = [1,1,1,1]T\n");
    Printf("vec * 0.1 =\n");
    v1.SetFill(1);
    (v1 * 0.1f).Print();
    Println;

    Printf("/* Operator: /(scalar) */\n");
    Printf("vec = [1,1,1,1]T\n");
    Printf("vec / 10 =\n");
    v1.SetFill(1);
    (v1 / 10).Print();
    Println;

    Printf("/* Operator: *(scalar, CdtVector) */\n");
    Printf("vec = [1,1,1,1]T\n");
    Printf("0.1 * vec =\n");
    v1.SetFill(1);
    (0.1f * v1).Print();
    Println;

    Printf("/* Operator: *(CdtMatrix<1, col>) */\n");
    Printf("vec = [1,1,1,1]T, mat14 = [1,1,1,1]\n");
    Printf("vec * mat14 =\n");
    v1.SetFill(1);
    m14.SetFill(1);
    (v1 * m14).Print();
    Println;

    Printf("/* Operator: dot(CdtVector4) */\n");
    Printf("vec1 = [1,1,1,1]T, vec2 = [1,1,1,1]T\n");
    v1.SetFill(1);
    v2.SetFill(1);
    Printf("vec1 dot vec2 = %f\n", v1.dot(v2));
    Println;

    Printf("/* Operator: dot(CdtVector) */\n");
    Printf("vec1 = [1,1,1,1]T, vec41 = [1,1,1,1]T\n");
    v1.SetFill(1);
    v41.SetFill(1);
    Printf("vec1 dot vec41 = %f\n", v1.dot(v41));
    Println;

    Printf("/* Operator: dot(CdtMatrix) */\n");
    Printf("vec1 = [1,1,1,1]T, mat41 = [1,1,1,1]T\n");
    v1.SetFill(1);
    m41.SetFill(1);
    Printf("vec1 dot mat41 = %f\n", v1.dot(m41));
    Println;
}

void Vec4ComparisonOperator()
{
    PrintHeading("Vector4 Comparison operators ");

    float v[4] = { 1,2,3,4 };
    CdtVector4<> v1(1, 2, 3, 4);
    CdtVector<4> v41(v, sizeof(v));
    CdtMatrix<4, 1> m41(v, sizeof(v));

    Printf("/* Operator: ==(CdtVector4 &v) and !=(CdtVector4 &v) */\n");
    if (v1 == v1) Printf("true");
    else Printf("false");
    Println;

    if (v1 != v1) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(CdtVector &v) and !=(CdtVector &v) */\n");
    if (v1 == v41) Printf("true");
    else Printf("false");
    Println;

    if (v1 != v41) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(CdtMatrix &v) and !=(CdtMatrix &v) */\n");
    if (v1 == m41) Printf("true");
    else Printf("false");
    Println;

    if (v1 != m41) Printf("false");
    else Printf("true");
    Println;
    Println;
}
