#include "testVector4.h"
#include "testPrint.h"

#include <dtMath/dtMath.h>

using dtMath::dtMatrix;
using dtMath::dtVector;
using dtMath::dtVector4;

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

    dtVector4<> vec;
    vec << 1, 2, 3, 4;

    vec.Print();
    Println;
}

void Vec4MemberFunc()
{
    PrintHeading("Vector4 Member Functions ");

    float v[4] = {1, 2, 3, 4};

    Printf("/* Class Create: () */\n");
    dtVector4<> v1;
    v1.Print();
    Println;

    Printf("/* Class Create: with array */\n");
    dtVector4<> v2(v, sizeof(v));
    Printf("array[4] = {1,2,3,4}\n");
    Printf("vec =\n");
    v2.Print();
    Println;

    Printf("/* Class Create: with element */\n");
    dtVector4<> v3(1, 2, 3, 4);
    Printf("elements are 1, 2, 3, 4\n");
    Printf("vec =\n");
    v3.Print();
    Println;

    Printf("/* Class Create: with dtVector4 */\n");
    dtVector4<> v4(v3);
    Printf("vec = [1,2,3,4]T\n");
    Printf("vec =\n");
    v4.Print();
    Println;

    Printf("/* Class Create: with dtVector */\n");
    dtVector<4> v41(v, sizeof(v));
    dtVector4<> v5(v41);
    Printf("vec = [1,2,3,4]T\n");
    Printf("vec =\n");
    v5.Print();
    Println;

    Printf("/* Class Create: with dtMatrix */\n");
    dtMatrix<4, 1> m41(v, sizeof(v));
    dtVector4<> v6(m41);
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

    Printf("/* Function: SetElement(dtVector4) */\n");
    Printf("vec = [1,2,3,4]T\n");
    Printf("vec =\n");
    v2.SetElement(v3);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(dtVector) */\n");
    Printf("vec = [1,2,3,4]T\n");
    Printf("vec =\n");
    v2.SetElement(v41);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(dtMatrix) */\n");
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
    dtVector4<> v1;

    PrintHeading("Vector4 Access operator ");
    Printf("/* Operator: () */\n");
    Printf("vec(i) = i; Printf("
           "%%f"
           ", vec(i));\n");

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

    float v[4] = {1, 2, 3, 4};
    dtVector4<> v1(1, 2, 3, 4);
    dtVector<4> v41(v, sizeof(v));
    dtMatrix<4, 1> m41(v, sizeof(v));

    Printf("vec =\n");
    v1.Print();
    Println;

    Printf("/* Operator: -() */\n");
    Printf("-vec =\n");
    (-v1).Print();
    Println;

    Printf("/* Operator: +(dtVector4) */\n");
    Printf("vec + (-vec) =\n");
    (v1 + (-v1)).Print();
    Println;

    Printf("/* Operator: -(dtVector4) */\n");
    Printf("vec - vec =\n");
    (v1 - v1).Print();
    Println;

    Printf("/* Operator: +(dtVector) */\n");
    Printf("vec = [1,2,3,4]T, vec41 = [1,2,3,4]T\n");
    Printf("vec + (-vec41) =\n");
    (v1 + (-v41)).Print();
    Println;

    Printf("/* Operator: -(dtVector) */\n");
    Printf("vec = [1,2,3,4]T, vec41 = [1,2,3,4]T\n");
    Printf("vec - vec41) =\n");
    (v1 - v41).Print();
    Println;

    Printf("/* Operator: +(dtMatrix) */\n");
    Printf("vec = [1,2,3,4]T, mat41 = [1,2,3,4]T\n");
    Printf("vec + (-mat41) =\n");
    (v1 + (-m41)).Print();
    Println;

    Printf("/* Operator: -(dtMatrix) */\n");
    Printf("vec = [1,2,3,4]T, mat41 = [1,2,3,4]T\n");
    Printf("vec - mat41 =\n");
    (v1 - m41).Print();
    Println;
}

void Vec4ArithmeticProduct()
{
    PrintHeading("Vector4 Arithmetic operators - Product ");

    float v[4] = {1, 2, 3, 4};
    dtVector4<> v1(1, 2, 3, 4);
    dtVector4<> v2(v, sizeof(v));
    dtVector<4> v41(v, sizeof(v));
    dtMatrix<4, 1> m41(v, sizeof(v));
    dtMatrix<1, 4> m14(v, sizeof(v));

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

    Printf("/* Operator: *(scalar, dtVector) */\n");
    Printf("vec = [1,1,1,1]T\n");
    Printf("0.1 * vec =\n");
    v1.SetFill(1);
    (0.1f * v1).Print();
    Println;

    Printf("/* Operator: *(dtMatrix<1, col>) */\n");
    Printf("vec = [1,1,1,1]T, mat14 = [1,1,1,1]\n");
    Printf("vec * mat14 =\n");
    v1.SetFill(1);
    m14.SetFill(1);
    (v1 * m14).Print();
    Println;

    Printf("/* Operator: dot(dtVector4) */\n");
    Printf("vec1 = [1,1,1,1]T, vec2 = [1,1,1,1]T\n");
    v1.SetFill(1);
    v2.SetFill(1);
    Printf("vec1 dot vec2 = %f\n", v1.dot(v2));
    Println;

    Printf("/* Operator: dot(dtVector) */\n");
    Printf("vec1 = [1,1,1,1]T, vec41 = [1,1,1,1]T\n");
    v1.SetFill(1);
    v41.SetFill(1);
    Printf("vec1 dot vec41 = %f\n", v1.dot(v41));
    Println;

    Printf("/* Operator: dot(dtMatrix) */\n");
    Printf("vec1 = [1,1,1,1]T, mat41 = [1,1,1,1]T\n");
    v1.SetFill(1);
    m41.SetFill(1);
    Printf("vec1 dot mat41 = %f\n", v1.dot(m41));
    Println;
}

void Vec4ComparisonOperator()
{
    PrintHeading("Vector4 Comparison operators ");

    float v[4] = {1, 2, 3, 4};
    dtVector4<> v1(1, 2, 3, 4);
    dtVector<4> v41(v, sizeof(v));
    dtMatrix<4, 1> m41(v, sizeof(v));

    Printf("/* Operator: ==(dtVector4 &v) and !=(dtVector4 &v) */\n");
    if (v1 == v1) Printf("true");
    else Printf("false");
    Println;

    if (v1 != v1) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(dtVector &v) and !=(dtVector &v) */\n");
    if (v1 == v41) Printf("true");
    else Printf("false");
    Println;

    if (v1 != v41) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(dtMatrix &v) and !=(dtMatrix &v) */\n");
    if (v1 == m41) Printf("true");
    else Printf("false");
    Println;

    if (v1 != m41) Printf("false");
    else Printf("true");
    Println;
    Println;
}
