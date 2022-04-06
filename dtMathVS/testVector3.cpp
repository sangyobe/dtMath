#include "testPrint.h"
#include "testVector3.h"
#include "./dtMath/dtMath.h"

void Test_Vector3()
{
    PrintTitle("Test Vector3");
    Vec3Init();
    Vec3MemberFunc();
    Vec3MemberAccessOperator();
    Vec3ArithmeticSum();
    Vec3ArithmeticProduct();
    Vec3ComparisonOperator();
}

void Vec3Init()
{
    PrintHeading("Initialize Vector3 ");
    Printf("CommaInit, vec << 1, 2, 3;\n");
    Printf("vec = \n");

    CdtVector3<> vec;
    vec << 1, 2, 3;

    vec.Print();
    Println;
}

void Vec3MemberFunc()
{
    PrintHeading("Vector3 Member Functions ");

    float v[3] = { 1,2,3 };

    Printf("/* Class Create: () */\n");
    CdtVector3<> v1;
    v1.Print();
    Println;

    Printf("/* Class Create: with array */\n");
    CdtVector3<> v2(v, sizeof(v));
    Printf("array[3] = {1,2,3}\n");
    Printf("vec =\n");
    v2.Print();
    Println;

    Printf("/* Class Create: with element */\n");
    CdtVector3<> v3(1, 2, 3);
    Printf("elements are 1, 2, 3\n");
    Printf("vec =\n");
    v3.Print();
    Println;

    Printf("/* Class Create: with CdtVector3 */\n");
    CdtVector3<> v4(v3);
    Printf("vec = [1,2,3]T\n");
    Printf("vec =\n");
    v4.Print();
    Println;

    Printf("/* Class Create: with CdtVector */\n");
    CdtVector<3> v31(v, sizeof(v));
    CdtVector3<> v5(v31);
    Printf("vec = [1,2,3]T\n");
    Printf("vec =\n");
    v5.Print();
    Println;

    Printf("/* Class Create: with CdtMatrix */\n");
    CdtMatrix<3, 1> m31(v, sizeof(v));
    CdtVector3<> v6(m31);
    Printf("mat31 = [1,2,3]T\n");
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
    Printf("element[3] = {1,2,3}\n");
    Printf("vec =\n");
    v2.SetElement(v, sizeof(v));
    v2.Print();
    Println;

    Printf("/* Function: SetElement(i, j, k) */\n");
    Printf("i=1, j=2, k=3\n");
    Printf("vec =\n");
    v2.SetElement(1, 1, 1);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(CdtVector3) */\n");
    Printf("vec = [1,2,3]T\n");
    Printf("vec =\n");
    v2.SetElement(v3);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(CdtVector) */\n");
    Printf("vec = [1,2,3]T\n");
    Printf("vec =\n");
    v2.SetElement(v31);
    v2.Print();
    Println;

    Printf("/* Function: SetElement(CdtMatrix) */\n");
    Printf("mat31 = [1,2,3]T\n");
    Printf("vec =\n");
    v2.SetElement(m31);
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

    Printf("/* Function: GetSkew() */\n");
    v2.SetElement(1, 2, 3);
    Printf("vec =\n");
    v2.Print();
    Printf("Skew-symmetric matrix is\n");
    v2.GetSkew().Print();
    Println;

    Printf("/* Function: Transpose() */\n");
    v2.SetElement(1, 2, 3);
    Printf("vec =\n");
    v2.Print();
    Printf("Transposed vec =\n");
    v2.Transpose().Print();
    Println;
}

void Vec3MemberAccessOperator()
{
    CdtVector3<> v1;

    PrintHeading("Vector3 Access operator ");
    Printf("/* Operator: () */\n");
    Printf("vec(i) = i; Printf(""%%f"", vec(i));\n");

    for (uint16_t i = 0; i < 3; i++)
    {
        v1(i) = i;
        Printf("%7.3f ", v1(i));
    }

    Println;
    Println;
}

void Vec3ArithmeticSum()
{
    PrintHeading("Vector3 Arithmetic operators - Sum ");

    float v[3] = { 1,2,3 };
    CdtVector3<> v1(1, 2, 3);
    CdtVector<3> v31(v, sizeof(v));
    CdtMatrix<3, 1> m31(v, sizeof(v));

    Printf("vec =\n");
    v1.Print();
    Println;

    Printf("/* Operator: -() */\n");
    Printf("-vec =\n");
    (-v1).Print();
    Println;

    Printf("/* Operator: +(CdtVector3) */\n");
    Printf("vec + (-vec) =\n");
    (v1 + (-v1)).Print();
    Println;

    Printf("/* Operator: -(CdtVector3) */\n");
    Printf("vec - vec =\n");
    (v1 - v1).Print();
    Println;

    Printf("/* Operator: +(CdtVector) */\n");
    Printf("vec = [1,2,3]T, vec31 = [1,2,3]T\n");
    Printf("vec + (-vec31) =\n");
    (v1 + (-v31)).Print();
    Println;

    Printf("/* Operator: -(CdtVector) */\n");
    Printf("vec = [1,2,3]T, vec31 = [1,2,3]T\n");
    Printf("vec - vec31 =\n");
    (v1 - v31).Print();
    Println;

    Printf("/* Operator: +(CdtMatrix) */\n");
    Printf("vec = [1,2,3]T, mat31 = [1,2,3]T\n");
    Printf("vec + (-mat31) =\n");
    (v1 + (-m31)).Print();
    Println;

    Printf("/* Operator: -(CdtMatrix) */\n");
    Printf("vec = [1,2,3]T, mat31 = [1,2,3]T\n");
    Printf("vec - mat31 =\n");
    (v1 - m31).Print();
    Println;
}

void Vec3ArithmeticProduct()
{
    PrintHeading("Vector3 Arithmetic operators - Product ");

    float v[3] = { 0,1,0 };
    CdtVector3<> v1(1, 0, 0);
    CdtVector3<> v2(v, sizeof(v));
    CdtVector<3> v31(v, sizeof(v));
    CdtMatrix<3, 1> m31(v, sizeof(v));
    CdtMatrix<1, 3> m13(v, sizeof(v));

    Printf("/* Operator: *(scalar) */\n");
    Printf("vec = [1,1,1]T\n");
    Printf("vec * 0.1 =\n");
    v1.SetFill(1);
    (v1 * 0.1f).Print();
    Println;

    Printf("/* Operator: /(scalar) */\n");
    Printf("vec = [1,1,1]T\n");
    Printf("vec / 10 =\n");
    v1.SetFill(1);
    (v1 / 10).Print();
    Println;

    Printf("/* Operator: *(scalar, CdtVector) */\n");
    Printf("vec = [1,1,1]T\n");
    Printf("0.1 * vec =\n");
    v1.SetFill(1);
    (0.1f * v1).Print();
    Println;

    Printf("/* Operator: &(CdtVector3) */\n");
    Printf("vec1 = [1,0,0]T, vec2 = [0,1,0]T\n");
    Printf("cross product: vec1 & vec2 =\n");
    v1.SetElement(1, 0, 0);
    (v1 & v2).Print();
    Println;

    Printf("/* Operator: &(CdtVector) */\n");
    Printf("vec1 = [1,0,0]T, vec31 = [0,1,0]T\n");
    Printf("cross product: vec1 & vec31 =\n");
    (v1 & v31).Print();
    Println;

    Printf("/* Operator: &(CdtMatrix) */\n");
    Printf("vec1 = [1,0,0]T, mat31 = [0,1,0]T\n");
    Printf("cross product: vec1 & mat31 =\n");
    (v1 & m31).Print();
    Println;

    Printf("/* Operator: *(CdtMatrix<1, col>) */\n");
    Printf("vec = [1,1,1]T, mat13 = [1,1,1]\n");
    Printf("vec * mat13 =\n");
    v1.SetFill(1);
    m13.SetFill(1);
    (v1 * m13).Print();
    Println;

    Printf("/* Operator: dot(CdtVector3) */\n");
    Printf("vec1 = [1,1,1]T, vec2 = [1,1,1]T\n");
    v1.SetFill(1);
    v2.SetFill(1);
    Printf("vec1 dot vec2 = %f\n", v1.dot(v2));
    Println;

    Printf("/* Operator: dot(CdtVector) */\n");
    Printf("vec1 = [1,1,1]T, vec31 = [1,1,1]T\n");
    v1.SetFill(1);
    v31.SetFill(1);
    Printf("vec1 dot vec31 = %f\n", v1.dot(v31));
    Println;

    Printf("/* Operator: dot(CdtMatrix) */\n");
    Printf("vec1 = [1,1,1]T, mat31 = [1,1,1]T\n");
    v1.SetFill(1);
    m31.SetFill(1);
    Printf("vec1 dot mat31 = %f\n", v1.dot(m31));
    Println;
}

void Vec3ComparisonOperator()
{
    PrintHeading("Vector3 Comparison operators ");

    float v[3] = { 1,2,3 };
    CdtVector3<> v1(1, 2, 3);
    CdtVector<3> v31(v, sizeof(v));
    CdtMatrix<3, 1> m31(v, sizeof(v));

    Printf("/* Operator: ==(CdtVector3 &v) and !=(CdtVector3 &v) */\n");
    if (v1 == v1) Printf("true");
    else Printf("false");
    Println;

    if (v1 != v1) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(CdtVector &v) and !=(CdtVector &v) */\n");
    if (v1 == v31) Printf("true");
    else Printf("false");
    Println;

    if (v1 != v31) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(CdtMatrix &v) and !=(CdtMatrix &v) */\n");
    if (v1 == m31) Printf("true");
    else Printf("false");
    Println;

    if (v1 != m31) Printf("false");
    else Printf("true");
    Println;
    Println;
}