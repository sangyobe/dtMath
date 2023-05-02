#include "testMatrix3.h"
#include "testPrint.h"

#include <dtMath/dtMath.h>

using dtMath::dtMatrix;
using dtMath::dtMatrix3;
using dtMath::dtRotation;
using dtMath::dtVector;
using dtMath::dtVector3;

void Test_Matrix3()
{
    PrintTitle("Test Matrix3");
    Mat3Init();
    Mat3MemberFunc();
    Mat3MemberAccessOperator();
    Mat3ArithmeticSum();
    Mat3ArithmeticProduct();
    Mat3ComparisonOperator();
}

void Mat3Init()
{
    PrintHeading("Initialize Matrix3 ");
    Printf("CommaInit, mat << 1, 2, 3, 4, 5, 6, 7, 8, 9;\n");
    Printf("mat = \n");

    dtMatrix3<> mat;
    mat << 1, 2, 3,
        4, 5, 6,
        7, 8, 9;

    mat.Print();
    Println;
}

void Mat3MemberFunc()
{
    PrintHeading("Matrix3 Member Functions ");

    Printf("/* Class Create: with array */\n");
    float a[3 * 3] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    dtMatrix3<> m1(a, sizeof(a));
    Printf("a[3*3] = {1,2,3,4,5,6,7,8,9}\n");
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Class Create: with diagonal */\n");
    Printf("d[3*3] = {1,2,3,4,5,6,7,8,9}\n");
    Printf("mat =\n");
    dtMatrix3<> m2('d', a, sizeof(a));
    m2.Print();
    Println;

    Printf("/* Class Create: with elements */\n");
    Printf("elements are 1,2,3,4,5,6,7,8,9\n");
    Printf("mat =\n");
    dtMatrix3<> m3(1, 2, 3, 4, 5, 6, 7, 8, 9);
    m3.Print();
    Println;

    Printf("/* Class Create: with Matrix3 */\n");
    Printf("dtMatrix3 is\n");
    m3.Print();
    Printf("mat =\n");
    dtMatrix3<> m4(m3);
    m4.Print();
    Println;

    Printf("/* Class Create: with Rotation matrix */\n");
    dtRotation<> r;
    dtMatrix3<> m5(r);
    Printf("rotation matrix is\n");
    r.Print();
    Printf("mat =\n");
    m5.Print();
    Println;

    Printf("/* Class Create: with Matrix */\n");
    dtMatrix<3, 3> m33;
    dtMatrix3<> m6(m33);
    Printf("dtMatrix is\n");
    m33.Print();
    Printf("mat =\n");
    m6.Print();
    Println;

    Printf("/* Function: SetZero() */\n");
    Printf("mat =\n");
    m1.SetZero();
    m1.Print();
    Println;

    Printf("/* Function: SetIdentity() */\n");
    Printf("mat =\n");
    m1.SetIdentity();
    m1.Print();
    Println;

    Printf("/* Function: SetDiagonal(d1, d2, d3) */\n");
    Printf("diagonal is 1,2,3\n");
    m1.SetDiagonal(1, 2, 3);
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Function: SetFill(value) */\n");
    Printf("Value is 1\n");
    m1.SetFill(1);
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Function: SetElement(*element) */\n");
    Printf("array[9] = {1,2,3,4,5,6,7,8,9}\n");
    m1.SetElement(a, sizeof(a));
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Function: SetElement(m00, ... , m22) */\n");
    Printf("elements are 1,2,3,4,5,6,7,8,9\n");
    m1.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Function: SetElement(dtMatrix3 &m) */\n");
    Printf("dtMatrix3 is\n");
    m2.Print();
    m1.SetElement(m2);
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Function: SetElement(dtRotation) */\n");
    Printf("dtRotation is\n");
    r.Print();
    m1.SetElement(r);
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Function: SetElement(dtMatrix) */\n");
    Printf("dtMatrix is\n");
    m1.SetElement(m33);
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Function: SetSwapRowVec(row1, row2) */\n");
    m1.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    Printf("mat =\n");
    m1.Print();
    m1.SetSwapRowVec(1, 2);
    Printf("Swap Row vector1, and 2\n");
    Printf("mat = \n");
    m1.Print();
    Println;

    Printf("/* Function: SetSwapColVec(col1, col2) */\n");
    m1.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    Printf("mat =\n");
    m1.Print();
    m1.SetSwapColVec(1, 2);
    Printf("Swap Col vector1, and 2\n");
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Function: GetRowVec(row) */\n");
    m1.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    Printf("mat =\n");
    m1.Print();
    Printf("Row-0 Vector is\n");
    m1.GetRowVec(0).Print();
    Println;

    Printf("/* Function: GetColVec(col) */\n");
    m1.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    Printf("mat =\n");
    m1.Print();
    Printf("Col-0 Vector is\n");
    m1.GetColVec(0).Print();
    Println;

    Printf("/* Function: Transpose() */\n");
    m1.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    Printf("mat =\n");
    m1.Print();
    Printf("Transposed mat =\n");
    m1.Transpose().Print();
    Println;
}

void Mat3MemberAccessOperator()
{
    PrintHeading("Matrix3 Member Access Operator ");
    dtMatrix3<> m1(1, 2, 3, 4, 5, 6, 7, 8, 9);

    Printf("/* Operator: () */\n");
    Printf("mat(i,j) = 100; Printf("
           "%%f"
           ", mat(i,j));\n");
    for (uint16_t i = 0; i < 3; i++)
    {
        for (uint16_t j = 0; j < 3; j++)
        {
            m1(i, j) = 100;
            Printf("%7.3f ", m1(i, j));
        }
        Printf("\n");
    }
    Println;
}

void Mat3ArithmeticSum()
{
    PrintHeading("Matrix3 Arithmetic operators - Sum ");
    float m[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    dtMatrix3<> m1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    dtRotation<> r1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    dtMatrix<3, 3> m33(m, sizeof(m));
    Printf("mat =\n");
    m1.Print();

    Printf("/* Operator: -() */\n");
    Printf("-mat =\n");
    (-m1).Print();
    Println;

    Printf("/* Operator: +(dtMatrix3) */\n");
    Printf("mat + (-mat) =\n");
    (m1 + (-m1)).Print();
    Println;

    Printf("/* Operator: -(dtMatrix3) */\n");
    Printf("mat - mat =\n");
    (m1 - m1).Print();
    Println;

    Printf("/* Operator: +(dtRotation) */\n");
    Printf("rot =\n");
    r1.Print();
    Printf("mat + (-rot) =\n");
    (m1 + (-r1)).Print();
    Println;

    Printf("/* Operator: -(dtRotation) */\n");
    Printf("mat - rot =\n");
    (m1 - r1).Print();
    Println;

    Printf("/* Operator: +(dtMatrix) */\n");
    Printf("dtMatrix =\n");
    m33.Print();
    Printf("mat + (-mat33) =\n");
    (m1 + (-m33)).Print();
    Println;

    Printf("/* Operator: -(dtMatrix) */\n");
    Printf("mat - mat33 =\n");
    (m1 - m33).Print();
    Println;
}

void Mat3ArithmeticProduct()
{
    PrintHeading("Matrix3 Arithmetic operators - Product ");
    dtMatrix3<> m1;
    dtMatrix3<> m2;
    dtRotation<> r1;
    dtMatrix<3, 6> m36;
    dtVector<3> v31;
    dtVector3<> v3;
    Printf("mat =\n");
    m1.SetFill(10);
    m1.Print();

    Printf("/* Operator: *(scalar) */\n");
    Printf("mat * 0.1 =\n");
    (m1 * 0.1f).Print();
    Println;

    Printf("/* Operator: /(scalar) */\n");
    Printf("mat / 10 =\n");
    (m1 / 10.0f).Print();
    Println;

    Printf("/* Operator: *(scalar, dtMatrix3) */\n");
    Printf("0.1 * mat =\n");
    (0.1f * m1).Print();
    Println;

    Printf("/* Operator: *(dtMatrix) */\n");
    m1.SetElement(-1, -2, 4, 4, -2, -1, -1, -2, 4);
    m36.SetFill(1);
    (m1 * m36).Print();
    Println;

    Printf("/* Operator: *(dtMatrix3) */\n");
    m1.SetElement(-1, -2, 4, 4, -2, -1, -1, -2, 4);
    m2.SetFill(1);
    Printf("matA =\n");
    m1.Print();
    Printf("matB =\n");
    m2.Print();
    Printf("matA * matB =\n");
    (m1 * m2).Print();
    Println;

    Printf("/* Operator: *(dtRotation) */\n");
    m1.SetElement(-1, -2, 4, 4, -2, -1, -1, -2, 4);
    r1.SetFill(1);
    Printf("matA =\n");
    m1.Print();
    Printf("rot =\n");
    r1.Print();
    Printf("matA * rot =\n");
    (m1 * r1).Print();
    Println;

    Printf("/* Operator: *(dtVector) */\n");
    m1.SetElement(-1, -2, 4, 4, -2, -1, -1, -2, 4);
    v31.SetFill(1);
    Printf("mat =\n");
    m1.Print();
    Printf("vec =\n");
    v31.Print();
    Printf("mat * vec =\n");
    (m1 * v31).Print();
    Println;

    Printf("/* Operator: *(dtVector3) */\n");
    m1.SetElement(-1, -2, 4, 4, -2, -1, -1, -2, 4);
    v3.SetFill(1);
    Printf("mat =\n");
    m1.Print();
    Printf("vec3 =\n");
    v3.Print();
    Printf("mat * vec3 =\n");
    (m1 * v3).Print();
    Println;
}

void Mat3ComparisonOperator()
{
    PrintHeading("Matrix3 Comparison operators ");
    float m[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    dtMatrix3<> m1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    dtRotation<> r1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    dtMatrix<3, 3> m33(m, sizeof(m));

    Printf("/* Operator: ==(dtMatrix3 &m) and !=(dtMatrix3 &m) */\n");
    if (m1 == m1) Printf("true");
    else Printf("false");
    Println;

    if (m1 != m1) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(dtRotation &m) and !=(dtRotation &m) */\n");
    if (m1 == r1) Printf("true");
    else Printf("false");
    Println;

    if (m1 != r1) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(dtMatrix &m) and !=(dtMatrix &m) */\n");
    if (m1 == m33) Printf("true");
    else Printf("false");
    Println;

    if (m1 != m33) Printf("false");
    else Printf("true");
    Println;
    Println;
}