#include "testPrint.h"
#include "testMatrix3.h"

#include <dtMath/dtMath.h>


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

    CdtMatrix3<> mat;
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
    float a[3 * 3] = { 1,2,3,4,5,6,7,8,9 };
    CdtMatrix3<> m1(a, sizeof(a));
    Printf("a[3*3] = {1,2,3,4,5,6,7,8,9}\n");
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Class Create: with diagonal */\n");
    Printf("d[3*3] = {1,2,3,4,5,6,7,8,9}\n");
    Printf("mat =\n");
    CdtMatrix3<> m2('d', a, sizeof(a));
    m2.Print();
    Println;

    Printf("/* Class Create: with elements */\n");
    Printf("elements are 1,2,3,4,5,6,7,8,9\n");
    Printf("mat =\n");
    CdtMatrix3<> m3(1, 2, 3, 4, 5, 6, 7, 8, 9);
    m3.Print();
    Println;

    Printf("/* Class Create: with Matrix3 */\n");
    Printf("CdtMatrix3 is\n");
    m3.Print();
    Printf("mat =\n");
    CdtMatrix3<> m4(m3);
    m4.Print();
    Println;

    Printf("/* Class Create: with Rotation matrix */\n");
    CdtRotation<> r;
    CdtMatrix3<> m5(r);
    Printf("rotation matrix is\n");
    r.Print();
    Printf("mat =\n");
    m5.Print();
    Println;

    Printf("/* Class Create: with Matrix */\n");
    CdtMatrix<3, 3> m33;
    CdtMatrix3<> m6(m33);
    Printf("CdtMatrix is\n");
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

    Printf("/* Function: SetElement(CdtMatrix3 &m) */\n");
    Printf("CdtMatrix3 is\n");
    m2.Print();
    m1.SetElement(m2);
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Function: SetElement(CdtRotation) */\n");
    Printf("CdtRotation is\n");
    r.Print();
    m1.SetElement(r);
    Printf("mat =\n");
    m1.Print();
    Println;

    Printf("/* Function: SetElement(CdtMatrix) */\n");
    Printf("CdtMatrix is\n");
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
    CdtMatrix3<> m1(1,2,3,4,5,6,7,8,9);
    
    Printf("/* Operator: () */\n");
    Printf("mat(i,j) = 100; Printf(""%%f"", mat(i,j));\n");
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
    float m[9] = { 1,2,3,4,5,6,7,8,9 };
    CdtMatrix3<> m1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    CdtRotation<> r1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    CdtMatrix<3, 3> m33(m, sizeof(m));
    Printf("mat =\n");
    m1.Print();

    Printf("/* Operator: -() */\n");
    Printf("-mat =\n");
    (-m1).Print();
    Println;

    Printf("/* Operator: +(CdtMatrix3) */\n");
    Printf("mat + (-mat) =\n");
    (m1 + (-m1)).Print();
    Println;

    Printf("/* Operator: -(CdtMatrix3) */\n");
    Printf("mat - mat =\n");
    (m1 - m1).Print();
    Println;

    Printf("/* Operator: +(CdtRotation) */\n");
    Printf("rot =\n");
    r1.Print();
    Printf("mat + (-rot) =\n");
    (m1 + (-r1)).Print();
    Println;

    Printf("/* Operator: -(CdtRotation) */\n");
    Printf("mat - rot =\n");
    (m1 - r1).Print();
    Println;

    Printf("/* Operator: +(CdtMatrix) */\n");
    Printf("CdtMatrix =\n");
    m33.Print();
    Printf("mat + (-mat33) =\n");
    (m1 + (-m33)).Print();
    Println;

    Printf("/* Operator: -(CdtMatrix) */\n");
    Printf("mat - mat33 =\n");
    (m1 - m33).Print();
    Println;
}

void Mat3ArithmeticProduct()
{
    PrintHeading("Matrix3 Arithmetic operators - Product ");
    CdtMatrix3<> m1;
    CdtMatrix3<> m2;
    CdtRotation<> r1;
    CdtMatrix<3, 6> m36;
    CdtVector<3> v31;
    CdtVector3<> v3;
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

    Printf("/* Operator: *(scalar, CdtMatrix3) */\n");
    Printf("0.1 * mat =\n");
    (0.1f * m1).Print();
    Println;

    Printf("/* Operator: *(CdtMatrix) */\n");
    m1.SetElement(-1, -2, 4, 4, -2, -1, -1, -2, 4);
    m36.SetFill(1);
    (m1 * m36).Print();
    Println;

    Printf("/* Operator: *(CdtMatrix3) */\n");
    m1.SetElement(-1, -2, 4, 4, -2, -1, -1, -2, 4);
    m2.SetFill(1);
    Printf("matA =\n");
    m1.Print();
    Printf("matB =\n");
    m2.Print();
    Printf("matA * matB =\n");
    (m1 * m2).Print();
    Println;

    Printf("/* Operator: *(CdtRotation) */\n");
    m1.SetElement(-1, -2, 4, 4, -2, -1, -1, -2, 4);
    r1.SetFill(1);
    Printf("matA =\n");
    m1.Print();
    Printf("rot =\n");
    r1.Print();
    Printf("matA * rot =\n");
    (m1 * r1).Print();
    Println;

    Printf("/* Operator: *(CdtVector) */\n");
    m1.SetElement(-1, -2, 4, 4, -2, -1, -1, -2, 4);
    v31.SetFill(1);
    Printf("mat =\n");
    m1.Print();
    Printf("vec =\n");
    v31.Print();
    Printf("mat * vec =\n");
    (m1 * v31).Print();
    Println;

    Printf("/* Operator: *(CdtVector3) */\n");
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
    float m[9] = { 1,2,3,4,5,6,7,8,9 };
    CdtMatrix3<> m1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    CdtRotation<> r1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    CdtMatrix<3, 3> m33(m, sizeof(m));

    Printf("/* Operator: ==(CdtMatrix3 &m) and !=(CdtMatrix3 &m) */\n");
    if (m1 == m1) Printf("true");
    else Printf("false");
    Println;

    if (m1 != m1) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(CdtRotation &m) and !=(CdtRotation &m) */\n");
    if (m1 == r1) Printf("true");
    else Printf("false");
    Println;

    if (m1 != r1) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(CdtMatrix &m) and !=(CdtMatrix &m) */\n");
    if (m1 == m33) Printf("true");
    else Printf("false");
    Println;

    if (m1 != m33) Printf("false");
    else Printf("true");
    Println;
    Println;
}