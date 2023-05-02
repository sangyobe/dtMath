#include "testRotation.h"
#include "testPrint.h"

#include <dtMath/dtMath.h>

using dtMath::dtMatrix;
using dtMath::dtMatrix3;
using dtMath::dtQuaternion;
using dtMath::dtRotation;
using dtMath::dtVector;
using dtMath::dtVector3;

void Test_Rotation()
{
    PrintTitle("Test Rotation matrix");
    RotInit();
    RotMemberFunc();
    RotMemberAccessOperator();
    RotArithmetic();
    RotComparisonOperator();
}

void RotInit()
{
    PrintHeading("Initialize Rot Matrix ");
    Printf("CommaInit, rot << 1, 2, 3, 4, 5, 6, 7, 8, 9;\n");
    Printf("rot = \n");

    dtRotation<> rot;
    rot << 1, 2, 3,
        4, 5, 6,
        7, 8, 9;

    rot.Print();
    Println;
}

void RotMemberFunc()
{
    PrintHeading("Rotation Member Functions ");

    Printf("/* Class Create: () */\n");
    dtRotation<> r1;
    r1.Print();
    Println;

    Printf("/* Class Create: with array */\n");
    Printf("array = {1,2,3,4,5,6,7,8,9}\n");
    float a[3 * 3] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    dtRotation<> r2(a, sizeof(a));
    r2.Print();
    Println;

    Printf("/* Class Create: with array, diagonal */\n");
    Printf("array = {1,2,3}\n");
    dtRotation<> r3('d', a, sizeof(a));
    r3.Print();
    Println;

    Printf("/* Class Create: with elements */\n");
    Printf("elements is 1,2,3,4,5,6,7,8,9\n");
    dtRotation<> r4(1, 2, 3, 4, 5, 6, 7, 8, 9);
    r4.Print();
    Println;

    Printf("/* Class Create: with angle */\n");
    Printf("angle is X-Axis 30deg\n");
    dtRotation<> r5(AXIS1(0), 30 * DEG2RADf);
    r5.Print();
    Println;

    Printf("/* Class Create: with angle1 and angle2 */\n");
    Printf("angle1 is Y-Axis 10deg, angle2 is X-Axis 20deg\n");
    dtRotation<> r6(AXIS2(1, 0), 10 * DEG2RADf, 20 * DEG2RADf);
    r6.Print();
    Printf("In terms of YXZ Euler angles again:\n");
    (r6.GetEulerAngles(AXIS3(1, 0, 2)) * RAD2DEGf).Print();
    Println;

    Printf("/* Class Create: with angle1, angle2 and angle3 */\n");
    Printf("angle1 is Z(10), angle2 is Y(20), angle3 is X(30)\n");
    dtRotation<> r7(AXIS3(2, 1, 0), 10 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf);
    r7.Print();
    Printf("In terms of ZYX Euler angles again:\n");
    (r7.GetEulerAngles(AXIS3(2, 1, 0)) * RAD2DEGf).Print();
    Println;

    Printf("/* Class Create: with dtRotation */\n");
    dtRotation<> r8(r7);
    r8.Print();
    Println;

    Printf("/* Class Create: with dtMatrix3 */\n");
    Printf("Matrix3 is [1,2,3;4,5,6;7,8,9]");
    dtMatrix3<> m3(1, 2, 3, 4, 5, 6, 7, 8, 9);
    dtRotation<> r9(m3);
    r9.Print();
    Println;

    Printf("/* Class Create: with dtMatrix */\n");
    Printf("Matrix is [1,2,3;4,5,6;7,8,9]");
    float m[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    dtMatrix<3, 3> m33(m, sizeof(m));
    dtRotation<> r10(m33);
    r10.Print();
    Println;

    Printf("/* Class Create: with Euler Angle Order and dtVector3 */\n");
    Printf("Vector is [10, 20, 30]T deg\n");
    dtVector3<> v3(10 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf);
    dtRotation<> r11(AXIS3(2, 1, 0), v3);
    r11.Print();
    Printf("In terms of ZYX Euler angles again:\n");
    (r11.GetEulerAngles(AXIS3(2, 1, 0)) * RAD2DEGf).Print();
    Println;

    Printf("/* Class Create: with Euler Angle Order and dtVector */\n");
    float v[3] = {10 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf};
    Printf("Vector is [10, 20, 30]T deg\n");
    dtVector<3> v31(v, sizeof(v));
    dtRotation<> r12(AXIS3(2, 1, 0), v31);
    r12.Print();
    Printf("In terms of ZYX Euler angles again:\n");
    (r12.GetEulerAngles(AXIS3(2, 1, 0)) * RAD2DEGf).Print();
    Println;

    Printf("/* Class Create: with Quaternion */\n");
    Printf("Quaternion is same ZYX Euler [10, 20, 30]T deg\n");
    dtQuaternion<> q(AXIS3(2, 1, 0), 10 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf);
    dtRotation<> r13(q);
    r13.Print();
    Printf("In terms of ZYX Euler angles again:\n");
    (r13.GetEulerAngles(AXIS3(2, 1, 0)) * RAD2DEGf).Print();
    Println;

    Printf("/* Function: SetZero() */\n");
    r1.SetZero();
    r1.Print();
    Println;

    Printf("/* Function: SetIdentity() */\n");
    r1.SetIdentity();
    r1.Print();
    Println;

    Printf("/* Function: SetDiagonal(d1, d2, d3) */\n");
    Printf("d1:1, d2:2, d3:3\n");
    r1.SetDiagonal(1, 2, 3);
    r1.Print();
    Println;

    Printf("/* Function: SetFill(value) */\n");
    Printf("value is 1\n");
    r1.SetFill(1);
    r1.Print();
    Println;

    Printf("/* Function: SetElement(*element) */\n");
    Printf("float a = {1,2,3,4,5,6,7,8,9};\n");
    r1.SetElement(a, sizeof(a));
    r1.Print();
    Println;

    Printf("/* Function: SetElement(m00, ... , r22) */\n");
    r1.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    r1.Print();
    Println;

    Printf("/* Function: SetElement(order, angle) */\n");
    Printf("order is X-axis, angle is 10deg\n");
    r1.SetElement(AXIS1(0), 10 * DEG2RADf);
    r1.Print();
    Printf("In terms of XYZ Euler angles again:\n");
    (r1.GetEulerAngles(AXIS3(0, 1, 2)) * RAD2DEGf).Print();
    Println;

    Printf("/* Function: SetElement(order, angle, angle) */\n");
    Printf("order is YX-axis, angle1 is 10, angle2 is 20deg\n");
    r1.SetElement(AXIS2(1, 0), 10 * DEG2RADf, 20 * DEG2RADf);
    r1.Print();
    Printf("In terms of YXZ Euler angles again:\n");
    (r1.GetEulerAngles(AXIS3(1, 0, 2)) * RAD2DEGf).Print();
    Println;

    Printf("/* Function: SetElement(order, angle, angle, angle) */\n");
    Printf("order is ZYX-axis, angle1 is 10, angle2 is 20deg, angle3 is 30deg\n");
    r1.SetElement(AXIS3(2, 1, 0), 10 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf);
    r1.Print();
    Printf("In terms of ZYX Euler angles again:\n");
    (r1.GetEulerAngles(AXIS3(2, 1, 0)) * RAD2DEGf).Print();
    Println;

    Printf("/* Function: SetElement(dtRotation) */\n");
    r1.SetElement(r2);
    r1.Print();
    Println;

    Printf("/* Function: SetElement(dtMatrix3) */\n");
    Printf("Matrix3 is {1,2,3;4,5,6;7,8,9}\n");
    r1.SetElement(m3);
    r1.Print();
    Println;

    Printf("/* Function: SetElement(dtMatrix) */\n");
    Printf("Matrix is {1,2,3;4,5,6;7,8,9}\n");
    r1.SetElement(m33);
    r1.Print();
    Println;

    Printf("/* Function: SetElement(order, dtVector3 euler) */\n");
    Printf("order is ZYX-axis, vec is [10,20,30]T deg\n");
    r1.SetElement(AXIS3(2, 1, 0), v3);
    r1.Print();
    Println;

    Printf("/* Function: SetElement(order, dtVector euler) */\n");
    Printf("order is ZYX-axis, vec is [10,20,30]T deg\n");
    r1.SetElement(AXIS3(2, 1, 0), v31);
    r1.Print();
    Println;

    Printf("/* Function: SetElement(dtQuaternion) */\n");
    Printf("Quaternion is in the same orientation with ZYX Euler [10,20,30]T deg\n");
    r1.SetElement(q);
    r1.Print();
    Println;

    Printf("/* Function: SetSwapRowVec */\n");
    r1.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    Printf("R =\n");
    r1.Print();
    r1.SetSwapRowVec(0, 2);
    Printf("SetSwapRowVec(0,2);\n");
    Printf("R =\n");
    r1.Print();
    Println;

    Printf("/* Function: SetSwapColVec */\n");
    r1.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    Printf("R =\n");
    r1.Print();
    r1.SetSwapColVec(0, 2);
    Printf("SetSwapColVec(0,2);\n");
    Printf("R =\n");
    r1.Print();
    Println;

    Printf("/* Function: GetRowVec */\n");
    r1.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    Printf("R =\n");
    r1.Print();
    Printf("1st Row Vector is\n");
    r1.GetRowVec(0).Print();
    Println;

    Printf("/* Function: GetColVec */\n");
    r1.SetElement(1, 2, 3, 4, 5, 6, 7, 8, 9);
    Printf("R =\n");
    r1.Print();
    Printf("3rd Col Vector is\n");
    r1.GetColVec(2).Print();
    Println;

    Printf("/* Function: Transpose */\n");
    Printf("R =\n");
    r1.Print();
    Printf("Transposed R =\n");
    r1.Transpose().Print();
    Println;

    Printf("/* Function: Inverse */\n");
    Printf("R =\n");
    r1.Print();
    Printf("Inverse of R =\n");
    r1.Inv().Print();
    Println;
}

void RotMemberAccessOperator()
{
    PrintHeading("Rotation Member Access Operator ");
    dtRotation<> r1(1, 2, 3, 4, 5, 6, 7, 8, 9);

    Printf("/* Operator: () */\n");
    Printf("R(i,j) = 100; & Printf("
           "%%f"
           ", R(i,j));\n");
    for (uint16_t i = 0; i < 3; i++)
    {
        for (uint16_t j = 0; j < 3; j++)
        {
            r1(i, j) = 100;
            Printf("%7.3f ", r1(i, j));
        }
        Printf("\n");
    }
    Println;
}

void RotArithmetic()
{
    PrintHeading("Rotation Arithmetic operators ");
    float m[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    dtRotation<> r1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    dtRotation<> r2;
    dtMatrix3<> m3(1, 2, 3, 4, 5, 6, 7, 8, 9);
    dtMatrix<3, 3> m33(m, sizeof(m));
    dtVector<3> v31;
    dtVector3<> v3;

    Printf("R = \n");
    r1.Print();
    Println;

    Printf("/* Operator: -() */\n");
    Printf("-R =\n");
    (-r1).Print();
    Println;

    Printf("/* Operator: +(const dtRotation) */\n");
    Printf("R + (-R) =\n");
    (r1 + (-r1)).Print();
    Println;

    Printf("/* Operator: -(const dtRotation) */\n");
    Printf("R - R =\n");
    (r1 - r1).Print();
    Println;

    Printf("/* Operator: +(const dtMatrix3) */\n");
    Printf("mat3 = \n");
    m3.Print();
    Printf("R + (-mat3) =\n");
    (r1 + (-m3)).Print();
    Println;

    Printf("/* Operator: -(const dtMatrix3) */\n");
    Printf("mat3 = \n");
    m3.Print();
    Printf("R - mat3 =\n");
    (r1 - m3).Print();
    Println;

    Printf("/* Operator: +(const dtMatrix) */\n");
    Printf("mat33 = \n");
    m33.Print();
    Printf("R + (-mat33) =\n");
    (r1 + (-m33)).Print();
    Println;

    Printf("/* Operator: -(const dtMatrix) */\n");
    Printf("mat33 = \n");
    m33.Print();
    Printf("R - mat33 =\n");
    (r1 - m33).Print();
    Println;

    Printf("/* Operator: *(scalar) */\n");
    r1.SetFill(1);
    Printf("R =\n");
    r1.Print();
    Printf("R * 0.1 =\n");
    (r1 * 0.1f).Print();
    Println;

    Printf("/* Operator: /(scalar) */\n");
    r1.SetFill(1);
    Printf("R =\n");
    r1.Print();
    Printf("R / 10 =\n");
    (r1 / 10).Print();
    Println;

    Printf("/* Operator: *(scalar, dtRotation) */\n");
    r1.SetFill(1);
    Printf("R =\n");
    r1.Print();
    Printf("0.1 * R =\n");
    (0.1f * r1).Print();
    Println;

    Printf("/* Operator: *(dtMatrix) */\n");
    m33.SetFill(1);
    Printf("R =\n");
    r1.Print();
    Printf("mat33 =\n");
    m33.Print();
    Printf("R * mat33 = \n");
    (r1 * m33).Print();
    Println;

    Printf("/* Operator: *(dtMatrix3) */\n");
    m3.SetFill(1);
    Printf("R =\n");
    r1.Print();
    Printf("mat3 =\n");
    m3.Print();
    Printf("R * mat3 = \n");
    (r1 * m3).Print();
    Println;

    Printf("/* Operator: *(dtRotation) */\n");
    r2.SetFill(1);
    Printf("R =\n");
    r1.Print();
    Printf("R2 =\n");
    r2.Print();
    Printf("R * R2 = \n");
    (r1 * r2).Print();
    Println;

    Printf("/* Operator: *(dtVector) */\n");
    v31.SetFill(1);
    Printf("R =\n");
    r1.Print();
    Printf("vec31 =\n");
    v31.Print();
    Printf("R * vec31 = \n");
    (r1 * v31).Print();
    Println;

    Printf("/* Operator: *(dtVector3) */\n");
    v3.SetFill(1);
    Printf("R =\n");
    r1.Print();
    Printf("vec3 =\n");
    v3.Print();
    Printf("R * vec3 = \n");
    (r1 * v3).Print();
    Println;
}

void RotComparisonOperator()
{
    PrintHeading("Rotation Comparison operators ");
    float m[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    dtRotation<> r1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    dtMatrix3<> m3(1, 2, 3, 4, 5, 6, 7, 8, 9);
    dtMatrix<3, 3> m33(m, sizeof(m));

    Printf("/* Operator: ==(dtRotation &m) and !=(dtRotation &m) */\n");
    if (r1 == r1) Printf("true");
    else Printf("false");
    Println;

    if (r1 != r1) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(dtMatrix3 &m) and !=(dtMatrix3 &m) */\n");
    if (r1 == m3) Printf("true");
    else Printf("false");
    Println;

    if (r1 != m3) Printf("false");
    else Printf("true");
    Println;

    Printf("/* Operator: ==(dtMatrix &m) and !=(dtMatrix &m) */\n");
    if (r1 == m33) Printf("true");
    else Printf("false");
    Println;

    if (r1 != m33) Printf("false");
    else Printf("true");
    Println;
    Println;
}