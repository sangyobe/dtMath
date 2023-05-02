#include "testQuaternion.h"
#include "testPrint.h"

#include <dtMath/dtMath.h>

using dtMath::dtQuaternion;
using dtMath::dtRotation;
using dtMath::dtVector;
using dtMath::dtVector3;

void Test_Quaternion()
{
    PrintTitle("Test Quaternion");
    QuatInit();
    QuatMemberFunc();
    QuatMemberAccessOperator();
    QuatArithmetic();
    QuatComparisonOperator();
}

void QuatInit()
{
    PrintHeading("Initialize Quaternion ");
    Printf("CommaInit, quat << 1, 2, 3, 4;\n");
    Printf("quat = \n");

    dtQuaternion<> quat;
    quat << 1, 2, 3, 4;

    quat.Print();
    Println;
}

void QuatMemberFunc()
{
    PrintHeading("Quaternion Member Functions ");

    float q[4] = {1, 2, 3, 4};

    Printf("/* Class Create */\n");
    dtQuaternion<> q1;
    q1.Print();
    Println;

    Printf("/* Class Create: with array */\n");
    Printf("arrary is q[4] = { 1, 2, 3, 4 }\n");
    dtQuaternion<> q2(q);
    q2.Print();
    Println;

    Printf("/* Class Create: with element arguments */\n");
    Printf("dtQuaternion<> q(1, 2, 3, 4)\n");
    dtQuaternion<> q3(1.0f, 2.0f, 3.0f, 4.0f);
    q3.Print();
    Println;

    // QuatX
    Printf("/* Class Create: with 1-axis angle */\n");
    Printf("x-axis, 10[deg]\n");
    dtQuaternion<> q4(AXIS1(0), 10.0f * DEG2RADf);
    q4.Print();
    Println;

    // QuatY * QuatX
    Printf("/* Class Create: with 2-axis angles */\n");
    Printf("Qy(10deg) * Qx(20deg)\n");
    dtQuaternion<> q5(AXIS2(1, 0), 10.0f * DEG2RADf, 20 * DEG2RADf);
    q5.Print();
    Printf("=> ZYX Euler Angle:\n");
    (RAD2DEGf * q5.GetEulerAngles(AXIS3(2, 1, 0))).Print();
    Println;

    // QuatZ * QuatY * QuatX
    Printf("/* Class Create: with 3-axis angle */\n");
    Printf("Qz(10deg) * Qy(20deg) * Qx(30deg)\n");
    dtQuaternion<> q6(AXIS3(2, 1, 0), 10.0f * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf);
    q6.Print();
    Printf("=> ZYX Euler Angle:\n");
    (RAD2DEGf * q6.GetEulerAngles(AXIS3(2, 1, 0))).Print();
    Println;

    Printf("/* Class Create: with dtQuaternion */\n");
    Printf("dtQuaternion<> q(q)\n");
    dtQuaternion<> q7(q6);
    q7.Print();
    Println;

    // euler ZYX
    Printf("/* Class Create: with euler angle (ZYX 10, 20, 30 deg) */\n");
    float zyx[3] = {10 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf};
    dtQuaternion<> q8(AXIS3(2, 1, 0), dtVector3<>(zyx, sizeof(zyx)));
    q8.Print();
    //(RAD2DEGf * q8.GetEulerAngles(AXIS3(2, 1, 0))).Print();
    Println;

    Printf("/* Class Create: with rotation matrix */\n");
    Printf("Rotation matrix is equal to ZYX Euler(10, 20, 30 deg)\n");
    dtRotation<> rot(AXIS3(2, 1, 0), 10.0f * DEG2RADf, 20.0f * DEG2RADf, 30.0f * DEG2RADf);
    dtQuaternion<> q9(rot);
    q9.Print();
    Printf("=> ZYX Euler Angle:\n");
    (RAD2DEGf * q9.GetEulerAngles(AXIS3(2, 1, 0))).Print();
    Println;

    Printf("/* Function: SetZero() */\n");
    Printf("q.SetZero()\n");
    q2.SetZero();
    q2.Print();
    Println;

    Printf("/* Function: SetFill() */\n");
    Printf("q.SetFill(1)\n");
    q2.SetFill(1);
    q2.Print();
    Println;

    Printf("/* Function: SetElement(*element) */\n");
    Printf("q.SetElement(array)");
    q2.SetElement(q);
    q2.Print();
    Println;

    Printf("/* Function: SetElement(w(1), x(2), y(3), z(4)) */\n");
    Printf("q.SetElement(1, 2, 3, 4)\n");
    q2.SetElement(1.0f, 2.0f, 3.0f, 4.0f);
    q2.Print();
    Println;

    Printf("/* Function: SetElement(order, ang) */\n");
    Printf("q.SetElement(AXIS1(0), 10 * DEG2RAD\n");
    q2.SetElement(AXIS1(0), 10 * DEG2RADf);
    q2.Print();
    Printf("=> XYZ Euler Angle:\n");
    (RAD2DEGf * q2.GetEulerAngles(AXIS3(0, 1, 2))).Print();
    Println;

    Printf("/* Function: SetElement(order, ang1, ang2) */\n");
    Printf("q.SetElement(AXIS2(1, 0), 10 * DEG2RAD, 20 * DEG2RAD\n");
    q2.SetElement(AXIS2(1, 0), 10 * DEG2RADf, 20 * DEG2RADf);
    q2.Print();
    Printf("=> YXZ Euler Angle:\n");
    (RAD2DEGf * q2.GetEulerAngles(AXIS3(1, 0, 2))).Print();
    Println;

    Printf("/* Function: SetElement(order, ang1, ang2, ang3) */\n");
    Printf("q.SetElement(AXIS2(2, 1, 0), 10 * DEG2RAD, 20 * DEG2RAD, 30 * DEG2RAD\n");
    q2.SetElement(AXIS3(2, 1, 0), 10 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf);
    q2.Print();
    Printf("=> ZYX Euler Angle:\n");
    (RAD2DEGf * q2.GetEulerAngles(AXIS3(2, 1, 0))).Print();
    Println;

    Printf("/* Function: SetElement(dtQuaternion) */\n");
    q2.SetElement(q3);
    q2.Print();
    Println;

    // euler ZYX
    Printf("/* Function: SetElement(order, dtVector3) */\n");
    q2.SetElement(AXIS3(2, 1, 0), dtVector3<>(zyx, sizeof(zyx)));
    q2.Print();
    Printf("=> ZYX Euler Angle:\n");
    (RAD2DEGf * q2.GetEulerAngles(AXIS3(2, 1, 0))).Print();
    Println;

    Printf("/* Function: SetElement(order, dtVector) */\n");
    q2.SetElement(AXIS3(2, 1, 0), dtVector<3>(zyx, sizeof(zyx)));
    q2.Print();
    Printf("=> ZYX Euler Angle:\n");
    (RAD2DEGf * q2.GetEulerAngles(AXIS3(2, 1, 0))).Print();
    Println;

    Printf("/* Function: SetElement(dtRotation) */\n");
    q2.SetElement(rot);
    q2.Print();
    Printf("=> ZYX Euler Angle:\n");
    (RAD2DEGf * q2.GetEulerAngles(AXIS3(2, 1, 0))).Print();
    Println;

    Printf("/* Function: SetSwap(i, j) of q(1, 2, 3, 4) and i=1, j=3 */\n");
    q2.SetElement(1.0f, 2.0f, 3.0f, 4.0f);
    q2.SetSwap(1, 3);
    q2.Print();
    Println;

    Printf("/* Function: SetNormalize() of q(1,1,1,1) */\n");
    q3 << 1, 1, 1, 1;
    q3.SetNormalize();
    q3.Print();
    Println;

    Printf("/* Function: GetNorm() of q(1,2,3,4) */\n");
    q3.SetElement(1.0f, 2.0f, 3.0f, 4.0f);
    Printf("%f\n", q3.GetNorm());
    Println;

    Printf("/* Function: GetSqareNorm() of q(1,2,3,4) */\n");
    Printf("%f\n", q3.GetSqNorm());
    Println;

    Printf("/* Function: GetSum() of q(1,2,3,4) */\n");
    Printf("%f\n", q3.GetSum());
    Println;

    Printf("/* Function: GetNormalized() of q(1,2,3,4) */\n");
    q3.GetNormalized().Print();
    Println;

    Printf("/* Function: GetConj() of q(1,2,3,4) */\n");
    q3.GetConj().Print();
    Println;

    Printf("/* Function: GetEulerAngles(order) */\n");
    Printf("qz(10deg), qy(20deg), qx(30deg)\n");
    Printf("Target quaternion = qz * qy * qx\n");
    dtQuaternion<> qz, qy, qx;
    qz(0) = std::cos(10 * DEG2RADf / 2);
    qz(3) = std::sin(10 * DEG2RADf / 2);
    qy(0) = std::cos(20 * DEG2RADf / 2);
    qy(2) = std::sin(20 * DEG2RADf / 2);
    qx(0) = std::cos(30 * DEG2RADf / 2);
    qx(1) = std::sin(30 * DEG2RADf / 2);
    q1 = qz * qy * qx;
    q1.Print();
    Printf("GetEulerAngles ZYX:\n");
    (RAD2DEGf * q1.GetEulerAngles(AXIS3(2, 1, 0))).Print();
    Println;
}

void QuatMemberAccessOperator()
{
    dtQuaternion<> q1;

    PrintHeading("Quaternion Access operators ");
    Printf("/* Operator: () */\n");

    for (uint16_t i = 0; i < 4; i++)
    {
        q1(i) = i;
        Printf("q(%d) = %7.3f\n", i, q1(i));
    }

    Println;
}

void QuatArithmetic()
{
    dtQuaternion<> q1(1.0f, 2.0f, 3.0f, 4.0f);
    dtQuaternion<> q2(1.0f, 2.0f, 3.0f, 4.0f);

    PrintHeading("Quaternion Arithmetic operators ");
    Printf("q = [w,x,y,z]T, q1 = [1,2,3,4]T, q2 = [1,2,3,4]T\n");
    Printf("/* Operator: -(), -q1 */\n");
    (-q1).Print();
    Println;

    Printf("/* Operator: +(dtQuaternion), q1 + (-q1) */\n");
    (q1 + (-q1)).Print();
    Println;

    Printf("/* Operator: -(dtQuaternion), q1 - q1 */\n");
    (q1 - q1).Print();
    Println;

    Printf("/* Operator: *(scalar), q1 * 0.1 */\n");
    (q1 * 0.1f).Print();
    Println;

    Printf("/* Operator: /(scalar), q1 / 10 */\n");
    (q1 / 10).Print();
    Println;

    Printf("/* Operator: *(scalar, dtQuaternion), 0.1 * q1 */\n");
    (0.1f * q1).Print();
    Println;

    Printf("/* Operator: *(dtQuaternion), q1 * q2 */\n");
    Printf("q1 is equal to ZYX Euler [10, 20, 30]\n");
    Printf("q2 is eual to XYZ Euler [-30, -20, -10]\n");
    q1.SetElement(AXIS3(2, 1, 0), 10 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf);
    q2.SetElement(AXIS3(0, 1, 2), -30 * DEG2RADf, -20 * DEG2RADf, -10 * DEG2RADf);
    (q1 * q2).Print();
    Println;
}

void QuatComparisonOperator()
{
    dtQuaternion<> q1(1.0f, 2.0f, 3.0f, 4.0f);
    dtQuaternion<> q2(1.0f, 2.0f, 3.0f, 4.0f);

    PrintHeading("Quaternion Comparison operators ");
    Printf("q = [w,x,y,z]T, q1 = [1,2,3,4]T, q2 = [1,2,3,4]T\n");
    Printf("/* Operator: ==(dtQuaternion &q), if q1 == q2 */\n");
    if (q1 == q2) Printf("true");
    else Printf("false");
    Println;

    Printf("/* Operator: !=(dtQuaternion &q), if q1 != q2 */\n");
    if (q1 != q2) Printf("true");
    else Printf("false");
    Println;
    Println;
}