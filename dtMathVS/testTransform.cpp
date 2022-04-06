#include "testPrint.h"
#include "testTransform.h"
#include "./dtMath/dtMath.h"

void Test_Transform()
{
    PrintTitle("Test Homogeneous transformation matrix");
    TrfMemberFunc();
    TrfMemberAccessOperator();
    TrfArithmetic();
    TrfComparisonOperator();
}

void TrfMemberFunc()
{
    PrintHeading("Transformation Member Functions ");

    Printf("/* Class Create: () */\n");
    CdtTransform<> t1;
    t1.Print();
    Println;

    Printf("/* Class Create: with Rotation */\n");
    Printf("Rot = Rz(10)Ry(20)Rx(30) deg\n");
    Printf("Pos = [1, 2, 3]T\n");
    CdtRotation<> r1(AXIS3(2, 1, 0), 10 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf);
    CdtVector3<> v1(1, 2, 3);
    CdtTransform<> t2(r1, v1);
    t2.Print();
    Println;

    Printf("/* Class Create: with Quaternion */\n");
    Printf("Quat = ZYX Euler [10, 20, 30]T deg\n");
    Printf("Pos = [1, 2, 3]T\n");
    CdtQuaternion<> q1(AXIS3(2, 1, 0), 10 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf);
    CdtTransform<> t3(q1, v1);
    t3.Print();
    Println;

    Printf("/* Class Create: with Euler */\n");
    Printf("ZYX Euler = [10, 20, 30]T deg\n");
    Printf("Pos = [1, 2, 3]T\n");
    CdtVector3<> euler(10 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf);
    CdtTransform<> t4(AXIS3(2,1,0), euler, v1);
    t4.Print();
    Println;

    Printf("/* Class Create: CdtTransform */\n");
    CdtTransform<> t5(t4);
    t5.Print();
    Println;
    
    Printf("/* Function: SetZero() */\n");
    t1.SetZero();
    t1.Print();
    Println;

    Printf("/* Function: SetIdentity() */\n");
    t1.SetIdentity();
    t1.Print();
    Println;

    Printf("/* Function: SetElement(CdtVector3 pos) */\n");
    Printf("Pos - [1,2,3]T\n");
    v1.SetElement(1, 2, 3);
    t1.SetElement(v1);
    t1.Print();
    Println;

    Printf("/* Function: SetElement(CdtRotation) */\n");
    Printf("Z-Axis 10deg Rotatin matrix\n");
    r1.SetElement(AXIS1(0), 10 * DEG2RADf);
    t1.SetElement(r1);
    t1.Print();
    Println;

    Printf("/* Function: SetElement(CdtQuaternion) */\n");
    Printf("X-Axis 10deg Quaternion\n");
    q1.SetElement(AXIS1(0), 10 * DEG2RADf);
    t1.SetElement(r1);
    t1.Print();
    Println;

    Printf("/* Function: SetElement(euler angle order, CdtVector3 euler angles) */\n");
    Printf("XYZ Euler Angle - [10,0,0]T deg\n");
    euler.SetElement(10 * DEG2RADf, 0, 0);
    t1.SetElement(AXIS3(0, 1, 2), euler);
    t1.Print();
    Println;

    Printf("/* Function: SetElement(CdtRotation, CdtVector3 Pos) */\n");
    Printf("Z-Axis 10deg Quaternion\n");
    Printf("Position - [1,0,0]T\n");
    r1.SetElement(AXIS1(2), 10 * DEG2RADf);
    v1.SetElement(1, 0, 0);
    t1.SetElement(r1, v1);
    t1.Print();
    Println;

    Printf("/* Function: SetElement(CdtQuaternion, CdtVector3 Pos) */\n");
    Printf("Z-Axis 10deg Quaternion\n");
    Printf("Position - [0,1,0]T\n");
    q1.SetElement(AXIS1(2), 10 * DEG2RADf);
    v1.SetElement(0, 1, 0);
    t1.SetElement(q1, v1);
    t1.Print();
    Println;

    Printf("/* Function: SetElement(euler angle order, CdtVector3 Euler, CdtVector3 Pos) */\n");
    Printf("XYZ Euler Angle - [0,0,10]T deg\n");
    Printf("Position - [0,0,1]T\n");
    euler.SetElement(0, 0, 10 * DEG2RADf);
    v1.SetElement(0, 0, 1);
    t1.SetElement(AXIS3(0, 1, 2), euler, v1);
    t1.Print();
    Println;

    Printf("/* Function: q() - quaternion */\n");
    Printf("Target Matrix T:\n");
    t1.Print();
    Printf("Quaternion of T:\n");
    t1.q().Print();
    Println;

    Printf("/* Function: R() - rotation matrix */\n");
    Printf("Target Matrix T:\n");
    t1.Print();
    Printf("Rotation matrix of T:\n");
    t1.R().Print();
    Println;

    Printf("/* Function: e() - euler angles */\n");
    Printf("Target Matrix T:\n");
    t1.Print();
    Printf("XYZ Euler Angle of T:\n");
    t1.e(AXIS3(0, 1, 2)).Print();
    Println;

    Printf("/* Function: p() - position */\n");
    Printf("Target Matrix T:\n");
    t1.Print();
    Printf("Position of T:\n");
    t1.p().Print();
    Println;

    Printf("/* Function: GetError(CdtTransform) */\n");
    Printf("T1_ori: ZYX Euler [30,0,0]T, T1_pos:[10,20,30]T\n");
    r1.SetElement(AXIS1(2), 30 * DEG2RADf);
    v1.SetElement(10, 20, 30);
    t1.SetElement(r1, v1);
    (t1.R().GetEulerAngles(AXIS3(2, 1, 0)) * RAD2DEGf).Print();
    Println;

    Printf("T2_ori: ZYX Euler [0,0,0]T, T2_pos:[0,0,0]T\n");
    r1.SetElement(AXIS1(2), 0 * DEG2RADf);
    v1.SetElement(0, 0, 0);
    t2.SetElement(r1, v1);
    (t2.R().GetEulerAngles(AXIS3(2, 1, 0))* RAD2DEGf).Print();
    Println;

    Printf("T1 and T2 position error:\n");
    t1.GetError(t2).GetPos().Print();
    Printf("T1 and T2 orientation error( = rotation vector or angular vector or angle-axis vector):\n");
    (t1.GetError(t2).GetOri() * RAD2DEGf).Print();
    Println;

    Printf("/* Function: Transpose() */\n");
    Printf("Target Matrix:\n");
    t1.Print();
    Printf("Transposed Matrix:\n");
    t1.Transpose().Print();
    Println;

    Printf("/* Function: Inv() */\n");
    Printf("Target Matrix:\n");
    t1.Print();
    Printf("Inverse:\n");
    t1.Inv().Print();
    Println;

    Printf("Target Matrix x Inverse Matrix:\n");
    (t1* t1.Inv()).Print();
    Println;

    Printf("Inverse Matrix x Inverse Matrix:\n");
    (t1.Inv() * t1).Print();
    Println;
}

void TrfMemberAccessOperator()
{
    PrintHeading("Transformation Member Access Operator ");
    CdtRotation<> r1(AXIS1(2), 30 * DEG2RADf);
    CdtVector3<> p1(1, 2, 3);
    CdtTransform<> t1(r1, p1);

    Printf("/* Operator: () */\n");
    Printf("T(i,j) = 100; & Printf(""%%f"", T(i,j));\n");
    for (uint16_t i = 0; i < 4; i++)
    {
        for (uint16_t j = 0; j < 4; j++)
        {
            t1(i, j) = 100;
            Printf("%7.3f ", t1(i, j));
        }
        Printf("\n");
    }
    Println;
}

void TrfArithmetic()
{
    PrintHeading("Transformation Arithmetic operators ");
    float v[3] = { 1,2,3 };
    CdtTransform<> t1;
    CdtTransform<> t2;
    CdtMatrix<4, 4> m44;
    CdtMatrix<4, 6> m46;
    CdtVector<3> v31(v, sizeof(v));
    CdtVector3<> v3(1, 2, 3);

    Printf("T =\n");
    t1.Print();
    Println;

    Printf("/* Operator: +(const CdtMatrix) */\n");
    m44.SetIdentity();
    Printf("mat44 =\n");
    m44.Print();
    Printf("T + -mat44 =\n");
    (t1 + -m44).Print();
    Println;

    Printf("/* Operator: -(const CdtMatrix) */\n");
    m44.SetIdentity();
    Printf("mat44 =\n");
    m44.Print();
    Printf("T - mat44 =\n");
    (t1 - m44).Print();
    Println;

    Printf("/* Operator: *(const CdtMatrix) */\n");
    m44.SetFill(1);
    Printf("mat44 =\n");
    m44.Print();
    Printf("T * mat44 =\n");
    (t1 * m44).Print();
    Println;

    Printf("/* Operator: *(const CdtMatrix) */\n");
    m46.SetFill(1);
    Printf("mat46 =\n");
    m44.Print();
    Printf("T * mat46 =\n");
    (t1 * m46).Print();
    Println;

    Printf("/* Operator: *(const CdtTransform) */\n");
    t2.SetElement(CdtVector3<>(1, 2, 3));
    Printf("T2 =\n");
    t2.Print();
    Printf("T * T2 =\n");
    (t1 * t2).Print();
    Println;

    Printf("/* Operator: *(const CdtVector) */\n");
    Printf("vec31 =\n");
    v31.Print();
    Printf("T * vec31 =\n");
    (t1 * v31).Print();
    Println;

    Printf("/* Operator: *(const CdtVector3) */\n");
    Printf("vec3 =\n");
    v3.Print();
    Printf("T * vec3 =\n");
    (t1 * v3).Print();
    Println;
}

void TrfComparisonOperator()
{
    PrintHeading("Transformation Comparison operators ");
    CdtTransform<> t1;

    Printf("/* Operator: ==(CdtTransform &m) and !=(CdtTransform &m) */\n");
    if (t1 == t1) Printf("true");
    else Printf("false");
    Println;

    if (t1 != t1) Printf("false");
    else Printf("true");
    Println;
    Println;
}
