#if defined(ARDUINO)
#include <Arduino.h>
#else
#include <string.h>
#include "dhTerm.h"
using namespace dhTerm;
#endif

//#include <dtMath/dtMath.h>

#include "testMatrix.h"
#include "testCscMatrix.h"
#include "testVector.h"
#include "testMatrix3.h"
#include "testVector3.h"
#include "testVector4.h"
#include "testVector6.h"
#include "testRotation.h"
#include "testTransform.h"
#include "testOriErr.h"
#include "testQuaternion.h"
#include "testEulerAngles.h"
#include "testDecomposition.h"
#include "testCompareClass.h"
#include "testMatrixInverse.h"
#include "testQuadProg.h"
#include "testMPC.h"
#include "testContiAlgebraicRiccatiEq.h"
#include "testGnuPlot.h"

#include <robotics/djBezier.h>
#include <robotics/djBezierOri.h>
#include <robotics/dtSwingBezier.h>
#include <robotics/dtStanceBezier.h>
#include <robotics/dtGaitPlanner.h>

/* Choose test functions */
int8_t testTB[22] = {
    1, // Test_dtMatrix
    0, // Test_dtCscMatrix
    0, // Test_dtVector
    0, // Test_dtMatrix3
    0, // Test_dtVector3
    0, // Test_dtVector4
    0, // Test_dtVector6
    0, // Test_dtRotation
    0, // Test_dtTransform
    0, // Test_OrientationError
    0, // Test_dtQuaternion
    0, // Test_EulerAngles
    0, // Test_Decomposition
    0, // Test_MatrixInverse
    0, // Test_CompareClass
    0, // Test_QuadProg
    0, // Test_QuadrupedRobotQP
    0, // Test_ContiAlgebraicRiccatiEq
    0, // Test_MPC
    0, // Test_dtGnuPlot
};

enum Test
{
    IDX_dtMatrix = 0,
    IDX_dtCscMatrix,
    IDX_dtVector,
    IDX_dtMatrix3,
    IDX_dtVector3,
    IDX_dtVector4,
    IDX_dtVector6,
    IDX_dtRotation,
    IDX_dtTransform,
    IDX_dtOriErr,
    IDX_dtQuaternion,
    IDX_dtEulerAngles,
    IDX_dtDecomposition,
    IDX_dtMatrixInverse,
    IDX_dtCompareClass,
    IDX_dtQuadProg,
    IDX_dtQuadProgExample,
    IDX_dtCARE,
    IDX_dtMPC,
    IDX_dtGnuPlot,
};

#if defined(ARDUINO)
void setup()
{
    Serial.begin(115200);
    delay(2000);
}

void loop()
{
    if (testTB[IDX_dtMatrix]) Test_Matrix();
    if (testTB[IDX_dtVector]) Test_Vector();
    if (testTB[IDX_dtMatrix3]) Test_Matrix3();
    if (testTB[IDX_dtVector3]) Test_Vector3();
    if (testTB[IDX_dtVector4]) Test_Vector4();
    if (testTB[IDX_dtVector6]) Test_Vector6();
    if (testTB[IDX_dtRotation]) Test_Rotation();
    if (testTB[IDX_dtTransform]) Test_Transform();
    if (testTB[IDX_dtOriErr]) Test_OriErr();
    if (testTB[IDX_dtQuaternion]) Test_Quaternion();
    if (testTB[IDX_dtEulerAngles]) Test_EulerAngles();
    if (testTB[IDX_dtDecomposition]) Test_Decomposition();
    if (testTB[IDX_dtMatrixInverse]) Test_MatrixInverse();
    if (testTB[IDX_dtCompareClass]) Test_CompareClass();
    if (testTB[IDX_dtQuadProg]) Test_QuadProg();
    if (testTB[IDX_dtQuadProgExample]) Test_QuadrupedRobotQP();
    if (testTB[IDX_dtCARE]) Test_ContiAlgebraicRiccatiEq();

    memset(testTB, 0, sizeof(testTB));
}
#else

void Test_Bezier()
{
    CdtMatrix<3, 1> ctrlPt;
    CdjBezier<3, 1> bezier;
    //Bezier<3, 1> bzOri;
    output<1> bzOriOut;

    //BezierOut<3> bezierOut;
    float period = 2;
    float elapseTime = 2;

    //for (int i = 0; i < 10; i++)
    //{
    //    bezier.BezierInterp(ctrlPt, elapseTime / period, period);
    //    elapseTime += 0.002f;
    //    bezier.pos.Print();
    //    if (elapseTime > period) elapseTime = period;
    //}
    elapseTime = 0;
    ctrlPt << 0, 60, 0;
    ctrlPt *= DEG2RADf;

    while (1)
    {
        bezier.BezierInterp(ctrlPt, elapseTime/period, period);
        elapseTime += 0.002f;
        bezier.pos.Print();
        if (elapseTime > period) break;
    }

    //bezier.BezierInterp(ctrlPt, elapseTime/period, period);
    //bezier.pos.Print('\n');
    //bzOriOut = bzOri.BezierInterp(ctrlPt, elapseTime/period, period);
    //(bzOriOut.pos - bezier.pos).Print();
    //bezierOut.pos.Print();
}

void Test_SwingBezier()
{
    CdtVector3<> startPos;
    CdtVector3<> endPos(0.3f, 0, 0);
    CdtVector3<> startVel;
    CdtVector3<> endVel;
    float h = 0.3f;
    float period_ms = 300.0f;

    CdtSwingBezier<> swingBezier;

    swingBezier.Init(1000, period_ms); // sampling time is 1ms
    swingBezier.SetTarget(startPos, endPos, startVel, endVel, h, period_ms);
    
    printf("Swing bezier control points\n");
    for (int i = 0; i < 16; i++)
    {
        swingBezier.CtrlPts(i).Transpose().Print();
    }

    printf("\nSwing bezier position\n");
    for (int i = 0; i < 300; i++)
    {
        swingBezier.Compute();
        //swingBezier.pos_x.Transpose().Print();
        swingBezier.vel_xps.Transpose().Print();
    }
}

void Test_GaitPlaner()
{
    CdtGaitPlaner<4> gait;
    CdtVector3<> startPos[4];
    CdtVector3<> startVel[4];
    CdtVector3<> step(0.3f, 0, 0);
    float swingHeight = 0.3f;

    float phaseDelay_pu[4] = { 0, 0.25f, 0.5f, 0.75f };
    gait.Init(1000, 300, 700, phaseDelay_pu);
    gait.BindData(0, &startPos[0], &startVel[0]);
    gait.BindData(1, &startPos[1], &startVel[1]);
    gait.BindData(2, &startPos[2], &startVel[2]);
    gait.BindData(3, &startPos[3], &startVel[3]);
    gait.SetStep(step, swingHeight);
    gait.Start();

    for (int i = 0; i < 1000; i++)
    {
        gait.Compute();
        startPos[0] = gait.Pos(0);
        startVel[0] = gait.Vel(0);
        startPos[1] = gait.Pos(1);
        startVel[1] = gait.Vel(1);
        startPos[2] = gait.Pos(2);
        startVel[2] = gait.Vel(2);
        startPos[3] = gait.Pos(3);
        startVel[3] = gait.Vel(3);
    }
}

int main()
{
    SetupTerminal(false);

    if (testTB[IDX_dtMatrix]) Test_Matrix();
    if (testTB[IDX_dtCscMatrix]) Test_CscMatrix();
    if (testTB[IDX_dtVector]) Test_Vector();
    if (testTB[IDX_dtMatrix3]) Test_Matrix3();
    if (testTB[IDX_dtVector3]) Test_Vector3();
    if (testTB[IDX_dtVector4]) Test_Vector4();
    if (testTB[IDX_dtVector6]) Test_Vector6();
    if (testTB[IDX_dtRotation]) Test_Rotation();
    if (testTB[IDX_dtTransform]) Test_Transform();
    if (testTB[IDX_dtOriErr]) Test_OrientationError();
    if (testTB[IDX_dtQuaternion]) Test_Quaternion();
    if (testTB[IDX_dtEulerAngles]) Test_EulerAngles();
    if (testTB[IDX_dtDecomposition]) Test_Decomposition();
    if (testTB[IDX_dtMatrixInverse]) Test_MatrixInverse();
    if (testTB[IDX_dtCompareClass]) Test_CompareClass();
    if (testTB[IDX_dtQuadProg]) Test_QuadProg();
    if (testTB[IDX_dtQuadProgExample]) Test_QuadrupedRobotQP();
    if (testTB[IDX_dtCARE]) Test_ContiAlgebraicRiccatiEq();
    if (testTB[IDX_dtMPC]) Test_MPC();
    if (testTB[IDX_dtGnuPlot]) Test_GnuPlot();
    
    memset(testTB, 0, sizeof(testTB));

    //Test_Bezier();
    //Test_SwingCtrlPts();
    Test_GaitPlaner();

    RestoreTerminal();

    return 0;
}
#endif