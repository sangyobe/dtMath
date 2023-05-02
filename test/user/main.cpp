#if defined(ARDUINO)
#include <Arduino.h>
#else
#include "dhTerm.h"
#include <string.h>
using namespace dhTerm;
#endif

// #include <dtMath/dtMath.h>

#include "testCompareClass.h"
#include "testCscMatrix.h"
#include "testDecomposition.h"
#include "testEulerAngles.h"
#include "testGnuPlot.h"
#include "testMPC.h"
#include "testMatrix.h"
#include "testMatrix3.h"
#include "testMatrixInverse.h"
#include "testOriErr.h"
#include "testQuadProg.h"
#include "testQuaternion.h"
#include "testRotation.h"
#include "testTransform.h"
#include "testVector.h"
#include "testVector3.h"
#include "testVector4.h"
#include "testVector6.h"

/* Choose test functions */
int8_t testTB[22] = {
    1, // Test_dtMatrix
    1, // Test_dtCscMatrix
    1, // Test_dtVector
    1, // Test_dtMatrix3
    1, // Test_dtVector3
    1, // Test_dtVector4
    1, // Test_dtVector6
    1, // Test_dtRotation
    1, // Test_dtTransform
    1, // Test_OrientationError
    1, // Test_dtQuaternion
    1, // Test_EulerAngles
    1, // Test_Decomposition
    1, // Test_MatrixInverse
    1, // Test_CompareClass
    1, // Test_QuadProg
    0, // Test_QuadrupedRobotQP
    1, // Test_ContiAlgebraicRiccatiEq
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

    memset(testTB, 0, sizeof(testTB));
}
#else

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
    if (testTB[IDX_dtMPC]) Test_MPC();
    if (testTB[IDX_dtGnuPlot]) Test_GnuPlot();

    memset(testTB, 0, sizeof(testTB));

    // Test_Bezier();
    // Test_SwingCtrlPts();
    // Test_GaitPlaner();

    RestoreTerminal();

    return 0;
}
#endif