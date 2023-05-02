#include "testEulerAngles.h"
#include "testPrint.h"

#include <dtMath/dtMath.h>

using dtMath::dtMatrix;
using dtMath::dtQuaternion;
using dtMath::dtRotMat;
using dtMath::dtVec3;
using dtMath::dtVector;
using dtMath::dtVector3;

void Test_EulerAngles()
{
    PrintTitle("Test Euler Angles");
    EulerToQuaternion();
    EulerToRotation();
}

void EulerToQuaternion()
{
    // euler ZYX to Quaternion
    PrintHeading("Euler Angles ZYX to Quaternion ");
    double zyx[3] = {12 * DEG2RADd, 67 * DEG2RADd, 32 * DEG2RADd};
    dtVector3<double> euler(zyx, sizeof(zyx));
    dtQuaternion<double> quat(AXIS3(2, 1, 0), euler); // AXIS(2,1,0) means the ZYX Euler order

    Printf("Euler ZYX[deg] = \n");
    (RAD2DEGd * euler).Print();
    Printf("Euler ZYX to Quaternion =\n");
    quat.Print();
    Printf("Quaternion to Euler ZYX[deg] =\n");
    (RAD2DEGd * quat.GetEulerAngles(AXIS3(Z_AXIS, Y_AXIS, X_AXIS))).Print();
    Printf("\n");
}

void EulerToRotation()
{
    // euler ZYX to Rotation
    PrintHeading("Euler Angles ZYX to Rotation ");
    float zyx[3] = {12 * DEG2RADf, 67 * DEG2RADf, 32 * DEG2RADf};
    dtVec3 euler(zyx, sizeof(zyx));                     // is equal to dtVector3<float> euler(zyx);
    dtRotMat rot(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), euler); // is equal to dtRotation<float> rot(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), euler);

    Printf("Euler ZYX[deg] = \n");
    (RAD2DEGf * euler).Print();
    Printf("Euler ZYX to Rotation matrix =\n");
    rot.Print();
    Printf("Rotation matrix to Euler ZYX[deg] =\n");
    (RAD2DEGf * rot.GetEulerAngles(AXIS3(Z_AXIS, Y_AXIS, X_AXIS))).Print();
    Printf("\n");
}