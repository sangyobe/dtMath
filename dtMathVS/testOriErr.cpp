#include "testPrint.h"
#include "testOriErr.h"
#include "./dtMath/dtMath.h"

void Test_OrientationError()
{
    PrintTitle("Test Orientation Error");
    OrientationErr();
    AngularVeloErr();
}

void OrientationErr()
{
    PrintHeading("Orientation Error Methods ");
    CdtRotation<> dR, R, Rerr;
    CdtVec3 Vr;

    CdtQuaternion<> dQ, Q, Qerr;
    CdtVec3 Vq;

    CdtRotation<> dHR, HR;
    CdtVec3 dHp, Hp;
    CdtTransform<> dH, H;
    CdtVec3 Vh;

    Printf("/* Method1 - Rotation Matrix: Log(dR*R^T): Log map from SO(3) to so(3) */\n");
    //R1.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 30 * DEG2RADf, 20 * DEG2RADf, -15 * DEG2RADf); // Desired
    //R2.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), -15 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf); // current
    dR.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 0 * DEG2RADf, 0 * DEG2RADf, 0 * DEG2RADf); // Desired
    R.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 0 * DEG2RADf, 0 * DEG2RADf, 30 * DEG2RADf); // current
    Printf("desired Rotation matrix wrt world frame, dR = I(0[deg])\n");
    dR.Print('\n');
    Printf("current Rotation matrix wrt world frame, R = Rz(30[deg])\n");
    R.Print('\n');
    Rerr = dR * R.Transpose(); // Rerr with respect to the fixed frame
    Vr = Rerr.Log()*RAD2DEGf;   // u*phi, u is unit vector of a rotation vector, Axis/Angle
    Printf("Orientation Error wrt world frame, Log(Rerr) = \n");
    Vr.Print('\n');
    //Printf("Normalized Orientation Error wrt world frame, Norm(Log(Rerr)) = \n");
    //Vr.SetNormalize();
    //Vr.Print('\n');

    Printf("/* Method2 - Quaternion: Log(dQ*inv(Q)): Log map from so(4) to so(3) */\n");
    //Q1.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 30 * DEG2RADf, 20 * DEG2RADf, -15 * DEG2RADf); // Desired
    //Q2.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), -15 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf); // current
    dQ.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 0 * DEG2RADf, 0 * DEG2RADf, 0 * DEG2RADf); // Desired
    Q.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 0 * DEG2RADf, 0 * DEG2RADf, 30 * DEG2RADf); // current
    Printf("desired Quaternion wrt world frame, dQ = 0[deg]\n");
    dQ.Print('\n');
    Printf("current Quaternion wrt world frame, Q = Rz(30[deg])\n");
    Q.Print('\n');
    Qerr = dQ * (Q.Inv());   // Qerr with respect to the fixed frame
    Vq = Qerr.Log()*RAD2DEGf; // u*phi, u is unit vector of a rotation vector, Axis/Angle representation
    Printf("Orientation Error wrt world frame, Log(Qerr) = \n");
    Vq.Print('\n');
    //Printf("Normalized Orientation Error wrt world frame, Norm(Log(Qerr)) = \n");
    //Vq.SetNormalize();
    //Vq.Print('\n');

    Printf("/* Method3 - Quaternion: KIST Method(ref:  S.-K. Kim's paper */\n");
    //r1.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 30 * DEG2RADf, 20 * DEG2RADf, -15 * DEG2RADf); // Desired
    //r2.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), -15 * DEG2RADf, 20 * DEG2RADf, 30 * DEG2RADf); // current
    dHR.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 0 * DEG2RADf, 0 * DEG2RADf, 0 * DEG2RADf);
    HR.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 0 * DEG2RADf, 0 * DEG2RADf, 30 * DEG2RADf);
    dHp.SetElement(0, 0, 0);
    Hp.SetElement(0, 0, 0);

    dH.SetElement(dHR, dHp);
    H.SetElement(HR, Hp);

    Printf("desired Rotation matrix wrt world frame, dHR = I(0[deg])\n");
    dHR.Print('\n');
    Printf("current Rotation matrix wrt world frame, HR = Rz(30[deg])\n");
    HR.Print('\n');
    Printf("Orientation Error for geometiric jacobian wrt world frame = \n");
    Vh = dH.GetError(H).GetOri()*RAD2DEGf;
    Vh.Print('\n');
    //Printf("Normalized Orientation Error wrt world frame = \n");
    //Vh.SetNormalize();
    //Vh.Print('\n');
}

void AngularVeloErr()
{
    PrintHeading("Angular Velocity Error Method ");
       
    // W^{g}_{err} = dR*W^{d}_{d} - R*W^{c}_{c}
    // W^{b}_{err} = R^T*Rd*W^{d}_{d} - W^{c}_{c}
    // W^{d}_{err} = W^{d}_{d} - Rd^T*R*W^{c}_{c}

    // Werr = Omega(Angular velocity) error wrt global
    // dR = global to desired, R^{g}_{d}
    // Wdd = desired omega wrt desired frame
    // R = global to current, R^{g}_{c}
    // Wcc = current omega wrt current frame

    CdtRotation<> dR, R;
    CdtVector3<> Wcc, Wdd;
    CdtVector3<> Werr;

    dR.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 0 * DEG2RADf, 0 * DEG2RADf,30 * DEG2RADf);
    R.SetElement(AXIS3(Z_AXIS, Y_AXIS, X_AXIS), 30 * DEG2RADf, 0 * DEG2RADf, 0 * DEG2RADf);

    Printf("desired Rotation(30deg @Xaxis) global to desired frame = \n");
    dR.Print('\n');

    Printf("current Rotation(30deg @Zaxis) global to current frame = \n");
    R.Print('\n');

    Printf("desired Angular Velocity(10dps @Xaxis) wrt desired frame = \n");
    Wdd << 10, 0, 0;
    Wdd.Print('\n');
    Wdd *= DEG2RADf;

    Printf("current Angular Velocity(10dps @Zaxis) wrt current frame = \n");
    Wcc << 0, 0, 10; // 10 deg/sec @ Z-axis
    Wcc.Print('\n');
    Wcc *= DEG2RADf; 

    Printf("Angular Velocity Error(so3) wrt global frame = \n");
    Printf("Werr = dR*Wdd - R*Wcc =\n");
    Werr = dR * Wdd - R * Wcc;
    (Werr*RAD2DEGf).Print('\n');
}