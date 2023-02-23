#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
#include <utility/dhTimeCheck.h>

#include <utility/dtGnuPlot.h>
#endif

#include <cmath>
#include <limits>
#include "testPrint.h"
#include "testMPC.h"

#include <dtMath/dtMath.h>


#define Nstep 5
static const int simulLen = 2000;
//static float smpT = 0.02f; // 50Hz(20ms)
static float smpT = 0.25f; // 50Hz(20ms)
static float h = 0.5f;     // height of CoM
static float g = 9.806f;   // gravity
static float h_g = h / g;
static float R = 1.0f;      // weight of jerk
static float Q = 100000.0f; // weight of (z-zref)

CdtMatrix<Nstep, 3> Pps, Pvs, Pzs;
CdtMatrix<Nstep, Nstep> Ppu, Pvu, Pzu;
CdtMatrix<3, 3> A;
CdtVector<3> B;
CdtVector<3> C;

void MakeRelationMat(float T);
void MakeTimeData(float *timeData);
void MakeRefXzmp(CdtVector<simulLen> &xZmpRef);
void MakeRefYzmp(CdtVector<simulLen> &yZmpRef);

void ADH_MPC_Test();
void ADH_MPC_Test2();

void Test_MPC()
{
    //AnalyticalPreviewCtrl();
    MPC_Method1();  // dtMath Style
    MPC_Method2();  // OSQP Style

    //ADH_MPC_Test();
    //ADH_MPC_Test2();
}

void MakeRelationMat(float T)
{
    // x(k+1) = A*x(k) + B*x_jerk(k)
    A << 1, T, T*T / 2.0f, 0, 1, T, 0, 0, 1;
    B << T * T*T / 6.0f, T*T / 2.0f, T;
    C << 1, 0, -h_g;

    // Generate Pps, Pvs, Pzs Matrix (Nx3)
    float rowVec[3];
    float Ppu_col[Nstep];
    float Pvu_col[Nstep];
    float Pzu_col[Nstep];

    for (int i = 0; i < Nstep; i++)
    {
        rowVec[0] = 0;
        rowVec[1] = 1;
        rowVec[2] = (i + 1) * T;
        Pvs.SetRowVec(i, rowVec, sizeof(float) * 3);

        rowVec[0] = 1;
        rowVec[1] = (i + 1) * T;
        rowVec[2] = (i + 1)*(i + 1) * T*T / 2.0f;
        Pps.SetRowVec(i, rowVec, sizeof(float) * 3);

        rowVec[2] -= h_g;
        Pzs.SetRowVec(i, rowVec, sizeof(float) * 3);

        // for Pxu matrix
        Ppu_col[i] = (1 + 3 * i + 3 * i*i) * T*T*T / 6.0f;
        Pvu_col[i] = (1 + 2 * i) * T*T / 2.0f;
        Pzu_col[i] = (1 + 3 * i + 3 * i*i) * T*T*T / 6.0f - (T * h_g);
    }

    // Generate Ppu, Pvu, Pzu Matrix (NxN)
    CdtVector<Nstep> vec;
    for (int i = 0; i < Nstep; i++)
    {
        vec.SetZero();
        vec.SetBlock(i, Ppu_col, sizeof(float)*Nstep);
        Ppu.SetColVec(i, vec);

        vec.SetZero();
        vec.SetBlock(i, Pvu_col, sizeof(float)*Nstep);
        Pvu.SetColVec(i, vec);

        vec.SetZero();
        vec.SetBlock(i, Pzu_col, sizeof(float)*Nstep);
        Pzu.SetColVec(i, vec);
    }
}

void MakeTimeData(float *timeData)
{
    float time = 0;
    for (int i = 0; i < simulLen; i++)
    {
        timeData[i] = time;
        time += smpT;
    }
}

void MakeRefXzmp(CdtVector<simulLen> &xZmpRef)
{
    // make reference zmp
    float time = 0;
    for (int i = 0; i < simulLen; i++)
    {
        if (time < 1) xZmpRef(i) = 0.0f;
        else if (time >= 1 && time < 2) xZmpRef(i) = 0.2f;
        else if (time >= 2 && time < 3) xZmpRef(i) = 0.4f;
        else if (time >= 3 && time < 4) xZmpRef(i) = 0.6f;
        else if (time >= 4 && time < 5) xZmpRef(i) = 0.8f;
        else if (time >= 5 && time < 6) xZmpRef(i) = 1.0f;
        else if (time >= 6 && time < 7) xZmpRef(i) = 1.2f;
        else if (time >= 7 && time < 8) xZmpRef(i) = 1.4f;
        else if (time >= 8 && time < 9) xZmpRef(i) = 1.6f;
        else if (time >= 9 && time < 10) xZmpRef(i) = 1.8f;
        else xZmpRef(i) = 1.8f;

        time += smpT;
    }
}

void MakeRefYzmp(CdtVector<simulLen> &yZmpRef)
{
    // make reference zmp
    float time = 0;
    for (int i = 0; i < simulLen; i++)
    {
        if (time < 1) yZmpRef(i) = 0.0f;
        else if (time >= 1 && time < 2) yZmpRef(i) = -0.2f;
        else if (time >= 2 && time < 3) yZmpRef(i) = 0.2f;
        else if (time >= 3 && time < 4) yZmpRef(i) = -0.2f;
        else if (time >= 4 && time < 5) yZmpRef(i) = 0.2f;
        else if (time >= 5 && time < 6) yZmpRef(i) = -0.2f;
        else if (time >= 6 && time < 7) yZmpRef(i) = 0.2f;
        else if (time >= 7 && time < 8) yZmpRef(i) = -0.2f;
        else if (time >= 8 && time < 9) yZmpRef(i) = 0.2f;
        else if (time >= 9 && time < 10) yZmpRef(i) = -0.2f;
        else yZmpRef(i) = 0.0f;

        time += smpT;
    }
}

void AnalyticalPreviewCtrl()
{// ZMP Preview controller
    PrintTitle("Analytical ZMP Preview Ctrl ");

    CdhTimeCheck clock;
    float timeData[simulLen] = { 0, };
    float time = 0;

    // graph variables
#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
    CdtGnuPlot<> graphSolvTime;
    CdtGnuPlot<> graphX;
    CdtGnuPlot<> graphY;
#endif
    float yZmpData[simulLen] = { 0, };
    float yComData[simulLen] = { 0, };
    float ySolvTimeData[simulLen] = { 0, };

    int8_t ok;
    CdtMatrix<Nstep, Nstep> Inn;
    CdtMatrix<Nstep, Nstep> Inv;

    CdtVector<3> xState;
    CdtVector<simulLen> xZmpRef;
    CdtVector<Nstep> xZmpRefK;
    CdtVector<Nstep> xJerk;
    float xZmpData[simulLen] = { 0, };
    float xComData[simulLen] = { 0, };

    CdtVector<3> yState;
    CdtVector<simulLen> yZmpRef;
    CdtVector<Nstep> yZmpRefK;
    CdtVector<Nstep> yJerk;

    float R_Q = R / Q; // R/Q

    Inn.SetIdentity();
    MakeRelationMat(smpT);
    MakeTimeData(timeData);
    MakeRefXzmp(xZmpRef);
    MakeRefYzmp(yZmpRef);

    time = 0;
    Inv = (Pzu.Transpose() * Pzu + R_Q * Inn).Inv(&ok);
    if (!ok)
    {
        Printf("Inverse error\n");
        return;
    }

    for (int i = 0; i < (simulLen - Nstep); i++)
    {
        clock.Start();
        /* update Zmp_ref */
        xZmpRefK.SetZero();
        yZmpRefK.SetZero();

        xZmpRef.GetBlock(i, xZmpRefK);
        yZmpRef.GetBlock(i, yZmpRefK);

        /* update state */
        // In practically, state is updated by the sensor and kinematics algo.
        // In simulations, state is updated with equations (idealy, there are no disturbunce)

        /* evaluate the jerk */
        xJerk = -Inv * Pzu.Transpose()*(Pzs*xState - xZmpRefK);
        yJerk = -Inv * Pzu.Transpose()*(Pzs*yState - yZmpRefK);

        /* desired CoM position, velocity and acceleration */
        xState = A * xState + B * xJerk(0);   // x(k+1) = A*x(k) + B*x_jerk(k)
        yState = A * yState + B * yJerk(0);   // y(k+1) = A*y(k) + B*y_jerk(k)

        clock.Stop();

        /* save graph data */
        ySolvTimeData[i] = (float)clock.GetElapsedTime_msec();
        xComData[i] = xState(0);
        xZmpData[i] = C.dot(xState);
        yComData[i] = yState(0);
        yZmpData[i] = C.dot(yState);
    }

#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
    graphSolvTime.SetTitle("PreviewCtrl Solver Time", 15, "arial");
    graphSolvTime.SetXlabel("time[sec]");
    graphSolvTime.SetYlabel("load[ms]");
    graphSolvTime.SetXtics(1);
    graphSolvTime.Point(timeData, ySolvTimeData, simulLen, "load");
    graphSolvTime.Draw();

    graphX.SetTitle("PreviewCtrl X ZMP \\\\& CoM", 15, "arial");
    graphX.SetXlabel("time[sec]");
    graphX.SetYlabel("X[m]");
    graphX.SetXrange(0, 20);
    graphX.SetYrange(0.f, 2.f);
    graphX.Dash(timeData, xZmpRef.GetElementsAddr(), simulLen, "zmp ref", "-", "dark-green");
    graphX.Line(timeData, xComData, simulLen, "com result", "red");
    graphX.Line(timeData, xZmpData, simulLen, "zmp result", "blue");
    graphX.Draw();

    graphY.SetTitle("PreviewCtrl Y ZMP \\\\& CoM", 15, "arial");
    graphY.SetXlabel("time[sec]");
    graphY.SetYlabel("Y[m]");
    graphY.SetXrange(0, 20);
    graphY.SetYrange(-0.3f, 0.3f);
    graphY.Dash(timeData, yZmpRef.GetElementsAddr(), simulLen, "zmp ref", "-", "dark-green");
    graphY.Line(timeData, yComData, simulLen, "com result", "red");
    graphY.Line(timeData, yZmpData, simulLen, "zmp result", "blue");
    graphY.Draw();
#endif
}

void MPC_Method1()
{
    PrintTitle("ZMP MPC-Method1 ");
    CdhTimeCheck clock;
    float timeData[simulLen] = { 0, };
    float time = 0;

    // graph vairables
#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
    CdtGnuPlot<> graphSolvTime;
    CdtGnuPlot<> graphCost;
    CdtGnuPlot<> graphY;
#endif
    float yZmpData[simulLen] = { 0, };
    float yPosComData[simulLen] = { 0, };
    float yVelComData[simulLen] = { 0, };
    float yCostData[simulLen] = { 0, };
    float ySolvTimeData[simulLen] = { 0, };

    // min a/2*jerk^2 + b/2*(Z-Zref)^2 wrt jerk
    CdtMatrix<2 * Nstep, 2 * Nstep> mGy;
    float w[2 * Nstep];
    CdtVector<2 * Nstep> vGy;

    CdtMatrix<2 * Nstep, 2 * Nstep> mCIy;
    CdtVector<2 * Nstep> vCIy;

    CdtMatrix<2 * Nstep, Nstep> mCEy;
    CdtVector<Nstep> vCEy;

    CdtMatrix<Nstep, Nstep> PzuT;
    CdtMatrix<Nstep, Nstep> Inn;
    CdtVector<2 * Nstep> yX;          // yX = [yJerk; (yZ - yZref)]
    CdtVector<3> yState;

    CdtVector<simulLen> yZmpRef;
    CdtVector<Nstep> yZmpRefK;
    CdtVector<simulLen> yZmpMax;
    CdtVector<Nstep> yZmpMaxK;
    CdtVector<simulLen> yZmpMin;
    CdtVector<Nstep> yZmpMinK;

    CdtQuadProg<2 * Nstep, 2 * Nstep, Nstep> yMPC;
    float yCost = 0;

    MakeRelationMat(smpT);
    Inn.SetIdentity();
    PzuT = Pzu.Transpose();

    /* Weight */
    for (int i = 0; i < Nstep; i++)
    {
        w[i] = R;
    }
    for (int i = Nstep; i < 2 * Nstep; i++)
    {
        w[i] = Q;
    }
    mGy.SetDiagonal(w, sizeof(w));
    vGy.SetZero();

    /* Equality */
    mCEy.SetBlock(0, 0, PzuT);
    mCEy.SetBlock(Nstep, 0, -Inn);
    vCEy = Pzs * yState - yZmpRefK;

    /* Inequality */
    mCIy.SetBlock(0, 0, PzuT);
    mCIy.SetBlock(0, Nstep, -PzuT);
    vCIy.SetBlock(0, Pzs*yState - yZmpMinK);
    vCIy.SetBlock(Nstep, yZmpMaxK - Pzs * yState);

    // make timedata & ref zmp
    MakeTimeData(timeData);
    MakeRefYzmp(yZmpRef);
    yZmpMax = yZmpRef + 0.05f;
    yZmpMin = yZmpRef - 0.05f;

    yMPC.SetObjectFunc(mGy, vGy);

    for (int i = 0; i < (simulLen - Nstep); i++)
    {
        clock.Start();
        /* update Zmp_ref */
        yZmpRefK.SetZero();
        yZmpMinK.SetZero();
        yZmpMaxK.SetZero();

        yZmpRef.GetBlock(i, yZmpRefK);
        yZmpMin.GetBlock(i, yZmpMinK);
        yZmpMax.GetBlock(i, yZmpMaxK);

        /* update State */
        vCEy = Pzs * yState - yZmpRefK;
        vCIy.SetBlock(0, Pzs*yState - yZmpMinK);
        vCIy.SetBlock(Nstep, yZmpMaxK - Pzs * yState);

        /* Solve */
        if (yMPC.Solve(mCEy, vCEy, mCIy, vCIy, yX))
        {
            Printf("QP Solution Fail!\n");
            return;
        }

        /* desired CoM position, velocity and acceleration */
        yState = A * yState + B * yX(0);   // y(k+1) = A*y(k) + B*y_jerk(k)

        clock.Stop();

        /* save graph data */
        yPosComData[i] = yState(0);
        yVelComData[i] = yState(1);
        yZmpData[i] = C.dot(yState);
        yCostData[i] = yMPC.GetObjectValue();
        ySolvTimeData[i] = (float)clock.GetElapsedTime_msec();
    }
#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
    graphSolvTime.SetTitle("Method1 MPC Solver Time", 15, "arial");
    graphSolvTime.SetXlabel("time[sec]");
    graphSolvTime.SetYlabel("load[ms]");
    graphSolvTime.SetXtics(1);
    graphSolvTime.Point(timeData, ySolvTimeData, simulLen, "load");
    graphSolvTime.Draw();

    graphCost.SetTitle("Method1 MPC Cost", 15, "arial");
    graphCost.SetXlabel("time[sec]");
    graphCost.SetYlabel("cost");
    graphCost.SetXrange(0, 20);
    graphCost.SetXtics(1);
    graphCost.Line(timeData, yCostData, simulLen, "cost");
    graphCost.Draw();

    graphY.SetTitle("Method1 MPC Y ZMP \\\\& CoM", 15, "arial");
    graphY.SetXlabel("time[sec]");
    graphY.SetYlabel("Y[m]");
    graphY.SetXrange(0, 20);
    graphY.SetYrange(-0.3f, 0.3f);
    graphY.SetXtics(1);
    graphY.Dash(timeData, yZmpRef.GetElementsAddr(), simulLen, "zmp ref", "-", "dark-green");
    graphY.Dash(timeData, yZmpMax.GetElementsAddr(), simulLen, "zmp max", "-", "dark-red");
    graphY.Dash(timeData, yZmpMin.GetElementsAddr(), simulLen, "zmp min", "-", "dark-blue");
    graphY.Line(timeData, yPosComData, simulLen, "com pos", "red");
    graphY.Line(timeData, yVelComData, simulLen, "com vel", "dark-violet");
    graphY.Line(timeData, yZmpData, simulLen, "zmp", "blue");
    graphY.Draw();
#endif
}

void MPC_Method2()
{
    PrintTitle("ZMP MPC-Method2 ");
    CdhTimeCheck clock;
    float timeData[simulLen] = { 0, };
    float time = 0;

    // graph vairables
#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
    CdtGnuPlot<> graphSolvTime;
    CdtGnuPlot<> graphCost;
    CdtGnuPlot<> graphY;
#endif
    float yZmpData[simulLen] = { 0, };
    float yComData[simulLen] = { 0, };
    float yCostData[simulLen] = { 0, };
    float ySolvTimeData[simulLen] = { 0, };

    // min a/2*jerk^2 + b/2*(Z-Zref)^2 wrt jerk
    CdtMatrix<2 * Nstep, 2 * Nstep> mGy;
    float w[2 * Nstep];
    CdtVector<2 * Nstep> vGy;

    CdtMatrix<2 * Nstep, 2 * Nstep> mCIy;
    CdtVector<2 * Nstep> vCIy;

    CdtMatrix<2 * Nstep, Nstep> mCEy;
    CdtVector<Nstep> vCEy;

    CdtMatrix<Nstep, Nstep> PzuT;
    CdtMatrix<Nstep, Nstep> Inn;
    CdtVector<2 * Nstep> yX;          // yX = [yJerk; yZ]
    CdtVector<3> yState;

    CdtVector<simulLen> yZmpRef;
    CdtVector<Nstep> yZmpRefK;
    CdtVector<simulLen> yZmpMax;
    CdtVector<Nstep> yZmpMaxK;
    CdtVector<simulLen> yZmpMin;
    CdtVector<Nstep> yZmpMinK;

    CdtQuadProg<2 * Nstep, 2 * Nstep, Nstep> yMPC;
    float yCost = 0;

    MakeRelationMat(smpT);
    Inn.SetIdentity();
    PzuT = Pzu.Transpose();

    /* Weight */
    for (int i = 0; i < Nstep; i++)
    {
        w[i] = R;
    }
    for (int i = Nstep; i < 2 * Nstep; i++)
    {
        w[i] = Q;
    }
    mGy.SetDiagonal(w, sizeof(w));
    vGy.SetBlock(Nstep, yZmpRefK * (-Q));
    yMPC.SetObjectFunc(mGy, vGy);

    /* Equality */
    mCEy.SetBlock(0, 0, PzuT);
    mCEy.SetBlock(Nstep, 0, -Inn);
    vCEy = Pzs * yState;

    /* Inequality */
    mCIy.SetBlock(Nstep, 0, Inn);
    mCIy.SetBlock(Nstep, Nstep, -Inn);
    vCIy.SetBlock(0, -yZmpMinK);
    vCIy.SetBlock(Nstep, yZmpMaxK);

    // make timedata & ref zmp
    MakeTimeData(timeData);
    MakeRefYzmp(yZmpRef);
    yZmpMax = yZmpRef + 0.05f;
    yZmpMin = yZmpRef - 0.05f;

    for (int i = 0; i < (simulLen - Nstep); i++)
    {
        clock.Start();
        /* update Zmp_ref */
        yZmpRefK.SetZero();
        yZmpMinK.SetZero();
        yZmpMaxK.SetZero();

        yZmpRef.GetBlock(i, yZmpRefK);
        yZmpMin.GetBlock(i, yZmpMinK);
        yZmpMax.GetBlock(i, yZmpMaxK);

        /* Update State */
        vGy.SetBlock(Nstep, yZmpRefK * (-Q));
        vCEy = Pzs * yState;
        vCIy.SetBlock(0, -yZmpMinK);
        vCIy.SetBlock(Nstep, yZmpMaxK);

        /* Solve */
        yMPC.UpdateVectorG(vGy);
        if (yMPC.Solve(mCEy, vCEy, mCIy, vCIy, yX))
        {
            Printf("QP Solution Fail!\n");
            return;
        }

        /* desired CoM position, velocity and acceleration */
        yState = A * yState + B * yX(0);   // y(k+1) = A*y(k) + B*y_jerk(k)

        clock.Stop();

        /* save graph data */
        yComData[i] = yState(0);
        yZmpData[i] = C.dot(yState);
        yCostData[i] = yMPC.GetObjectValue();
        ySolvTimeData[i] = (float)clock.GetElapsedTime_msec();
    }
#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
    graphSolvTime.SetTitle("Method2 MPC Solver Time", 15, "arial");
    graphSolvTime.SetXlabel("time[sec]");
    graphSolvTime.SetYlabel("load[ms]");
    graphSolvTime.SetXtics(1);
    graphSolvTime.Point(timeData, ySolvTimeData, simulLen, "load");
    graphSolvTime.Draw();

    graphCost.SetTitle("Method2 MPC Cost", 15, "arial");
    graphCost.SetXlabel("time[sec]");
    graphCost.SetYlabel("cost");
    graphCost.SetXrange(0, 20);
    graphCost.SetXtics(1);
    graphCost.Line(timeData, yCostData, simulLen, "cost");
    graphCost.Draw();

    graphY.SetTitle("Method2 MPC Y ZMP \\\\& CoM", 15, "arial");
    graphY.SetXlabel("time[sec]");
    graphY.SetYlabel("Y[m]");
    graphY.SetXrange(0, 20);
    graphY.SetYrange(-0.3f, 0.3f);
    graphY.SetXtics(1);
    graphY.Dash(timeData, yZmpRef.GetElementsAddr(), simulLen, "zmp ref", "-", "dark-green");
    graphY.Dash(timeData, yZmpMax.GetElementsAddr(), simulLen, "zmp max", "-", "dark-red");
    graphY.Dash(timeData, yZmpMin.GetElementsAddr(), simulLen, "zmp min", "-", "dark-blue");
    graphY.Line(timeData, yComData, simulLen, "com result", "red");
    graphY.Line(timeData, yZmpData, simulLen, "zmp result", "blue");
    graphY.Draw();
#endif
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Dong-hyun Ann, MIT Cheetah3 MPC
#define FORMAT  double
#define IDX_FL  0   // Front Left Leg
#define IDX_BL  1   // Back Left Leg
#define IDX_FR  2   // Front Right Leg
#define IDX_BR  3   // Back Right Leg

const uint8_t mpc_Ns = 13;  // state number
const uint8_t mpc_Nuf = 12; // grf number
const uint8_t mpc_Nu = 25;  // input number
//FORMAT mpc_Ts = 0.05f;        // sampling time
FORMAT mpc_Ts = 0.04f;        // sampling time
const uint8_t mpc_N = 8;    // horizon number

const uint16_t mpc_n = mpc_Nu * mpc_N; // decision variable number
const uint16_t mpc_m = 24 * mpc_N;     // inequality constraint number
const uint16_t mpc_p = mpc_Ns * mpc_N; // equality constaint number

FORMAT mass = 18.0f; //17.141;
FORMAT Ixx = 0.4f;
FORMAT Iyy = 1.041f;
FORMAT Izz = 0.954f;

FORMAT com_height = 0.50f;
FORMAT foot_height = 0.04f;

CdtVector<mpc_Ns, FORMAT> act_state;
CdtVector<mpc_Ns, FORMAT> des_state;
CdtVector<mpc_Ns * mpc_N, FORMAT> des_state_tilde;
CdtVector3<FORMAT> com2foot_des_pos[4];

CdtMatrix<mpc_Ns, mpc_Ns, FORMAT> Ad;
CdtMatrix<mpc_Ns * mpc_N, mpc_Ns, FORMAT> Ad_tilde;

CdtMatrix<mpc_Ns, mpc_Nuf, FORMAT> Bd;
CdtMatrix<mpc_Ns * mpc_N, mpc_Nuf * mpc_N, FORMAT> Bd_tilde;

CdtMatrix3<FORMAT> Ig, inv_Ig;

//FORMAT Q_vec[mpc_Ns] = { 500.0f, 500.0f, 0.01f,  100.0f, 100.0f, 100.0f,  5.0f, 5.0f, 5.0f,  10.0f, 10.0f, 10.0f,  1.0f };
FORMAT Q_vec[mpc_Ns] = { 500.0f, 500.0f, 0.01f,  100.0f, 100.0f, 100.0f,  5.0f, 5.0f, 5.0f,  1.0f, 1.0f, 1.0f,  1.0f };
CdtMatrix<mpc_Ns, mpc_Ns, FORMAT> Q_mat;
CdtMatrix<mpc_Ns * mpc_N, mpc_Ns * mpc_N, FORMAT> Q_tilde;

//FORMAT alpha_mpc_uf = 0.0001f; //0.000001;
FORMAT alpha_mpc_uf = 0.00001;
CdtMatrix<mpc_Nuf, mpc_Nuf, FORMAT> R_mat;
CdtMatrix<mpc_Nuf * mpc_N, mpc_Nuf * mpc_N, FORMAT> R_tilde;

FORMAT mpc_cost = 0;
CdtMatrix<mpc_n, mpc_n, FORMAT> H_mat;
CdtVector<mpc_n, FORMAT> b_vec;
CdtMatrix<mpc_n, mpc_p, FORMAT> mpc_Ce;
CdtVector<mpc_p, FORMAT> mpc_Ce0;
CdtMatrix<mpc_n, mpc_m, FORMAT> mpc_Ci;
CdtVector<mpc_m, FORMAT> mpc_Ci0;
CdtVector<mpc_n, FORMAT> mpc_X;
FORMAT mu = 0.7f;

void ADH_MPC_Init()
{
    CdhTimeCheck clock;
    CdtMatrix<mpc_Ns, mpc_Ns, FORMAT> tmp_Ad;
    CdtMatrix<mpc_Ns, mpc_Nuf, FORMAT> tmp_Bd;
    CdtMatrix<mpc_Ns * mpc_N, mpc_Ns * mpc_N, FORMAT> tmp_x;
    CdtMatrix<6, 3, FORMAT> tmp_Ci;
    CdtMatrix<mpc_Nuf, 24, FORMAT> tmp_Ci_2;
    CdtVector<6, FORMAT> tmp_mCi0[4];
    CdtVector<24, FORMAT> tmp_Ci0_2;

    CdtQuadProg<mpc_n, mpc_m, mpc_p, FORMAT> mpc_qp; //decision variables, inequality constraint, equality constraint

    // Init
    Ig(0, 0) = Ixx;
    Ig(1, 1) = Iyy;
    Ig(2, 2) = Izz;

    inv_Ig = Ig.Inv();

    // destination Foot pos
    //com2foot_des_pos[IDX_FL] << 0.26316f, 0.056f + 0.09041f, -0.477f;
    //com2foot_des_pos[IDX_BL] << -0.26316f, 0.056f + 0.09041f, -com_height;
    //com2foot_des_pos[IDX_FR] << 0.26316f, -(0.056f + 0.09041f), -com_height;
    //com2foot_des_pos[IDX_BR] << -0.26316f, -(0.056f + 0.09041f), -0.477f;

    com2foot_des_pos[IDX_FL] << 0.26316f, 0.056f + 0.09041f, -com_height;
    com2foot_des_pos[IDX_BL] << -0.26316f, 0.056f + 0.09041f, -com_height;
    com2foot_des_pos[IDX_FR] << 0.26316f, -(0.056f + 0.09041f), -com_height;
    com2foot_des_pos[IDX_BR] << -0.26316f, -(0.056f + 0.09041f), -com_height;

    // define Ad
    Ad << 1, 0, 0, 0, 0, 0, mpc_Ts, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, mpc_Ts, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, mpc_Ts, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, mpc_Ts, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, mpc_Ts, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, mpc_Ts, -mpc_Ts * mpc_Ts / 2.0f,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -mpc_Ts,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;

    // Define Ad tilde
    for (int i = 0; i < mpc_N; ++i)
    {
        tmp_Ad = Ad;

        for (int j = 0; j < i; ++j)
        {
            tmp_Ad = tmp_Ad * Ad;
        }

        Ad_tilde.SetBlock(mpc_Ns * i, 0, tmp_Ad);
    }

    // Define Bd
    Bd << -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(0, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(0, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(0, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(0, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(0, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(0, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(0, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(0, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(0, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(0, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(0, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(0, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(1, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(1, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(1, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(1, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(1, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(1, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(1, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(1, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(1, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(1, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(1, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(1, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(2, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(2, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(2, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(2, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(2, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(2, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(2, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(2, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(2, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(2, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(2, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(2, 0))) / 2.0f,
        mpc_Ts * mpc_Ts / (2.0f * mass), 0, 0,
        mpc_Ts * mpc_Ts / (2.0f * mass), 0, 0,
        mpc_Ts * mpc_Ts / (2.0f * mass), 0, 0,
        mpc_Ts * mpc_Ts / (2.0f * mass), 0, 0,
        0, mpc_Ts * mpc_Ts / (2.0f * mass), 0,
        0, mpc_Ts * mpc_Ts / (2.0f * mass), 0,
        0, mpc_Ts * mpc_Ts / (2.0f * mass), 0,
        0, mpc_Ts * mpc_Ts / (2.0f * mass), 0,
        0, 0, mpc_Ts * mpc_Ts / (2.0f * mass),
        0, 0, mpc_Ts * mpc_Ts / (2.0f * mass),
        0, 0, mpc_Ts * mpc_Ts / (2.0f * mass),
        0, 0, mpc_Ts * mpc_Ts / (2.0f * mass),
        -mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(0, 1)), mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(0, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(0, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(0, 1)), mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(0, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(0, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(0, 1)), mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(0, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(0, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(0, 1)), mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(0, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(0, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(1, 1)), mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(1, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(1, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(1, 1)), mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(1, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(1, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(1, 1)), mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(1, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(1, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(1, 1)), mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(1, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(1, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(2, 1)), mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(2, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(2, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(2, 1)), mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(2, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(2, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(2, 1)), mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(2, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(2, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(2, 1)), mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(2, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(2, 0)),
        mpc_Ts / mass, 0, 0,
        mpc_Ts / mass, 0, 0,
        mpc_Ts / mass, 0, 0,
        mpc_Ts / mass, 0, 0,
        0, mpc_Ts / mass, 0,
        0, mpc_Ts / mass, 0,
        0, mpc_Ts / mass, 0,
        0, mpc_Ts / mass, 0,
        0, 0, mpc_Ts / mass,
        0, 0, mpc_Ts / mass,
        0, 0, mpc_Ts / mass,
        0, 0, mpc_Ts / mass,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    // Define Bd tilde
    for (int i = 0; i < mpc_N; ++i)
    {
        tmp_Ad.SetIdentity();
        for (int j = 0; j < i; ++j)
        {
            tmp_Ad = tmp_Ad * Ad;
        }
        tmp_Bd = tmp_Ad * Bd;

        for (int j = i; j < mpc_N; ++j)
        {
            Bd_tilde.SetBlock(mpc_Ns * j, mpc_Nuf * (j - i), tmp_Bd);
        }
    }

    // Define matQ
    for (int i = 0; i < mpc_Ns; ++i)
    {
        Q_mat(i, i) = Q_vec[i];
    }

    // Define Q tilde
    for (int i = 0; i < mpc_N; ++i)
    {
        Q_tilde.SetBlock(mpc_Ns * i, mpc_Ns * i, Q_mat);
    }

    // Define matR
    for (int i = 0; i < mpc_Nuf; ++i)
    {
        R_mat(i, i) = alpha_mpc_uf;
    }

    // Define R tilde
    for (int i = 0; i < mpc_N; ++i)
    {
        R_tilde.SetBlock(mpc_Nuf * i, mpc_Nuf * i, R_mat);
    }

    // Define matH
    H_mat.SetBlock(0, 0, R_tilde);
    H_mat.SetBlock(mpc_Nuf * mpc_N, mpc_Nuf * mpc_N, Q_tilde);

    // Define Desired state
    des_state.SetZero();
    des_state(5) = com_height;
    des_state(12) = g;

    // Define Actual State
    act_state = des_state;

    // Define Desired state tilde
    for (int i = 0; i < mpc_N; ++i)
    {
        des_state_tilde.SetBlock(mpc_Ns * i, des_state);
    }

    // Define vec b
    b_vec.SetBlock(mpc_Nuf * mpc_N, -Q_tilde * des_state_tilde);

    // Equality Constraints
    tmp_x.SetIdentity();
    mpc_Ce.SetBlock(mpc_Nuf * mpc_N, 0, tmp_x);
    mpc_Ce.SetBlock(0, 0, -Bd_tilde.Transpose());

    mpc_Ce0 = -Ad_tilde * act_state;

    // Inequality Constraints
    tmp_Ci << 1, 0, mu,
        -1, 0, mu,
        0, 1, mu,
        0, -1, mu,
        0, 0, 1,
        0, 0, -1;

    tmp_Ci_2.SetBlock(0, 0, tmp_Ci.Transpose());
    tmp_Ci_2.SetBlock(3, 6, tmp_Ci.Transpose());
    tmp_Ci_2.SetBlock(6, 12, tmp_Ci.Transpose());
    tmp_Ci_2.SetBlock(9, 18, tmp_Ci.Transpose());

    for (int i = 0; i < mpc_N; ++i)
    {
        mpc_Ci.SetBlock(mpc_Nuf * i, 24 * i, tmp_Ci_2);
    }

    // tmp
    tmp_mCi0[0] << 0, 0, 0, 0, 0, 0;
    tmp_mCi0[1] << 0, 0, 0, 0, 0, 250;
    tmp_mCi0[2] << 0, 0, 0, 0, 0, 250;
    tmp_mCi0[3] << 0, 0, 0, 0, 0, 0;

    for (int i = 0; i < 4; ++i)
    {
        tmp_Ci0_2.SetBlock(6 * i, tmp_mCi0[i]);
    }

    for (int i = 0; i < mpc_N; ++i)
    {
        mpc_Ci0.SetBlock(24 * i, tmp_Ci0_2);
    }

    mpc_qp.SetObjectFunc(H_mat, b_vec);

    clock.Start();
    mpc_qp.UpdateVectorG(b_vec);
    int rtn = mpc_qp.Solve(mpc_Ce, mpc_Ce0, mpc_Ci, mpc_Ci0, mpc_X);
    clock.Stop();

    printf("MPC_Init\n");
    printf("mcp_X=\n");
    mpc_X.Print('\n');
    printf("load[ms]: %f\n", clock.GetElapsedTime_msec());
    printf("Iteration: %d\n", mpc_qp.GetIteration());
    printf("Cost: %f\n", mpc_qp.GetObjectValue());
    printf("rtn = %d\n", rtn);
}

void ADH_MPC_Process()
{
    FORMAT desState[] = { 0.000000f, 0.000000f, 0.000000f,-0.046316f, 0.000000f, 0.458223f, 0.000000f, 0.000000f, 0.000000f, 0.020518f, 0.000000f, 0.018508f, 9.806000f };
    FORMAT actState[] = { 0.172336f, 0.177815f, 0.000000f,-0.008411f,-0.011410f, 0.510281f, 1.417905f,-1.398757f, 0.476276f, 3.246905f, 0.042101f,-0.159293f, 9.806000f };

    CdhTimeCheck clock;
    CdtMatrix<mpc_Ns, mpc_Ns, FORMAT> tmp_Ad;
    CdtMatrix<mpc_Ns, mpc_Nuf, FORMAT> tmp_Bd;
    CdtVector<6, FORMAT> tmp_mCi0[4];
    CdtVector<24, FORMAT> tmp_Ci0_2;
    CdtVector<4, FORMAT> contact;

    FORMAT f_ext_max = 250;
    FORMAT F_ext_max = 250;

    CdtQuadProg<mpc_n, mpc_m, mpc_p, FORMAT> mpc_qp2;

    //com2foot_des_pos[IDX_FL] << 0.309476f, 0.146410f, -0.458223f;
    //com2foot_des_pos[IDX_BL] <<-0.309476f, 0.146410f, -0.458223f;
    //com2foot_des_pos[IDX_FR] << 0.309476f,-0.146410f, -0.458223f;
    //com2foot_des_pos[IDX_BR] <<-0.309476f,-0.146410f, -0.458223f;

    com2foot_des_pos[IDX_FL] << 0.263160f, 0.146410f, -0.470911f;
    com2foot_des_pos[IDX_BL] << -0.263160f, 0.146410f, -0.5f;
    com2foot_des_pos[IDX_FR] << 0.263160f, -0.146410f, -0.5f;
    com2foot_des_pos[IDX_BR] << -0.263160f, -0.146410f, -0.470911f;

    Bd << -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(0, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(0, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(0, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(0, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(0, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(0, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(0, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(0, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(0, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(0, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(0, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(0, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(1, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(1, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(1, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(1, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(1, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(1, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(1, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(1, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(1, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(1, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(1, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(1, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(2, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(2, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(2, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(2, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(2, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(2, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(2, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(2, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(2, 0))) / 2.0f,
        -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(2, 1))) / 2.0f, (mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(2, 0))) / 2.0f, -(mpc_Ts * mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(2, 0))) / 2.0f,
        mpc_Ts * mpc_Ts / (2.0f * mass), 0, 0,
        mpc_Ts * mpc_Ts / (2.0f * mass), 0, 0,
        mpc_Ts * mpc_Ts / (2.0f * mass), 0, 0,
        mpc_Ts * mpc_Ts / (2.0f * mass), 0, 0,
        0, mpc_Ts * mpc_Ts / (2.0f * mass), 0,
        0, mpc_Ts * mpc_Ts / (2.0f * mass), 0,
        0, mpc_Ts * mpc_Ts / (2.0f * mass), 0,
        0, mpc_Ts * mpc_Ts / (2.0f * mass), 0,
        0, 0, mpc_Ts * mpc_Ts / (2.0f * mass),
        0, 0, mpc_Ts * mpc_Ts / (2.0f * mass),
        0, 0, mpc_Ts * mpc_Ts / (2.0f * mass),
        0, 0, mpc_Ts * mpc_Ts / (2.0f * mass),
        -mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(0, 1)), mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(0, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(0, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(0, 1)), mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(0, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(0, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(0, 1)), mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(0, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(0, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(0, 1)), mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(0, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(0, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(0, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(0, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(1, 1)), mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(1, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(1, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(1, 1)), mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(1, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(1, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(1, 1)), mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(1, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(1, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(1, 1)), mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(1, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(1, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(1, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(1, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_FL](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(2, 1)), mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FL](2) * inv_Ig(2, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FL](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_FL](1) * inv_Ig(2, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BL](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(2, 1)), mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BL](2) * inv_Ig(2, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BL](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_BL](1) * inv_Ig(2, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_FR](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(2, 1)), mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_FR](2) * inv_Ig(2, 0)), -mpc_Ts * (com2foot_des_pos[IDX_FR](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_FR](1) * inv_Ig(2, 0)),
        -mpc_Ts * (com2foot_des_pos[IDX_BR](1) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(2, 1)), mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(2, 2) - com2foot_des_pos[IDX_BR](2) * inv_Ig(2, 0)), -mpc_Ts * (com2foot_des_pos[IDX_BR](0) * inv_Ig(2, 1) - com2foot_des_pos[IDX_BR](1) * inv_Ig(2, 0)),
        mpc_Ts / mass, 0, 0,
        mpc_Ts / mass, 0, 0,
        mpc_Ts / mass, 0, 0,
        mpc_Ts / mass, 0, 0,
        0, mpc_Ts / mass, 0,
        0, mpc_Ts / mass, 0,
        0, mpc_Ts / mass, 0,
        0, mpc_Ts / mass, 0,
        0, 0, mpc_Ts / mass,
        0, 0, mpc_Ts / mass,
        0, 0, mpc_Ts / mass,
        0, 0, mpc_Ts / mass,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    // Define Bd tilde
    for (int i = 0; i < mpc_N; ++i)
    {
        tmp_Ad.SetIdentity();
        for (int j = 0; j < i; ++j)
        {
            tmp_Ad = tmp_Ad * Ad;
        }
        tmp_Bd = tmp_Ad * Bd;

        for (int j = i; j < mpc_N; ++j)
        {
            Bd_tilde.SetBlock(mpc_Ns * j, mpc_Nuf * (j - i), tmp_Bd);
        }
    }

    des_state.SetElement(desState, sizeof(desState));
    act_state.SetElement(actState, sizeof(actState));

    for (int i = 0; i < mpc_N; ++i)
    {
        des_state_tilde.SetBlock(mpc_Ns * i, des_state);
    }

    // Define vec b
    b_vec.SetBlock(mpc_Nuf * mpc_N, -Q_tilde * des_state_tilde);

    // Equality constraints
    mpc_Ce.SetBlock(0, 0, -Bd_tilde.Transpose());
    mpc_Ce0 = -Ad_tilde * act_state;

    // Inequality constraints
    //contact.SetFill(1.0f);
    contact(IDX_FL) = 0;
    contact(IDX_BL) = 1;
    contact(IDX_FR) = 1;
    contact(IDX_BR) = 0;

    for (int i = 0; i < 4; ++i)
    {
        if (contact(i) > 0)
        {
            F_ext_max = f_ext_max;
        }
        else
        {
            F_ext_max = 0;
        }

        tmp_mCi0[i] << 0, 0, 0, 0, 0, F_ext_max;
    }

    for (int i = 0; i < 4; ++i)
    {
        tmp_Ci0_2.SetBlock(6 * i, tmp_mCi0[i]);
    }

    for (int i = 0; i < mpc_N; ++i)
    {
        mpc_Ci0.SetBlock(24 * i, tmp_Ci0_2);
    }

    mpc_qp2.SetObjectFunc(H_mat, b_vec);

    clock.Start();
    mpc_qp2.UpdateVectorG(b_vec);
    int rtn = mpc_qp2.Solve(mpc_Ce, mpc_Ce0, mpc_Ci, mpc_Ci0, mpc_X);
    clock.Stop();
    
    printf("MPC_Process\n");
    printf("mpc_X=\n");
    mpc_X.Print('\n');
    printf("load[ms]: %f\n", clock.GetElapsedTime_msec());
    printf("Iteration: %d\n", mpc_qp2.GetIteration());
    printf("Cost: %f\n", mpc_qp2.GetObjectValue());
    printf("rtn = %d\n", rtn);
}

void ADH_MPC_Test2()
{
    ADH_MPC_Init();

    //for (int i=0; i<50; i++)
        ADH_MPC_Process();
}