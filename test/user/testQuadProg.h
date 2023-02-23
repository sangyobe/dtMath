#pragma once

void Test_QuadProg();
void EigenQuadProg();
void dtMathQuadProg();
#if defined(_WIN32) || defined(__linux__)
void Test_QuadrupedRobotQP();
void EigenQuadrupedRobotQP();
void dtMathQuadrupedRobotQP();
void dtMathQuadrupedRobotQP_DH();
void dtMathQuadrupedRobotQP_DJ();
void dtMathMoBedQP();
#endif /* defined(_WIN32) || defined(__linux__) */