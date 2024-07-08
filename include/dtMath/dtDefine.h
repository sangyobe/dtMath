/*!
\file       dtDefine.h
\brief      dtMath, definition of utils
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#pragma once

#if defined(_WIN32) || defined(__linux__)
#define _USE_MATH_DEFINES
#include <cmath>
#include <stdint.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

namespace dtMath
{

#ifndef M_PI
#define M_PI (3.14159265358979323846264338328) /* pi */
#define M_PId (3.141592653589793)
#define M_PIf (3.14159274)
#endif

#define DEG2RADd (double)(0.017453292519943295) // M_PI/180.0
#define RAD2DEGd (double)(57.295779513082323)   // 180.0/M_PI
#define DEG2RADf (float)(0.0174532924)          // M_PI/180.0f
#define RAD2DEGf (float)(57.2957802)            // 180.0f/M_PI

template <typename T>
inline T dtPI()
{
    return T(3.14159265358979323846264338328);
}

template <typename T>
inline T dt2PI() { return T(6.28318530717958647692528676656); } // 2 * PI

template <typename T>
inline T dtPI2() { return T(1.57079632679489661923132169164); } // PI / 2

template <typename T>
inline T dtDEG2RAD() { return T(0.01745329251994329576923690768489); }

template <typename T>
inline T dtRAD2DEG() { return T(57.295779513082320876798154814105); } // 180.0/M_PI

#define X_AXIS (uint16_t)(0)
#define Y_AXIS (uint16_t)(1)
#define Z_AXIS (uint16_t)(2)

#define AXIS3(A1, A2, A3) (uint16_t)(A3 << 8 | A2 << 4 | A1)
#define AXIS2(A1, A2) (uint16_t)(A2 << 4 | A1)
#define AXIS1(A1) (uint16_t)(A1)

} // namespace dtMath


/**
 * dt_assert - macro with one argument that is used inside this for assertions. By default, it is basically defined to be assert, which aborts the program if the assertion is violated. Redefine this macro if you want to do something else, like throwing an exception.
 */
// #include <cstdlib>     // for abort
// #include <iostream>    // for std::cerr
// #include <stdlib.h>    // for free
// #include <execinfo.h>  // for backtrace
// #define dt_assert(var) \
//     do { \
//         if (!(bool)(var)) { \
//             std::cerr << "assertion failed: " << #var << " in function " << __PRETTY_FUNCTION__ << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
//             void *callstack[128]; \
//             int i, nr_frames; \
//             char **strs; \
//             nr_frames = backtrace(callstack, sizeof(callstack)/sizeof(void *)); \
//             strs = backtrace_symbols(callstack, nr_frames); \
//             std::cerr << "[call stack trace]: \n"; \
//             for (i = 0; i < nr_frames; i++) { \
//                 std::cerr << "\t(" << i << ") " << strs[i] << std::endl; \
//             } \
//             free(strs); \
//             abort(); \
//         } \
//     } while(false)

#include "dtAssert.h"