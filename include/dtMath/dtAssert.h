/*!
\file       dtAssert.h
\brief      art framework assertion routine
\author     Sean Yi, sean_yi@hyundai.com
\date       Last modified on 2024. 04. 01
\version    1.0.0
*/

#ifndef __DTCORE_DTASSERT_H__
#define __DTCORE_DTASSERT_H__

#ifdef NDEBUG
#define DT_NDEBUG
#endif

#ifdef DT_NDEBUG
#define dt_plain_assert(x)
#else
#if defined(_WIN32) || defined(__linux__)
#include <assert.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif
#define dt_plain_assert(x) assert(x)
#endif

#ifndef dt_assert
#define dt_assert(x) dt_plain_assert(x)
#endif

namespace dtCore
{
} // namespace dtMath

#endif // __DTCORE_DTASSERT_H__