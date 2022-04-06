/*!
\file       dhTimeCheck.hpp
\brief      Elapsed time Checker for Algorithm or loop, uint sec/msec/usec
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2018. 05. 28
\version    1.0.0
\copyright  (c) Dong-hyun Lee All rights reserved.
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DH_TIME_CHECK_H
#define DH_TIME_CHECK_H

#if defined(_WIN32)
#if defined(_WIN64)
//define something for Windows (64-bit only)
#else
//define something for Windows (32-bit only)
#endif /* defined(_WIN64) */
//define something for Windows (32-bit or 64-bit)
#define NOMINMAX
#include <Windows.h>
#elif defined(__linux__)
#include <time.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

#if defined(_WIN32)
class CdhTimeCheck
{
public:
    CdhTimeCheck()
    {
        bStart = false;
        // Gets the frequency of the performance counter.
        // The frequency of the performance counter is fixed at system boot
        QueryPerformanceFrequency(&freq);
    }
    ~CdhTimeCheck() {}

    /*! Start Time Check */
    inline void Start(void)
    {
        // Gets start point value of the performance counter
        QueryPerformanceCounter(&startTime);
        bStart = true;
    }

    /*! Stop Time Check */
    inline int Stop()
    {
        // Gets end point value of the performance counter
        QueryPerformanceCounter(&endTime);

        if (!bStart) return -1;

        if (freq.QuadPart != 0)
        {
            elapsedTime_msec = (endTime.QuadPart - startTime.QuadPart) * 1E3
                / static_cast<double>(freq.QuadPart);
            elapsedTime_usec = (endTime.QuadPart - startTime.QuadPart) * 1E6
                / static_cast<double>(freq.QuadPart);
        }

        bStart = false;
        return 0;
    }

    void Reset(void)
    {
        bStart = false;
    }

    inline double GetElapsedTime_sec(void)
    {
        return elapsedTime_msec / 1E3;
    }

    inline double GetElapsedTime_msec(void)
    {
        return elapsedTime_msec;
    }

    inline double GetElapsedTime_usec()
    {
        return elapsedTime_usec;
    }

private:
    bool bStart;
    double elapsedTime_msec;
    double elapsedTime_usec;
    LARGE_INTEGER freq;
    LARGE_INTEGER startTime;
    LARGE_INTEGER endTime;
};
#elif defined(__linux__)
class CdhTimeCheck
{
public:
    CdhTimeCheck()
    {
        bStart = false;
    }
    ~CdhTimeCheck() {}

    /*! Start Time Check */
    inline void Start(void)
    {
        clock_gettime(CLOCK_MONOTONIC, &startTime);
        bStart = true;
    }

    /*! Stop Time Check */
    inline int Stop()
    {
        clock_gettime(CLOCK_MONOTONIC, &endTime);

        if (!bStart) return -1;

        elapsedTime_msec = (endTime.tv_sec - startTime.tv_sec) * 1E3
            + (endTime.tv_nsec - startTime.tv_nsec) / 1E6;
        elapsedTime_usec = (endTime.tv_sec - startTime.tv_sec) * 1E6
            + (endTime.tv_nsec - startTime.tv_nsec) / 1E3;
        
        bStart = false;
        return 0;
    }

    void Reset(void)
    {
        bStart = false;
    }

    inline double GetElapsedTime_sec(void)
    {
        return elapsedTime_msec / 1E3;
    }

    inline double GetElapsedTime_msec(void)
    {
        return elapsedTime_msec;
    }

    inline double GetElapsedTime_usec()
    {
        return elapsedTime_usec;
    }

private:
    bool bStart;
    double elapsedTime_msec;
    double elapsedTime_usec;
    struct timespec startTime;
    struct timespec endTime;
};
#elif defined(ARDUINO)
class CdhTimeCheck
{
public:
    CdhTimeCheck()
    {
        bStart = false;
    }
    ~CdhTimeCheck() {}

    /*! Start Time Check */
    inline void Start(void)
    {
        startTime = micros();
        bStart = true;
    }

    /*! Stop Time Check */
    inline int Stop()
    {
        endTime = micros();

        if (!bStart)
        {
            elapsedTime_usec = -1;
            bStart = false;
            return -1;
        }

        elapsedTime_usec = (double)(endTime - startTime);

        bStart = false;
        return 0;
    }

    void Reset(void)
    {
        bStart = false;
    }

    inline double GetElapsedTime_sec(void)
    {
        return (elapsedTime_usec) / 1E6;
    }

    inline double GetElapsedTime_msec(void)
    {
        return (elapsedTime_usec) / 1E3;
    }

    inline double GetElapsedTime_usec()
    {
        return elapsedTime_usec;
    }

private:
    bool bStart;
    double elapsedTime_msec;
    double elapsedTime_usec;
    uint32_t startTime;
    uint32_t endTime;
};
#endif
#endif /* DH_TIME_CHECK_H */
