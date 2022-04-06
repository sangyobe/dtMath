/*!
\file       dtCommaInit.hpp
\brief      dtMath, Matrix & Vector comma initializer class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       2020. 10. 21
\version    1.0.0
\warning    Do Not delete this comment for document history!
*/

#pragma once

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

template <uint16_t m_size, typename m_type = float>
class CdtCommaInit
{
private:
    m_type *m_elem;
    int m_idx = 1;

public:
    CdtCommaInit(m_type *elem) : m_elem(elem) {}

    CdtCommaInit& operator ,(const m_type s)
    {
        if (m_idx < m_size) m_elem[m_idx++] = s;
        return *this;
    }
};

#include "dtCommaInit.ipp"