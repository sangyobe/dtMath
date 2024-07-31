/*!
\file       dtCommaInit.h
\brief      dtMath, Matrix & Vector comma initializer class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history!
*/

#ifndef DTMATH_DTCOMMA_INIT_H_
#define DTMATH_DTCOMMA_INIT_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__) || defined(__APPLE__)
#include <assert.h>
#include <stdint.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

namespace dt
{
namespace Math
{

template <uint16_t t_size, typename t_type = float>
class CommaInit
{
private:
    t_type *m_elem;
    int m_idx = 1;

public:
    CommaInit(t_type *elem) : m_elem(elem) {}

    CommaInit &operator,(const t_type s)
    {
        assert(m_idx < t_size && "Index out of range");

        if (m_idx < t_size) m_elem[m_idx++] = s;
        return *this;
    }
};

// Template Partial Specializtion Declare
template <typename t_type>
class CommaInit<0, t_type>
{
private:
    t_type *m_elem;
    uint16_t t_size;
    int m_idx = 1;

public:
    CommaInit(t_type *elem, uint16_t size) : m_elem(elem), t_size(size) {}

    CommaInit &operator,(const t_type s)
    {
        assert(m_idx < t_size && "Index out of range");

        if (m_idx < t_size) m_elem[m_idx++] = s;
        return *this;
    }
};

} // namespace Math
} // namespace dt

#include "dtCommaInit.tpp"

#endif // DTMATH_DTCOMMA_INIT_H_