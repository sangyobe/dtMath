/*!
\file       dtStorage.tpp
\brief      dtMath, Memory Storage Management class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 04. 05
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTSTORAGE_TPP_
#define DTMATH_DTSTORAGE_TPP_

#include "dtMathMem.h"
#include "dtStorage.h"

#include <assert.h>

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t m_col, typename t_type>
inline Storage<t_row, m_col, t_type>::Storage() : m_elem()
{
}

// Template Partial Specialization Definition
// Case1: row is dynamic, col is static allocation
template <uint16_t t_col, typename t_type>
inline Storage<0, t_col, t_type>::Storage() : m_elem(nullptr), m_row(0)
{
}

template <uint16_t t_col, typename t_type>
inline Storage<0, t_col, t_type>::Storage(uint16_t row)
{
    m_row = row;
    m_elem = MemAllocZeroInit<t_type>(m_row * t_col);
}

template <uint16_t t_col, typename t_type>
inline Storage<0, t_col, t_type>::~Storage()
{
    if (m_elem)
    {
        MemFree<t_type>(m_elem);
        m_elem = nullptr;
        m_row = 0;
    }
}

template <uint16_t t_col, typename t_type>
inline t_type *Storage<0, t_col, t_type>::NewSize(const uint16_t row)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    m_row = row;
    m_elem = MemAllocZeroInit<t_type>(m_row * t_col);
    return m_elem;
}

template <uint16_t t_col, typename t_type>
inline t_type *Storage<0, t_col, t_type>::ReSize(const uint16_t row)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    m_row = row;
    m_elem = MemReAlloc<t_type>(m_elem, m_row * t_col);
    return m_elem;
}

template <uint16_t t_col, typename t_type>
inline void Storage<0, t_col, t_type>::Release()
{
    if (m_elem)
    {
        MemFree<t_type>(m_elem);
        m_elem = nullptr;
        m_row = 0;
    }
}

// Case2: row is static, col is dynamic allocation

// Case3: row and col are dynamic allocation
template <typename t_type>
inline Storage<0, 0, t_type>::Storage() : m_elem(nullptr), m_row(0), m_col(0)
{
}

template <typename t_type>
inline Storage<0, 0, t_type>::Storage(uint16_t row, uint16_t col)
{
    m_row = row;
    m_col = col;
    m_elem = MemAllocZeroInit<t_type>(m_row * m_col);
}

template <typename t_type>
inline Storage<0, 0, t_type>::~Storage()
{
    if (m_elem)
    {
        MemFree<t_type>(m_elem);
        m_elem = nullptr;
    }
}

template <typename t_type>
inline t_type *Storage<0, 0, t_type>::NewSize(const uint16_t row, const uint16_t col)
{
    assert(m_elem == nullptr && "Memory has been allocated");

    m_row = row;
    m_col = col;
    m_elem = MemAllocZeroInit<t_type>(m_row * m_col);
    return m_elem;
}

template <typename t_type>
inline t_type *Storage<0, 0, t_type>::ReSize(const uint16_t row, const uint16_t col)
{
    assert(m_elem != nullptr && "Memory has not been allocated");
    m_row = row;
    m_col = col;
    m_elem = MemReAlloc<t_type>(m_elem, m_row * m_col);
    return m_elem;
}

template <typename t_type>
inline void Storage<0, 0, t_type>::Release()
{
    if (m_elem)
    {
        MemFree<t_type>(m_elem);
        m_elem = nullptr;
        m_row = 0;
        m_col = 0;
    }
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTSTORAGE_TPP_