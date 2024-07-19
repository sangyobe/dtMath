/*!
\file       dtStorage.h
\brief      dtMath, Memory Storage Management class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 04. 05
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTSTORAGE_H_
#define DTMATH_DTSTORAGE_H_

#include "dtDefine.h"

#if defined(_WIN32) || defined(__linux__)
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#elif defined(ARDUINO)
#include <Arduino.h>
#endif

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type = float>
class Storage
{
private:
    t_type m_elem[t_row * t_col];

public:
    Storage();
    ~Storage() {}

    t_type *NewSize(const uint16_t row, const uint16_t col) { return m_elem; }
    t_type *ReSize(const uint16_t row, const uint16_t col) { return m_elem; }
    t_type *GetAddr() { return m_elem; }
    static constexpr uint16_t GetRow() noexcept { return t_row; }
    static constexpr uint16_t GetCol() noexcept { return t_col; }
};

// Template Partial Specialization
// Case1: row is dynamic, col is static allocation
template <uint16_t t_col, typename t_type>
class Storage<0, t_col, t_type>
{
private:
    t_type *m_elem;
    uint16_t m_row;

public:
    Storage();
    Storage(const uint16_t row);
    ~Storage();

    t_type *GetAddr() { return m_elem; }
    uint16_t GetRow() const noexcept { return m_row; }
    static constexpr uint16_t GetCol() noexcept { return t_col; }
    t_type *NewSize(const uint16_t row);
    t_type *ReSize(const uint16_t row);
    void Release();
};

// template <uint16_t t_row, typename t_type>
// class Storage<t_row, 0, t_type>
// {
// private:
//     t_type *m_elem;
//     uint16_t m_col;

// public:
//     Storage();
//     Storage(const uint16_t row, const uint16_t col);
//     ~Storage();

//     t_type *GetAddr() { return m_elem; }
//     uint16_t GetRow() const noexcept { return t_row; }
//     uint16_t GetCol() const noexcept { return m_col; }
//     t_type *NewSize(const uint16_t row, const uint16_t col);
//     t_type *ReSize(const uint16_t row, const uint16_t col);
//     void Release();
// };

template <typename t_type>
class Storage<0, 0, t_type>
{
private:
    t_type *m_elem;
    uint16_t m_row, m_col;

public:
    Storage();
    Storage(const uint16_t row, const uint16_t col);
    ~Storage();

    t_type *GetAddr() { return m_elem; }
    uint16_t GetRow() const noexcept { return m_row; }
    uint16_t GetCol() const noexcept { return m_col; }
    t_type *NewSize(const uint16_t row, const uint16_t col);
    t_type *ReSize(const uint16_t row, const uint16_t col);
    void Release();
};

} // namespace Math
} // namespace dt

#include "dtStorage.tpp"

#endif // DTMATH_DTSTORAGE_H_