/*!
\file       dtMatrix.tpp
\brief      dtMath, General Matrix(m x n) class
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Joonhee Jo, allusivejune@gmail.com
\author     Who is next author?
\date       Last modified on 2023. 05. 02
\version    1.1.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATRIX_TPP_
#define DTMATH_DTMATRIX_TPP_

#include "dtMatrix.h"

#include <cassert>

namespace dt
{
namespace Math
{

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type>::Matrix() : m_elem()
{
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type>::Matrix(const t_type *element)
{
    memcpy(m_elem, element, sizeof(t_type) * t_row * t_col);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type>::Matrix(const t_type *element, const size_t n_byte)
{
    size_t matSz = sizeof(t_type) * t_row * t_col;

    if (matSz <= n_byte)
    {
        memcpy(m_elem, element, matSz);
    }
    else
    {
        memset(m_elem, 0, matSz);
        memcpy(m_elem, element, n_byte);
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type>::Matrix(const char c, const t_type *element, const size_t n_byte)
{
    if (c == 'a')
    {
        size_t matSz = sizeof(t_type) * t_row * t_col;

        if (matSz <= n_byte)
        {
            memcpy(m_elem, element, matSz);
        }
        else
        {
            memset(m_elem, 0, matSz);
            memcpy(m_elem, element, n_byte);
        }
    }

    else if (c == 'd')
    {
        memset(m_elem, 0, sizeof(t_type) * t_row * t_col);

        uint16_t num = (t_row > t_col) ? t_col : t_row;
        uint16_t offset = t_col + 1;
        uint16_t elemNum = (uint16_t)(n_byte / sizeof(t_type));

        num = (num > elemNum) ? elemNum : num;

        for (uint16_t i = 0; i < num; i++)
            m_elem[i * offset] = element[i];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type>::Matrix(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetZero()
{
    memset(m_elem, 0, sizeof(t_type) * t_row * t_col);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetIdentity()
{
    memset(m_elem, 0, sizeof(t_type) * t_row * t_col);
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = num >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * offset] = 1;
        m_elem[(irow + 1) * offset] = 1;
        m_elem[(irow + 2) * offset] = 1;
        m_elem[(irow + 3) * offset] = 1;
    }

    for (cnt = num % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * offset] = 1;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetDiagonal(const t_type *element, const size_t n_byte)
{
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;
    uint16_t elemNum = (uint16_t)(n_byte / sizeof(t_type));
    uint16_t cnt;
    uint16_t irow = 0;

    num = (num > elemNum) ? elemNum : num;

    for (uint16_t i = 0; i < num; i++)
        m_elem[i * offset] = element[i];

    for (cnt = num >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * offset] = element[irow];
        m_elem[(irow + 1) * offset] = element[irow + 1];
        m_elem[(irow + 2) * offset] = element[irow + 2];
        m_elem[(irow + 3) * offset] = element[irow + 3];
    }

    for (cnt = num % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * offset] = element[irow];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetDiagonal(const t_type value)
{
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = num >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * offset] = value;
        m_elem[(irow + 1) * offset] = value;
        m_elem[(irow + 2) * offset] = value;
        m_elem[(irow + 3) * offset] = value;
    }

    for (cnt = num % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * offset] = value;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t row>
inline void Matrix<t_row, t_col, t_type>::SetDiagonal(const Vector<row, t_type> &v)
{
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;
    uint16_t cnt;
    uint16_t irow = 0;

    num = (num > row) ? row : num;

    for (cnt = num >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * offset] = v.m_elem[irow];
        m_elem[(irow + 1) * offset] = v.m_elem[irow + 1];
        m_elem[(irow + 2) * offset] = v.m_elem[irow + 2];
        m_elem[(irow + 3) * offset] = v.m_elem[irow + 3];
    }

    for (cnt = num % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * offset] = v.m_elem[irow];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetDiagonal(const Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;
    uint16_t cnt;
    uint16_t irow = 0;

    num = (num > v.m_row) ? v.m_row : num;

    for (cnt = num >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * offset] = v.m_elem[irow];
        m_elem[(irow + 1) * offset] = v.m_elem[irow + 1];
        m_elem[(irow + 2) * offset] = v.m_elem[irow + 2];
        m_elem[(irow + 3) * offset] = v.m_elem[irow + 3];
    }

    for (cnt = num % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * offset] = v.m_elem[irow];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetDiagonal(const Vector3<t_type, 3> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 3) ? 3 : num;

    switch (num)
    {
    case 1:
        m_elem[0] = v.m_elem[0];
        break;
    case 2:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        break;
    case 3:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        m_elem[2 * offset] = v.m_elem[2];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetDiagonal(const Vector4<t_type, 4> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 4) ? 4 : num;

    switch (num)
    {
    case 1:
        m_elem[0] = v.m_elem[0];
        break;
    case 2:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        break;
    case 3:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        m_elem[2 * offset] = v.m_elem[2];
        break;
    case 4:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        m_elem[2 * offset] = v.m_elem[2];
        m_elem[3 * offset] = v.m_elem[3];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetDiagonal(const Vector6<t_type, 6> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 6) ? 6 : num;

    switch (num)
    {
    case 1:
        m_elem[0] = v.m_elem[0];
        break;
    case 2:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        break;
    case 3:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        m_elem[2 * offset] = v.m_elem[2];
        break;
    case 4:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        m_elem[2 * offset] = v.m_elem[2];
        m_elem[3 * offset] = v.m_elem[3];
        break;
    case 5:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        m_elem[2 * offset] = v.m_elem[2];
        m_elem[3 * offset] = v.m_elem[3];
        m_elem[4 * offset] = v.m_elem[4];
        break;
    case 6:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        m_elem[2 * offset] = v.m_elem[2];
        m_elem[3 * offset] = v.m_elem[3];
        m_elem[4 * offset] = v.m_elem[4];
        m_elem[5 * offset] = v.m_elem[5];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetDiagonal(const Quaternion<t_type, 4> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 4) ? 4 : num;

    switch (num)
    {
    case 1:
        m_elem[0] = v.m_elem[0];
        break;
    case 2:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        break;
    case 3:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        m_elem[2 * offset] = v.m_elem[2];
        break;
    case 4:
        m_elem[0] = v.m_elem[0];
        m_elem[offset] = v.m_elem[1];
        m_elem[2 * offset] = v.m_elem[2];
        m_elem[3 * offset] = v.m_elem[3];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetFill(const t_type value)
{
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 2u; cnt > 0u; cnt--, i += 4)
    {
        m_elem[i] = value;
        m_elem[i + 1] = value;
        m_elem[i + 2] = value;
        m_elem[i + 3] = value;
    }

    for (cnt = (t_row * t_col) % 4u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] = value;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetFillRow(const uint16_t idxRow, const t_type value)
{
    assert(t_row > idxRow && "Index out of range");

    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = t_col >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        m_elem[idxRow * t_col + icol] = value;
        m_elem[idxRow * t_col + icol + 1] = value;
        m_elem[idxRow * t_col + icol + 2] = value;
        m_elem[idxRow * t_col + icol + 3] = value;
    }

    for (cnt = t_col % 4u; cnt > 0u; cnt--, icol++)
    {
        m_elem[idxRow * t_col + icol] = value;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetFillCol(const uint16_t idxCol, const t_type value)
{
    assert(t_col > idxCol && "Index out of range");

    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[irow * t_col + idxCol] = value;
        m_elem[(irow + 1) * t_col + idxCol] = value;
        m_elem[(irow + 2) * t_col + idxCol] = value;
        m_elem[(irow + 3) * t_col + idxCol] = value;
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[irow * t_col + idxCol] = value;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetElement(const t_type *element, const size_t n_byte)
{
    size_t matSz = sizeof(t_type) * t_row * t_col;

    if (matSz <= n_byte)
    {
        memcpy(m_elem, element, matSz);
    }
    else
    {
        memset(m_elem, 0, matSz);
        memcpy(m_elem, element, n_byte);
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t row, uint16_t col>
inline void Matrix<t_row, t_col, t_type>::SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix<row, col, t_type> &m, const uint16_t jdxRow, const uint16_t jdxCol, const int16_t rSize, const int16_t cSize)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(row > jdxRow && "Index out of range");
    assert(col > jdxCol && "Index out of range");
    assert((row - jdxRow) >= rSize && "Over size");
    assert((col - jdxCol) >= cSize && "Over size");

    uint16_t irow, icol, cnt;
    uint16_t rowSize = t_row - idxRow;
    uint16_t colSize = t_col - idxCol;

    if (rowSize > rSize) rowSize = rSize;
    if (colSize > cSize) colSize = cSize;

    for (uint16_t irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            m_elem[(irow + idxRow) * t_col + idxCol + icol] = m.m_elem[(jdxRow + irow) * col + jdxCol + icol];
            m_elem[(irow + idxRow) * t_col + idxCol + icol + 1] = m.m_elem[(jdxRow + irow) * col + jdxCol + icol + 1];
            m_elem[(irow + idxRow) * t_col + idxCol + icol + 2] = m.m_elem[(jdxRow + irow) * col + jdxCol + icol + 2];
            m_elem[(irow + idxRow) * t_col + idxCol + icol + 3] = m.m_elem[(jdxRow + irow) * col + jdxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            m_elem[(irow + idxRow) * t_col + idxCol + icol] = m.m_elem[(jdxRow + irow) * col + jdxCol + icol];
        }
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix<0, 0, t_type> &m, const uint16_t jdxRow, const uint16_t jdxCol, int16_t rSize, int16_t cSize)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(m.m_row > jdxRow && "Index out of range");
    assert(m.m_col > jdxCol && "Index out of range");
    assert((m.m_row - jdxRow) >= rSize && "Over size");
    assert((m.m_col - jdxCol) >= cSize && "Over size");

    uint16_t irow, icol, cnt;
    uint16_t rowSize = t_row - idxRow;
    uint16_t colSize = t_col - idxCol;

    if (rSize < 0) rSize = m.m_row;
    if (cSize < 0) cSize = m.m_col;
    if (rowSize > rSize) rowSize = rSize;
    if (colSize > cSize) colSize = cSize;

    for (uint16_t irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            m_elem[(irow + idxRow) * t_col + idxCol + icol] = m.m_elem[(jdxRow + irow) * m.m_col + jdxCol + icol];
            m_elem[(irow + idxRow) * t_col + idxCol + icol + 1] = m.m_elem[(jdxRow + irow) * m.m_col + jdxCol + icol + 1];
            m_elem[(irow + idxRow) * t_col + idxCol + icol + 2] = m.m_elem[(jdxRow + irow) * m.m_col + jdxCol + icol + 2];
            m_elem[(irow + idxRow) * t_col + idxCol + icol + 3] = m.m_elem[(jdxRow + irow) * m.m_col + jdxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            m_elem[(irow + idxRow) * t_col + idxCol + icol] = m.m_elem[(jdxRow + irow) * m.m_col + jdxCol + icol];
        }
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Matrix3<t_type, 3, 3> &m)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");

    m_elem[idxRow * t_col + idxCol + 0] = m.m_elem[0];
    m_elem[idxRow * t_col + idxCol + 1] = m.m_elem[1];
    m_elem[idxRow * t_col + idxCol + 2] = m.m_elem[2];

    m_elem[(1 + idxRow) * t_col + idxCol + 0] = m.m_elem[3];
    m_elem[(1 + idxRow) * t_col + idxCol + 1] = m.m_elem[4];
    m_elem[(1 + idxRow) * t_col + idxCol + 2] = m.m_elem[5];

    m_elem[(2 + idxRow) * t_col + idxCol + 0] = m.m_elem[6];
    m_elem[(2 + idxRow) * t_col + idxCol + 1] = m.m_elem[7];
    m_elem[(2 + idxRow) * t_col + idxCol + 2] = m.m_elem[8];
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetBlock(const uint16_t idxRow, const uint16_t idxCol, const Rotation<t_type, 3, 3> &m)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");

    m_elem[idxRow * t_col + idxCol + 0] = m.m_elem[0];
    m_elem[idxRow * t_col + idxCol + 1] = m.m_elem[1];
    m_elem[idxRow * t_col + idxCol + 2] = m.m_elem[2];

    m_elem[(1 + idxRow) * t_col + idxCol + 0] = m.m_elem[3];
    m_elem[(1 + idxRow) * t_col + idxCol + 1] = m.m_elem[4];
    m_elem[(1 + idxRow) * t_col + idxCol + 2] = m.m_elem[5];

    m_elem[(2 + idxRow) * t_col + idxCol + 0] = m.m_elem[6];
    m_elem[(2 + idxRow) * t_col + idxCol + 1] = m.m_elem[7];
    m_elem[(2 + idxRow) * t_col + idxCol + 2] = m.m_elem[8];
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t col>
inline void Matrix<t_row, t_col, t_type>::SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Vector<col, t_type> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(col > jdx && "Index out of range");
    assert((col - jdx) >= size && "Over size");

    uint16_t colSize = t_col - idxCol;
    uint16_t cnt;
    uint16_t icol = 0;

    if (colSize > size) colSize = size;

    for (cnt = colSize >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        m_elem[idxRow * t_col + idxCol + icol] = v.m_elem[jdx + icol];
        m_elem[idxRow * t_col + idxCol + icol + 1] = v.m_elem[jdx + icol + 1];
        m_elem[idxRow * t_col + idxCol + icol + 2] = v.m_elem[jdx + icol + 2];
        m_elem[idxRow * t_col + idxCol + icol + 3] = v.m_elem[jdx + icol + 3];
    }

    for (cnt = colSize % 4u; cnt > 0u; cnt--, icol++)
    {
        m_elem[idxRow * t_col + idxCol + icol] = v.m_elem[jdx + icol];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Vector<0, t_type> &v, const uint16_t jdx, int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(v.m_row > jdx && "Index out of range");
    assert((v.m_row - jdx) >= size && "Over size");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t colSize = t_col - idxCol;
    uint16_t cnt;
    uint16_t icol = 0;

    if (size < 0) size = v.m_row;
    if (colSize > size) colSize = size;

    for (cnt = colSize >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        m_elem[idxRow * t_col + idxCol + icol] = v.m_elem[jdx + icol];
        m_elem[idxRow * t_col + idxCol + icol + 1] = v.m_elem[jdx + icol + 1];
        m_elem[idxRow * t_col + idxCol + icol + 2] = v.m_elem[jdx + icol + 2];
        m_elem[idxRow * t_col + idxCol + icol + 3] = v.m_elem[jdx + icol + 3];
    }

    for (cnt = colSize % 4u; cnt > 0u; cnt--, icol++)
    {
        m_elem[idxRow * t_col + idxCol + icol] = v.m_elem[jdx + icol];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Vector3<t_type, 3> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(3 > jdx && "Index out of range");
    assert((3 - jdx) >= size && "Over size");

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[jdx + 1];
        break;
    case 3:
    default:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[0];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[1];
        m_elem[idxRow * t_col + idxCol + 2] = v.m_elem[2];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Vector4<t_type, 4> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "Over size");

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[jdx + 1];
        break;
    case 3:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[jdx + 1];
        m_elem[idxRow * t_col + idxCol + 2] = v.m_elem[jdx + 2];
        break;
    case 4:
    default:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[0];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[1];
        m_elem[idxRow * t_col + idxCol + 2] = v.m_elem[2];
        m_elem[idxRow * t_col + idxCol + 3] = v.m_elem[3];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Quaternion<t_type, 4> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "Over size");

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[jdx + 1];
        break;
    case 3:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[jdx + 1];
        m_elem[idxRow * t_col + idxCol + 2] = v.m_elem[jdx + 2];
        break;
    case 4:
    default:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[0];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[1];
        m_elem[idxRow * t_col + idxCol + 2] = v.m_elem[2];
        m_elem[idxRow * t_col + idxCol + 3] = v.m_elem[3];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const Vector6<t_type, 6> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(6 > jdx && "Index out of range");
    assert((6 - jdx) >= size && "Over size");

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[jdx + 1];
        break;
    case 3:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[jdx + 1];
        m_elem[idxRow * t_col + idxCol + 2] = v.m_elem[jdx + 2];
        break;
    case 4:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[jdx + 1];
        m_elem[idxRow * t_col + idxCol + 2] = v.m_elem[jdx + 2];
        m_elem[idxRow * t_col + idxCol + 3] = v.m_elem[jdx + 3];
        break;
    case 5:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[jdx + 1];
        m_elem[idxRow * t_col + idxCol + 2] = v.m_elem[jdx + 2];
        m_elem[idxRow * t_col + idxCol + 3] = v.m_elem[jdx + 3];
        m_elem[idxRow * t_col + idxCol + 4] = v.m_elem[jdx + 4];
        break;
    case 6:
    default:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[0];
        m_elem[idxRow * t_col + idxCol + 1] = v.m_elem[1];
        m_elem[idxRow * t_col + idxCol + 2] = v.m_elem[2];
        m_elem[idxRow * t_col + idxCol + 3] = v.m_elem[3];
        m_elem[idxRow * t_col + idxCol + 4] = v.m_elem[4];
        m_elem[idxRow * t_col + idxCol + 5] = v.m_elem[5];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetRowVec(const uint16_t idxRow, const uint16_t idxCol, const t_type *v, const size_t n_byte)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");

    uint16_t col = (uint16_t)(n_byte / sizeof(t_type));
    uint16_t colSize = t_col - idxCol;
    uint16_t cnt;
    uint16_t icol = 0;

    if (colSize > col) colSize = col;

    for (cnt = colSize >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        m_elem[idxRow * t_col + idxCol + icol] = v[icol];
        m_elem[idxRow * t_col + idxCol + icol + 1] = v[icol + 1];
        m_elem[idxRow * t_col + idxCol + icol + 2] = v[icol + 2];
        m_elem[idxRow * t_col + idxCol + icol + 3] = v[icol + 3];
    }

    for (cnt = colSize % 4u; cnt > 0u; cnt--, icol++)
    {
        m_elem[idxRow * t_col + idxCol + icol] = v[icol];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t row>
inline void Matrix<t_row, t_col, t_type>::SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Vector<row, t_type> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(row > jdx && "Index out of range");
    assert((row - jdx) >= size && "Over size");

    uint16_t rowSize = t_row - idxRow;
    uint16_t cnt;
    uint16_t irow = 0;

    if (rowSize > size) rowSize = size;

    for (cnt = rowSize >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[(idxRow + irow) * t_col + idxCol] = v.m_elem[jdx + irow];
        m_elem[(idxRow + irow + 1) * t_col + idxCol] = v.m_elem[jdx + irow + 1];
        m_elem[(idxRow + irow + 2) * t_col + idxCol] = v.m_elem[jdx + irow + 2];
        m_elem[(idxRow + irow + 3) * t_col + idxCol] = v.m_elem[jdx + irow + 3];
    }

    for (cnt = rowSize % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[(idxRow + irow) * t_col + idxCol] = v.m_elem[jdx + irow];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Vector<0, t_type> &v, const uint16_t jdx, int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(v.m_row > jdx && "Index out of range");
    assert((v.m_row - jdx) >= size && "Over size");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t rowSize = t_row - idxRow;
    uint16_t cnt;
    uint16_t irow = 0;

    if (size < 0) size = v.m_row;
    if (rowSize > size) rowSize = size;

    for (cnt = rowSize >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[(idxRow + irow) * t_col + idxCol] = v.m_elem[jdx + irow];
        m_elem[(idxRow + irow + 1) * t_col + idxCol] = v.m_elem[jdx + irow + 1];
        m_elem[(idxRow + irow + 2) * t_col + idxCol] = v.m_elem[jdx + irow + 2];
        m_elem[(idxRow + irow + 3) * t_col + idxCol] = v.m_elem[jdx + irow + 3];
    }

    for (cnt = rowSize % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[(idxRow + irow) * t_col + idxCol] = v.m_elem[jdx + irow];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Vector3<t_type, 3> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(3 > jdx && "Index out of range");
    assert((3 - jdx) >= size && "Over size");

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[jdx + 1];
        break;
    case 3:
    default:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[0];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[1];
        m_elem[(idxRow + 2) * t_col + idxCol] = v.m_elem[2];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Vector4<t_type, 4> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "Over size");

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[jdx + 1];
        break;
    case 3:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[jdx + 1];
        m_elem[(idxRow + 2) * t_col + idxCol] = v.m_elem[jdx + 2];
        break;
    case 4:
    default:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[0];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[1];
        m_elem[(idxRow + 2) * t_col + idxCol] = v.m_elem[2];
        m_elem[(idxRow + 3) * t_col + idxCol] = v.m_elem[3];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Quaternion<t_type, 4> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "Over size");

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[jdx + 1];
        break;
    case 3:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[jdx + 1];
        m_elem[(idxRow + 2) * t_col + idxCol] = v.m_elem[jdx + 2];
        break;
    case 4:
    default:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[0];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[1];
        m_elem[(idxRow + 2) * t_col + idxCol] = v.m_elem[2];
        m_elem[(idxRow + 3) * t_col + idxCol] = v.m_elem[3];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetColVec(const uint16_t idxRow, const uint16_t idxCol, const Vector6<t_type, 6> &v, const uint16_t jdx, const int16_t size)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(6 > jdx && "Index out of range");
    assert((6 - jdx) >= size && "Over size");

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        break;
    case 2:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[jdx + 1];
        break;
    case 3:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[jdx + 1];
        m_elem[(idxRow + 2) * t_col + idxCol] = v.m_elem[jdx + 2];
        break;
    casejdx4:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[jdx + 1];
        m_elem[(idxRow + 2) * t_col + idxCol] = v.m_elem[jdx + 2];
        m_elem[(idxRow + 3) * t_col + idxCol] = v.m_elem[jdx + 3];
        break;
    case 5:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[jdx];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[jdx + 1];
        m_elem[(idxRow + 2) * t_col + idxCol] = v.m_elem[jdx + 2];
        m_elem[(idxRow + 3) * t_col + idxCol] = v.m_elem[jdx + 3];
        m_elem[(idxRow + 4) * t_col + idxCol] = v.m_elem[jdx + 4];
        break;
    case 6:
    default:
        m_elem[idxRow * t_col + idxCol] = v.m_elem[0];
        m_elem[(idxRow + 1) * t_col + idxCol] = v.m_elem[1];
        m_elem[(idxRow + 2) * t_col + idxCol] = v.m_elem[2];
        m_elem[(idxRow + 3) * t_col + idxCol] = v.m_elem[3];
        m_elem[(idxRow + 4) * t_col + idxCol] = v.m_elem[4];
        m_elem[(idxRow + 5) * t_col + idxCol] = v.m_elem[5];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetColVec(const uint16_t idxRow, const uint16_t idxCol, const t_type *v, const size_t n_byte)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");

    uint16_t row = n_byte / sizeof(t_type);
    uint16_t rowSize = t_row - idxRow;
    uint16_t cnt;
    uint16_t irow = 0;

    if (rowSize > row) rowSize = row;

    for (cnt = rowSize >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        m_elem[(idxRow + irow) * t_col + idxCol] = v[irow];
        m_elem[(idxRow + irow + 1) * t_col + idxCol] = v[irow + 1];
        m_elem[(idxRow + irow + 2) * t_col + idxCol] = v[irow + 2];
        m_elem[(idxRow + irow + 3) * t_col + idxCol] = v[irow + 3];
    }

    for (cnt = rowSize % 4u; cnt > 0u; cnt--, irow++)
    {
        m_elem[(idxRow + irow) * t_col + idxCol] = v[irow];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetSwapRowVec(const uint16_t idxRow1, const uint16_t idxRow2)
{
    assert(t_row > idxRow1 && "Index out of range");
    assert(t_row > idxRow2 && "Index out of range");

    t_type tmpVec[t_col];
    uint16_t cnt;
    uint16_t icol = 0;

    for (cnt = t_col >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        tmpVec[icol] = m_elem[idxRow1 * t_col + icol];
        m_elem[idxRow1 * t_col + icol] = m_elem[idxRow2 * t_col + icol];
        m_elem[idxRow2 * t_col + icol] = tmpVec[icol];

        tmpVec[icol + 1] = m_elem[idxRow1 * t_col + icol + 1];
        m_elem[idxRow1 * t_col + icol + 1] = m_elem[idxRow2 * t_col + icol + 1];
        m_elem[idxRow2 * t_col + icol + 1] = tmpVec[icol + 1];

        tmpVec[icol + 2] = m_elem[idxRow1 * t_col + icol + 2];
        m_elem[idxRow1 * t_col + icol + 2] = m_elem[idxRow2 * t_col + icol + 2];
        m_elem[idxRow2 * t_col + icol + 2] = tmpVec[icol + 2];

        tmpVec[icol + 3] = m_elem[idxRow1 * t_col + icol + 3];
        m_elem[idxRow1 * t_col + icol + 3] = m_elem[idxRow2 * t_col + icol + 3];
        m_elem[idxRow2 * t_col + icol + 3] = tmpVec[icol + 3];
    }

    for (cnt = t_col % 4u; cnt > 0u; cnt--, icol++)
    {
        tmpVec[icol] = m_elem[idxRow1 * t_col + icol];
        m_elem[idxRow1 * t_col + icol] = m_elem[idxRow2 * t_col + icol];
        m_elem[idxRow2 * t_col + icol] = tmpVec[icol];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::SetSwapColVec(const uint16_t idxCol1, const uint16_t idxCol2)
{
    assert(t_col > idxCol1 && "Index out of range");
    assert(t_col > idxCol2 && "Index out of range");

    t_type tmpVec[t_row];
    uint16_t cnt;
    uint16_t irow = 0;

    for (cnt = t_row >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        tmpVec[irow] = m_elem[irow * t_col + idxCol1];
        m_elem[irow * t_col + idxCol1] = m_elem[irow * t_col + idxCol2];
        m_elem[irow * t_col + idxCol2] = tmpVec[irow];

        tmpVec[irow + 1] = m_elem[(irow + 1) * t_col + idxCol1];
        m_elem[(irow + 1) * t_col + idxCol1] = m_elem[(irow + 1) * t_col + idxCol2];
        m_elem[(irow + 1) * t_col + idxCol2] = tmpVec[irow + 1];

        tmpVec[irow + 2] = m_elem[(irow + 2) * t_col + idxCol1];
        m_elem[(irow + 2) * t_col + idxCol1] = m_elem[(irow + 2) * t_col + idxCol2];
        m_elem[(irow + 2) * t_col + idxCol2] = tmpVec[irow + 2];

        tmpVec[irow + 3] = m_elem[(irow + 3) * t_col + idxCol1];
        m_elem[(irow + 3) * t_col + idxCol1] = m_elem[(irow + 3) * t_col + idxCol2];
        m_elem[(irow + 3) * t_col + idxCol2] = tmpVec[irow + 3];
    }

    for (cnt = t_row % 4u; cnt > 0u; cnt--, irow++)
    {
        tmpVec[irow] = m_elem[irow * t_col + idxCol1];
        m_elem[irow * t_col + idxCol1] = m_elem[irow * t_col + idxCol2];
        m_elem[irow * t_col + idxCol2] = tmpVec[irow];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline const t_type *const Matrix<t_row, t_col, t_type>::GetElementsAddr() const
{
    return m_elem;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t row>
inline void Matrix<t_row, t_col, t_type>::GetDiagonal(Vector<row, t_type> &v)
{
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;
    uint16_t cnt;
    uint16_t irow = 0;

    num = (num > row) ? row : num;

    for (cnt = num >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        v.m_elem[irow] = m_elem[irow * offset];
        v.m_elem[irow + 1] = m_elem[(irow + 1) * offset];
        v.m_elem[irow + 2] = m_elem[(irow + 2) * offset];
        v.m_elem[irow + 3] = m_elem[(irow + 3) * offset];
    }

    for (cnt = num % 4u; cnt > 0u; cnt--, irow++)
    {
        v.m_elem[irow] = m_elem[irow * offset];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::GetDiagonal(Vector<0, t_type> &v)
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;
    uint16_t cnt;
    uint16_t irow = 0;

    num = (num > v.m_row) ? v.m_row : num;

    for (cnt = num >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        v.m_elem[irow] = m_elem[irow * offset];
        v.m_elem[irow + 1] = m_elem[(irow + 1) * offset];
        v.m_elem[irow + 2] = m_elem[(irow + 2) * offset];
        v.m_elem[irow + 3] = m_elem[(irow + 3) * offset];
    }

    for (cnt = num % 4u; cnt > 0u; cnt--, irow++)
    {
        v.m_elem[irow] = m_elem[irow * offset];
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::GetDiagonal(Vector3<t_type, 3> &v)
{
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 3) ? 3 : num;

    switch (num)
    {
    case 1:
        v.m_elem[0] = m_elem[0];
        break;
    case 2:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        break;
    case 3:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        v.m_elem[2] = m_elem[2 * offset];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::GetDiagonal(Vector4<t_type, 4> &v)
{
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 4) ? 4 : num;

    switch (num)
    {
    case 1:
        v.m_elem[0] = m_elem[0];
        break;
    case 2:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        break;
    case 3:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        v.m_elem[2] = m_elem[2 * offset];
        break;
    case 4:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        v.m_elem[2] = m_elem[2 * offset];
        v.m_elem[3] = m_elem[3 * offset];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::GetDiagonal(Vector6<t_type, 6> &v)
{
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 6) ? 6 : num;

    switch (num)
    {
    case 1:
        v.m_elem[0] = m_elem[0];
        break;
    case 2:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        break;
    case 3:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        v.m_elem[2] = m_elem[2 * offset];
        break;
    case 4:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        v.m_elem[2] = m_elem[2 * offset];
        v.m_elem[3] = m_elem[3 * offset];
        break;
    case 5:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        v.m_elem[2] = m_elem[2 * offset];
        v.m_elem[3] = m_elem[3 * offset];
        v.m_elem[4] = m_elem[4 * offset];
        break;
    case 6:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        v.m_elem[2] = m_elem[2 * offset];
        v.m_elem[3] = m_elem[3 * offset];
        v.m_elem[4] = m_elem[4 * offset];
        v.m_elem[5] = m_elem[5 * offset];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::GetDiagonal(Quaternion<t_type, 4> &v)
{
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 4) ? 4 : num;

    switch (num)
    {
    case 1:
        v.m_elem[0] = m_elem[0];
        break;
    case 2:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        break;
    case 3:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        v.m_elem[2] = m_elem[2 * offset];
        break;
    case 4:
        v.m_elem[0] = m_elem[0];
        v.m_elem[1] = m_elem[offset];
        v.m_elem[2] = m_elem[2 * offset];
        v.m_elem[3] = m_elem[3 * offset];
        break;
    }
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t row>
inline Vector<row, t_type> Matrix<t_row, t_col, t_type>::GetDiagonal()
{
    t_type elem[row]{0};
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;
    uint16_t cnt;
    uint16_t irow = 0;

    num = (num > row) ? row : num;

    for (cnt = num >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        elem[irow] = m_elem[irow * offset];
        elem[irow + 1] = m_elem[(irow + 1) * offset];
        elem[irow + 2] = m_elem[(irow + 2) * offset];
        elem[irow + 3] = m_elem[(irow + 3) * offset];
    }

    for (cnt = num % 4u; cnt > 0u; cnt--, irow++)
    {
        elem[irow] = m_elem[irow * offset];
    }

    return Vector<row, t_type>(elem);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<0, t_type> Matrix<t_row, t_col, t_type>::GetDiagonalVec0()
{
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;
    uint16_t cnt;
    uint16_t irow = 0;
    Vector<0, t_type> v(num);

    for (cnt = num >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        v.m_elem[irow] = m_elem[irow * offset];
        v.m_elem[irow + 1] = m_elem[(irow + 1) * offset];
        v.m_elem[irow + 2] = m_elem[(irow + 2) * offset];
        v.m_elem[irow + 3] = m_elem[(irow + 3) * offset];
    }

    for (cnt = num % 4u; cnt > 0u; cnt--, irow++)
    {
        v.m_elem[irow] = m_elem[irow * offset];
    }

    return v;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector3<t_type, 3> Matrix<t_row, t_col, t_type>::GetDiagonalVec3()
{
    t_type elem[3]{0};
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 3) ? 3 : num;

    switch (num)
    {
    case 1:
        elem[0] = m_elem[0];
        break;
    case 2:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        break;
    case 3:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        elem[2] = m_elem[2 * offset];
        break;
    }

    return Vector3<t_type, 3>(elem);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector4<t_type, 4> Matrix<t_row, t_col, t_type>::GetDiagonalVec4()
{
    t_type elem[4]{0};
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 4) ? 4 : num;

    switch (num)
    {
    case 1:
        elem[0] = m_elem[0];
        break;
    case 2:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        break;
    case 3:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        elem[2] = m_elem[2 * offset];
        break;
    case 4:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        elem[2] = m_elem[2 * offset];
        elem[3] = m_elem[3 * offset];
        break;
    }

    return Vector4<t_type, 4>(elem);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector6<t_type, 6> Matrix<t_row, t_col, t_type>::GetDiagonalVec6()
{
    t_type elem[6]{0};
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 6) ? 6 : num;

    switch (num)
    {
    case 1:
        elem[0] = m_elem[0];
        break;
    case 2:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        break;
    case 3:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        elem[2] = m_elem[2 * offset];
        break;
    case 4:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        elem[2] = m_elem[2 * offset];
        elem[3] = m_elem[3 * offset];
        break;
    case 5:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        elem[2] = m_elem[2 * offset];
        elem[3] = m_elem[3 * offset];
        elem[4] = m_elem[4 * offset];
        break;
    case 6:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        elem[2] = m_elem[2 * offset];
        elem[3] = m_elem[3 * offset];
        elem[4] = m_elem[4 * offset];
        elem[5] = m_elem[5 * offset];
        break;
    }

    return Vector6<t_type, 6>(elem);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Quaternion<t_type, 4> Matrix<t_row, t_col, t_type>::GetDiagonalQuat()
{
    t_type elem[4]{0};
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;

    num = (num > 4) ? 4 : num;

    switch (num)
    {
    case 1:
        elem[0] = m_elem[0];
        break;
    case 2:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        break;
    case 3:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        elem[2] = m_elem[2 * offset];
        break;
    case 4:
        elem[0] = m_elem[0];
        elem[1] = m_elem[offset];
        elem[2] = m_elem[2 * offset];
        elem[3] = m_elem[3 * offset];
        break;
    }

    return Quaternion<t_type, 4>(elem);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t row, uint16_t col>
inline Matrix<row, col, t_type> Matrix<t_row, t_col, t_type>::GetBlock(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdxRow, const uint16_t jdxCol, const int16_t rSize, const int16_t cSize)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(row > jdxRow && "Index out of range");
    assert(col > jdxCol && "Index out of range");
    assert((row - jdxRow) >= rSize && "Over size");
    assert((col - jdxCol) >= cSize && "Over size");

    uint16_t irow, icol, cnt;
    t_type elem[row * col]{0};
    uint16_t rowSize = t_row - idxRow;
    uint16_t colSize = t_col - idxCol;

    if (idxRow >= t_row) return Matrix<row, col, t_type>(elem);
    if (idxCol >= t_col) return Matrix<row, col, t_type>(elem);
    if (rowSize > rSize) rowSize = rSize;
    if (colSize > cSize) colSize = cSize;

    for (irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            elem[(jdxRow + irow) * col + jdxCol + icol] = m_elem[(irow + idxRow) * t_col + idxCol + icol];
            elem[(jdxRow + irow) * col + jdxCol + icol + 1] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 1];
            elem[(jdxRow + irow) * col + jdxCol + icol + 2] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 2];
            elem[(jdxRow + irow) * col + jdxCol + icol + 3] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            elem[(jdxRow + irow) * col + jdxCol + icol] = m_elem[(irow + idxRow) * t_col + idxCol + icol];
        }
    }

    return Matrix<row, col, t_type>(elem);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<0, 0, t_type> Matrix<t_row, t_col, t_type>::GetBlock(const uint16_t idxRow, const uint16_t idxCol, const uint16_t row, const uint16_t col, const uint16_t jdxRow, const uint16_t jdxCol, int16_t rSize, int16_t cSize)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(row > jdxRow && "Index out of range");
    assert(col > jdxCol && "Index out of range");
    assert((row - jdxRow) >= rSize && "Over size");
    assert((col - jdxCol) >= cSize && "Over size");

    uint16_t irow, icol, cnt;
    uint16_t rowSize = t_row - idxRow;
    uint16_t colSize = t_col - idxCol;
    Matrix<0, 0, t_type> mat(row, col);

    if (rSize < 0) rSize = row;
    if (cSize < 0) cSize = col;
    if (rowSize > rSize) rowSize = rSize;
    if (colSize > cSize) colSize = cSize;

    for (irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            mat.m_elem[(jdxRow + irow) * col + jdxCol + icol] = m_elem[(irow + idxRow) * t_col + idxCol + icol];
            mat.m_elem[(jdxRow + irow) * col + jdxCol + icol + 1] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 1];
            mat.m_elem[(jdxRow + irow) * col + jdxCol + icol + 2] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 2];
            mat.m_elem[(jdxRow + irow) * col + jdxCol + icol + 3] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            mat.m_elem[(jdxRow + irow) * col + jdxCol + icol] = m_elem[(irow + idxRow) * t_col + idxCol + icol];
        }
    }

    return mat;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t row, uint16_t col>
inline int8_t Matrix<t_row, t_col, t_type>::GetBlock(const uint16_t idxRow, const uint16_t idxCol, Matrix<row, col, t_type> &m, const uint16_t jdxRow, const uint16_t jdxCol, const int16_t rSize, const int16_t cSize)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(row > jdxRow && "Index out of range");
    assert(col > jdxCol && "Index out of range");
    assert((row - jdxRow) >= rSize && "Over size");
    assert((col - jdxCol) >= cSize && "Over size");

    uint16_t irow, icol, cnt;
    uint16_t rowSize = t_row - idxRow;
    uint16_t colSize = t_col - idxCol;

    if (idxRow >= t_row) return -1;
    if (idxCol >= t_col) return -1;
    if (rowSize > rSize) rowSize = rSize;
    if (colSize > cSize) colSize = cSize;

    for (irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            m.m_elem[(jdxRow + irow) * col + jdxCol + icol] = m_elem[(irow + idxRow) * t_col + idxCol + icol];
            m.m_elem[(jdxRow + irow) * col + jdxCol + icol + 1] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 1];
            m.m_elem[(jdxRow + irow) * col + jdxCol + icol + 2] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 2];
            m.m_elem[(jdxRow + irow) * col + jdxCol + icol + 3] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            m.m_elem[(jdxRow + irow) * col + jdxCol + icol] = m_elem[(irow + idxRow) * t_col + idxCol + icol];
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetBlock(const uint16_t idxRow, const uint16_t idxCol, Matrix3<t_type, 3, 3> &m, const uint16_t jdxRow, const uint16_t jdxCol, const int16_t rSize, const int16_t cSize)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(3 > jdxRow && "Index out of range");
    assert(3 > jdxCol && "Index out of range");
    assert((3 - jdxRow) >= rSize && "Over size");
    assert((3 - jdxCol) >= cSize && "Over size");

    uint16_t irow, icol, cnt;
    uint16_t rowSize = t_row - idxRow;
    uint16_t colSize = t_col - idxCol;

    if (idxRow >= t_row) return -1;
    if (idxCol >= t_col) return -1;
    if (rowSize > rSize) rowSize = rSize;
    if (colSize > cSize) colSize = cSize;

    for (irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            m.m_elem[(jdxRow + irow) * 3 + jdxCol + icol] = m_elem[(irow + idxRow) * t_col + idxCol + icol];
            m.m_elem[(jdxRow + irow) * 3 + jdxCol + icol + 1] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 1];
            m.m_elem[(jdxRow + irow) * 3 + jdxCol + icol + 2] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 2];
            m.m_elem[(jdxRow + irow) * 3 + jdxCol + icol + 3] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            m.m_elem[(jdxRow + irow) * 3 + jdxCol + icol] = m_elem[(irow + idxRow) * t_col + idxCol + icol];
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetBlock(const uint16_t idxRow, const uint16_t idxCol, Matrix<0, 0, t_type> &m, const uint16_t jdxRow, const uint16_t jdxCol, int16_t rSize, int16_t cSize)
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(m.m_row > jdxRow && "Index out of range");
    assert(m.m_col > jdxCol && "Index out of range");
    assert((m.m_row - jdxRow) >= rSize && "Over size");
    assert((m.m_col - jdxCol) >= cSize && "Over size");
    assert(m.m_elem != nullptr && "Memory has not been allocated");

    uint16_t irow, icol, cnt;
    uint16_t rowSize = t_row - idxRow;
    uint16_t colSize = t_col - idxCol;

    if (idxRow >= t_row) return -1;
    if (idxCol >= t_col) return -1;
    if (rSize < 0) rSize = m.m_row;
    if (cSize < 0) cSize = m.m_col;
    if (rowSize > rSize) rowSize = rSize;
    if (colSize > cSize) colSize = cSize;

    for (irow = 0; irow < rowSize; ++irow)
    {
        for (cnt = colSize >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4)
        {
            m.m_elem[(jdxRow + irow) * m.m_col + jdxCol + icol] = m_elem[(irow + idxRow) * t_col + idxCol + icol];
            m.m_elem[(jdxRow + irow) * m.m_col + jdxCol + icol + 1] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 1];
            m.m_elem[(jdxRow + irow) * m.m_col + jdxCol + icol + 2] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 2];
            m.m_elem[(jdxRow + irow) * m.m_col + jdxCol + icol + 3] = m_elem[(irow + idxRow) * t_col + idxCol + icol + 3];
        }
        for (cnt = colSize % 4; cnt > 0u; cnt--, icol++)
        {
            m.m_elem[(jdxRow + irow) * m.m_col + jdxCol + icol] = m_elem[(irow + idxRow) * t_col + idxCol + icol];
        }
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t col>
inline int8_t Matrix<t_row, t_col, t_type>::GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Vector<col, t_type> &v, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(col > jdx && "Index out of range");
    assert(((int16_t)col - jdx) >= size && "over size");

    uint16_t cnt;
    uint16_t icol = 0;
    uint16_t colSize = t_col - idxCol;

    if (colSize > size) colSize = size;

    for (cnt = colSize >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        v.m_elem[jdx + icol] = m_elem[idxRow * t_col + idxCol + icol];
        v.m_elem[jdx + icol + 1] = m_elem[idxRow * t_col + idxCol + icol + 1];
        v.m_elem[jdx + icol + 2] = m_elem[idxRow * t_col + idxCol + icol + 2];
        v.m_elem[jdx + icol + 3] = m_elem[idxRow * t_col + idxCol + icol + 3];
    }

    for (cnt = colSize % 4u; cnt > 0u; cnt--, icol++)
    {
        v.m_elem[jdx + icol] = m_elem[idxRow * t_col + idxCol + icol];
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Vector<0, t_type> &v, const uint16_t jdx, int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(v.m_row > jdx && "Index out of range");
    assert(((int16_t)v.m_row - jdx) >= size && "over size");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt;
    uint16_t icol = 0;
    uint16_t colSize = t_col - idxCol;

    if (size < 0) size = v.m_row;
    if (colSize > size) colSize = size;

    for (cnt = colSize >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        v.m_elem[jdx + icol] = m_elem[idxRow * t_col + idxCol + icol];
        v.m_elem[jdx + icol + 1] = m_elem[idxRow * t_col + idxCol + icol + 1];
        v.m_elem[jdx + icol + 2] = m_elem[idxRow * t_col + idxCol + icol + 2];
        v.m_elem[jdx + icol + 3] = m_elem[idxRow * t_col + idxCol + icol + 3];
    }

    for (cnt = colSize % 4u; cnt > 0u; cnt--, icol++)
    {
        v.m_elem[jdx + icol] = m_elem[idxRow * t_col + idxCol + icol];
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Vector3<t_type, 3> &v, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(3 > jdx && "Index out of range");
    assert((3 - jdx) >= size && "over size");

    if (idxRow >= t_row) return -1;
    if (idxCol >= t_col) return -1;

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        break;
    case 3:
    default:
        v.m_elem[0] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[1] = m_elem[idxRow * t_col + idxCol + 1];
        v.m_elem[2] = m_elem[idxRow * t_col + idxCol + 2];
        break;
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Vector4<t_type, 4> &v, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "over size");

    if (idxRow >= t_row) return -1;
    if (idxCol >= t_col) return -1;

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        break;
    case 3:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        v.m_elem[jdx + 2] = m_elem[idxRow * t_col + idxCol + 2];
        break;
    case 4:
    default:
        v.m_elem[0] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[1] = m_elem[idxRow * t_col + idxCol + 1];
        v.m_elem[2] = m_elem[idxRow * t_col + idxCol + 2];
        v.m_elem[3] = m_elem[idxRow * t_col + idxCol + 3];
        break;
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Vector6<t_type, 6> &v, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(6 > jdx && "Index out of range");
    assert((6 - jdx) >= size && "over size");

    if (idxRow >= t_row) return -1;
    if (idxCol >= t_col) return -1;

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        break;
    case 3:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        v.m_elem[jdx + 2] = m_elem[idxRow * t_col + idxCol + 2];
        break;
    case 4:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        v.m_elem[jdx + 2] = m_elem[idxRow * t_col + idxCol + 2];
        v.m_elem[jdx + 3] = m_elem[idxRow * t_col + idxCol + 3];
        break;
    case 5:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        v.m_elem[jdx + 2] = m_elem[idxRow * t_col + idxCol + 2];
        v.m_elem[jdx + 3] = m_elem[idxRow * t_col + idxCol + 3];
        v.m_elem[jdx + 4] = m_elem[idxRow * t_col + idxCol + 4];
        break;
    case 6:
    default:
        v.m_elem[0] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[1] = m_elem[idxRow * t_col + idxCol + 1];
        v.m_elem[2] = m_elem[idxRow * t_col + idxCol + 2];
        v.m_elem[3] = m_elem[idxRow * t_col + idxCol + 3];
        v.m_elem[4] = m_elem[idxRow * t_col + idxCol + 4];
        v.m_elem[5] = m_elem[idxRow * t_col + idxCol + 5];
        break;
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetRowVec(const uint16_t idxRow, const uint16_t idxCol, Quaternion<t_type, 4> &v, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "Over size");

    if (idxRow >= t_row) return -1;
    if (idxCol >= t_col) return -1;

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        break;
    case 3:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        v.m_elem[jdx + 2] = m_elem[idxRow * t_col + idxCol + 2];
        break;
    case 4:
    default:
        v.m_elem[0] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[1] = m_elem[idxRow * t_col + idxCol + 1];
        v.m_elem[2] = m_elem[idxRow * t_col + idxCol + 2];
        v.m_elem[3] = m_elem[idxRow * t_col + idxCol + 3];
        break;
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t row>
inline int8_t Matrix<t_row, t_col, t_type>::GetColVec(const uint16_t idxRow, const uint16_t idxCol, Vector<row, t_type> &v, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(row > jdx && "Index out of range");
    assert(((int16_t)row - jdx) >= size && "over size");

    uint16_t cnt;
    uint16_t irow = 0;
    uint16_t rowSize = t_row - idxRow;

    if (rowSize > size) rowSize = size;

    for (cnt = rowSize >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        v.m_elem[jdx + irow] = m_elem[(idxRow + irow) * t_col + idxCol];
        v.m_elem[jdx + irow + 1] = m_elem[(idxRow + irow + 1) * t_col + idxCol];
        v.m_elem[jdx + irow + 2] = m_elem[(idxRow + irow + 2) * t_col + idxCol];
        v.m_elem[jdx + irow + 3] = m_elem[(idxRow + irow + 3) * t_col + idxCol];
    }

    for (cnt = rowSize % 4u; cnt > 0u; cnt--, irow++)
    {
        v.m_elem[jdx + irow] = m_elem[(idxRow + irow) * t_col + idxCol];
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetColVec(const uint16_t idxRow, const uint16_t idxCol, Vector<0, t_type> &v, const uint16_t jdx, int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(v.m_row > jdx && "Index out of range");
    assert(((int16_t)v.m_row - jdx) >= size && "Over size");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    uint16_t cnt;
    uint16_t irow = 0;
    uint16_t rowSize = t_row - idxRow;

    if (size < 0) size = v.m_row;
    if (rowSize > size) rowSize = size;

    for (cnt = rowSize >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        v.m_elem[jdx + irow] = m_elem[(idxRow + irow) * t_col + idxCol];
        v.m_elem[jdx + irow + 1] = m_elem[(idxRow + irow + 1) * t_col + idxCol];
        v.m_elem[jdx + irow + 2] = m_elem[(idxRow + irow + 2) * t_col + idxCol];
        v.m_elem[jdx + irow + 3] = m_elem[(idxRow + irow + 3) * t_col + idxCol];
    }

    for (cnt = rowSize % 4u; cnt > 0u; cnt--, irow++)
    {
        v.m_elem[jdx + irow] = m_elem[(idxRow + irow) * t_col + idxCol];
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetColVec(const uint16_t idxRow, const uint16_t idxCol, Vector3<t_type, 3> &v, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(3 > jdx && "Index out of range");
    assert((3 - jdx) >= size && "Over size");

    if (idxRow >= t_row) return -1;
    if (idxCol >= t_col) return -1;

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        break;
    case 3:
    default:
        v.m_elem[0] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[1] = m_elem[(idxRow + 1) * t_col + idxCol];
        v.m_elem[2] = m_elem[(idxRow + 2) * t_col + idxCol];
        break;
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetColVec(const uint16_t idxRow, const uint16_t idxCol, Vector4<t_type, 4> &v, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "Over size");

    if (idxRow >= t_row) return -1;
    if (idxCol >= t_col) return -1;

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        break;
    case 3:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        v.m_elem[jdx + 2] = m_elem[(idxRow + 2) * t_col + idxCol];
        break;
    case 4:
    default:
        v.m_elem[0] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[1] = m_elem[(idxRow + 1) * t_col + idxCol];
        v.m_elem[2] = m_elem[(idxRow + 2) * t_col + idxCol];
        v.m_elem[3] = m_elem[(idxRow + 3) * t_col + idxCol];
        break;
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetColVec(const uint16_t idxRow, const uint16_t idxCol, Vector6<t_type, 6> &v, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(6 > jdx && "Index out of range");
    assert((6 - jdx) >= size && "Over size");

    if (idxRow >= t_row) return -1;
    if (idxCol >= t_col) return -1;

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        break;
    case 3:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        v.m_elem[jdx + 2] = m_elem[(idxRow + 2) * t_col + idxCol];
        break;
    case 4:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        v.m_elem[jdx + 2] = m_elem[(idxRow + 2) * t_col + idxCol];
        v.m_elem[jdx + 3] = m_elem[(idxRow + 3) * t_col + idxCol];
        break;
    case 5:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        v.m_elem[jdx + 2] = m_elem[(idxRow + 2) * t_col + idxCol];
        v.m_elem[jdx + 3] = m_elem[(idxRow + 3) * t_col + idxCol];
        v.m_elem[jdx + 4] = m_elem[(idxRow + 4) * t_col + idxCol];
        break;
    case 6:
    default:
        v.m_elem[0] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[1] = m_elem[(idxRow + 1) * t_col + idxCol];
        v.m_elem[2] = m_elem[(idxRow + 2) * t_col + idxCol];
        v.m_elem[3] = m_elem[(idxRow + 3) * t_col + idxCol];
        v.m_elem[4] = m_elem[(idxRow + 4) * t_col + idxCol];
        v.m_elem[5] = m_elem[(idxRow + 5) * t_col + idxCol];
        break;
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline int8_t Matrix<t_row, t_col, t_type>::GetColVec(const uint16_t idxRow, const uint16_t idxCol, Quaternion<t_type, 4> &v, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "Over size");

    if (idxRow >= t_row) return -1;
    if (idxCol >= t_col) return -1;

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        break;
    case 3:
        v.m_elem[jdx] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        v.m_elem[jdx + 2] = m_elem[(idxRow + 2) * t_col + idxCol];
        break;
    case 4:
    default:
        v.m_elem[0] = m_elem[idxRow * t_col + idxCol];
        v.m_elem[1] = m_elem[(idxRow + 1) * t_col + idxCol];
        v.m_elem[2] = m_elem[(idxRow + 2) * t_col + idxCol];
        v.m_elem[3] = m_elem[(idxRow + 3) * t_col + idxCol];
        break;
    }

    return 0;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t col>
inline Vector<col, t_type> Matrix<t_row, t_col, t_type>::GetRowVec(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(col > jdx && "Index out of range");
    assert(((int16_t)col - jdx) >= size && "Over size");

    t_type vec[col];
    uint16_t cnt;
    uint16_t icol = 0;
    uint16_t colSize = t_col - idxCol;

    if (colSize > size) colSize = size;

    for (cnt = colSize >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        vec[jdx + icol] = m_elem[idxRow * t_col + idxCol + icol];
        vec[jdx + icol + 1] = m_elem[idxRow * t_col + idxCol + icol + 1];
        vec[jdx + icol + 2] = m_elem[idxRow * t_col + idxCol + icol + 2];
        vec[jdx + icol + 3] = m_elem[idxRow * t_col + idxCol + icol + 3];
    }

    for (cnt = colSize % 4u; cnt > 0u; cnt--, icol++)
    {
        vec[jdx + icol] = m_elem[idxRow * t_col + idxCol + icol];
    }

    return Vector<col, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<0, t_type> Matrix<t_row, t_col, t_type>::GetRowVec0(const uint16_t idxRow, const uint16_t idxCol, const uint16_t col, const uint16_t jdx, int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(col > jdx && "Index out of range");
    assert(((int16_t)col - jdx) >= size && "Over size");

    Vector<0, t_type> vec(col);
    uint16_t cnt;
    uint16_t icol = 0;
    uint16_t colSize = t_col - idxCol;

    if (size < 0) size = col;
    if (colSize > size) colSize = size;

    for (cnt = colSize >> 2u; cnt > 0u; cnt--, icol += 4)
    {
        vec.m_elem[jdx + icol] = m_elem[idxRow * t_col + idxCol + icol];
        vec.m_elem[jdx + icol + 1] = m_elem[idxRow * t_col + idxCol + icol + 1];
        vec.m_elem[jdx + icol + 2] = m_elem[idxRow * t_col + idxCol + icol + 2];
        vec.m_elem[jdx + icol + 3] = m_elem[idxRow * t_col + idxCol + icol + 3];
    }

    for (cnt = colSize % 4u; cnt > 0u; cnt--, icol++)
    {
        vec.m_elem[jdx + icol] = m_elem[idxRow * t_col + idxCol + icol];
    }

    return vec;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector3<t_type, 3> Matrix<t_row, t_col, t_type>::GetRowVec3(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(3 > jdx && "Index out of range");
    assert((3 - jdx) >= size && "over size");

    t_type vec[3]{0};

    if (idxRow >= t_row) return Vector3<t_type, 3>();
    if (idxCol >= t_col) return Vector3<t_type, 3>();

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        break;
    case 3:
    default:
        vec[0] = m_elem[idxRow * t_col + idxCol];
        vec[1] = m_elem[idxRow * t_col + idxCol + 1];
        vec[2] = m_elem[idxRow * t_col + idxCol + 2];
        break;
    }

    return Vector3<t_type, 3>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector4<t_type, 4> Matrix<t_row, t_col, t_type>::GetRowVec4(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "over size");

    t_type vec[4]{0};

    if (idxRow >= t_row) return Vector4<t_type, 4>();
    if (idxCol >= t_col) return Vector4<t_type, 4>();

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        break;
    case 3:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        vec[jdx + 2] = m_elem[idxRow * t_col + idxCol + 2];
        break;
    case 4:
    default:
        vec[0] = m_elem[idxRow * t_col + idxCol];
        vec[1] = m_elem[idxRow * t_col + idxCol + 1];
        vec[2] = m_elem[idxRow * t_col + idxCol + 2];
        vec[3] = m_elem[idxRow * t_col + idxCol + 3];
        break;
    }

    return Vector4<t_type, 4>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector6<t_type, 6> Matrix<t_row, t_col, t_type>::GetRowVec6(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(6 > jdx && "Index out of range");
    assert((6 - jdx) >= size && "over size");

    t_type vec[6]{0};

    if (idxRow >= t_row) return Vector6<t_type, 6>();
    if (idxCol >= t_col) return Vector6<t_type, 6>();

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        break;
    case 3:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        vec[jdx + 2] = m_elem[idxRow * t_col + idxCol + 2];
        break;
    case 4:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        vec[jdx + 2] = m_elem[idxRow * t_col + idxCol + 2];
        vec[jdx + 3] = m_elem[idxRow * t_col + idxCol + 3];
        break;
    case 5:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        vec[jdx + 2] = m_elem[idxRow * t_col + idxCol + 2];
        vec[jdx + 3] = m_elem[idxRow * t_col + idxCol + 3];
        vec[jdx + 4] = m_elem[idxRow * t_col + idxCol + 4];
        break;
    case 6:
    default:
        vec[0] = m_elem[idxRow * t_col + idxCol];
        vec[1] = m_elem[idxRow * t_col + idxCol + 1];
        vec[2] = m_elem[idxRow * t_col + idxCol + 2];
        vec[3] = m_elem[idxRow * t_col + idxCol + 3];
        vec[4] = m_elem[idxRow * t_col + idxCol + 4];
        vec[5] = m_elem[idxRow * t_col + idxCol + 5];
        break;
    }

    return Vector6<t_type, 6>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Quaternion<t_type, 4> Matrix<t_row, t_col, t_type>::GetRowQuat(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "over size");

    t_type vec[4]{0};

    if (idxRow >= t_row) return Vector4<t_type, 4>();
    if (idxCol >= t_col) return Vector4<t_type, 4>();

    uint16_t colSize = t_col - idxCol;
    if (colSize > size) colSize = size;

    switch (colSize)
    {
    case 1:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        break;
    case 3:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[idxRow * t_col + idxCol + 1];
        vec[jdx + 2] = m_elem[idxRow * t_col + idxCol + 2];
        break;
    case 4:
    default:
        vec[0] = m_elem[idxRow * t_col + idxCol];
        vec[1] = m_elem[idxRow * t_col + idxCol + 1];
        vec[2] = m_elem[idxRow * t_col + idxCol + 2];
        vec[3] = m_elem[idxRow * t_col + idxCol + 3];
        break;
    }

    return Quaternion<t_type, 4>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t row>
inline Vector<row, t_type> Matrix<t_row, t_col, t_type>::GetColVec(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(row > jdx && "Index out of range");
    assert(((int16_t)row - jdx) >= size && "over size");

    t_type vec[row];
    uint16_t cnt;
    uint16_t irow = 0;
    uint16_t rowSize = t_row - idxRow;

    if (rowSize > size) rowSize = size;

    for (cnt = rowSize >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec[jdx + irow] = m_elem[(idxRow + irow) * t_col + idxCol];
        vec[jdx + irow + 1] = m_elem[(idxRow + irow + 1) * t_col + idxCol];
        vec[jdx + irow + 2] = m_elem[(idxRow + irow + 2) * t_col + idxCol];
        vec[jdx + irow + 3] = m_elem[(idxRow + irow + 3) * t_col + idxCol];
    }

    for (cnt = rowSize % 4u; cnt > 0u; cnt--, irow++)
    {
        vec[jdx + irow] = m_elem[(idxRow + irow) * t_col + idxCol];
    }

    return Vector<row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<0, t_type> Matrix<t_row, t_col, t_type>::GetColVec0(const uint16_t idxRow, const uint16_t idxCol, const uint16_t row, const uint16_t jdx, int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(row > jdx && "Index out of range");
    assert(((int16_t)row - jdx) >= size && "over size");

    Vector<0, t_type> vec(row);
    uint16_t cnt;
    uint16_t irow = 0;
    uint16_t rowSize = t_row - idxRow;

    if (size < 0) size = row;
    if (rowSize > size) rowSize = size;

    for (cnt = rowSize >> 2u; cnt > 0u; cnt--, irow += 4)
    {
        vec.m_elem[jdx + irow] = m_elem[(idxRow + irow) * t_col + idxCol];
        vec.m_elem[jdx + irow + 1] = m_elem[(idxRow + irow + 1) * t_col + idxCol];
        vec.m_elem[jdx + irow + 2] = m_elem[(idxRow + irow + 2) * t_col + idxCol];
        vec.m_elem[jdx + irow + 3] = m_elem[(idxRow + irow + 3) * t_col + idxCol];
    }

    for (cnt = rowSize % 4u; cnt > 0u; cnt--, irow++)
    {
        vec.m_elem[jdx + irow] = m_elem[(idxRow + irow) * t_col + idxCol];
    }

    return vec;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector3<t_type, 3> Matrix<t_row, t_col, t_type>::GetColVec3(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(3 > jdx && "Index out of range");
    assert((3 - jdx) >= size && "over size");

    t_type vec[3]{0};

    if (idxRow >= t_row) return Vector3<t_type, 3>();
    if (idxCol >= t_col) return Vector3<t_type, 3>();

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        break;
    case 3:
    default:
        vec[0] = m_elem[idxRow * t_col + idxCol];
        vec[1] = m_elem[(idxRow + 1) * t_col + idxCol];
        vec[2] = m_elem[(idxRow + 2) * t_col + idxCol];
        break;
    }

    return Vector3<t_type, 3>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector4<t_type, 4> Matrix<t_row, t_col, t_type>::GetColVec4(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "over size");

    t_type vec[4]{0};

    if (idxRow >= t_row) return Vector4<t_type, 4>();
    if (idxCol >= t_col) return Vector4<t_type, 4>();

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        break;
    case 3:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        vec[jdx + 2] = m_elem[(idxRow + 2) * t_col + idxCol];
        break;
    case 4:
    default:
        vec[0] = m_elem[idxRow * t_col + idxCol];
        vec[1] = m_elem[(idxRow + 1) * t_col + idxCol];
        vec[2] = m_elem[(idxRow + 2) * t_col + idxCol];
        vec[3] = m_elem[(idxRow + 3) * t_col + idxCol];
    }

    return Vector4<t_type, 4>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector6<t_type, 6> Matrix<t_row, t_col, t_type>::GetColVec6(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(6 > jdx && "Index out of range");
    assert((6 - jdx) >= size && "over size");

    t_type vec[6]{0};

    if (idxRow >= t_row) return Vector6<t_type, 6>();
    if (idxCol >= t_col) return Vector6<t_type, 6>();

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        break;
    case 3:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        vec[jdx + 2] = m_elem[(idxRow + 2) * t_col + idxCol];
        break;
    case 4:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        vec[jdx + 2] = m_elem[(idxRow + 2) * t_col + idxCol];
        vec[jdx + 3] = m_elem[(idxRow + 3) * t_col + idxCol];
        break;
    case 5:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        vec[jdx + 2] = m_elem[(idxRow + 2) * t_col + idxCol];
        vec[jdx + 3] = m_elem[(idxRow + 3) * t_col + idxCol];
        vec[jdx + 4] = m_elem[(idxRow + 4) * t_col + idxCol];
        break;
    case 6:
    default:
        vec[0] = m_elem[idxRow * t_col + idxCol];
        vec[1] = m_elem[(idxRow + 1) * t_col + idxCol];
        vec[2] = m_elem[(idxRow + 2) * t_col + idxCol];
        vec[3] = m_elem[(idxRow + 3) * t_col + idxCol];
        vec[4] = m_elem[(idxRow + 4) * t_col + idxCol];
        vec[5] = m_elem[(idxRow + 5) * t_col + idxCol];
    }

    return Vector6<t_type, 6>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Quaternion<t_type, 4> Matrix<t_row, t_col, t_type>::GetColQuat(const uint16_t idxRow, const uint16_t idxCol, const uint16_t jdx, const int16_t size) const
{
    assert(t_row > idxRow && "Index out of range");
    assert(t_col > idxCol && "Index out of range");
    assert(4 > jdx && "Index out of range");
    assert((4 - jdx) >= size && "over size");

    t_type vec[4]{0};

    if (idxRow >= t_row) return Vector4<t_type, 4>();
    if (idxCol >= t_col) return Vector4<t_type, 4>();

    uint16_t rowSize = t_row - idxRow;
    if (rowSize > size) rowSize = size;

    switch (rowSize)
    {
    case 1:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        break;
    case 2:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        break;
    case 3:
        vec[jdx] = m_elem[idxRow * t_col + idxCol];
        vec[jdx + 1] = m_elem[(idxRow + 1) * t_col + idxCol];
        vec[jdx + 2] = m_elem[(idxRow + 2) * t_col + idxCol];
        break;
    case 4:
    default:
        vec[0] = m_elem[idxRow * t_col + idxCol];
        vec[1] = m_elem[(idxRow + 1) * t_col + idxCol];
        vec[2] = m_elem[(idxRow + 2) * t_col + idxCol];
        vec[3] = m_elem[(idxRow + 3) * t_col + idxCol];
    }

    return Quaternion<t_type, 4>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CscMatrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::GetCscMat() const
{
    return CscMatrix<t_row, t_col, t_type>(*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_col, t_row, t_type> Matrix<t_row, t_col, t_type>::Transpose() const
{
    t_type mat[t_row * t_col];
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; ++irow)
    {
        for (cnt = t_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
        {
            mat[irow + t_row * icol] = m_elem[irow * t_col + icol];
            mat[irow + t_row * (icol + 1)] = m_elem[irow * t_col + icol + 1];
            mat[irow + t_row * (icol + 2)] = m_elem[irow * t_col + icol + 2];
            mat[irow + t_row * (icol + 3)] = m_elem[irow * t_col + icol + 3];
        }

        for (cnt = t_col % 4u; cnt > 0u; cnt--, icol++)
        {
            mat[irow + t_row * icol] = m_elem[irow * t_col + icol];
        }
    }

    return Matrix<t_col, t_row, t_type>(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type Matrix<t_row, t_col, t_type>::Trace() const
{
    int cnt, i;
    uint16_t num = (t_row > t_col) ? t_col : t_row;
    uint16_t offset = t_col + 1;
    t_type sum = 0;

    // for (uint16_t i = 0; i < num; i++)
    //     sum += m_elem[i * offset];
    for (cnt = num >> 2, i = 0; cnt > 0; cnt--, i += 4)
    {
        sum += m_elem[i * offset];
        sum += m_elem[(i + 1) * offset];
        sum += m_elem[(i + 2) * offset];
        sum += m_elem[(i + 3) * offset];
    }
    for (cnt = num % 4; cnt > 0; cnt--, i++)
    {
        sum += m_elem[i * offset];
    }

    return sum;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type Matrix<t_row, t_col, t_type>::GetNorm() const
{
    uint16_t cnt, i = 0;
    t_type sqSum = 0;

    for (cnt = (t_row * t_col) >> 2u; cnt > 0u; cnt--, i += 4)
    {
        sqSum += m_elem[i] * m_elem[i];
        sqSum += m_elem[i + 1] * m_elem[i + 1];
        sqSum += m_elem[i + 2] * m_elem[i + 2];
        sqSum += m_elem[i + 3] * m_elem[i + 3];
    }

    for (cnt = (t_row * t_col) % 4u; cnt > 0u; cnt--, i++)
    {
        sqSum += m_elem[i] * m_elem[i];
    }

    return std::sqrt(sqSum);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type Matrix<t_row, t_col, t_type>::GetSqNorm() const
{
    uint16_t cnt, i = 0;
    t_type sqSum = 0;

    for (cnt = (t_row * t_col) >> 2u; cnt > 0u; cnt--, i += 4)
    {
        sqSum += m_elem[i] * m_elem[i];
        sqSum += m_elem[i + 1] * m_elem[i + 1];
        sqSum += m_elem[i + 2] * m_elem[i + 2];
        sqSum += m_elem[i + 3] * m_elem[i + 3];
    }

    for (cnt = (t_row * t_col) % 4u; cnt > 0u; cnt--, i++)
    {
        sqSum += m_elem[i] * m_elem[i];
    }

    return sqSum;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type Matrix<t_row, t_col, t_type>::GetLpNorm(const int p) const
{
    uint16_t cnt, i = 0;
    t_type powSum = 0;

    for (cnt = (t_row * t_col) >> 2u; cnt > 0u; cnt--, i += 4)
    {
        powSum += std::pow(std::abs(m_elem[i]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 1]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 2]), (t_type)p);
        powSum += std::pow(std::abs(m_elem[i + 3]), (t_type)p);
    }

    for (cnt = (t_row * t_col) % 4u; cnt > 0u; cnt--, i++)
    {
        powSum += std::pow(std::abs(m_elem[i]), (t_type)p);
    }

    return std::pow(powSum, (t_type)1 / p);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type Matrix<t_row, t_col, t_type>::Determinant() const
{
    return PartialPivLU<t_row, t_col, t_type>(*this).Determinant();
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline NoPivLU<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::GetNoPivLU() const
{
    return NoPivLU<t_row, t_col, t_type>(*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline PartialPivLU<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::GetPartialPivLU() const
{
    return PartialPivLU<t_row, t_col, t_type>(*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline FullPivLU<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::GetFullPivLU() const
{
    return FullPivLU<t_row, t_col, t_type>(*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline LLT<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::GetLLT() const
{
    return LLT<t_row, t_col, t_type>(*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline LDLT<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::GetLDLT() const
{
    return LDLT<t_row, t_col, t_type>(*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline QR<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::GetQR() const
{
    return QR<t_row, t_col, t_type>(*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline SVD<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::GetSVD() const
{
    return SVD<t_row, t_col, t_type>(*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::Inv(int8_t *isOk) const
{
    return PartialPivLU<t_row, t_col, t_type>(*this).Inverse(isOk);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::FInv(int8_t *isOk) const
{
    return FullPivLU<t_row, t_col, t_type>(*this).Inverse(isOk);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_col, t_row, t_type> Matrix<t_row, t_col, t_type>::PInv(int8_t *isOk, t_type tolerance) const
{
    return SVD<t_row, t_col, t_type>(*this).Inverse(isOk, tolerance);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type Matrix<t_row, t_col, t_type>::Max(uint16_t *irow, uint16_t *icol)
{
    t_type *pElem = m_elem;
    t_type max, value;     // max value, elem value
    uint16_t outIdx = 0;   // max value index
    uint16_t cnt = 0, idx; // loop count, index

    max = *pElem++; // max value initialization

    cnt = (t_row * t_col - 1) >> 2u;
    while (cnt > 0u)
    {
        value = *pElem++;
        if (max < value)
        {
            max = value;
            outIdx = idx + 1u;
        }

        value = *pElem++;
        if (max < value)
        {
            max = value;
            outIdx = idx + 2u;
        }

        value = *pElem++;
        if (max < value)
        {
            max = value;
            outIdx = idx + 3u;
        }

        value = *pElem++;
        if (max < value)
        {
            max = value;
            outIdx = idx + 4u;
        }

        idx += 4u;
        cnt--;
    }

    cnt = (t_row * t_col - 1) % 4u;
    while (cnt > 0u)
    {
        value = *pElem++;
        if (max < value)
        {
            max = value;
            outIdx = t_row * t_col - cnt;
        }

        cnt--;
    }

    if (irow != nullptr) *irow = outIdx / t_col;
    if (icol != nullptr) *icol = outIdx % t_col;

    return max;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type Matrix<t_row, t_col, t_type>::Min(uint16_t *irow, uint16_t *icol)
{
    t_type *pElem = m_elem;
    t_type min, value;     // min value, elem value
    uint16_t outIdx = 0;   // min value index
    uint16_t cnt = 0, idx; // loop count, index

    min = *pElem++; // min value initialization

    cnt = (t_row * t_col - 1) >> 2u;
    while (cnt > 0u)
    {
        value = *pElem++;
        if (min > value)
        {
            min = value;
            outIdx = idx + 1u;
        }

        value = *pElem++;
        if (min > value)
        {
            min = value;
            outIdx = idx + 1u;
        }

        value = *pElem++;
        if (min > value)
        {
            min = value;
            outIdx = idx + 1u;
        }

        value = *pElem++;
        if (min > value)
        {
            min = value;
            outIdx = idx + 1u;
        }

        idx += 4u;
        cnt--;
    }

    cnt = (t_row * t_col - 1) % 4u;
    while (cnt > 0u)
    {
        value = *pElem++;
        if (min > value)
        {
            min = value;
            outIdx = t_row * t_col - cnt;
        }

        cnt--;
    }

    if (irow != nullptr) *irow = outIdx / t_col;
    if (icol != nullptr) *icol = outIdx % t_col;

    return min;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type Matrix<t_row, t_col, t_type>::AbsMax(uint16_t *irow, uint16_t *icol)
{
    t_type *pElem = m_elem;
    t_type max, value;     // max value, elem value
    uint16_t outIdx = 0;   // max value index
    uint16_t cnt = 0, idx; // loop count, index

    max = std::abs(*pElem++); // max value initialization

    cnt = (t_row * t_col - 1) >> 2u;
    while (cnt > 0u)
    {
        value = std::abs(*pElem++);
        if (max < value)
        {
            max = value;
            outIdx = idx + 1u;
        }

        value = std::abs(*pElem++);
        if (max < value)
        {
            max = value;
            outIdx = idx + 2u;
        }

        value = std::abs(*pElem++);
        if (max < value)
        {
            max = value;
            outIdx = idx + 3u;
        }

        value = std::abs(*pElem++);
        if (max < value)
        {
            max = value;
            outIdx = idx + 4u;
        }

        idx += 4u;
        cnt--;
    }

    cnt = (t_row * t_col - 1) % 4u;
    while (cnt > 0u)
    {
        value = std::abs(*pElem++);
        if (max < value)
        {
            max = value;
            outIdx = t_row * t_col - cnt;
        }

        cnt--;
    }

    if (irow != nullptr) *irow = outIdx / t_col;
    if (icol != nullptr) *icol = outIdx % t_col;

    return max;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type Matrix<t_row, t_col, t_type>::AbsMin(uint16_t *irow, uint16_t *icol)
{
    t_type *pElem = m_elem;
    t_type min, value;     // min value, elem value
    uint16_t outIdx = 0;   // min value index
    uint16_t cnt = 0, idx; // loop count, index

    min = std::abs(*pElem++); // min value initialization

    cnt = (t_row * t_col - 1) >> 2u;
    while (cnt > 0u)
    {
        value = std::abs(*pElem++);
        if (min > value)
        {
            min = value;
            outIdx = idx + 1u;
        }

        value = std::abs(*pElem++);
        if (min > value)
        {
            min = value;
            outIdx = idx + 1u;
        }

        value = std::abs(*pElem++);
        if (min > value)
        {
            min = value;
            outIdx = idx + 1u;
        }

        value = std::abs(*pElem++);
        if (min > value)
        {
            min = value;
            outIdx = idx + 1u;
        }

        idx += 4u;
        cnt--;
    }

    cnt = (t_row * t_col - 1) % 4u;
    while (cnt > 0u)
    {
        value = *pElem++;
        if (min > value)
        {
            min = value;
            outIdx = t_row * t_col - cnt;
        }

        cnt--;
    }

    if (irow != nullptr) *irow = outIdx / t_col;
    if (icol != nullptr) *icol = outIdx % t_col;

    return min;
}

/* Member access operators */
template <uint16_t t_row, uint16_t t_col, typename t_type>
inline t_type &Matrix<t_row, t_col, t_type>::operator()(uint16_t irow, uint16_t icol)
{
    assert(irow < t_row && "Index out of range");
    assert(icol < t_col && "Index out of range");

    return m_elem[irow * t_col + icol];
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline const t_type &Matrix<t_row, t_col, t_type>::operator()(uint16_t irow, uint16_t icol) const
{
    assert(irow < t_row && "Index out of range");
    assert(icol < t_col && "Index out of range");

    return m_elem[irow * t_col + icol];
}

/* Assignment operators */
template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator=(const Matrix<t_row, t_col, t_type> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator+=(const Matrix<t_row, t_col, t_type> &m)
{
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] += m.m_elem[i];
        m_elem[i + 2] += m.m_elem[i + 2];
        m_elem[i + 4] += m.m_elem[i + 4];
        m_elem[i + 6] += m.m_elem[i + 6];
        m_elem[i + 1] += m.m_elem[i + 1];
        m_elem[i + 3] += m.m_elem[i + 3];
        m_elem[i + 5] += m.m_elem[i + 5];
        m_elem[i + 7] += m.m_elem[i + 7];
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] += m.m_elem[i];
    }

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator-=(const Matrix<t_row, t_col, t_type> &m)
{
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] -= m.m_elem[i];
        m_elem[i + 2] -= m.m_elem[i + 2];
        m_elem[i + 4] -= m.m_elem[i + 4];
        m_elem[i + 6] -= m.m_elem[i + 6];
        m_elem[i + 1] -= m.m_elem[i + 1];
        m_elem[i + 3] -= m.m_elem[i + 3];
        m_elem[i + 5] -= m.m_elem[i + 5];
        m_elem[i + 7] -= m.m_elem[i + 7];
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] -= m.m_elem[i];
    }

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator=(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == m.m_row && "Row dimensions do not matched");
    assert(t_col == m.m_col && "Col dimensions do not matched");

    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator+=(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == m.m_row && "Row dimensions do not matched");
    assert(t_col == m.m_col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] += m.m_elem[i];
        m_elem[i + 2] += m.m_elem[i + 2];
        m_elem[i + 4] += m.m_elem[i + 4];
        m_elem[i + 6] += m.m_elem[i + 6];
        m_elem[i + 1] += m.m_elem[i + 1];
        m_elem[i + 3] += m.m_elem[i + 3];
        m_elem[i + 5] += m.m_elem[i + 5];
        m_elem[i + 7] += m.m_elem[i + 7];
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] += m.m_elem[i];
    }

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator-=(const Matrix<0, 0, t_type> &m)
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == m.m_row && "Row dimensions do not matched");
    assert(t_col == m.m_col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] -= m.m_elem[i];
        m_elem[i + 2] -= m.m_elem[i + 2];
        m_elem[i + 4] -= m.m_elem[i + 4];
        m_elem[i + 6] -= m.m_elem[i + 6];
        m_elem[i + 1] -= m.m_elem[i + 1];
        m_elem[i + 3] -= m.m_elem[i + 3];
        m_elem[i + 5] -= m.m_elem[i + 5];
        m_elem[i + 7] -= m.m_elem[i + 7];
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] -= m.m_elem[i];
    }

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator=(const Matrix3<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator+=(const Matrix3<t_type, t_row, t_col> &m)
{
    m_elem[0] += m.m_elem[0];
    m_elem[1] += m.m_elem[1];
    m_elem[2] += m.m_elem[2];
    m_elem[3] += m.m_elem[3];
    m_elem[4] += m.m_elem[4];
    m_elem[5] += m.m_elem[5];
    m_elem[6] += m.m_elem[6];
    m_elem[7] += m.m_elem[7];
    m_elem[8] += m.m_elem[8];

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator-=(const Matrix3<t_type, t_row, t_col> &m)
{
    m_elem[0] -= m.m_elem[0];
    m_elem[1] -= m.m_elem[1];
    m_elem[2] -= m.m_elem[2];
    m_elem[3] -= m.m_elem[3];
    m_elem[4] -= m.m_elem[4];
    m_elem[5] -= m.m_elem[5];
    m_elem[6] -= m.m_elem[6];
    m_elem[7] -= m.m_elem[7];
    m_elem[8] -= m.m_elem[8];

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator=(const Rotation<t_type, t_row, t_col> &m)
{
    memcpy(m_elem, m.m_elem, sizeof(t_type) * t_row * t_col);

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator+=(const Rotation<t_type, t_row, t_col> &m)
{
    m_elem[0] += m.m_elem[0];
    m_elem[1] += m.m_elem[1];
    m_elem[2] += m.m_elem[2];
    m_elem[3] += m.m_elem[3];
    m_elem[4] += m.m_elem[4];
    m_elem[5] += m.m_elem[5];
    m_elem[6] += m.m_elem[6];
    m_elem[7] += m.m_elem[7];
    m_elem[8] += m.m_elem[8];

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator-=(const Rotation<t_type, t_row, t_col> &m)
{
    m_elem[0] -= m.m_elem[0];
    m_elem[1] -= m.m_elem[1];
    m_elem[2] -= m.m_elem[2];
    m_elem[3] -= m.m_elem[3];
    m_elem[4] -= m.m_elem[4];
    m_elem[5] -= m.m_elem[5];
    m_elem[6] -= m.m_elem[6];
    m_elem[7] -= m.m_elem[7];
    m_elem[8] -= m.m_elem[8];

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator=(const Transform<t_type, t_row, t_col> &m)
{
    m_elem[0] = m.m_R.m_elem[0];
    m_elem[1] = m.m_R.m_elem[1];
    m_elem[2] = m.m_R.m_elem[2];
    m_elem[3] = m.m_p.m_elem[0];
    m_elem[4] = m.m_R.m_elem[3];
    m_elem[5] = m.m_R.m_elem[4];
    m_elem[6] = m.m_R.m_elem[5];
    m_elem[7] = m.m_p.m_elem[1];
    m_elem[8] = m.m_R.m_elem[6];
    m_elem[9] = m.m_R.m_elem[7];
    m_elem[10] = m.m_R.m_elem[8];
    m_elem[11] = m.m_p.m_elem[2];
    m_elem[12] = 0;
    m_elem[13] = 0;
    m_elem[14] = 0;
    m_elem[15] = 1;

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator+=(const Transform<t_type, t_row, t_col> &m)
{
    m_elem[0] += m.m_R.m_elem[0];
    m_elem[1] += m.m_R.m_elem[1];
    m_elem[2] += m.m_R.m_elem[2];
    m_elem[3] += m.m_p.m_elem[0];
    m_elem[4] += m.m_R.m_elem[3];
    m_elem[5] += m.m_R.m_elem[4];
    m_elem[6] += m.m_R.m_elem[5];
    m_elem[7] += m.m_p.m_elem[1];
    m_elem[8] += m.m_R.m_elem[6];
    m_elem[9] += m.m_R.m_elem[7];
    m_elem[10] += m.m_R.m_elem[8];
    m_elem[11] += m.m_p.m_elem[2];
    m_elem[15] += 1;

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator-=(const Transform<t_type, t_row, t_col> &m)
{
    m_elem[0] -= m.m_R.m_elem[0];
    m_elem[1] -= m.m_R.m_elem[1];
    m_elem[2] -= m.m_R.m_elem[2];
    m_elem[3] -= m.m_p.m_elem[0];
    m_elem[4] -= m.m_R.m_elem[3];
    m_elem[5] -= m.m_R.m_elem[4];
    m_elem[6] -= m.m_R.m_elem[5];
    m_elem[7] -= m.m_p.m_elem[1];
    m_elem[8] -= m.m_R.m_elem[6];
    m_elem[9] -= m.m_R.m_elem[7];
    m_elem[10] -= m.m_R.m_elem[8];
    m_elem[11] -= m.m_p.m_elem[2];
    m_elem[15] -= 1;

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator=(const t_type s)
{
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] = s;
        m_elem[i + 2] = s;
        m_elem[i + 4] = s;
        m_elem[i + 6] = s;
        m_elem[i + 1] = s;
        m_elem[i + 3] = s;
        m_elem[i + 5] = s;
        m_elem[i + 7] = s;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] = s;
    }

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator+=(const t_type s)
{
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] += s;
        m_elem[i + 2] += s;
        m_elem[i + 4] += s;
        m_elem[i + 6] += s;
        m_elem[i + 1] += s;
        m_elem[i + 3] += s;
        m_elem[i + 5] += s;
        m_elem[i + 7] += s;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] += s;
    }

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator-=(const t_type s)
{
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] -= s;
        m_elem[i + 2] -= s;
        m_elem[i + 4] -= s;
        m_elem[i + 6] -= s;
        m_elem[i + 1] -= s;
        m_elem[i + 3] -= s;
        m_elem[i + 5] -= s;
        m_elem[i + 7] -= s;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] -= s;
    }

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator*=(const t_type s)
{
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] *= s;
        m_elem[i + 2] *= s;
        m_elem[i + 4] *= s;
        m_elem[i + 6] *= s;
        m_elem[i + 1] *= s;
        m_elem[i + 3] *= s;
        m_elem[i + 5] *= s;
        m_elem[i + 7] *= s;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] *= s;
    }

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> &Matrix<t_row, t_col, t_type>::operator/=(const t_type s)
{
    t_type scalar = s;

    if (std::abs(scalar) < std::numeric_limits<t_type>::epsilon())
    {
        if (scalar < 0) scalar = -std::numeric_limits<t_type>::epsilon();
        else scalar = std::numeric_limits<t_type>::epsilon();
    }

    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        m_elem[i] /= scalar;
        m_elem[i + 2] /= scalar;
        m_elem[i + 4] /= scalar;
        m_elem[i + 6] /= scalar;
        m_elem[i + 1] /= scalar;
        m_elem[i + 3] /= scalar;
        m_elem[i + 5] /= scalar;
        m_elem[i + 7] /= scalar;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        m_elem[i] /= scalar;
    }

    return (*this);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline CommaInit<t_row * t_col, t_type> Matrix<t_row, t_col, t_type>::operator<<(const t_type s)
{
    m_elem[0] = s;
    return CommaInit<t_row * t_col, t_type>(m_elem);
}

/* Arithmetic operators */
template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator-() const
{
    uint16_t cnt, i = 0;
    t_type mat[t_row * t_col];

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = -m_elem[i];
        mat[i + 2] = -m_elem[i + 2];
        mat[i + 4] = -m_elem[i + 4];
        mat[i + 6] = -m_elem[i + 6];
        mat[i + 1] = -m_elem[i + 1];
        mat[i + 3] = -m_elem[i + 3];
        mat[i + 5] = -m_elem[i + 5];
        mat[i + 7] = -m_elem[i + 7];
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = -m_elem[i];
    }

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator+(const Matrix<t_row, t_col, t_type> &m) const
{
    uint16_t cnt, i = 0;
    t_type mat[t_row * t_col];

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] + m.m_elem[i];
        mat[i + 2] = m_elem[i + 2] + m.m_elem[i + 2];
        mat[i + 4] = m_elem[i + 4] + m.m_elem[i + 4];
        mat[i + 6] = m_elem[i + 6] + m.m_elem[i + 6];
        mat[i + 1] = m_elem[i + 1] + m.m_elem[i + 1];
        mat[i + 3] = m_elem[i + 3] + m.m_elem[i + 3];
        mat[i + 5] = m_elem[i + 5] + m.m_elem[i + 5];
        mat[i + 7] = m_elem[i + 7] + m.m_elem[i + 7];
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] + m.m_elem[i];
    }

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator-(const Matrix<t_row, t_col, t_type> &m) const
{
    uint16_t cnt, i = 0;
    t_type mat[t_row * t_col];

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] - m.m_elem[i];
        mat[i + 2] = m_elem[i + 2] - m.m_elem[i + 2];
        mat[i + 4] = m_elem[i + 4] - m.m_elem[i + 4];
        mat[i + 6] = m_elem[i + 6] - m.m_elem[i + 6];
        mat[i + 1] = m_elem[i + 1] - m.m_elem[i + 1];
        mat[i + 3] = m_elem[i + 3] - m.m_elem[i + 3];
        mat[i + 5] = m_elem[i + 5] - m.m_elem[i + 5];
        mat[i + 7] = m_elem[i + 7] - m.m_elem[i + 7];
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] - m.m_elem[i];
    }

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator+(const Matrix<0, 0, t_type> &m) const
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == m.m_row && "Row dimensions do not matched");
    assert(t_col == m.m_col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;
    t_type mat[t_row * t_col];

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] + m.m_elem[i];
        mat[i + 2] = m_elem[i + 2] + m.m_elem[i + 2];
        mat[i + 4] = m_elem[i + 4] + m.m_elem[i + 4];
        mat[i + 6] = m_elem[i + 6] + m.m_elem[i + 6];
        mat[i + 1] = m_elem[i + 1] + m.m_elem[i + 1];
        mat[i + 3] = m_elem[i + 3] + m.m_elem[i + 3];
        mat[i + 5] = m_elem[i + 5] + m.m_elem[i + 5];
        mat[i + 7] = m_elem[i + 7] + m.m_elem[i + 7];
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] + m.m_elem[i];
    }

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator-(const Matrix<0, 0, t_type> &m) const
{
    assert(m.m_elem != nullptr && "Memory has not been allocated");
    assert(t_row == m.m_row && "Row dimensions do not matched");
    assert(t_col == m.m_col && "Col dimensions do not matched");

    uint16_t cnt, i = 0;
    t_type mat[t_row * t_col];

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] - m.m_elem[i];
        mat[i + 2] = m_elem[i + 2] - m.m_elem[i + 2];
        mat[i + 4] = m_elem[i + 4] - m.m_elem[i + 4];
        mat[i + 6] = m_elem[i + 6] - m.m_elem[i + 6];
        mat[i + 1] = m_elem[i + 1] - m.m_elem[i + 1];
        mat[i + 3] = m_elem[i + 3] - m.m_elem[i + 3];
        mat[i + 5] = m_elem[i + 5] - m.m_elem[i + 5];
        mat[i + 7] = m_elem[i + 7] - m.m_elem[i + 7];
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] - m.m_elem[i];
    }

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator+(const Matrix3<t_type, t_row, t_col> &m) const
{
    t_type mat[t_row * t_col];

    mat[0] = m_elem[0] + m.m_elem[0];
    mat[1] = m_elem[1] + m.m_elem[1];
    mat[2] = m_elem[2] + m.m_elem[2];
    mat[3] = m_elem[3] + m.m_elem[3];
    mat[4] = m_elem[4] + m.m_elem[4];
    mat[5] = m_elem[5] + m.m_elem[5];
    mat[6] = m_elem[6] + m.m_elem[6];
    mat[7] = m_elem[7] + m.m_elem[7];
    mat[8] = m_elem[8] + m.m_elem[8];

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator-(const Matrix3<t_type, t_row, t_col> &m) const
{
    t_type mat[t_row * t_col];

    mat[0] = m_elem[0] - m.m_elem[0];
    mat[1] = m_elem[1] - m.m_elem[1];
    mat[2] = m_elem[2] - m.m_elem[2];
    mat[3] = m_elem[3] - m.m_elem[3];
    mat[4] = m_elem[4] - m.m_elem[4];
    mat[5] = m_elem[5] - m.m_elem[5];
    mat[6] = m_elem[6] - m.m_elem[6];
    mat[7] = m_elem[7] - m.m_elem[7];
    mat[8] = m_elem[8] - m.m_elem[8];

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator+(const Rotation<t_type, t_row, t_col> &m) const
{
    t_type mat[t_row * t_col];

    mat[0] = m_elem[0] + m.m_elem[0];
    mat[1] = m_elem[1] + m.m_elem[1];
    mat[2] = m_elem[2] + m.m_elem[2];
    mat[3] = m_elem[3] + m.m_elem[3];
    mat[4] = m_elem[4] + m.m_elem[4];
    mat[5] = m_elem[5] + m.m_elem[5];
    mat[6] = m_elem[6] + m.m_elem[6];
    mat[7] = m_elem[7] + m.m_elem[7];
    mat[8] = m_elem[8] + m.m_elem[8];

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator-(const Rotation<t_type, t_row, t_col> &m) const
{
    t_type mat[t_row * t_col];

    mat[0] = m_elem[0] - m.m_elem[0];
    mat[1] = m_elem[1] - m.m_elem[1];
    mat[2] = m_elem[2] - m.m_elem[2];
    mat[3] = m_elem[3] - m.m_elem[3];
    mat[4] = m_elem[4] - m.m_elem[4];
    mat[5] = m_elem[5] - m.m_elem[5];
    mat[6] = m_elem[6] - m.m_elem[6];
    mat[7] = m_elem[7] - m.m_elem[7];
    mat[8] = m_elem[8] - m.m_elem[8];

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator+(const Transform<t_type, t_row, t_col> &m) const
{
    t_type mat[t_row * t_col];

    mat[0] = m_elem[0] + m.m_R.m_elem[0];
    mat[1] = m_elem[1] + m.m_R.m_elem[1];
    mat[2] = m_elem[2] + m.m_R.m_elem[2];
    mat[3] = m_elem[3] + m.m_p.m_elem[0];
    mat[4] = m_elem[4] + m.m_R.m_elem[3];
    mat[5] = m_elem[5] + m.m_R.m_elem[4];
    mat[6] = m_elem[6] + m.m_R.m_elem[5];
    mat[7] = m_elem[7] + m.m_p.m_elem[1];
    mat[8] = m_elem[8] + m.m_R.m_elem[6];
    mat[9] = m_elem[9] + m.m_R.m_elem[7];
    mat[10] = m_elem[10] + m.m_R.m_elem[8];
    mat[11] = m_elem[11] + m.m_p.m_elem[2];
    mat[12] = m_elem[12];
    mat[13] = m_elem[13];
    mat[14] = m_elem[14];
    mat[15] = m_elem[15] + 1;

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator-(const Transform<t_type, t_row, t_col> &m) const
{
    t_type mat[t_row * t_col];

    mat[0] = m_elem[0] - m.m_R.m_elem[0];
    mat[1] = m_elem[1] - m.m_R.m_elem[1];
    mat[2] = m_elem[2] - m.m_R.m_elem[2];
    mat[3] = m_elem[3] - m.m_p.m_elem[0];
    mat[4] = m_elem[4] - m.m_R.m_elem[3];
    mat[5] = m_elem[5] - m.m_R.m_elem[4];
    mat[6] = m_elem[6] - m.m_R.m_elem[5];
    mat[7] = m_elem[7] - m.m_p.m_elem[1];
    mat[8] = m_elem[8] - m.m_R.m_elem[6];
    mat[9] = m_elem[9] - m.m_R.m_elem[7];
    mat[10] = m_elem[10] - m.m_R.m_elem[8];
    mat[11] = m_elem[11] - m.m_p.m_elem[2];
    mat[12] = m_elem[12];
    mat[13] = m_elem[13];
    mat[14] = m_elem[14];
    mat[15] = m_elem[15] - 1;

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator+(const t_type s) const
{
    t_type mat[t_row * t_col];
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] + s;
        mat[i + 2] = m_elem[i + 2] + s;
        mat[i + 4] = m_elem[i + 4] + s;
        mat[i + 6] = m_elem[i + 6] + s;
        mat[i + 1] = m_elem[i + 1] + s;
        mat[i + 3] = m_elem[i + 3] + s;
        mat[i + 5] = m_elem[i + 5] + s;
        mat[i + 7] = m_elem[i + 7] + s;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] + s;
    }

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator-(const t_type s) const
{
    t_type mat[t_row * t_col];
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] - s;
        mat[i + 2] = m_elem[i + 2] - s;
        mat[i + 4] = m_elem[i + 4] - s;
        mat[i + 6] = m_elem[i + 6] - s;
        mat[i + 1] = m_elem[i + 1] - s;
        mat[i + 3] = m_elem[i + 3] - s;
        mat[i + 5] = m_elem[i + 5] - s;
        mat[i + 7] = m_elem[i + 7] - s;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] - s;
    }

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator*(const t_type s) const
{
    t_type mat[t_row * t_col];
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] * s;
        mat[i + 2] = m_elem[i + 2] * s;
        mat[i + 4] = m_elem[i + 4] * s;
        mat[i + 6] = m_elem[i + 6] * s;
        mat[i + 1] = m_elem[i + 1] * s;
        mat[i + 3] = m_elem[i + 3] * s;
        mat[i + 5] = m_elem[i + 5] * s;
        mat[i + 7] = m_elem[i + 7] * s;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] * s;
    }

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator/(const t_type s) const
{
    t_type mat[t_row * t_col];
    uint16_t cnt, i = 0;
    t_type scalar = s;

    if (std::abs(scalar) < std::numeric_limits<t_type>::epsilon())
    {
        if (scalar < 0) scalar = -std::numeric_limits<t_type>::epsilon();
        else scalar = std::numeric_limits<t_type>::epsilon();
    }

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m_elem[i] / scalar;
        mat[i + 2] = m_elem[i + 2] / scalar;
        mat[i + 4] = m_elem[i + 4] / scalar;
        mat[i + 6] = m_elem[i + 6] / scalar;
        mat[i + 1] = m_elem[i + 1] / scalar;
        mat[i + 3] = m_elem[i + 3] / scalar;
        mat[i + 5] = m_elem[i + 5] / scalar;
        mat[i + 7] = m_elem[i + 7] / scalar;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m_elem[i] / scalar;
    }

    return Matrix(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
template <uint16_t col>
inline Matrix<t_row, col, t_type> Matrix<t_row, t_col, t_type>::operator*(const Matrix<t_col, col, t_type> &m) const
{
    t_type mat[t_row * col]{0};
    uint16_t cnt;
    uint16_t irow, icol, i;

    for (irow = 0; irow < t_row; ++irow)
    {
        for (icol = 0; icol < col; ++icol)
        {
            for (cnt = t_col >> 2u, i = 0; cnt > 0u; cnt--, i += 4u)
            {
                mat[irow * col + icol] += m_elem[irow * t_col + i] * m.m_elem[i * col + icol];
                mat[irow * col + icol] += m_elem[irow * t_col + i + 1] * m.m_elem[(i + 1) * col + icol];
                mat[irow * col + icol] += m_elem[irow * t_col + i + 2] * m.m_elem[(i + 2) * col + icol];
                mat[irow * col + icol] += m_elem[irow * t_col + i + 3] * m.m_elem[(i + 3) * col + icol];
            }

            for (cnt = t_col % 4u; cnt > 0u; cnt--, i++)
            {
                mat[irow * col + icol] += m_elem[irow * t_col + i] * m.m_elem[i * col + icol];
            }
        }
    }

    return Matrix<t_row, col, t_type>(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<0, 0, t_type> Matrix<t_row, t_col, t_type>::operator*(const Matrix<0, 0, t_type> &m) const
{
    assert(t_col == m.m_row && "Dimensions do not matched");

    uint16_t cnt;
    uint16_t irow, icol, i;
    Matrix<0, 0, t_type> mat(t_row, m.m_col);

    for (irow = 0; irow < t_row; ++irow)
    {
        for (icol = 0; icol < m.m_col; ++icol)
        {
            for (cnt = t_col >> 2u, i = 0; cnt > 0u; cnt--, i += 4u)
            {
                mat[irow * m.m_col + icol] += m_elem[irow * t_col + i] * m.m_elem[i * m.m_col + icol];
                mat[irow * m.m_col + icol] += m_elem[irow * t_col + i + 1] * m.m_elem[(i + 1) * m.m_col + icol];
                mat[irow * m.m_col + icol] += m_elem[irow * t_col + i + 2] * m.m_elem[(i + 2) * m.m_col + icol];
                mat[irow * m.m_col + icol] += m_elem[irow * t_col + i + 3] * m.m_elem[(i + 3) * m.m_col + icol];
            }

            for (cnt = t_col % 4u; cnt > 0u; cnt--, i++)
            {
                mat[irow * m.m_col + icol] += m_elem[irow * t_col + i] * m.m_elem[i * m.m_col + icol];
            }
        }
    }

    return mat;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator*(const Matrix3<t_type, t_col, t_col> &m) const
{
    t_type mat[t_row * t_col]{0};

    for (uint16_t irow = 0; irow < t_row; ++irow)
    {
        for (uint16_t icol = 0; icol < t_col; ++icol)
        {
            mat[irow * t_col + icol] += m_elem[irow * t_col] * m.m_elem[icol];
            mat[irow * t_col + icol] += m_elem[irow * t_col + 1] * m.m_elem[3 + icol];
            mat[irow * t_col + icol] += m_elem[irow * t_col + 2] * m.m_elem[6 + icol];
        }
    }

    return Matrix<t_row, t_col, t_type>(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator*(const Rotation<t_type, t_col, t_col> &m) const
{
    t_type mat[t_row * t_col]{0};

    for (uint16_t irow = 0; irow < t_row; ++irow)
    {
        for (uint16_t icol = 0; icol < t_col; ++icol)
        {
            mat[irow * t_col + icol] += m_elem[irow * t_col] * m.m_elem[icol];
            mat[irow * t_col + icol] += m_elem[irow * t_col + 1] * m.m_elem[3 + icol];
            mat[irow * t_col + icol] += m_elem[irow * t_col + 2] * m.m_elem[6 + icol];
        }
    }

    return Matrix<t_row, t_col, t_type>(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator*(const Transform<t_type, t_col, t_col> &m) const
{
    t_type mat[t_row * t_col]{0};

    for (uint16_t irow = 0; irow < t_row * 4; irow += 4)
    {
        for (uint16_t icol = 0; icol < 3; ++icol)
        {
            mat[irow + icol] =
                m_elem[irow] * m.m_R.m_elem[icol] +
                m_elem[irow + 1] * m.m_R.m_elem[icol + 3] +
                m_elem[irow + 2] * m.m_R.m_elem[icol + 6];
        }
        mat[irow + 3] =
            m_elem[irow] * m.m_p.m_elem[0] +
            m_elem[irow + 1] * m.m_p.m_elem[1] +
            m_elem[irow + 2] * m.m_p.m_elem[2] +
            m_elem[irow + 3];
    }

    return Matrix<t_row, t_col, t_type>(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> Matrix<t_row, t_col, t_type>::operator*(const Vector<t_col, t_type> &v) const
{
    t_type vec[t_row]{0};
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; ++irow)
    {
        for (cnt = t_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
        {
            vec[irow] += m_elem[irow * t_col + icol] * v.m_elem[icol];
            vec[irow] += m_elem[irow * t_col + icol + 1] * v.m_elem[icol + 1];
            vec[irow] += m_elem[irow * t_col + icol + 2] * v.m_elem[icol + 2];
            vec[irow] += m_elem[irow * t_col + icol + 3] * v.m_elem[icol + 3];
        }

        for (cnt = t_col % 4u; cnt > 0u; cnt--, icol++)
            vec[irow] += m_elem[irow * t_col + icol] * v.m_elem[icol];
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> Matrix<t_row, t_col, t_type>::operator*(const Vector<0, t_type> &v) const
{
    assert(v.m_elem != nullptr && "Memory has not been allocated");
    assert(t_col == v.m_row && "Check dimensions");

    t_type vec[t_row]{0};
    uint16_t cnt;
    uint16_t irow, icol;

    for (irow = 0; irow < t_row; ++irow)
    {
        for (cnt = t_col >> 2u, icol = 0; cnt > 0u; cnt--, icol += 4u)
        {
            vec[irow] += m_elem[irow * t_col + icol] * v.m_elem[icol];
            vec[irow] += m_elem[irow * t_col + icol + 1] * v.m_elem[icol + 1];
            vec[irow] += m_elem[irow * t_col + icol + 2] * v.m_elem[icol + 2];
            vec[irow] += m_elem[irow * t_col + icol + 3] * v.m_elem[icol + 3];
        }

        for (cnt = t_col % 4u; cnt > 0u; cnt--, icol++)
            vec[irow] += m_elem[irow * t_col + icol] * v.m_elem[icol];
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> Matrix<t_row, t_col, t_type>::operator*(const Vector3<t_type, t_col> &v) const
{
    t_type vec[t_row];

    for (uint16_t irow = 0; irow < t_row; ++irow)
    {
        vec[irow] = m_elem[irow * 3] * v.m_elem[0];
        vec[irow] += m_elem[irow * 3 + 1] * v.m_elem[1];
        vec[irow] += m_elem[irow * 3 + 2] * v.m_elem[2];
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> Matrix<t_row, t_col, t_type>::operator*(const Vector4<t_type, t_col> &v) const
{
    t_type vec[t_row];

    for (uint16_t irow = 0; irow < t_row; ++irow)
    {
        vec[irow] = m_elem[irow * 4] * v.m_elem[0];
        vec[irow] += m_elem[irow * 4 + 1] * v.m_elem[1];
        vec[irow] += m_elem[irow * 4 + 2] * v.m_elem[2];
        vec[irow] += m_elem[irow * 4 + 3] * v.m_elem[3];
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> Matrix<t_row, t_col, t_type>::operator*(const Quaternion<t_type, t_col> &v) const
{
    t_type vec[t_row];

    for (uint16_t irow = 0; irow < t_row; ++irow)
    {
        vec[irow] = m_elem[irow * 4] * v.m_elem[0];
        vec[irow] += m_elem[irow * 4 + 1] * v.m_elem[1];
        vec[irow] += m_elem[irow * 4 + 2] * v.m_elem[2];
        vec[irow] += m_elem[irow * 4 + 3] * v.m_elem[3];
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Vector<t_row, t_type> Matrix<t_row, t_col, t_type>::operator*(const Vector6<t_type, t_col> &v) const
{
    t_type vec[t_row];

    for (uint16_t irow = 0; irow < t_row; ++irow)
    {
        vec[irow] = m_elem[irow * 6] * v.m_elem[0];
        vec[irow] += m_elem[irow * 6 + 1] * v.m_elem[1];
        vec[irow] += m_elem[irow * 6 + 2] * v.m_elem[2];
        vec[irow] += m_elem[irow * 6 + 3] * v.m_elem[3];
        vec[irow] += m_elem[irow * 6 + 4] * v.m_elem[4];
        vec[irow] += m_elem[irow * 6 + 5] * v.m_elem[5];
    }

    return Vector<t_row, t_type>(vec);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator&(const Vector<t_col, t_type> &v) const
{ // matrix3 * [v]x, []x is skew-symmetric matrix
    static_assert(t_row == 3, "This method is only for 3 x 3 matrix");
    static_assert(t_col == 3, "This method is only for 3 x 3 matrix");

    t_type mat[t_row * t_col];

    mat[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    mat[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    mat[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    mat[3] = m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1];
    mat[4] = m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2];
    mat[5] = m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0];

    mat[6] = m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1];
    mat[7] = m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2];
    mat[8] = m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0];

    return Matrix<t_row, t_col, t_type>(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator&(const Vector<0, t_type> &v) const
{ // matrix3 * [v]x, []x is skew-symmetric matrix
    static_assert(t_row == 3, "This method is only for 3 x 3 matrix");
    static_assert(t_col == 3, "This method is only for 3 x 3 matrix");
    assert(v.m_row == 3 && "This method is only for 3 x 3 matrix");
    assert(v.m_elem != nullptr && "Memory has not been allocated");

    t_type mat[t_row * t_col];

    mat[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    mat[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    mat[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    mat[3] = m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1];
    mat[4] = m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2];
    mat[5] = m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0];

    mat[6] = m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1];
    mat[7] = m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2];
    mat[8] = m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0];

    return Matrix<t_row, t_col, t_type>(mat);
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline Matrix<t_row, t_col, t_type> Matrix<t_row, t_col, t_type>::operator&(const Vector3<t_type, t_col> &v) const
{ // matrix3 * [v]x, []x is skew-symmetric matrix
    static_assert(t_row == 3, "This method is only for 3 x 3 matrix");
    static_assert(t_col == 3, "This method is only for 3 x 3 matrix");

    t_type mat[t_row * t_col];

    mat[0] = m_elem[1] * v.m_elem[2] - m_elem[2] * v.m_elem[1];
    mat[1] = m_elem[2] * v.m_elem[0] - m_elem[0] * v.m_elem[2];
    mat[2] = m_elem[0] * v.m_elem[1] - m_elem[1] * v.m_elem[0];

    mat[3] = m_elem[4] * v.m_elem[2] - m_elem[5] * v.m_elem[1];
    mat[4] = m_elem[5] * v.m_elem[0] - m_elem[3] * v.m_elem[2];
    mat[5] = m_elem[3] * v.m_elem[1] - m_elem[4] * v.m_elem[0];

    mat[6] = m_elem[7] * v.m_elem[2] - m_elem[8] * v.m_elem[1];
    mat[7] = m_elem[8] * v.m_elem[0] - m_elem[6] * v.m_elem[2];
    mat[8] = m_elem[6] * v.m_elem[1] - m_elem[7] * v.m_elem[0];

    return Matrix<t_row, t_col, t_type>(mat);
}

/* Comparison operators */
template <uint16_t t_row, uint16_t t_col, typename t_type>
inline bool Matrix<t_row, t_col, t_type>::operator==(const Matrix<t_row, t_col, t_type> &m) const
{
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - m.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 4] - m.m_elem[i + 4]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 6] - m.m_elem[i + 6]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - m.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - m.m_elem[i + 3]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 5] - m.m_elem[i + 5]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 7] - m.m_elem[i + 7]) > m_tolerance) return false;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline bool Matrix<t_row, t_col, t_type>::operator!=(const Matrix<t_row, t_col, t_type> &m) const
{
    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - m.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 4] - m.m_elem[i + 4]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 6] - m.m_elem[i + 6]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - m.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - m.m_elem[i + 3]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 5] - m.m_elem[i + 5]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 7] - m.m_elem[i + 7]) > m_tolerance) return true;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline bool Matrix<t_row, t_col, t_type>::operator==(const Matrix<0, 0, t_type> &m) const
{
    assert(t_row == m.m_row && "Row Dimensions do not matched");
    assert(t_col == m.m_col && "Col Dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 2] - m.m_elem[i + 2]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 4] - m.m_elem[i + 4]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 6] - m.m_elem[i + 6]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 1] - m.m_elem[i + 1]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 3] - m.m_elem[i + 3]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 5] - m.m_elem[i + 5]) > m_tolerance) return false;
        if (std::abs(m_elem[i + 7] - m.m_elem[i + 7]) > m_tolerance) return false;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return false;
    }

    return true;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline bool Matrix<t_row, t_col, t_type>::operator!=(const Matrix<0, 0, t_type> &m) const
{
    assert(t_row == m.m_row && "Row Dimensions do not matched");
    assert(t_col == m.m_col && "Col Dimensions do not matched");

    uint16_t cnt, i = 0;

    for (cnt = (t_row * t_col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 2] - m.m_elem[i + 2]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 4] - m.m_elem[i + 4]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 6] - m.m_elem[i + 6]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 1] - m.m_elem[i + 1]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 3] - m.m_elem[i + 3]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 5] - m.m_elem[i + 5]) > m_tolerance) return true;
        if (std::abs(m_elem[i + 7] - m.m_elem[i + 7]) > m_tolerance) return true;
    }

    for (cnt = (t_row * t_col) % 8u; cnt > 0u; cnt--, i++)
    {
        if (std::abs(m_elem[i] - m.m_elem[i]) > m_tolerance) return true;
    }

    return false;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline bool Matrix<t_row, t_col, t_type>::operator==(const Matrix3<t_type, t_row, t_col> &m) const
{
    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance) return false;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance) return false;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance) return false;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance) return false;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance) return false;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance) return false;

    return true;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline bool Matrix<t_row, t_col, t_type>::operator!=(const Matrix3<t_type, t_row, t_col> &m) const
{
    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance) return true;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance) return true;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance) return true;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance) return true;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance) return true;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance) return true;

    return false;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline bool Matrix<t_row, t_col, t_type>::operator==(const Rotation<t_type, t_row, t_col> &m) const
{
    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance) return false;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance) return false;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance) return false;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance) return false;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance) return false;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance) return false;

    return true;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline bool Matrix<t_row, t_col, t_type>::operator!=(const Rotation<t_type, t_row, t_col> &m) const
{
    if (std::abs(m_elem[0] - m.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - m.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - m.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - m.m_elem[3]) > m_tolerance) return true;
    if (std::abs(m_elem[4] - m.m_elem[4]) > m_tolerance) return true;
    if (std::abs(m_elem[5] - m.m_elem[5]) > m_tolerance) return true;
    if (std::abs(m_elem[6] - m.m_elem[6]) > m_tolerance) return true;
    if (std::abs(m_elem[7] - m.m_elem[7]) > m_tolerance) return true;
    if (std::abs(m_elem[8] - m.m_elem[8]) > m_tolerance) return true;

    return false;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline bool Matrix<t_row, t_col, t_type>::operator==(const Transform<t_type, t_row, t_col> &m) const
{
    if (std::abs(m_elem[0] - m.m_R.m_elem[0]) > m_tolerance) return false;
    if (std::abs(m_elem[1] - m.m_R.m_elem[1]) > m_tolerance) return false;
    if (std::abs(m_elem[2] - m.m_R.m_elem[2]) > m_tolerance) return false;
    if (std::abs(m_elem[3] - m.m_p.m_elem[0]) > m_tolerance) return false;

    if (std::abs(m_elem[4] - m.m_R.m_elem[3]) > m_tolerance) return false;
    if (std::abs(m_elem[5] - m.m_R.m_elem[4]) > m_tolerance) return false;
    if (std::abs(m_elem[6] - m.m_R.m_elem[5]) > m_tolerance) return false;
    if (std::abs(m_elem[7] - m.m_p.m_elem[1]) > m_tolerance) return false;

    if (std::abs(m_elem[8] - m.m_R.m_elem[6]) > m_tolerance) return false;
    if (std::abs(m_elem[9] - m.m_R.m_elem[7]) > m_tolerance) return false;
    if (std::abs(m_elem[10] - m.m_R.m_elem[8]) > m_tolerance) return false;
    if (std::abs(m_elem[11] - m.m_p.m_elem[2]) > m_tolerance) return false;

    if (std::abs(m_elem[12]) > m_tolerance) return false;
    if (std::abs(m_elem[13]) > m_tolerance) return false;
    if (std::abs(m_elem[14]) > m_tolerance) return false;
    if (std::abs(1 - m_elem[15]) > m_tolerance) return false;

    return true;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline bool Matrix<t_row, t_col, t_type>::operator!=(const Transform<t_type, t_row, t_col> &m) const
{
    if (std::abs(m_elem[0] - m.m_R.m_elem[0]) > m_tolerance) return true;
    if (std::abs(m_elem[1] - m.m_R.m_elem[1]) > m_tolerance) return true;
    if (std::abs(m_elem[2] - m.m_R.m_elem[2]) > m_tolerance) return true;
    if (std::abs(m_elem[3] - m.m_p.m_elem[0]) > m_tolerance) return true;

    if (std::abs(m_elem[4] - m.m_R.m_elem[3]) > m_tolerance) return true;
    if (std::abs(m_elem[5] - m.m_R.m_elem[4]) > m_tolerance) return true;
    if (std::abs(m_elem[6] - m.m_R.m_elem[5]) > m_tolerance) return true;
    if (std::abs(m_elem[7] - m.m_p.m_elem[1]) > m_tolerance) return true;

    if (std::abs(m_elem[8] - m.m_R.m_elem[6]) > m_tolerance) return true;
    if (std::abs(m_elem[9] - m.m_R.m_elem[7]) > m_tolerance) return true;
    if (std::abs(m_elem[10] - m.m_R.m_elem[8]) > m_tolerance) return true;
    if (std::abs(m_elem[11] - m.m_p.m_elem[2]) > m_tolerance) return true;

    if (std::abs(m_elem[12]) > m_tolerance) return true;
    if (std::abs(m_elem[13]) > m_tolerance) return true;
    if (std::abs(m_elem[14]) > m_tolerance) return true;
    if (std::abs(1 - m_elem[15]) > m_tolerance) return true;

    return false;
}

template <uint16_t t_row, uint16_t t_col, typename t_type>
inline void Matrix<t_row, t_col, t_type>::Print(const char endChar)
{
#if defined(ARDUINO)
    for (uint16_t irow = 0; irow < t_row; irow++)
    {
        for (uint16_t icol = 0; icol < t_col; icol++)
        {
            Serial.printf("%7.3f ", (t_type)(m_elem[irow * t_col + icol]));
        }
        Serial.write('\n');
    }
    Serial.write(endChar);
#else
    for (uint16_t irow = 0; irow < t_row; irow++)
    {
        for (uint16_t icol = 0; icol < t_col; icol++)
        {
            printf("%10.6f ", (t_type)(m_elem[irow * t_col + icol]));
        }
        printf("\n");
    }
    printf("%c", endChar);
#endif
}

//-- Template Function ------------------------------------------------------//
// scalar * matrix
template <uint16_t row, uint16_t col, typename type>
inline Matrix<row, col, type> operator*(const type s, const Matrix<row, col, type> &m)
{
    type mat[row * col];
    uint16_t cnt, i = 0;

    for (cnt = (row * col) >> 3u; cnt > 0u; cnt--, i += 8)
    {
        mat[i] = m.m_elem[i] * s;
        mat[i + 2] = m.m_elem[i + 2] * s;
        mat[i + 4] = m.m_elem[i + 4] * s;
        mat[i + 6] = m.m_elem[i + 6] * s;
        mat[i + 1] = m.m_elem[i + 1] * s;
        mat[i + 3] = m.m_elem[i + 3] * s;
        mat[i + 5] = m.m_elem[i + 5] * s;
        mat[i + 7] = m.m_elem[i + 7] * s;
    }

    for (cnt = (row * col) % 8u; cnt > 0u; cnt--, i++)
    {
        mat[i] = m.m_elem[i] * s;
    }

    return Matrix<row, col, type>(mat);
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTMATRIX_TPP_
