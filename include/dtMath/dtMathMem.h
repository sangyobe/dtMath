/*!
\file       dtMathMem.h
\brief      dtMath, Dynamic Allocation Management Functions
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 03. 29
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATHMEM_H_
#define DTMATH_DTMATHMEM_H_

#include <cstdlib>

namespace dt
{
namespace Math
{

template <typename t_type> inline void CheckSizeOverflow(std::size_t num);
template <typename t_type> inline t_type *MemAlloc(std::size_t num);
template <typename t_type> inline t_type *MemAllocZeroInit(std::size_t num);
template <typename t_type> inline t_type *MemReAlloc(t_type *ptr, std::size_t num);
template <typename t_type> inline void MemFree(t_type *ptr);
template <typename t_type> inline size_t MemSize(t_type *ptr);

} // namespace Math
} // namespace dt

#include "dtMathMem.ipp"

#endif // DTMATH_DTMATHMEM_H_