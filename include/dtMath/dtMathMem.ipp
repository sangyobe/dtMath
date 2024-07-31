/*!
\file       dtMathMem.h
\brief      dtMath, Dynamic Allocation Management Functions
\author     Dong-hyun Lee, phenom8305@gmail.com
\author     Who is next author?
\date       Last modified on 2024. 03. 29
\version    1.0.0
\warning    Do Not delete this comment for document history! This is minimal manners!
*/

#ifndef DTMATH_DTMATHMEM_IPP_
#define DTMATH_DTMATHMEM_IPP_

#include "dtMathMem.h"

#include <cstdlib>  // malloc, calloc, realloc, free
#if defined(__APPLE__)
#include <malloc/malloc.h>
#else
#include <malloc.h> // malloc_usable_size
#endif
#include <new>      // bad_alloc

namespace dt
{
namespace Math
{

template <typename t_type>
inline void CheckSizeOverflow(std::size_t num)
{
    if (num > std::size_t(-1) / sizeof(t_type))
        throw std::bad_alloc();
}

template <typename t_type>
inline t_type *MemAlloc(std::size_t num)
{
    if (num == 0) return 0;
    CheckSizeOverflow<t_type>(num);
    void *rtn = std::malloc(sizeof(t_type) * num);
    if (!rtn) throw std::bad_alloc();

    return static_cast<t_type *>(rtn);
}

template <typename t_type>
inline t_type *MemAllocZeroInit(std::size_t num)
{
    if (num == 0) return 0;
    CheckSizeOverflow<t_type>(num);
    void *rtn = std::calloc(num, sizeof(t_type));
    if (!rtn) throw std::bad_alloc();

    return static_cast<t_type *>(rtn);
}

template <typename t_type>
inline t_type *MemReAlloc(t_type *ptr, std::size_t num)
{
    if (num == 0) return 0;
    CheckSizeOverflow<t_type>(num);
    void *rtn = std::realloc(ptr, sizeof(t_type) * num);
    if (!rtn) throw std::bad_alloc();

    return static_cast<t_type *>(rtn);
}

template <typename t_type>
inline void MemFree(t_type *ptr)
{
    std::free(ptr);
}

template <typename t_type>
inline size_t MemSize(t_type *ptr)
{
#if defined(__APPLE__)
    return malloc_size((void *)ptr);
#else
    return malloc_usable_size((void *)ptr);
#endif
}

} // namespace Math
} // namespace dt

#endif // DTMATH_DTMATHMEM_IPP_