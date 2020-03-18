#pragma once
#ifndef T3_UTILITY_HH
#define T3_UTILITY_HH

#include "T3Types.h"

#include <type_traits>
#include <limits>

namespace T3Utility
{
  template<size_t size> T3double* bin_search(T3double* ptr, T3double val);
  template<size_t size> const T3double* bin_search(const T3double* ptr, T3double val);
  template<size_t size> inline T3double* bin_search(T3double* ptr, T3double val)
  {
    return val > ptr[size/2] ?  bin_search<size-size/2 - 1>(ptr+size/2 + 1, val) :
                                bin_search<size/2>(ptr, val);
  }

  template<size_t size> inline const T3double* bin_search(const T3double* ptr,
                                                          T3double val)
  {
    return val > ptr[size/2] ?  bin_search<size-size/2 - 1>(ptr+size/2 + 1, val) :
    bin_search<size/2>(ptr, val);
  }

  template<> inline const T3double* bin_search<0>(const T3double* ptr, T3double)
  {
    return ptr;
  }

  template<> inline T3double* bin_search<0>(T3double* ptr, T3double) {return ptr;}

}

#endif


