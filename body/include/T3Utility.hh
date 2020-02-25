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


  //this function should generate random x from the normalized (!!!)
  //Gaussian normal distribution (taken from random.pdf): 
  inline T3double RandomizeNormalDistribution(unsigned int & generator_seed,/*random number generator seed*/
                                       T3double mu,/*mean value*/
                                       T3double sigma/*root mean square - SQRT(dispersion)*/)
  {
    if(sigma<0.0) printf("***ERROR: T3Utility::RandomizeNormalDistribution(...): sigma=%f < 0!", sigma);
    T3double p, p1, p2;
    do {
       p1=2.0*RND01(generator_seed)-1.0;
       p2=2.0*RND01(generator_seed)-1.0;
       p=p1*p1+p2*p2;
    } while(p>=1.0);
    //std::cout<<"p1="<<p1<<" p2="<<p2<<" p="<<p<<std::endl;
    T3double res=mu + sigma * p1 * sqrt(-2.0 * log(p) / p);
    //added this not to make very small losses:
    if(res<-3*sigma) res=-3*sigma;
    return res;
  }
  

}

#endif


