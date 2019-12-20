#pragma once
#ifndef T3GLOBALS_H
#define T3GLOBALS_H

#include "T3Defs.h"

#ifdef OPENACC
#pragma acc routine seq
#endif
inline unsigned int Rand32(unsigned int & xn)
{
  u_quad_t a=0x5DEECE66D;
  u_quad_t c=0xB;
  return (unsigned int)((a*xn+c) & 0xFFFFFFFF);
}

#ifdef OPENACC
#pragma acc routine seq
#endif
inline double rndv(unsigned int xn)
{
  return (double) xn / (double) 0x100000000LL;
}

#ifdef OPENACC
#pragma acc routine seq
#endif
inline double RND01(unsigned int & xn)
{
  xn=Rand32(xn);
  return (double)(xn) / (double) 0x100000000LL;
}

typedef int PDG_t;
typedef int MatID_t;

#endif//T3GLOBALS_H
