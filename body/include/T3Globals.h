#pragma once
#ifndef T3GLOBALS_H
#define T3GLOBALS_H

#include "T3Defs.h"
#include "T3RNG.h"

#ifdef OPENACC
#include <accelmath.h>
#include <openacc.h>
#ifdef CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#else
#include <omp.h>
  struct float3{
    float x, y, z;
  };
  struct float4{
    float x, y, z, w;
  };
#endif

constexpr bool HISTOGRAM=true;

// flag to open debug of random number generator seed
// initialization in Inject() in T3DataHolder.h:
//
// #define DEBSEED

namespace t3{

  using FloatingType=double;
  using RandGen=t3::RNDGenerator;

  const unsigned int GL=1100000;
  constexpr bool report = true;
  constexpr long int G = 27;
  constexpr int N = 999999;
  constexpr int Np = N+1;
  constexpr int INJ = 10000;
  constexpr long int K = N+1;
  constexpr unsigned int max_loop = 10;
  constexpr int cuba = 16;
  constexpr int cubn = cuba + cuba;
  constexpr int cub2 = cubn * cubn;
  constexpr int cub3 = cub2 * cubn;
  constexpr FloatingType TLS=1.0*units::MeV;
  constexpr double fuse = .25;

  constexpr unsigned int Nbin = 8;
  constexpr int DN = (N+1) / Nbin +1;
  constexpr unsigned int BLt = GL/Nbin;
  constexpr unsigned int GL1 = BLt - 1;
  constexpr double cgam = 5.0;
  constexpr auto dcgam = cgam + cgam;

  constexpr int BinNumber1=100;
  constexpr int BinNumber2=200;
  constexpr FloatingType ag = 1.0e-5 * units::cm;
  constexpr FloatingType InitParticlex0 = 0.5;
  constexpr FloatingType InitParticley0 = 0.5;
  
  typedef int PDG_t;
  typedef int MatID_t;
  
}//end of namespace t3.

#endif//T3GLOBALS_H
