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

//KNL ip: 10.41.12.79.

//#define PROBLEM

constexpr bool HISTOGRAM=false;//true;//false;//true;

namespace t3{

  using FloatingType=double;
  using RandGen=t3::RNDGenerator;


  const unsigned int GL=10000000;//For CPU Intel Core i7.
  //const unsigned int GL=1000000;//For GPU GeForce 650 Ti.
  constexpr bool report = true;//print output on main() of nbody.cpp.
  constexpr long int G = 27;//2000;
  constexpr int N = 99999999;//1999999999;//4;//999999999;
  constexpr int Np = N+1;
  constexpr int INJ = 500000;//500000;//100000;//1000;
  constexpr long int K = N+1;//GL;
  constexpr unsigned int max_loop = 10;

  constexpr int NNN=1;
  
  constexpr int cuba = NNN;//16;
  constexpr int cubn = cuba + cuba;
  constexpr int cub2 = cubn * cubn;
  constexpr int cub3 = cub2 * cubn;
  constexpr FloatingType TLS=0.1*units::MeV;
  constexpr double fuse = .25;

  constexpr unsigned int Nbin = 64;//8;//64;//For GeForce GTX 650 Ti GPU.
  constexpr int DN = (N+1) / Nbin +1;
  constexpr unsigned int BLt = GL/Nbin;//200000000/Nbin;//GL/Nbin;//2000 / Nbin;
  constexpr unsigned int GL1 = BLt - 1;
  constexpr double cgam = 5.0;
  constexpr auto dcgam = cgam + cgam;

  constexpr int BinNumber1=100;
  constexpr int BinNumber2=200;

  const FloatingType TARGET_WIDTH = 1.0 * units::um;

  constexpr FloatingType ag = TARGET_WIDTH/2/NNN; //10.0 * units::um;//the width of the cell (voxel).
  constexpr FloatingType InitParticlex0 = 0.5;
  constexpr FloatingType InitParticley0 = 0.5;

  const FloatingType Tls_CUT=10.0 * units::keV;
  
  typedef int PDG_t;
  typedef int MatID_t;
  
}//end of namespace t3.

#endif//T3GLOBALS_H

//Beginning with PGI 16.1 pgcollect has been superseded. Its functionally is now part of PGPROF. 
//To profile an application, enter this command: 
//pgprof -o a.prof ./a.out 
//To read the profile data, enter this command: 
//pgprof -i a.prof

//pgprof --cpu-profiling-mode top-down -i test.prof

//Can not use pgprof GUI in the same GPU without sudo:
//==36714== Warning: The user does not have permission to profile on the target device. See the following link for instructions to enable permissions and get more information: https://developer.nvidia.com/NVSOLN1000 

