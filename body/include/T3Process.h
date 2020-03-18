#pragma once
#ifndef T3PROCESS_H
#define T3PROCESS_H

#include <array>
#include <cstdint>
#include <memory>
#include <tuple>
#include <utility>

#include "T3Defs.h"
#include "T3LorentzVector.h"
#include "T3Particle.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

#include "T3AllocateData.h"

using namespace data;

namespace t3 {

template <typename ProcessImpl> class Process {
public:
  using Base_t = ProcessImpl;

  template <typename... Args>
  Process(Args... args) : fProcessImpl(ProcessImpl(args...)) {}
  //Process &operator=(Process &) = default;
  
  template <typename Floating>
  void GetCS(MatID_t matID, int N,
             Floating * outputCSArray,
             ParticleTable & aParticleTable,
             MaterialTable & aMaterialTable) const;

  void GetFS(MatID_t matID, int N, int SHIFT,
             ParticleTable & aParticleTable,
             MaterialTable & aMaterialTable) const;

private:
  ProcessImpl fProcessImpl;
};

//______________________________________________________________________________
template <class ProcessImpl>
template <typename Floating>
void Process<ProcessImpl>::GetCS(MatID_t matID, int N,
                                 Floating * outputCSArray,
                                 ParticleTable & aParticleTable,
                                 MaterialTable & aMaterialTable) const
{         
//
#ifdef OPENACC
//without present(particles) fails at runtime with an error:
//Failing in Thread:1
//call to cuStreamSynchronize returned error 700: Illegal address
//during kernel execution.
#pragma acc parallel loop gang vector present(particles,outputCSArray,aParticleTable,aMaterialTable)
#else
#pragma omp parallel for simd
#endif
//
  for (int ind=0; ind < N; ++ind)
  {
    outputCSArray[ind] = fProcessImpl.GetCS(ind, matID, aParticleTable, aMaterialTable);
  }
}  

template <class ProcessImpl>
void Process<ProcessImpl>::GetFS(
    MatID_t matID, int N, int SHIFT,
    ParticleTable & aParticleTable,
    MaterialTable & aMaterialTable) const
{
  
//
#ifdef OPENACC
#pragma acc parallel loop gang vector present(particles,aParticleTable,aMaterialTable)
#else 
#pragma omp parallel for simd
#endif  
//
  for(int ind=0; ind < N; ++ind)
  {
    fProcessImpl.GetFS(ind, SHIFT, matID, aParticleTable, aMaterialTable);
  }
}

} // namespace t3

#endif // T3PROCESS_H
