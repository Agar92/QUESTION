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

namespace t3 {

template <typename ProcessImpl> class Process {
public:
  using Base_t = ProcessImpl;
  //Process(Process &&) = default;
  //Process(Process const &) = default;
  //Process(ProcessImpl &&aProcessImpl)
  //    : fProcessImpl(ProcessImpl(aProcessImpl)) {}

//since in the definition of multiplescatteringProcess
//anonymous temporary r-value objects are created,
//here: Process(Args... args)
//l-value-reference can not be used (l-value-reference to anb l-value is an error.)
//HERE WE CAN PASS THE ARGUMENTS ONLY EITHER BY VALUE:
// Process(Args... args)
//OR BY R-VALUE REFERENCE:
// Process(Args... args) : fProcessImpl(ProcessImpl(std::forward<Args>(args)...)) {}

  template <typename... Args>
  Process(Args... args) : fProcessImpl(ProcessImpl(args...)) {}
  //Process &operator=(Process &) = default;
  template <typename Floating>
  auto GetCS(Particle<Floating> * particles, MatID_t matID, uint64_t N,
             Floating * outputCSArray, ParticleTable & aParticleTable) const
      -> bool;
  template <typename Floating>
  auto GetCS(Particle<Floating> * particles, MatID_t matID, uint64_t N,
             Floating * outputCSArray,
             ParticleTable & aParticleTable,
             MaterialTable & aMaterialTable) const
      -> bool;
  //1 secondary particle:
  template <typename Floating>
  auto GetFS(Particle<Floating> * particles, MatID_t matID, uint64_t N,
             PDG_t * outPDG1, T3LorentzVector<Floating> * outP1,
             Floating * csBorderDataFS,
             ParticleTable & aParticleTable,
             /*CSBorderStep is for using csBorderDataCS/FS in nbody.cpp=max number of elements in materials*/
             MaterialTable & aMaterialTable, int CSBorderStep) const
      -> bool;
  //2 secondary particles:
  template <typename Floating>
  auto GetFS(Particle<Floating> * particles, MatID_t matID, uint64_t N,
             PDG_t * outPDG1, T3LorentzVector<Floating> * outP1, PDG_t * outPDG2,
             T3LorentzVector<Floating> * outP2, Floating * csBorderDataFS,
             ParticleTable & aParticleTable,
             /*CSBorderStep is for using csBorderDataCS/FS in nbody.cpp=max number of elements in materials*/
             MaterialTable & aMaterialTable, int CSBorderStep) const
      -> bool;

private:
  ProcessImpl fProcessImpl;
};

//______________________________________________________________________________
template <class ProcessImpl>
template <typename Floating>
auto Process<ProcessImpl>::GetCS(Particle<Floating> * particles,
                                 MatID_t matID, uint64_t N,
                                 Floating * outputCSArray,
                                 ParticleTable & aParticleTable) const
  -> bool {         
//
#ifdef OPENACC
//without present(particles) fails at runtime with an error:
//Failing in Thread:1
//call to cuStreamSynchronize returned error 700: Illegal address
//during kernel execution.
#pragma acc parallel loop  gang vector present(particles,outputCSArray,aParticleTable)
#else
#pragma omp parallel for simd
#endif
//
  for (auto ind = 0u; ind < N; ++ind) {
    outputCSArray[ind] = fProcessImpl.GetCS(particles[ind].p.E()-aParticleTable.GetMass(particles[ind].pdg), particles[ind].pdg, matID, aParticleTable);
  }
  return true;
}

template <class ProcessImpl>
template <typename Floating>
auto Process<ProcessImpl>::GetCS(Particle<Floating> * particles,
                                 MatID_t matID, uint64_t N,
                                 Floating * outputCSArray,
                                 ParticleTable & aParticleTable,
                                 MaterialTable & aMaterialTable) const
  -> bool {         
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
  for (int ind=0; ind < N; ++ind) {
    outputCSArray[ind] = fProcessImpl.GetCS(particles[ind].p.E()-aParticleTable.GetMass(particles[ind].pdg), particles[ind].pdg, matID,
                                            aParticleTable, aMaterialTable);
    //std::cout<<"^^^: cs="<<outputCSArray[ind]/mbarn
    //         <<" ir="<<particles[ind].ir<<std::endl;
  }
  return true;
}  

template <class ProcessImpl>
template <typename Floating>
auto Process<ProcessImpl>::GetFS(
    Particle<Floating> * particles, MatID_t matID, uint64_t N, PDG_t * outPDG1,
    T3LorentzVector<Floating> * outP1, Floating * csBorderDataFS,
    /*CSBorderStep is for using csBorderDataCS/FS in nbody.cpp=max number of elements in materials*/
    ParticleTable & aParticleTable, MaterialTable & aMaterialTable, int CSBorderStep) const
    -> bool {
  
//
#ifdef OPENACC
#pragma acc parallel loop gang vector present(particles,outPDG1,outP1,csBorderDataFS,aParticleTable,aMaterialTable) copy(CSBorderStep)
#else 
#pragma omp parallel for simd
#endif  
//
  for (auto ind = 0u; ind < N; ++ind) {
    auto const fProcessImplReturnTemporary=fProcessImpl.GetFS(particles[ind].p, particles[ind].pdg, matID, particles[ind].rs, csBorderDataFS,
                                                              particles[ind].tr, ind, aParticleTable, aMaterialTable, CSBorderStep);
    outPDG1[ind]=fProcessImplReturnTemporary.pdg1;
    outP1[ind]=fProcessImplReturnTemporary.P1;
  }
  return true;
}
  
template <class ProcessImpl>
template <typename Floating>
auto Process<ProcessImpl>::GetFS(
    Particle<Floating> * particles, MatID_t matID, uint64_t N, PDG_t * outPDG1,
    T3LorentzVector<Floating> * outP1, PDG_t * outPDG2,
    T3LorentzVector<Floating> * outP2, Floating * csBorderDataFS,
    ParticleTable & aParticleTable,
    /*CSBorderStep is for using csBorderDataCS/FS in nbody.cpp=max number of elements in materials*/
    MaterialTable & aMaterialTable, int CSBorderStep) const
    -> bool {
  
//
#ifdef OPENACC
#pragma acc parallel loop gang vector present(particles,outPDG1,outP1,outPDG2,outP2,csBorderDataFS,aParticleTable,aMaterialTable) copy(CSBorderStep)
#else 
#pragma omp parallel for simd
#endif  
// 
  for (auto ind = 0u; ind < N; ++ind) {
    auto const fProcessImplReturnTemporary=fProcessImpl.GetFS(particles[ind].p, particles[ind].pdg,
                                                              matID, particles[ind].rs, csBorderDataFS,
                                                              particles[ind].tr, ind, aParticleTable, aMaterialTable, CSBorderStep);
    outPDG1[ind]=fProcessImplReturnTemporary.pdg1;
    outPDG2[ind]=fProcessImplReturnTemporary.pdg2;
    outP1[ind]=fProcessImplReturnTemporary.P1;
    outP2[ind]=fProcessImplReturnTemporary.P2;
  }
  return true;
}
} // namespace t3

#endif // T3PROCESS_H
