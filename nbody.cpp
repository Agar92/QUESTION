#include <chrono>
#include <iostream>
#include "unistd.h"

#include "T3DataHolder.h"
#include "T3InelasticddCSImpl.h"
#include "T3InelasticddFSImpl.h"

//////////////////////////////////////////////////////////////////////////////
// Program main
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
#ifdef OPENACC
  #ifdef CUDA
    std::cout<<"OPENACC IS DEFINED, CUDA IS DEFINED"<<std::endl;
  #else
    std::cout<<"OPENACC IS DEFINED, CUDA IS NOT DEFINED"<<std::endl;
  #endif
#else
  #ifdef CUDA
    std::cout<<"OPENACC IS NOT DEFINED, CUDA IS DEFINED"<<std::endl;
  #else
    std::cout<<"OPENACC IS NOT DEFINED, CUDA IS NOT DEFINED"<<std::endl;
  #endif
#endif
    
  auto begin=std::chrono::steady_clock::now();

  DataHolder<FloatingType> d;
  
  if (report)
    std::cout << "size of DataHolder is ~"
              << sizeof(DataHolder<FloatingType>) / float(1ul << 30ul)
              << "GB, size of d is ~" << sizeof(d) / float(1ul << 30ul) << "GB"
              << std::endl;
  
  int count=0;
  LIFE=0;
  MAX_ELEMENT=0;
#ifdef OPENACC
#pragma acc data create(ind01,ind23,arr1,arr2,arr3,outPDG1Inelasticdd,outPDG2Inelasticdd,outP1Inelasticdd,outP2Inelasticdd,outPDG1ElasticEMIonIon,outP1ElasticEMIonIon,outPDG1ElasticStrongIonIon,outP1ElasticStrongIonIon,csBorderDataFS,csMultipleScattering,csInelasticdd,csElasticEMIonIon,csElasticStrongIonIon) \
  copyin(particles,d,inelasticddProcess,ElasticEMIonIonProcess,ElasticStrongIonIonProcess)
  {
#endif
    
    bool InitLoop=true;
    for(unsigned int step=1; GetNumOfAliveParticles()>0 || InitLoop==true || GetNumOfInjectedParticles()<Np; ++step)
    { 
      d.Inject();     
      d.Propagate();
      //d.Compress();
      d.React();
      if(report)
      {
        std::cout << step << "   " << GetNumOfAliveParticles() << "    "
                  << GetNoNew() <<  "   " << GetNumOfInjectedParticles() << " "
                  << GetSumDGam() <<std::endl;
      }
      if(InitLoop) InitLoop=false;
      ++count;
    }
#ifdef OPENACC
  }
#endif

  auto end=std::chrono::steady_clock::now();
	auto elapsed_ms=std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
	std::cout<<"time="<<elapsed_ms.count()<<" ms, G="<<G<<", K="<<K<<", Ntop="
           <<Ntop<<", SumDG="<<SumDGam<<std::endl;
  std::cout<<"Nbin="<<Nbin<<" FloatingType="<<typeid(FloatingType).name()<<std::endl;
}
