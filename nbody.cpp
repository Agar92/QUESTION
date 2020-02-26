#include <chrono>
#include <iostream>
#include "unistd.h"
#include <typeinfo>

#include "T3DataHolder.h"


#include "T3NSGangular_node.hh"
#include "T3Inelasticdd_DB.hh"
#include "T3InelasticddCSImpl.h"
#include "T3InelasticddFSImpl.h"
#include "T3Utility.hh"


#include "T3ElasticEMIonIonCSImpl.h"
#include "T3ElasticEMIonIonFSImpl.h"

#include <limits>

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
    

    ParticleTable cParticleTable;
    MaterialTable cMaterialTable;
    const PDG_t DeuteronPDG=cParticleTable.makePDGfromZandA(1,2);
    const PDG_t titanPDG=cParticleTable.makePDGfromZandA(22,48);
    const MatID_t material=MatID_t(2);
    
    T3Inelasticdd<double> inelasticprocess;
    auto begin=std::chrono::steady_clock::now();
    DataHolder<FloatingType> d;
    
    if(report)
      std::cout << "size of DataHolder is ~"
                << sizeof(DataHolder<FloatingType>) / float(1ul << 30ul)
                << "GB, size of d is ~" << sizeof(d) / float(1ul << 30ul) << "GB"
                << std::endl;
    int count=0;
    LIFE=0;
    MAX_ELEMENT=0;
    
#ifdef OPENACC
#pragma acc data create(ind01,ind23,arr1,arr2,arr3,outPDG1Inelasticdd,outPDG2Inelasticdd,outP1Inelasticdd,outP2Inelasticdd,outPDG1ElasticEMIonIon,outP1ElasticEMIonIon,outPDG1ElasticStrongIonIon,outP1ElasticStrongIonIon,csBorderDataFS,csMultipleScattering,csInelasticdd,csElasticEMIonIon,csElasticStrongIonIon) \
  copyin(particles,d,inelasticddProcess,ElasticEMIonIonProcess,ElasticStrongIonIonProcess/*,ElasticEMIonIonImpl,ElasticStrongIonIonImpl*/)
  {
#endif
    
    bool InitLoop=true;
    for(unsigned int step=1; GetNumOfAliveParticles()>0 || InitLoop==true || GetNumOfInjectedParticles()<Np; ++step)
    {
      d.Inject();
      d.Propagate();
      d.Compress();
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

    std::cout<<"NEUTRON_COUNT="<<NEUTRON_COUNT<<" He3_COUNT="<<He3_COUNT
             <<" PROTON_COUNT="<<PROTON_COUNT<<" TRITON_COUNT="<<TRITON_COUNT
             <<" Ti48_COUNT="<<Ti48_COUNT<<" KILLED_DEUTERON_COUNT="
             <<KILLED_DEUTERON_COUNT
             <<" COUNT_IR3_IN_PROPAGATE="<<COUNT_IR3_IN_PROPAGATE
             <<" COUNT_IR3_IN_REACT="<<COUNT_IR3_IN_REACT
             <<" IR1="<<IR1<<" IR2="<<IR2
             <<" IR1/IR2="<<static_cast<double>(IR1)/IR2
             <<std::endl;

    std::cout<<"deltaTls="<<deltaTls/keV<<std::endl;
    for(int m=0; m<NINT; ++m) ENINT[m]=TlsMin+(m+0.5)*deltaTls;
    for(int m=0; m<NINT; ++m) LF[m]/=Np;
    std::ofstream outf("fluence.dat");
    for(int m=0; m<NINT; ++m)
    {
      outf<<m<<"   ";
      outf<<ENINT[m]/keV<<"   "<<NEVENTS[m]<<"   "<<LF[m]/um<<std::endl;
    }
    outf.close();
    const MatID_t matTiD2=MatID_t(2);
    const double ntid2=cMaterialTable.GetConcentrations(matTiD2);
    double Integral=0.0;
    for(int m=0; m<NINT; ++m)
    {
      double Tlsi;
      if(m!=NINT-1) Tlsi=(ENINT[m]+ENINT[m+1])/2;
      else          Tlsi=(ENINT[m]+TlsMax)/2;
      
      const double InelasticDDCSi=inelasticprocess.GetCS(Tlsi, DeuteronPDG, material,
                                                        cParticleTable, cMaterialTable);
      const double dInt=LF[m]*ntid2*InelasticDDCSi;
      Integral+=dInt;
    }
    std::cout<<"Integral="<<Integral<<" ntid2="<<ntid2*cm*cm*cm<<std::endl;
    std::ofstream outn("probability_of_neutron_output.dat");
    for(int m=0; m<NINT; ++m)
    {
      double Tlsi;
      if(m!=NINT-1) Tlsi=(ENINT[m]+ENINT[m+1])/2;
      else          Tlsi=(ENINT[m]+TlsMax)/2;
      const double InelasticDDCSi=inelasticprocess.GetCS(Tlsi, DeuteronPDG, material,
                                                         cParticleTable, cMaterialTable);
    }
    outn.close();
#ifdef OPENACC
  }
#endif

  auto end=std::chrono::steady_clock::now();
	auto elapsed_ms=std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
	std::cout<<"time="<<elapsed_ms.count()<<" ms, G="<<G<<", K="<<K<<", Ntop="
           <<Ntop<<", SumDG="<<SumDGam<<std::endl;
  std::cout<<"Nbin="<<Nbin<<" FloatingType="<<typeid(FloatingType).name()<<std::endl;

  std::cout<<"NNN="<<NNN<<" cuba="<<cuba<<" INJ="<<INJ<<" ag="<<ag/um<<" um"
           <<" TARGET_WIDTH="<<TARGET_WIDTH/um<<" um"<<std::endl;
  
}
