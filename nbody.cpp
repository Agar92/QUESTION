#include <chrono>
#include <iostream>
#include "unistd.h"

#include "T3Globals.h"
#include "T3MSContiniousProcessImpl.h"

using namespace t3;
using namespace units;

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
    const PDG_t protonPDG=PDG_t(2212);
    const PDG_t DeuteronPDG=cParticleTable.makePDGfromZandA(1,2);
    const PDG_t titanPDG=cParticleTable.makePDGfromZandA(22,48);
    const MatID_t material=MatID_t(2);
    const double md=cParticleTable.GetMass(DeuteronPDG);
    const double mtitan=cParticleTable.GetMass(titanPDG);
    const double TTls= 10.0 * MeV;
    const double E=md+TTls;
    const double pls=std::sqrt(TTls*(2*md+TTls));
    t3::T3LorentzVector<double> pinc(0.0, 0.0, pls, E);

    T3MSContiniousProcess<double> mscontprocess;
    unsigned int rnd_seed=1;
    //*
    for(int i=0; i<100; ++i)
    {
      mscontprocess.MakeMSRandomAngleScatteringOnParticle4Momentum(DeuteronPDG, 1.0*units::um, pinc,
                                                                   rnd_seed, cParticleTable, material);
    }
    //*/
  
}
