#pragma once
#ifndef T3USERSTEPPINGACTIONININJECT_H
#define T3USERSTEPPINGACTIONININJECT_H

#include <iostream>
#include "T3Globals.h"
#include "T3ParticleTable.h"

using namespace t3;

//NUMBER OF BINS IN THE HISTOGRAM:
constexpr int HNbinInject=1024;

class UserSteppingActionInInject
{
public:
  UserSteppingActionInInject():
    deltacoscmInject(2.0/HNbinInject),
    deuteronPDG(t3::ParticleTable().makePDGfromZandA(1,2)),
    protonPDG(t3::PDG_t(2212)),neutronPDG(t3::PDG_t(2112))
  {
    _IsRegistered=false;
    //t3::ParticleTable aParticleTable=t3::ParticleTable();
    //deuteronPDG=aParticleTable.makePDGfromZandA(1,2);
    //std::cout<<"C: "<<std::endl;
    //std::cout<<"deuteronPDG="<<deuteronPDG<<std::endl;
    for(int m=0; m<HNbinInject; ++m) HistogramInelasticDD[m]=0;
  }
  void Register(){_IsRegistered=true;}
  bool IsRegistered(){return _IsRegistered;}
  void UserAction(Particle<FloatingType> & particlei)
  {
    if(particlei.pdg==neutronPDG || particlei.pdg==protonPDG)
      FillHistogramThetaInelasticDD(particlei);

    //the first user action is to kill all the deuterons after their first
    //inelastic d-d reaction:
    //std::cout<<"pdg="<<particlei.pdg<<std::endl;
    if(particlei.pdg==neutronPDG || particlei.pdg==protonPDG)
      particlei.ir=0;
  }
  void FillHistogramThetaInelasticDD(Particle<FloatingType> & particlei)
  {
    const double costcm=particlei.tr;
    const int bin=(costcm+1.0)/deltacoscmInject;
#ifdef OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif
    ++HistogramInelasticDD[bin];
  }
  void Histogram_theta_UserSteppingActionInInject()
  {
#ifdef OPENACC
#pragma acc update host(HistogramInelasticDD[0:HNbinInject])
#endif
    std::ofstream foutne_theta;
    foutne_theta.open("ne_thetaineldd.dat");
    for(int m=0; m<HNbinInject; ++m)
    {
      foutne_theta<<m<<"   ";
      const double x=-1.0+deltacoscmInject*(m+0.5);
      foutne_theta<<std::setw(8)<<x<<"   "
                  <<HistogramInelasticDD[m]<<"   "
                  <<static_cast<double>(HistogramInelasticDD[m])/Np/deltacoscmInject<<"   ";
      foutne_theta<<std::endl;
    }
    foutne_theta.close();
  }

private:
  const double deltacoscmInject;
  int HistogramInelasticDD[HNbinInject];
  bool _IsRegistered;
  const t3::PDG_t deuteronPDG;
  const t3::PDG_t protonPDG;
  const t3::PDG_t neutronPDG;
};


#endif//T3USERSTEPPINGACTIONININJECT_H
