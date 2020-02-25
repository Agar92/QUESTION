#pragma once
#ifndef T3ELASTICEMIONIONFSIMPL_H
#define T3ELASTICEMIONIONFSIMPL_H

#include <cmath>
#include <random>
#include <typeinfo>
#include "T3Defs.h"
#include "T3ThreeVector.h"
#include "T3LorentzVector.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

#include "T3NSGangular_RW.hh"
#include "T3R_DDCS.hh"
#include "T3ElasticEMIonIonCSImpl.h"

namespace t3 {

template <typename Floating = double>
class T3ElasticEMIonIonFS
{
public:
  T3ElasticEMIonIonFS()
  {
    ParticleTable aParticleTable;
    deuteronPDG = aParticleTable.makePDGfromZandA(1, 2);
    titanPDG = aParticleTable.makePDGfromZandA(22, 48);
  }

  inline Floating GetRNDCosThetaCMInEqualRutherford(
                  PDG_t incPDG,
                  Floating mr, int z,
                  Floating pcm2,
                  Floating tmin,
                  Floating tmax,
                  Floating tup,
                  Floating R1,
                  Floating R2,
                  ParticleTable & aParticleTable) const;
  
  inline auto GetFS(T3LorentzVector<Floating> const & p, PDG_t incPDG, MatID_t matID, unsigned int & generator_seed,
                    Floating * csBorderDataFS, Floating & tr,
                    int ind, ParticleTable & aParticleTable, MaterialTable & aMaterialTable,
                    int csBorderSTEP) const;
private:
  //******************************************************************************//
  PDG_t deuteronPDG;
  PDG_t titanPDG;
  const Floating alpha=1.0/137;
  const Floating hc=200.0 * MeV * fm;
  const Floating kb=1.0e-3;
  const Floating Edisplace_deuteron = 10 * eV;
  const Floating Edisplace_titan    = 25 * eV;
  const T3ElasticEMIonIonCS<Floating> aCS;
};

template <typename Floating>
Floating T3ElasticEMIonIonFS<Floating>::GetRNDCosThetaCMInEqualRutherford(
                              PDG_t incPDG,
                              Floating mr, int z,
                              Floating pcm2,
                              Floating tmin, Floating tmax,
                              Floating tup,
                              Floating R1,
                              Floating R2,
                              ParticleTable & aParticleTable) const                                                             
{
  const Floating coef=2*mr*alpha*z*z;
  const Floating COEF = M_PI/pcm2 * coef * coef;
  const Floating CS1=(1.0/tmin-1.0/tup) * COEF*hc*hc;
  const Floating CS2=(1.0/(tmax-tup)-1.0/(tmax-tmin)) * COEF*hc*hc;
  const Floating spin=aParticleTable.GetSpin(incPDG);
  const bool cond=(std::fmod(spin,1.0))<1.0e-7;
  const int SIGN=(cond)?1:-1;
  const Floating beta_r_cm = 1.0/std::sqrt(1.0+mr*mr/pcm2);
  const Floating cs3factor1=std::sin( alpha/beta_r_cm * std::log(tup/(tmax-tup)) );
  const Floating cs3factor2=std::sin( alpha/beta_r_cm * std::log(tmin/(tmax-tmin)) );
  const Floating CS3=SIGN*2.0/(2*spin+1)/tmax*beta_r_cm/alpha*
                                        (cs3factor1-cs3factor2) * COEF*hc*hc;
  int channel=-1;
  const Floating CSSum=CS1+CS2+CS3;
  const Floating R1ar=R1*CSSum;
  Floating t=-1.0;
  for(int i=0; i<3; ++i)
  {
    if(R1ar>=0.0) channel=1;
    if(R1ar>=CS1) channel=2;
    if(R1ar>=CS2) channel=3;
  }
  if(channel==1)
  {
    const Floating tR=(1.0-R2)/tmin+R2/tup;
    t=1.0/tR;
  }
  else if(channel==2)
  {
    const Floating tR=(1.0-R2)/(tmax-tmin)+R2/(tmax-tup);
    t=tmax-1.0/tR;
  }
  else if(channel==3)
  {
    const Floating Argument1=std::sin( alpha/beta_r_cm * std::log( tmin/(tmax-tmin) ) );
    const Floating Argument2=std::sin( alpha/beta_r_cm * std::log( tup/(tmax-tup) ) );
    const Floating A=beta_r_cm/alpha*std::asin((1.0-R2)*Argument1+Argument2);
    const Floating eA=std::exp(A);
    t=tmax*eA/(1.0+eA);
  }
  const Floating tm=2*pcm2;
  const Floating CosThetaCM=1.0-t/tm;
  return CosThetaCM;
}

template <typename Floating>
auto T3ElasticEMIonIonFS<Floating>::GetFS(T3LorentzVector<Floating> const & p, PDG_t incPDG, MatID_t matID,
                                          unsigned int & generator_seed,
                                          Floating * csBorderDataFS, Floating & tr,
                                          int ind, ParticleTable & aParticleTable, MaterialTable & aMaterialTable,
                                          int csBorderSTEP) const
{
  const Floating m = aParticleTable.GetMass(incPDG);
  const Floating E = p.E();
  const Floating Tls = E - m;
  const int numi = aMaterialTable.GetNumberOfIsotopes(matID);
  const size_t init=ind * csBorderSTEP;
  const size_t fin =init+numi;
  csBorderDataFS[init] = 0.0;
  Floating sigma_average = 0.0;
  for(size_t i=0; i<numi; ++i)
  {
    sigma_average += aMaterialTable.GetFractions(matID,i)
      * aCS.GetCS(Tls, incPDG, aMaterialTable.GetIsotopes(matID,i), aParticleTable);
    csBorderDataFS[init+i+1] = sigma_average;
  }
  const Floating aR = RND01(generator_seed);
  const Floating raR = aR * sigma_average;
  PDG_t isotopePDG=-1;
  for(size_t i=0; i<numi; ++i)
  {
    if(csBorderDataFS[init+i]<=raR) isotopePDG = aMaterialTable.GetIsotopes(matID,i);
  }
  const Floating M=aParticleTable.GetMass(isotopePDG);
  const int z=aParticleTable.GetZ(incPDG);            //z
  const int Z=aParticleTable.GetZ(isotopePDG);        //Z
  Floating tmin=0.0;
  if(isotopePDG == deuteronPDG)   tmin=2 * Edisplace_deuteron * M;
  else if(isotopePDG == titanPDG) tmin=2 * Edisplace_titan * M;
  T3ThreeVector<Floating> InitMomentum = p.Vect();
  const Floating plsinc=InitMomentum.R();
  const Floating Elsinc=sqrt(m*m+plsinc*plsinc);
  const Floating s=m*m+2*Elsinc*M+M*M;
  const Floating mr=m*M/sqrt(s);
  const Floating pls2=Tls*(2*m+Tls);
  const Floating rat=m/M;
  const Floating rat2=rat*rat;
  const Floating pcm2=pls2/(1.0+rat2+2.0*sqrt(rat2+pls2/M/M));
  const Floating tm  =2*pcm2;
  const Floating tmax=4*pcm2;
  const Floating coef=M/sqrt(s);
  T3ThreeVector<Floating> Pcm = InitMomentum*coef;
  const Floating xy=Pcm.x()*Pcm.y(), xz=Pcm.x()*Pcm.z(), yz=Pcm.y()*Pcm.z();
  const Floating x2=Pcm.x()*Pcm.x(), y2=Pcm.y()*Pcm.y(), z2=Pcm.z()*Pcm.z();
  T3ThreeVector<Floating> e1, e2, e3;
  if(Pcm.x() < Pcm.y())
  {
    if(Pcm.x() < Pcm.z())
    {
      e2={0., Pcm.z(), -Pcm.y()};
      e3={y2+z2, -xy, -xz};
    }
    else
    {
      e2={Pcm.y(), -Pcm.x(), 0.};
      e3={-xz, -yz, y2+x2};
    }
  }
  else
  {
    if(Pcm.y() < Pcm.z())
    {
      e2={Pcm.z(), 0., -Pcm.x()};
      e3={xy, -x2-z2, yz};
    }
    else
    {
      e2={Pcm.y(), -Pcm.x(), 0.};
      e3={-xz, -yz, y2+x2};
    }
  }
  e1=Pcm;
  e1.Unit();
  e2.Unit();
  e3.Unit();
  const Floating R1=RND01(generator_seed);
  const Floating R2=RND01(generator_seed);
  Floating CosThetaCM=0.0;
  Floating t=-1.0;
  if(incPDG == isotopePDG)
  {
    if(incPDG == deuteronPDG && Tls >= 30*keV)
      CosThetaCM=GetRNDCosThetaCMInEqualRutherford(incPDG,mr,z,pcm2,tmin,tmax,
                                                   kb*tm,R1,R2,aParticleTable);
    else
      CosThetaCM=GetRNDCosThetaCMInEqualRutherford(incPDG,mr,z,pcm2,tmin,tmax,
                                                   tmax-tmin,R1,R2,aParticleTable);
  }
  else
  {
    t=tmin/(1.0-R1*(1.0-tmin/tmax));
    CosThetaCM=1.0-t/tm;
  }
  tr=CosThetaCM;      
  const Floating phi=2*M_PI*RND01(generator_seed);
  const Floating cosPhi=cos(phi), sinPhi=sin(phi);
  const Floating SinThetaCM=sqrt(1.0-CosThetaCM*CosThetaCM);
  const Floating pcminc=Pcm.R();
  PDG_t outPDG1=incPDG;
  PDG_t outPDG2=isotopePDG;
  const Floating ss=pcminc*SinThetaCM*sinPhi, sc=pcminc*SinThetaCM*cosPhi;
  T3ThreeVector<Floating> Pcmfin;
  Pcmfin.SetX(e1.x()*pcminc*CosThetaCM+e2.x()*ss+e3.x()*sc);
  Pcmfin.SetY(e1.y()*pcminc*CosThetaCM+e2.y()*ss+e3.y()*sc);
  Pcmfin.SetZ(e1.z()*pcminc*CosThetaCM+e2.z()*ss+e3.z()*sc);
  T3ThreeVector<Floating> Vcm=InitMomentum/(Elsinc+M);
  const Floating E1cmfin=sqrt(m*m+pcminc*pcminc);
  const Floating AbsVcm=Vcm.R();
  const Floating gamma=1.0/std::sqrt(1-AbsVcm*AbsVcm);
  T3ThreeVector<Floating> ParallelProjection=(e1*Pcmfin)*e1;
  T3ThreeVector<Floating> PerpendicularProjection=Pcmfin-ParallelProjection;
  T3ThreeVector<Floating> plsfin1=gamma*(ParallelProjection+Vcm*E1cmfin);
  plsfin1+=PerpendicularProjection;
  T3ThreeVector<Floating> plsfin2=InitMomentum-plsfin1;
  const Floating Absplsfin1=plsfin1.R();
  const Floating Elsfin1=std::sqrt(m*m+Absplsfin1*Absplsfin1);
  const Floating Absplsfin2=plsfin2.R();
  const Floating Elsfin2=std::sqrt(M*M+Absplsfin2*Absplsfin2);
  T3LorentzVector<Floating> outP1 = T3LorentzVector<Floating>(plsfin1.x(), plsfin1.y(), plsfin1.z(), Elsfin1);
  T3LorentzVector<Floating> outP2 = T3LorentzVector<Floating>(plsfin2.x(), plsfin2.y(), plsfin2.z(), Elsfin2);

  const Floating Efin=Elsfin1+Elsfin2;
  const T3LorentzVector<Floating> tot = outP1 + outP2;
  
  return Four<Floating>(outPDG1, outPDG2, outP1, outP2);
}
  
}//namespace t3.
#endif//T3ELASTICEMIONIONFSIMPL_H
