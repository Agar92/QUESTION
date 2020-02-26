#pragma once
#ifndef T3ELASTICEMIONIONIMPL_H
#define T3ELASTICEMIONIONIMPL_H

#include <cmath>
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

#include "T3Defs.h"
#include "T3ThreeVector.h"
#include "T3LorentzVector.h"
#include "T3NSGangular_RW.hh"
#include "T3R_DDCS.hh"

namespace t3 {

using namespace units;  
template <typename Floating>
class T3ElasticEMIonIon
{
public:
  T3ElasticEMIonIon()
  {
    ParticleTable aParticleTable;
    deuteronPDG=aParticleTable.makePDGfromZandA(1,2);
    titanPDG=aParticleTable.makePDGfromZandA(22,48);
  }

  inline Floating
  GetCS(Floating Tls, PDG_t incPDG, PDG_t targetPDG,
        ParticleTable & aParticleTable) const;

  inline Floating
  GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
        ParticleTable & aParticleTable, MaterialTable & aMaterialTable) const;


  inline Floating GetRNDCosThetaCMInEqualRutherford(
                  PDG_t incPDG,
                  Floating Tls/*LS kinetic energy of the inc particle*/,
                  Floating mr/*reduced mass*/, int z/*charge of inc particle*/,
                  Floating pcm2/*pcm^2*/,
                  Floating tmin/*|t|min=2*M*Ede*/,
                  Floating tmax/*|t|max=4*pcm^2*/,
                  Floating tup/*=|t|max-|t|min or kb*|t|m*/,
                  Floating R1/*random number 1*/,
                  Floating R2/*random number 2*/,
                  ParticleTable & aParticleTable) const;


  inline auto GetFS(T3LorentzVector<Floating> const & p, PDG_t incPDG, MatID_t matID, unsigned int & generator_seed,
                    Floating * csBorderDataFS, Floating & tr/*temporary register of Particle struct*/,
                    int ind, ParticleTable & aParticleTable, MaterialTable & aMaterialTable,
                    int csBorderSTEP/*maximim number of isotopes in any material*/) const;
  
private:
  PDG_t deuteronPDG;
  PDG_t titanPDG;
  const Floating alpha=1.0/137;
  const Floating hc=200.0 * MeV * fm;
  const Floating kb=1.0e-3;
  const Floating Edisplace_deuteron = 10 * eV;
  const Floating Edisplace_titan    = 25 * eV;
  //-----------------------------------------------//

};

template <typename Floating>
Floating T3ElasticEMIonIon<Floating>::GetCS(Floating Tls, PDG_t incPDG, PDG_t targetPDG,
                                            ParticleTable & aParticleTable) const
{
  if(!aParticleTable.IsNucleus(incPDG))
  {
    return 0.0;
  }

  const Floating m=aParticleTable.GetMass(incPDG);   // m
  const Floating M=aParticleTable.GetMass(targetPDG);// M
  const int z=aParticleTable.GetZ(incPDG);           // z
  const int Z=aParticleTable.GetZ(targetPDG);        // Z
  const Floating Elsinc = m+Tls;//LS full energy of the inc particle.
  const Floating s      = m*m+2*Elsinc*M+M*M;//s mandelstahm variable.
  const Floating mr     = m*M/sqrt(s);//invariant mass.
  const Floating pls2   = Tls*(2*m+Tls);
  const Floating rat    = m/M;
  const Floating rat2   = rat*rat;
  const Floating pcm2   = pls2/(1.0+rat2+2*sqrt(rat2+pls2/M/M));
  const Floating tm     = 2*pcm2;
  const Floating tmax   = 4*pcm2;
  const Floating COEF1  = M_PI/pcm2;
  const Floating coef2  = 2*mr*alpha*z*Z;
  const Floating COEF2  = coef2*coef2;
  const Floating COEF   = COEF1*COEF2;
  Floating tmin=0.0;
  if(targetPDG == deuteronPDG)   tmin=2*Edisplace_deuteron*M;
  else if(targetPDG == titanPDG) tmin=2*Edisplace_titan*M;
  if(tmax<=tmin)
  { 
    return 0.0;
  }
  
  int channel=-1;
  Floating cs=0.0;
  if(incPDG == targetPDG)
  {
    const Floating beta_r_cm = 1.0/std::sqrt(1.0+mr*mr/pcm2);
    if(incPDG == deuteronPDG)
    {
      if(Tls < 10 * keV)
      {
        channel=1;
        const Floating term1=1.0/tmin - 1.0/(tmax-tmin);
        const Floating term21 = 2.0/3.0/tmax*beta_r_cm/alpha;
        const Floating term22 = alpha/beta_r_cm * std::log( tmax/tmin-1.0 );
        const Floating term2  = term21*std::sin(term22);
        cs = 2*COEF*(term1+term2)*hc*hc;
      }
      else
      {
        channel=2;

        if(1.0e-3*tm<=tmin)
        {          
          return 0.0;
        }
                
        const Floating term1=1.0/tmin - 1.0/kb/tm;
        const Floating term2=1.0/(tmax-kb*tm)-1.0/(tmax-tmin);
        const Floating term3c = 2.0/3.0/tmax*beta_r_cm/alpha;
        const Floating term31 = alpha/beta_r_cm * std::log( kb*tm/(tmax-kb*tm) );
        const Floating term32 = alpha/beta_r_cm * std::log( tmin/(tmax-tmin) );
        const Floating term3  = term3c*(std::sin(term31)-std::sin(term32));
        cs = 2*COEF*(term1+term2+term3)*hc*hc;
      }
    }
    else
    {
      channel=3; 
      const Floating term1 = 1.0/tmin - 1.0/(tmax-tmin);
      const Floating spin=aParticleTable.GetSpin(incPDG);
      const bool cond=(fmod(spin,1.0))<1.0e-7;
      const int SIGN=cond?1:-1;
      const Floating coef21 = SIGN * 2.0/(2*spin+1)/tmax*beta_r_cm/alpha;
      const Floating term21 = alpha/beta_r_cm * std::log( tmax/tmin-1.0 );
      const Floating term2  = coef21*std::sin(term21);
      cs = 2*COEF*(term1+term2)*hc*hc;
    }
  }
  else
  {
    channel=4; 
    const Floating term = 1.0/tmin - 1.0/tmax;
    cs = COEF*term*hc*hc;
  }
  return cs;
}

template <typename Floating>
Floating T3ElasticEMIonIon<Floating>::GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
                                            ParticleTable & aParticleTable,
                                            MaterialTable & aMaterialTable) const
{
  Floating cs=0.0;
  for(int i=0; i<aMaterialTable.GetNumberOfIsotopes(matID); ++i)
  {
    PDG_t isotope = aMaterialTable.GetIsotopes(matID, i);
    cs += aMaterialTable.GetFractions(matID, i) *
      GetCS(Tls, incPDG, isotope, aParticleTable);
  }  
  return cs;
}



template <typename Floating>
Floating T3ElasticEMIonIon<Floating>::GetRNDCosThetaCMInEqualRutherford(
                           PDG_t incPDG,
                           Floating Tls/*LS kinetic energy of the incident particle*/,
                           Floating mr/*reduced mass*/,
                           int z/*the particles are equal=>z=Z*/,
                           Floating pcm2/*pcm^2*/,
                           /*lower bound:*/
                           Floating tmin/*|t|min*/, Floating tmax/*|t|max*/,
                           /*upper bound:*/
                           Floating tM/*this is tup, =|t|m=2*pcm^2 or kb*|t|m for D-D*/,
                           Floating R1/*random number 1*/,
                           Floating R2/*random number 2*/,
                           ParticleTable & aParticleTable) const
{
  const Floating coef=2*mr*alpha*z*z;
  const Floating COEF = M_PI/pcm2 * coef * coef;
  Floating CS1=(1.0/tmin-1.0/tM);
  Floating CS2=(1.0/(tmax-tM)-1.0/(tmax-tmin));
  const Floating spin=aParticleTable.GetSpin(incPDG);
  const bool cond=(fmod(spin,1.0))<1.0e-7;
  const int SIGN=(cond)?1:-1;
  
  if(tM<=tmin)
  { 
    printf("*** ERROR: T3ElasticEMIonIon::GetRNDCosThetaCMInEqualRutherford(): tM=%f <= tmin=%f", tM, tmin);
    return 1.0;//cos(theta_cm)=1 => theta_cm=0.
  }
  if(incPDG == deuteronPDG)
  {
    if(Tls >= 10 * keV && fabs(tM-kb*tmax/2)>1.0e-7)
      printf("***ERROR: Tls=%f >= 10 keV and tM=%f != 1.0e-3*tm=%f", Tls, tM, kb*tmax/2);
    if(Tls < 10 * keV && fabs(tM-tmax/2)>1.0e-7)
      printf("***ERROR: Tls=%f < 10 keV and tM=%f != tm=%f", Tls, tM, tmax/2);
  }
  
  if(CS1<0.0 || CS2<0.0) printf("***ERROR: CS1=%f CS2=%f\n",CS1,CS2);
  int channel=-1;
  Floating CSSum=CS1+CS2;
  Floating R1ar=R1*CSSum;
  if(R1ar<=CS1)             channel=1;
  else if(R1ar<=(CS1+CS2))  channel=2;
  Floating t=-1.0;
  if(channel==1)
  {
     const Floating tR=(1.0-R2)/tmin+R2/tM;
      t=1.0/tR;
  }
  else if(channel==2)
  {
      const Floating tR=(1.0-R2)/(tmax-tmin)+R2/(tmax-tM);
      t=tmax-1.0/tR;
  }
  const Floating tm=2*pcm2;
  const Floating CosThetaCM=1.0-t/tm;
  return CosThetaCM;
}


template <typename Floating>
auto T3ElasticEMIonIon<Floating>::GetFS(T3LorentzVector<Floating> const & p, PDG_t incPDG, MatID_t matID,
                                        unsigned int & generator_seed/*random number seed*/,
                                        Floating * csBorderDataFS, Floating & tr/*temporary register of Particle struct*/,
                                        int ind, ParticleTable & aParticleTable, MaterialTable & aMaterialTable,
                                        int csBorderSTEP/*maximum number of isotopes in any material*/) const
{  
  if(!aParticleTable.IsNucleus(incPDG))
  {
    printf("***ERROR: T3ElasticEMIonIon::GetFS(): incPDG=%d - not a nucleus!", incPDG);
    return Four<Floating>(PDG_t(0), PDG_t(0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0));
  }
  const Floating m = aParticleTable.GetMass(incPDG);//mass of the incident ion
  const Floating E = p.E();//LS full energy of the incident ion
  const Floating Tls = E - m;//LS kinetic energy of the incident ion
  const int numi = aMaterialTable.GetNumberOfIsotopes(matID);
  const size_t init=ind * csBorderSTEP;
  const size_t fin =init+numi;
  csBorderDataFS[init] = 0.0;
  Floating sigma_average = 0.0;
  for(size_t i=0; i<numi; ++i)
  {
    sigma_average += aMaterialTable.GetFractions(matID,i) *
      GetCS(Tls, incPDG, aMaterialTable.GetIsotopes(matID,i), aParticleTable);
    csBorderDataFS[init+i+1] = sigma_average;
  }
  const Floating aR = RND01(generator_seed);
  const Floating raR = aR * sigma_average;
  PDG_t isotopePDG=-1;
  for(size_t i=0; i<numi; ++i)
  {
    if(csBorderDataFS[init+i]<=raR) isotopePDG = aMaterialTable.GetIsotopes(matID,i);
  }
  const Floating M=aParticleTable.GetMass(isotopePDG);//mass of the target ion.
  const int z=aParticleTable.GetZ(incPDG);            //z
  const int Z=aParticleTable.GetZ(isotopePDG);        //Z
  Floating tmin=0.0;
  if(isotopePDG == deuteronPDG)   tmin=2 * Edisplace_deuteron * M;
  else if(isotopePDG == titanPDG) tmin=2 * Edisplace_titan * M;
  T3ThreeVector<Floating> InitMomentum = p.Vect();
  const Floating plsinc=InitMomentum.R();
  const Floating Elsinc=sqrt(m*m+plsinc*plsinc);
  const Floating s=m*m+2*Elsinc*M+M*M;
  const Floating mr=m*M/sqrt(s);//invariant mass.
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
  int channel=-1;
  if(incPDG == isotopePDG)
  {
    if(incPDG == deuteronPDG)
    {
      channel=1;
      if(Tls >= 10*keV)
      {
        CosThetaCM=GetRNDCosThetaCMInEqualRutherford(incPDG,Tls,mr,z,pcm2,tmin,tmax,
                                                     kb*tm,R1,R2,aParticleTable);
      }
      else
      {
        CosThetaCM=GetRNDCosThetaCMInEqualRutherford(incPDG,Tls,mr,z,pcm2,tmin,tmax,
                                                     tm,R1,R2,aParticleTable);
      }
    }
    else
    {
      channel=2;
      CosThetaCM=GetRNDCosThetaCMInEqualRutherford(incPDG,Tls,mr,z,pcm2,tmin,tmax,
                                                   tm,R1,R2,aParticleTable);
    }
  }
  else
  { 
    Floating t=tmin/(1.0-R1*(1.0-tmin/tmax));
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

  if(outPDG1<0 || outPDG2<0) printf("***ERROR: T3ElasticEMIonIon::GetFS(): Wrong outPDG!!!");
  if(isnan(Pcmfin.x()) || isnan(Pcmfin.y()) || isnan(Pcmfin.z()))
  {
    printf("***ERROR2!!!\n");
    printf("\nx=%f y=%f z=%f\n", Pcmfin.x(), Pcmfin.y(), Pcmfin.z());
    printf("Pinc: %f %f %f\n", InitMomentum.x(), InitMomentum.y(), InitMomentum.z());
    printf("Pcminc: %f %f %f\n", Pcm.x(), Pcm.y(), Pcm.z());
    printf("e1: %f %f %f\n", e1.x(), e1.y(), e1.z());
    printf("e2: %f %f %f\n", e2.x(), e2.y(), e2.z());
    printf("e3: %f %f %f\n", e3.x(), e3.y(), e3.z());
    printf("CosThetaCM=%f\n", CosThetaCM);
    printf("%d %d\n", outPDG1, outPDG2);
    printf("%i\n", channel);
  }
  
  return Four<Floating>(outPDG1, outPDG2, outP1, outP2);
}
  

}
#endif//T3ELASTICEMIONIONIMPL_H
