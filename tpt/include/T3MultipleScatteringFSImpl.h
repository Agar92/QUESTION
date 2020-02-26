#pragma once
#ifndef T3MULTIPLESCATTERINGFSIMPL_H
#define T3MULTIPLESCATTERINGFSIMPL_H

#include "T3Defs.h"
#include "T3LorentzVector.h"
#include "T3ParticleTable.h"
//#include "T3MaterialTable.h"
#include "T3MultipleScatteringCSImpl.h"

#include "T3Utility.hh"
#include "T3NSGangular_node.hh"

namespace t3 {
using namespace units;
//<typename Floating = double> this means that the default template type is double.
//If to write RutherfordCS<>(), the type of Floating will be double.
template <bool generateRecoil = true, typename Floating = double> class MultipleScatteringFS {
public:
  MultipleScatteringFS()
  {
    Edisplace_deuteron=aCS.Edisplace_deuteron;
    Edisplace_titan   =aCS.Edisplace_titan;
    ParticleTable aParticleTable=ParticleTable();
    deuteronPDG=aParticleTable.makePDGfromZandA(1, 2);
    titanPDG   =aParticleTable.makePDGfromZandA(22, 48);
    mdeuteron  =aParticleTable.GetMass(deuteronPDG);
    mtitan     =aParticleTable.GetMass(titanPDG);
    me         =aParticleTable.GetMass(PDG_t(11));
    CTF        =1.0/2.0*pow(3*M_PI/4,2.0/3.0);
    c          =1.0;
    alpha      =1.0/137;
    ec         =sqrt(alpha);
    h          =1.0;
    a0         =h*h/me/ec/ec;
    Zdeuteron  =aParticleTable.GetZ(deuteronPDG);
    Ztitan     =aParticleTable.GetZ(titanPDG);
  }
  inline auto GetFS(T3LorentzVector<Floating> const &p/*4-momentum of the inc particle*/, PDG_t incPDG,
                    MatID_t matID, unsigned int & generator_seed,
                    Floating * csBorderDataFS, Floating & tr/*temporary register of Particle struct*/,
                    int ind, const ParticleTable & aParticleTable,
                    /*csBorderSTEP is for using csBorderDataCS/FS in nbody.cpp=max number of elements in materials*/
                    const MaterialTable & aMaterialTable, int csBorderSTEP) const;
private:
  MultipleScatteringCS<Floating> aCS;
  Floating Edisplace_deuteron;
  Floating Edisplace_titan;
  PDG_t deuteronPDG; //deuteron PDG
  PDG_t titanPDG;    //titan PDG
  Floating mdeuteron;//deuteron mass
  Floating mtitan;   //titan mass
  int Zdeuteron;
  int Ztitan;
  Floating me;
  Floating CTF;
  Floating c;
  Floating alpha;
  Floating ec;
  Floating h;
  Floating a0;
};

template <bool generateRecoil, typename Floating>
auto MultipleScatteringFS<generateRecoil, Floating>::GetFS(
    T3LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID,
    unsigned int & generator_seed,
    Floating * csBorderDataFS, Floating & tr/*temporary register of Particle struct*/,
    int ind, const ParticleTable & aParticleTable,
    /*csBorderSTEP is for using csBorderDataCS/FS in nbody.cpp=max number of elements in materials*/
    const MaterialTable & aMaterialTable, int csBorderSTEP) const {
  //generateSubCanonical() returns random number from 0.0 to 1.0.

  auto const generateSubCanonical = [&generator_seed]() {
    return RND01(generator_seed);
  };
  
  //if the particle is a gamma or a neutron, the Rutherford cross section is 0.0.
  /*if (incPDG == PDG_t(22) || incPDG == PDG_t(2112))
  {
    if constexpr (!generateRecoil) return Pair<Floating>(PDG_t(0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0));
    else return Four<Floating>(PDG_t(0), PDG_t(0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0));
  }*/  

  ///\\\///const PDG_t deuteronPDG = aParticleTable.makePDGfromZandA(1, 2);
  ///\\\///const PDG_t titanPDG    = aParticleTable.makePDGfromZandA(22, 48);
  auto const elsfull=p.energy();
  
  ///\\\///auto const m=aParticleTable.GetMass(incPDG);
  Floating m;
  if(incPDG==deuteronPDG)   m=mdeuteron;
  else if(incPDG==titanPDG) m=mtitan;
  
  auto const Tls=elsfull-m;
  auto const numi = aMaterialTable.GetNumberOfIsotopes(matID);
  auto sigma_average = 0.0;
  const size_t init=ind * csBorderSTEP;
  const size_t fin =init+numi;
  csBorderDataFS[init] = 0.0;
  for(size_t i=0; i<numi; ++i)
  {
    sigma_average += aMaterialTable.GetFractions(matID,i)
      * aCS.GetCS(Tls, incPDG, aMaterialTable.GetIsotopes(matID,i), aParticleTable);
    csBorderDataFS[init+i+1] = sigma_average;
  }
  auto const aR = RND01(generator_seed);
  auto const raR = aR * sigma_average;
  PDG_t isotopePDG=PDG_t(-1);
  for(size_t i=0; i<numi; ++i)
  {
    //HERE WAS AN ERROR:
    //if(csBorderDataFS[init+numi-1-i]<raR && raR<=csBorderDataFS[init+numi-i])
    //THIS DID NOT WORK ON GPU WHEN WAS aR=0.000000,
    //SO I CORRECTED:
    //if(csBorderDataFS[init+numi-1-i]<=raR && raR<=csBorderDataFS[init+numi-i])
    //THE CORRECTED VARIANT:
    if(csBorderDataFS[init+numi-1-i]<=raR && raR<=csBorderDataFS[init+numi-i])
      isotopePDG = aMaterialTable.GetIsotopes(matID,numi-1-i);
  }
  auto tmin=0.0;
  ///\\\//Floating const M=aParticleTable.GetMass(isotopePDG);
  Floating M;
  if(isotopePDG==deuteronPDG)   M=mdeuteron;
  else if(isotopePDG==titanPDG) M=mtitan;
  

  ///\\\///if(isotopePDG == deuteronPDG)   tmin=2 * aCS.Edisplace_deuteron * M;
  ///\\\///else if(isotopePDG == titanPDG) tmin=2 * aCS.Edisplace_titan * M;

  if(isotopePDG == deuteronPDG)   tmin=2 * Edisplace_deuteron * M;
  else if(isotopePDG == titanPDG) tmin=2 * Edisplace_titan * M;

  ///\\\///auto const c=1.0;
  ///\\\///auto const alpha=1.0/137;
  ///\\\///auto const ec=sqrt(alpha);
  ///\\\///auto const me=aParticleTable.GetMass(PDG_t(11));
  ///\\\///auto const h=1.0;
  ///\\\///auto const a0=h*h/me/ec/ec;
  ///\\\///auto const hc=200.0 * MeV * fm;
  ///\\\///auto const CTF = 1.0/2.0 * pow(3*M_PI/4,2.0/3.0);
  
  ///\\\///auto const z=aParticleTable.GetZ(incPDG);
  int z;
  if(incPDG==deuteronPDG)   z=Zdeuteron;
  else if(incPDG==titanPDG) z=Ztitan;
  
  ///\\\///auto const Z=aParticleTable.GetZ(isotopePDG);
  int Z;
  if(isotopePDG==deuteronPDG)   Z=Zdeuteron;
  else if(isotopePDG==titanPDG) Z=Ztitan;


  
  auto const s=m*m+2*elsfull*M+M*M;
  auto const mr=m*M/sqrt(s);
  auto const pls2=Tls*(2*m+Tls);
  auto const rat=m/M;
  auto const rat2=rat*rat;
  auto const pcm2=pls2/(1.0+rat2+2.0*sqrt(rat2+pls2/M/M));
  auto const pcm=sqrt(pcm2);
  auto const aI=CTF*a0/sqrt(pow(z, 2.0/3.0) + pow(Z, 2.0/3.0));
//!!!!!!!!!!!!!!!!!!!!!!!!//  
  auto const beta_r=1.0;
//!!!!!!!!!!!!!!!!!!!!!!!!//    
  auto const c12=alpha*z*Z/beta_r;
  auto const c122=c12*c12;
  auto const const1=h/2/aI/pcm;
  auto const const2=const1*const1;
  auto const AS=const2*(1.13+3.76*c122);
  T3ThreeVector<Floating> InitMomentum = p.Vect();
  auto const plsinc=InitMomentum.R();
  auto const Elsinc=sqrt(m*m+plsinc*plsinc);
  Floating const coef=M/sqrt(s);
  T3ThreeVector<Floating> Pcm = InitMomentum*coef;
  auto costhetacm=0.0;
  auto const R=RND01(generator_seed);
  auto const tmax_real = 4 * pcm2;
  auto tmax = 0.0;
  if(tmax_real<=tmin) tmax = tmax_real;
  else                tmax = tmin;
  auto const rc1=4*pcm2*AS;
  auto const rc2=tmax/(tmax+rc1);
  auto const t=rc1*R*rc2/(1.0-R*rc2);
  auto sinhalfthetacm = sqrt(t/4/pcm2);
  auto thetacm = 2 * asin(sinhalfthetacm);
  costhetacm = cos(thetacm);
  auto const xy=Pcm.x()*Pcm.y(), xz=Pcm.x()*Pcm.z(), yz=Pcm.y()*Pcm.z();
  auto const x2=Pcm.x()*Pcm.x(), y2=Pcm.y()*Pcm.y(), z2=Pcm.z()*Pcm.z();
  T3ThreeVector<Floating> e1, e2;
  if(Pcm.x() < Pcm.y())
  {
    if(Pcm.x() < Pcm.z())
    {
        e1={0., Pcm.z(), -Pcm.y()};
        e2={y2+z2, -xy, -xz};
    }
    else
    {
        e1={Pcm.y(), -Pcm.x(), 0.};
        e2={-xz, -yz, y2+x2};
    }
  }
  else
  {
    if(Pcm.y() < Pcm.z())
    {
        e1={Pcm.z(), 0., -Pcm.x()};
        e2={xy, -x2-z2, yz};
    }
    else
    {
        e1={Pcm.y(), -Pcm.x(), 0.};
        e2={-xz, -yz, y2+x2};
    }
  }
  e1.Unit();
  e2.Unit();
  auto const phi=RND01(generator_seed) * TWOPI;
  auto const sinthetacm = sin(thetacm);
  auto const cosPhi = cos(phi);
  auto const sinPhi = sin(phi);
  auto const pcminc=Pcm.R();
  auto const ss=pcminc*sinthetacm*sinPhi, sc=pcminc*sinthetacm*cosPhi;
  T3ThreeVector<Floating> Pcmfin;
  Pcmfin.SetX(Pcm.x()*costhetacm+e1.x()*ss+e2.x()*sc);
  Pcmfin.SetY(Pcm.y()*costhetacm+e1.y()*ss+e2.y()*sc);
  Pcmfin.SetZ(Pcm.z()*costhetacm+e1.z()*ss+e2.z()*sc);
  Floating const Ecminc=sqrt(m*m+pcminc*pcminc);
  auto const Vcm=InitMomentum/(Elsinc+M);
  auto const AbsVcm=Vcm.R();
  Floating const gamma=1.0/sqrt(1-AbsVcm*AbsVcm);
  T3ThreeVector<Floating> Plsfin1=gamma*(Pcmfin+Vcm*Ecminc);
  T3ThreeVector<Floating> Plsfin2=InitMomentum-Plsfin1;
  auto const Absplsfin1=Plsfin1.R();
  auto const Elsfin1=sqrt(m*m+Absplsfin1*Absplsfin1);
  auto const Absplsfin2=Plsfin2.R();
  auto const Elsfin2=sqrt(M*M+Absplsfin2*Absplsfin2);
  auto const outPinc    = T3LorentzVector<Floating>(Plsfin1.x(), Plsfin1.y(), Plsfin1.z(), Elsfin1);
  auto const outPtarget = T3LorentzVector<Floating>(Plsfin2.x(), Plsfin2.y(), Plsfin2.z(), Elsfin2);
  auto const outPDGinc = incPDG;
  auto const outPDGtarget = isotopePDG;
  auto const Els1=p.E();
  auto const de=Elsfin1-Els1;



#ifdef DEBUG
  if(isnan(incPDG)) printf("GetFS(): i=%d incPDG=%d\n", ind, incPDG);
  if(isnan(isotopePDG)) printf("GetFS(): i=%d isotopePDG=%d\n", ind, isotopePDG);
  if(isnan(Elsfin1) || isnan(Plsfin1.x()) || isnan(Plsfin1.y()) || isnan(Plsfin1.z()))
  {
    printf("Elsfin1:\n");
    printf("GetFS(): i=%d px=%f py=%f pz=%f E=%f\n", ind, p.x(), p.y(), p.z());
    printf("GetFS(): i=%d p1x=%f p1y=%f p1z=%f E1=%f incPDG=%d isotopePDG=%d matID=%d\n",
           ind, Plsfin1.x(), Plsfin1.y(), Plsfin1.z(), Elsfin1, matID);
    printf("aR=%f phi=%f\n", aR, phi);
  }
  if(isnan(Elsfin2) || isnan(Plsfin2.x()) || isnan(Plsfin2.y()) || isnan(Plsfin2.z()))
  {
    printf("Elsfin2:\n");
    printf("GetFS(): i=%d px=%f py=%f pz=%f E=%f\n", ind, p.x(), p.y(), p.z());
    printf("GetFS(): i=%d p2x=%f p2y=%f p2z=%f E2=%f incPDG=%d isotopePDG=%d matID=%d\n",
           ind, Plsfin2.x(), Plsfin2.y(), Plsfin2.z(), Elsfin2, incPDG, isotopePDG, matID);
    printf("aR=%f phi=%f\n", aR, phi);
  }
#endif  
  
  

  /*
  std::cout<<"MultipleScatteringFS::GetFS():"<<std::endl;
  std::cout<<"i="<<ind<<std::endl;
  std::cout<<"incPDG="<<incPDG<<" isotopePDG="<<isotopePDG
           <<" E1="<<Elsfin1<<" E2="<<Elsfin2
           <<" P1="<<Absplsfin1<<" P2="<<Absplsfin2<<std::endl;
  */

  //\\//if constexpr (!generateRecoil)
  //\\//        return std::make_tuple(outPDGinc, outPinc);
  //\\//else
  //\\//        return std::make_tuple(outPDGinc, outPDGtarget, outPinc, outPtarget);

  //g++ is old on i7 and it does not know constexpr
  //with constexpr the compiler can deduce the return type
  //(marked as auto) at compile time.
  //But without constexpr it can not deduce the return type at compile time
  //and gives the error:
  // error: inconsistent deduction for ‘auto’: ‘t3::Pair<double>’ and then ‘t3::Four<double>’
  //       return Four<Floating>(outPDGinc, outPDGtarget, outPinc, outPtarget);
//!!!That is why i cooment the 1st if:!!!
  //if /*constexpr*/ (!generateRecoil)
  //      return Pair<Floating>(outPDGinc, outPinc);
  //else
  //      return Four<Floating>(outPDGinc, outPDGtarget, outPinc, outPtarget);

  return Four<Floating>(outPDGinc, outPDGtarget, outPinc, outPtarget);
 
  }
//end of new code.

}// namespace t3
#endif // T3MULTIPLESCATTERINGFSIMPL_H
