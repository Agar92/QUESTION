#pragma once
#ifndef T3ELASTICSTRONGIONIONFSIMPL_H
#define T3ELASTICSTRONGIONIONFSIMPL_H

#include <random>
#include "T3Defs.h"
#include "T3ThreeVector.h"
#include "T3LorentzVector.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"
#include "T3InelasticddCSImpl.h"

#include "T3NSGangular_RW.hh"
#include "T3R_DDCS.hh"
#include "T3ElasticStrongIonIon_DB.hh"
#include "T3ElasticStrongIonIonCSImpl.h"

#include <typeinfo>

namespace t3 {

template <typename Floating = double>
class T3ElasticStrongIonIonFS
{
public:
  T3ElasticStrongIonIonFS()
  {
    ParticleTable aParticleTable;
    deuteronPDG = aParticleTable.makePDGfromZandA(1, 2);
    titanPDG    = aParticleTable.makePDGfromZandA(22, 48);
    md=aParticleTable.GetMass(deuteronPDG);
    //---------------------------------------------------------------------//
    const int tgZ = 1;//target deuteron Z=1
    const int tgA = 2;//target deuteron A=2
    //reaction product D+D. this is PDG of 1 of the deutrons in the reaction product.
    const t3::PDG_t sPDG = deuteronPDG;
    const int incZA = 1002;//1000*Z+A//inc deuteron
    //MT index of reaction in ENDF
    //For elastic scattering MT index of reaction in ENDF = 2. See ENDF6.
    const auto rid = "2";
// load the data if elastic D-D scattering approximation
// (Tls, normalized partial sums) from the binary file from
// ~/T3_DATA/angular/D/incD/2/T3DSGangular_100010020.bin .
    
    T3R_DDCS rdcs;
    rdcs.Load_from_PS(tgZ, tgA, rid, sPDG, incZA);
    //get vector of records:
    T3NSGangular_RW rw=rdcs.GetT3NSGangular_RW();
    //have only 1 record:
    T3NSGangular_RWrecord rwrec=rw.at(0);
    const Floating Emin = rwrec.front().Get_E();
    const Floating Emax = rwrec.back().Get_E();
    const Floating Ediff = Emax - Emin;
    const size_t nbins = 512 - 1;
    const Floating dE = Ediff / nbins;
    db=T3ElasticStrongIonIon_DB(-1.0, 1.0);
    for(size_t i = 0; i < 512; ++i)
    {
      const Floating Ei = Emin + i*dE;
      const T3NSGangular_RWnode rwnode = rwrec.interpolate(Ei);
      T3NSGangular_node nodes_i = T3NSGangular_node(rwnode);
      for(int j=0; j<127; ++j)
      {
        db.Set_V(i, j, nodes_i.Get_V(j));
        db.Set_a(i, j, nodes_i.Get_a(j));
        db.Set_b(i, j, nodes_i.Get_b(j));
        db.Set_c(i, j, nodes_i.Get_c(j));
      }
      db.Set_Einc(i, Ei);
    }

    //aCS=T3ElasticStrongIonIonCS<Floating>();
    
  }
  inline auto GetFS(T3LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID, unsigned int &generator_seed,
                    Floating * csBorderDataFS, Floating & tr/*temporary register of Particle struct*/,
                    int ind, ParticleTable & aParticleTable, MaterialTable & aMaterialTable,
                    int csBorderSTEP) const;
private:
  //******************************************************************************//
  const Floating alpha=1.0/137;
  PDG_t deuteronPDG;
  PDG_t titanPDG;
  Floating md;//deuteron mass.
  //displacement energies for D and Ti48:
  const Floating Edisplace_deuteron = 10 * eV;//displacement energy for deuteron is 10 eV.
  const Floating Edisplace_titan    = 25 * eV;//displacement energy for Ti48 is 25 eV.
  //T3ElasticStrongIonIonCS<Floating> aCS;
  T3ElasticStrongIonIon_DB db;//{-1.0, 1.0}
  //the Ox range for database with partial sums:
  //1-cos(theta_cm)_min=1.0e-3           1-cos(theta_cm)_max=1.0
  //ln(1-cos(theta_cm))_min=ln(1.0e-3)   ln(1-cos(theta_cm))_max=ln(1.0)=0.0:
  //Ox axis of the database, which we use for quadratic interpolation,
  //is in logarithmic scale
  const Floating Xmin=-1.0;
  const Floating Xmax= 1.0;
};

template <typename Floating>
auto T3ElasticStrongIonIonFS<Floating>::GetFS(T3LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID,
                                              unsigned int & generator_seed/*random number seed*/,
                                              Floating * csBorderDataFS, Floating & tr/*temporary register of Particle struct*/,
                                              int ind, ParticleTable & aParticleTable, MaterialTable & aMaterialTable,
                                              int csBorderSTEP) const
{

//1. This class is for D-D elastic scattering at Tls=30 keV - 250 MeV
//and at 10^(-3)*2*pcm^2 <= |t| <= 2*pcm^2.
//If the incident Particle is not a deuteron, do nothing: 
  if(incPDG != deuteronPDG) return Four<Floating>(PDG_t(0), PDG_t(0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0));
//2. HERE WE CHOOSE ON WHICH ISOTOPE OF MATERIAL TO SCATTER:
  const Floating m = aParticleTable.GetMass(incPDG);//mass of the incident ion
  const Floating E = p.E();//LS full enrgy of the incident ion
  const Floating Tls = E - m;//LS kinetic energy of the incident ion
//3. Check if Tls is in the necessary energy range for this proocess:
  const bool cond=(Tls >= 30 * keV) && (Tls <= 250 * MeV);
  if(!cond) return Four<Floating>(PDG_t(0), PDG_t(0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0));
//Number of isotopes in the material:
  const int numi = aMaterialTable.GetNumberOfIsotopes(matID);
//4. Now we should check if deuteron is among the isotopes
//of the material. If not found, do nothing.
  const bool isotopeFound=false;
  PDG_t isotopePDG=-1;//this is targetPDG.
  for(size_t i=0; i<numi; ++i)
  {
    if(deuteronPDG == aMaterialTable.GetIsotopes(matID,i))
    {
      isotopeFound=true;
      isotopePDG=deuteronPDG;
    }
  }
//5. This class is for D-D elastic scattering at Tls=30 keV - 250 MeV
//and at 10^(-3)*2*pcm^2 <= |t| <= 2*pcm^2.
//If the target Particle is not a deuteron, do nothing: 
  if(!isotopeFound) return Four<Floating>(PDG_t(0), PDG_t(0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0));
//6. Find |t|min:
  const Floating M=aParticleTable.GetMass(isotopePDG);//mass of the target ion.
  const int z=aParticleTable.GetZ(incPDG);           //z
  const int Z=aParticleTable.GetZ(isotopePDG);       //Z
  Floating tmin=0.0;
  if(isotopePDG == deuteronPDG)   tmin=2 * Edisplace_deuteron * M;
  else if(isotopePDG == titanPDG) tmin=2 * Edisplace_titan * M;
//7. Pls->Pcm:
  T3ThreeVector<Floating> InitMomentum = p.Vect();
  const Floating plsinc=InitMomentum.R();
  const Floating Elsinc=sqrt(m*m+plsinc*plsinc);
  const Floating s=m*m+2*Elsinc*M+M*M;
//8. Find |t|m and |t|max:
  const Floating pls2=Tls*(2*m+Tls);
  const Floating rat=m/M;
  const Floating rat2=rat*rat;
  const Floating pcm2=pls2/(1.0+rat2+2.0*sqrt(rat2+pls2/M/M));
  const Floating tm  =2*pcm2;
  const Floating tmax=4*pcm2;
  const Floating coef=M/sqrt(s);
  T3ThreeVector<Floating> Pcm = InitMomentum*coef;
//9. Find coordinate orts in CM:
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
//10. Get random cos(theta_cm):
  const Floating CosThetaCM=0.5;//db.RandomizeCost(generator_seed, Tls);
  tr=CosThetaCM;
//11. Random angle phi:
  const Floating phi=2*M_PI*RND01(generator_seed);
  const Floating cosPhi=cos(phi), sinPhi=sin(phi);
  const Floating SinThetaCM=sqrt(1.0-CosThetaCM*CosThetaCM);
  const Floating pcminc=Pcm.R();
  PDG_t outPDG1=incPDG;
  PDG_t outPDG2=isotopePDG;
  const Floating ss=pcminc*SinThetaCM*sinPhi, sc=pcminc*SinThetaCM*cosPhi;
//12. Rotate momentum vector in CM:
  T3ThreeVector<Floating> Pcmfin;
  Pcmfin.SetX(e1.x()*pcminc*CosThetaCM+e2.x()*ss+e3.x()*sc);
  Pcmfin.SetY(e1.y()*pcminc*CosThetaCM+e2.y()*ss+e3.y()*sc);
  Pcmfin.SetZ(e1.z()*pcminc*CosThetaCM+e2.z()*ss+e3.z()*sc);
//13. CM->LS:
  T3ThreeVector<Floating> Vcm=InitMomentum/(Elsinc+M);
  const Floating Ecminc=sqrt(m*m+pcminc*pcminc);
  const Floating E1cmfin=sqrt(m*m+pcminc*pcminc);
  ///\\\///const Floating E2cmfin=sqrt(M*M+pcminc*pcminc);
  const Floating AbsVcm=Vcm.R();
  const Floating gamma=1.0/std::sqrt(1-AbsVcm*AbsVcm);
//!!!
//14. Only the parallel to Vcm momentum component must
//change according to the Lorentz transformation.
//The perpendicular to Vcm component of the momentum vector
//is the same in LS and CM.
//!!!
  T3ThreeVector<Floating> ParallelProjection=(e1*Pcmfin)*e1;
  T3ThreeVector<Floating> PerpendicularProjection=Pcmfin-ParallelProjection;
  T3ThreeVector<Floating> plsfin1=gamma*(ParallelProjection+Vcm*E1cmfin);
  plsfin1+=PerpendicularProjection;
//15. The LS momentum of the target particle is found from the
//momentum conservation law.
  T3ThreeVector<Floating> plsfin2=InitMomentum-plsfin1;
  const Floating Absplsfin1=plsfin1.R();
  const Floating Elsfin1=std::sqrt(m*m+Absplsfin1*Absplsfin1);
  const Floating Absplsfin2=plsfin2.R();
  const Floating Elsfin2=std::sqrt(M*M+Absplsfin2*Absplsfin2);

  const Floating Efin=Elsfin1+Elsfin2;
  T3LorentzVector<Floating> outP1 = T3LorentzVector<Floating>(plsfin1.x(), plsfin1.y(), plsfin1.z(), Elsfin1);
  T3LorentzVector<Floating> outP2 = T3LorentzVector<Floating>(plsfin2.x(), plsfin2.y(), plsfin2.z(), Elsfin2);
  const T3LorentzVector<Floating> tot = outP1 + outP2;
  return Four<Floating>(outPDG1, outPDG2, outP1, outP2);  
}
  
}//namespace t3.
#endif//T3ELASTICSTRONGIONIONFSIMPL_H
