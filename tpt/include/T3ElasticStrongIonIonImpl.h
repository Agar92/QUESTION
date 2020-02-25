#pragma once
#ifndef T3ELASTICSTRONGIONIONIMPL_H
#define T3ELASTICSTRONGIONIONIMPL_H

#include <cmath>
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

#include "T3Defs.h"
#include "T3ThreeVector.h"
#include "T3LorentzVector.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

#include "T3NSGangular_RW.hh"
#include "T3R_DDCS.hh"
#include "T3ElasticStrongIonIon_DB.hh"

namespace t3 {

using namespace units;  
template <typename Floating>
class T3ElasticStrongIonIon
{
public:
  T3ElasticStrongIonIon()
  {
    ParticleTable aParticleTable;
    deuteronPDG=aParticleTable.makePDGfromZandA(1,2);
    titanPDG=aParticleTable.makePDGfromZandA(22,48);
//---------------------------------------------------------------------//
    lnEmin=std::log(Emin);
    lnEmax=std::log(Emax);
    DeltaLnE=(lnEmax-lnEmin)/(100-1);
    const int tgZ = 1;//target deuteron Z=1
    const int tgA = 2;//target deuteron A=2
    const t3::PDG_t sPDG = deuteronPDG;
    const int incZA = 1002;//1000*Z+A//inc deuteron
    const auto rid = "2";
    T3R_DDCS rdcs;
    std::cout<<"STEP #1"<<std::endl;
    rdcs.Load_from_CS();
    T3R_RW trw=rdcs.GetT3R_RW();
    std::cout<<"trw: "<<" size="<<trw.size()<<std::endl;
    for(int i=0; i<trw.size(); ++i)
    {
      E[i]=trw.NDI(i)->E();
      lnE[i]=std::log(E[i]);
      CS[i]=trw.NDI(i)->XS();
    }
    std::cout<<"Check cross sections in T3ElasticStrongIonIon():"<<std::endl;
    for(int i=0; i<100; ++i) std::cout<<"("<<E[i]<<","<<CS[i]/mbarn<<")   ";
    std::cout<<std::endl;

    md=aParticleTable.GetMass(deuteronPDG);
    //---------------------------------------------------------------------//
    rdcs.Load_from_PS(tgZ, tgA, rid, sPDG, incZA);
    T3NSGangular_RW rw=rdcs.GetT3NSGangular_RW();
    T3NSGangular_RWrecord rwrec=rw.at(0);
    const Floating Emin = rwrec.front().Get_E();
    const Floating Emax = rwrec.back().Get_E();
    const Floating Ediff = Emax - Emin;
    const size_t nbins = 512 - 1;
    const Floating dE = Ediff / nbins;
    db=T3ElasticStrongIonIon_DB(1.0e-3, 1.0);
    for(size_t i = 0; i < 512; ++i)
    {
      const Floating Ei = Emin + i*dE;
      const T3NSGangular_RWnode rwnode = rwrec.interpolate(Ei);
      T3NSGangular_node nodes_i = T3NSGangular_node(rwnode, Xmin, Xmax);
      for(int j=0; j<127; ++j)
      {
        db.Set_V(i, j, nodes_i.Get_V(j));
        db.Set_a(i, j, nodes_i.Get_a(j));
        db.Set_b(i, j, nodes_i.Get_b(j));
        db.Set_c(i, j, nodes_i.Get_c(j));
      }
      db.Set_Einc(i, Ei);
    }
    db.SetXmin(Xmin);
    db.SetXmax(Xmax);
  }
  
  inline Floating
  GetCS(Floating Tls, PDG_t incPDG, PDG_t targetPDG,
        ParticleTable & aParticleTable) const;
  
  inline Floating
  GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
        ParticleTable & aParticleTable, MaterialTable & aMaterialTable) const;

  inline auto GetFS(T3LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID, unsigned int &generator_seed,
                    Floating * csBorderDataFS, Floating & tr/*temporary register of Particle struct*/,
                    int ind, ParticleTable & aParticleTable, MaterialTable & aMaterialTable,
                    int csBorderSTEP) const;
  
private:
  PDG_t deuteronPDG;
  PDG_t titanPDG;
//-----------------------------------------------------------//
  const Floating Emin=10 * keV;
  const Floating Emax=250 * MeV;
  Floating lnEmin;
  Floating lnEmax;
  Floating DeltaLnE;
  Floating E[100]; //energies.
  Floating lnE[100];
  Floating CS[100];//cross sections.
  const Floating alpha=1.0/137;
  Floating md;//deuteron mass.
  const Floating Edisplace_deuteron = 10 * eV;//displacement energy for deuteron is 10 eV.
  const Floating Edisplace_titan    = 25 * eV;//displacement energy for Ti48 is 25 eV.
  T3ElasticStrongIonIon_DB db;//{-1.0, 1.0}
  const Floating Xmin= std::log(1.0e-3);
  const Floating Xmax= std::log(1.0);
  
};

template <typename Floating>
Floating T3ElasticStrongIonIon<Floating>::GetCS(Floating Tls, PDG_t incPDG, PDG_t targetPDG/*deuteronPDG*/,
                                                ParticleTable & aParticleTable) const
{
  const bool cond=(Tls >= Emin) && (Tls <= Emax);
  if(incPDG != deuteronPDG || !cond)
  {
    return 0.0;
  }
  const Floating lnTls=std::log(Tls);
  const int bin=(lnTls-lnEmin)/DeltaLnE;
  Floating cs=0.0;
  if(bin >= 0 && bin <= 99)
  {
    cs = CS[bin] + (CS[bin+1]-CS[bin]) * (Tls-E[bin]) / (E[bin+1]-E[bin]);
  }
  else
  {
    printf("***ERROR: T3ElasticStrongIonIon::GetCS() Out of energy range!!!");
  }
  return cs;
}

template <typename Floating>
Floating T3ElasticStrongIonIon<Floating>::GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
                                                ParticleTable & aParticleTable,
                                                MaterialTable & aMaterialTable) const
{
  const bool cond=(Tls >= Emin) && (Tls <= Emax); 
  if(incPDG != deuteronPDG || !cond)
  {
    return 0.0;
  }
  Floating cs=0.0;
  for(int i=0; i<aMaterialTable.GetNumberOfIsotopes(matID); ++i)
  {
    PDG_t isotopePDG = aMaterialTable.GetIsotopes(matID, i);
    if(isotopePDG == deuteronPDG)
    {
      cs += aMaterialTable.GetFractions(matID, i) *
        GetCS(Tls, incPDG, isotopePDG, aParticleTable);
    }
  }
  return cs;
}



template <typename Floating>
auto T3ElasticStrongIonIon<Floating>::GetFS(T3LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID,
                                            unsigned int & generator_seed/*random number seed*/,
                                            Floating * csBorderDataFS, Floating & tr/*temporary register of Particle struct*/,
                                            int ind, ParticleTable & aParticleTable, MaterialTable & aMaterialTable,
                                            int csBorderSTEP) const
{
  if(incPDG != deuteronPDG)
  {
    printf("***ERROR: T3ElasticStrongIonIon::GetFS(): incPDG is not deuteron!");
    return Four<Floating>(PDG_t(0), PDG_t(0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0));
  }
  const Floating m = aParticleTable.GetMass(incPDG);//mass of the incident ion
  const Floating E = p.E();//LS full energy of the incident ion
  const Floating Tls = E - m;//LS kinetic energy of the incident ion
  const bool cond=(Tls >= Emin) && (Tls <= Emax);
  if(!cond)
  {
    printf("***ERROR: T3ElasticStrongIonIon::GetFS(): Tls=%f keV is out of energy range! incPDG=%d\n", Tls/keV, incPDG);
    return Four<Floating>(PDG_t(0), PDG_t(0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0));
  }
  const int numi = aMaterialTable.GetNumberOfIsotopes(matID);
  bool isotopeFound=false;
  PDG_t isotopePDG=-1;//this is targetPDG.
  for(size_t i=0; i<numi; ++i)
  {
    if(deuteronPDG == aMaterialTable.GetIsotopes(matID,i))
    {
      isotopeFound=true;
      isotopePDG=deuteronPDG;
    }
  }
  if(!isotopeFound)
  {
    printf("***ERROR: T3ElasticStrongIonIon::GetFS(): deuteron is not found in the material!\n");
    return Four<Floating>(PDG_t(0), PDG_t(0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0));
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
  const Floating Ln1MinusCosThetaCM=db.RandomizeCost(generator_seed, Tls);
  const Floating CosThetaCM=1.0-std::exp(Ln1MinusCosThetaCM);
  tr=Ln1MinusCosThetaCM;
  const Floating phi=2*M_PI*RND01(generator_seed);
  const Floating cosPhi=cos(phi), sinPhi=sin(phi);
  const Floating SinThetaCM=sqrt(1.0-CosThetaCM*CosThetaCM);
  const Floating pcminc=Pcm.R();
  const PDG_t outPDG1=incPDG;
  const PDG_t outPDG2=isotopePDG;
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

  const Floating Efin=Elsfin1+Elsfin2;
  
  T3LorentzVector<Floating> outP1 = T3LorentzVector<Floating>(plsfin1.x(), plsfin1.y(), plsfin1.z(), Elsfin1);
  T3LorentzVector<Floating> outP2 = T3LorentzVector<Floating>(plsfin2.x(), plsfin2.y(), plsfin2.z(), Elsfin2);
  
  const T3LorentzVector<Floating> tot = outP1 + outP2;

  if(outPDG1!=deuteronPDG || outPDG2!=deuteronPDG)
    printf("***ERROR: T3ElasticStrongIonIon::GetFS(): wrong return PDG's: outPDG1=%d  outPDG2=%d!\n", outPDG1, outPDG2);
  
  return Four<Floating>(outPDG1, outPDG2, outP1, outP2);  
}

  

}
#endif//T3ELASTICSTRONGIONIONIMPL_H
