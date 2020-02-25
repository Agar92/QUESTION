#pragma once
#ifndef T3INELASTICDDIMPL_H
#define T3INELASTICDDIMPL_H

#include <cmath>
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

#include "T3Defs.h"
#include "T3ThreeVector.h"
#include "T3LorentzVector.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

#include "T3Inelasticdd_DB.hh"
#include "T3NSGangular_RW.hh"
#include "T3NSGangular_node.hh"

namespace t3 {

using namespace units;
template <typename Floating>
class T3Inelasticdd
{
public:
  T3Inelasticdd()
  {
    ParticleTable aParticleTable;
    deuteronPDG=aParticleTable.makePDGfromZandA(1,2);
    neutronPDG  = PDG_t(2112);
    protonPDG   = PDG_t(2212);
    tritonPDG   = aParticleTable.makePDGfromZandA(1, 3);
    he3PDG      = aParticleTable.makePDGfromZandA(2, 3);
    mn          = aParticleTable.GetMass(neutronPDG);
    mp          = aParticleTable.GetMass(protonPDG);
    md          = aParticleTable.GetMass(deuteronPDG);
    mt          = aParticleTable.GetMass(tritonPDG);
    mhe3        = aParticleTable.GetMass(he3PDG);
    md2         = md * md;
    D2Plusn_he3  = (mn + mhe3)*(mn + mhe3);
    D2Minusn_he3 = (mn - mhe3)*(mn - mhe3);
    D2Plusp_t    = (mp + mt)*(mp + mt);
    D2Minusp_t   = (mp - mt)*(mp - mt);
    T3NSGangular_RW rw;  
    constexpr int tgZ = 1;//target deuteron Z=1
    constexpr int tgA = 2;//target deuteron A=2
    constexpr int sPDG = 2112;//reaction product=n+3He//this is a neutron from a reaction product.
    constexpr int incZA = 1002;//1000*Z+A//inc deuteron
    constexpr auto rid = "50";//MT index of inelastic reaction in ENDF
    rw.load_binary(tgZ, tgA, rid, sPDG, incZA);//load the data from the binary file.
    T3NSGangular_RWrecord rwrec=rw.at(0);
    for(int i=0; i<rwrec.size(); ++i)
    {
      for(int j=0; j<127; ++j) rwrec.at(i).Set_V(j+1, rwrec.at(i).Get_V(j+1)*rwrec.at(i).Get_V(j+1));
    }
    const Floating Emin = rwrec.front().Get_E();
    const Floating Emax = rwrec.back().Get_E();
    const Floating Ediff = Emax - Emin;
    const size_t nbins = 512 - 1;
    const Floating dE = Ediff / nbins;
    db=T3Inelasticdd_DB(-1.0, 1.0);
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
      const Floating pri=nodes_i.Get_pr();
      const Floating sli=nodes_i.Get_sl();
      db.Set_pr(i, pri);
      db.Set_sl(i, sli);
      db.Set_Einc(i, Ei);
    }        
  }

  inline Floating
  GetDDInelasticIntegralCS(Floating Tls) const;

  inline Floating
  GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
        ParticleTable & aParticleTable, MaterialTable & aMaterialTable) const;

  inline auto GetFS(T3LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID, unsigned int & generator_seed,
                    Floating * csBorderDataFS, Floating & tr,
                    int ind, ParticleTable & aParticleTable, MaterialTable & aMaterialTable,
                    int csBorderSTEP) const;

  inline Floating C2nhe3(Floating Tls) const
  {
    Floating t1=(1.0+pow(Tls/0.0585, 0.86));
    Floating t2=(1.0+pow(Tls/0.44, 1.46));
    Floating t3=(1.0+pow(Tls/1.35, 4.7));
    Floating res=0.142*t1/t2/t3;
    return res;
  }

  inline Floating C2pt(Floating Tls) const
  {
    Floating t1=(1.0+pow(Tls/0.047, 0.855));
    Floating t2=(1.0+pow(Tls/0.43, 1.8));
    Floating t3=(1.0+pow(Tls/0.95, 6.0));
    Floating t4=(1.0+pow(Tls/1.056, 27.5));
    Floating res=0.08*t1/t2/t3/t4;
    return res;
  }

  inline Floating C4nhe3(Floating Tls) const
  {
    Floating t1=0.414e-3*pow(Tls/0.11e-3, 0.0543);
    Floating t2=(1.0+Tls/0.011);
    Floating t3=(1.0+pow(Tls/0.0451, 0.5));
    Floating res=t1*t2*t3;
    return res;
  }

  inline Floating C4pt(Floating Tls) const
  {
    Floating t1=0.361e-3*pow(Tls/0.16e-3, 0.034);
    Floating t2=(1.+Tls/0.01);
    Floating t3=(1.+pow(Tls/0.013, 0.524));
    Floating res=t1*t2*t3;
    return res;
  }

  inline Floating C6nhe3(Floating Tls) const
  {
    Floating t1=0.15e-6*pow(Tls, 0.13);
    Floating t2=(1.0+pow(Tls/0.703e-2, 1.5));
    Floating t3=(1.0+pow(Tls/0.0359, 0.8));
    Floating t4=(1.0+pow(Tls/0.788, 0.84));
    Floating res=t1*t2*t3*t4;
    return res;
  }
  
  inline Floating C6pt(Floating Tls) const
  {
    Floating t1=0.185e-6*pow(Tls, 0.0785);
    Floating t2=(1.0+pow(Tls/0.538e-2, 1.46));
    Floating t3=(1.0+pow(Tls/0.045, 1.19));
    Floating res=t1*t2*t3;
    return res;
  }
  
  inline Floating C0nhe3(Floating Tls) const
  {
    Floating res=0.5-C2nhe3(Tls)/3-C4nhe3(Tls)/5-C6nhe3(Tls)/7;
    return res;
  }
  
  inline Floating C0pt(Floating Tls) const
  {
    Floating res=0.5-C2pt(Tls)/3-C4pt(Tls)/5-C6pt(Tls)/7;
    return res;
  }

  inline Floating GetCosCMFromDDInelasticChannelsApproximationFunctions(PDG_t secondaryPDG, unsigned int & generator_seed, Floating Tls) const;

private:
  PDG_t deuteronPDG;
  const Floating Eg = 0.986 * MeV;
  //---------------------------------------------//
  PDG_t neutronPDG;
  PDG_t protonPDG;
  PDG_t tritonPDG;
  PDG_t he3PDG;
  Floating mn;
  Floating mp;
  Floating md;
  Floating mt;
  Floating mhe3;
  Floating md2;
  Floating D2Plusn_he3;
  Floating D2Minusn_he3;
  Floating D2Plusp_t;
  Floating D2Minusp_t;
  T3Inelasticdd_DB db;
  const Floating Xmin=-1.0;
  const Floating Xmax= 1.0;
  
};

template <typename Floating>
Floating T3Inelasticdd<Floating>::GetDDInelasticIntegralCS(Floating Tls) const
{
  const Floating Tcm = Tls / 2;
  const Floating rat = std::sqrt(Eg / Tcm);
  const Floating g = std::exp(-0.5*rat);
  const Floating nhe31 = 1.0 + std::pow(2.52e-2/Tcm, 1.5);
  const Floating nhe32 = 1.0 + std::pow(4.7e-3/Tcm, 4.5);
  const Floating csnhe3 = 0.21 / nhe31 / nhe32 / (1.+0.23*Tcm) * g * units::barn;
  const Floating pt1 = 1.0 + std::pow(2.0e-2/Tcm, 1.7);
  const Floating pt2 = 1.0 + std::pow(4.3e-3/Tcm, 4.4);
  const Floating pt3 = 1.0 + std::pow(0.17*Tcm, 0.85);
  const Floating cspt = 0.176/ pt1 / pt2 / pt3 * g * units::barn;
  const Floating cs = csnhe3 + cspt;
  return cs;
}

template <typename Floating>
Floating T3Inelasticdd<Floating>::GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
                                        ParticleTable & aParticleTable,
                                        MaterialTable & aMaterialTable) const
{
  if(incPDG != deuteronPDG) return 0.0;
  Floating cs=0.0;
  for(int i=0; i<aMaterialTable.GetNumberOfIsotopes(matID); ++i)
  {
    PDG_t isotope=aMaterialTable.GetIsotopes(matID, i);
    if(isotope == deuteronPDG)
    {
      cs += aMaterialTable.GetFractions(matID, i) *
        GetDDInelasticIntegralCS(Tls);
    }
  }
  return cs;
}

template <typename Floating>
Floating T3Inelasticdd<Floating>::GetCosCMFromDDInelasticChannelsApproximationFunctions(PDG_t secondaryPDG, unsigned int &generator_seed/*the particle random generator seed*/,
                                                                                        Floating Tls/*LS kinetic energy of the incident deuteron*/) const
{
  Floating b0, b2, b4, b6;
  if(secondaryPDG==neutronPDG)
  {
    const Floating C21=C2nhe3(Tls);
    const Floating C41=C4nhe3(Tls);
    const Floating C61=C6nhe3(Tls);
    const Floating C01=0.5-C21/3-C41/5-C61/7;
    const Floating Ctotnhe3=C01+C21/3+C41/5+C61/7;
    b0=C01;
    b2=b0+C21/3;
    b4=b2+C41/5;
    b6=b4+C61/7;
    b0/=Ctotnhe3;
    b2/=Ctotnhe3;
    b4/=Ctotnhe3;
    b6/=Ctotnhe3;
  }
  else if(secondaryPDG==protonPDG)
  {
    const Floating C22=C2pt(Tls);
    const Floating C42=C4pt(Tls);
    const Floating C62=C6pt(Tls);
    const Floating C02=0.5-C22/3-C42/5-C62/7;
    const Floating Ctotpt=C02+C22/3+C42/5+C62/7;
    b0=C02;
    b2=b0+C22/3;
    b4=b2+C42/5;
    b6=b4+C62/7;
    b0/=Ctotpt;
    b2/=Ctotpt;
    b4/=Ctotpt;
    b6/=Ctotpt;
  }  
  const Floating R=RND01(generator_seed);
  Floating coscm=0.0;
  const Floating r=RND01(generator_seed);
  if     (R<b0) coscm=r;
  else if(R<b2) coscm=pow(r, 1.0/3.0);
  else if(R<b4) coscm=pow(r, 1.0/5.0);
  else if(R<b6) coscm=pow(r, 1.0/7.0);
  const Floating RSIGN=RND01(generator_seed);
  const int SIGN=((RSIGN>0.5)?1:-1);
  coscm=coscm*SIGN;
  return coscm;
}

  
template <typename Floating>
auto T3Inelasticdd<Floating>::GetFS(T3LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID,
                                    unsigned int & generator_seed,
                                    Floating * csBorderDataFS, Floating & tr,
                                    int ind, ParticleTable & aParticleTable, MaterialTable & aMaterialTable,
                                    int csBorderSTEP) const
{
  auto const generateSubCanonical = [&generator_seed]() {
    generator_seed=Rand32(generator_seed);
    const Floating result=rndv(generator_seed);
    return result;
  };
  if(incPDG != deuteronPDG)
  {
    printf("***ERROR: T3Inelasticdd::GetFS() incPDG=%d != deuteronPDG=%d", incPDG, deuteronPDG);
    return Four<Floating>(PDG_t(0), PDG_t(0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0));
  }
  const Floating m = aParticleTable.GetMass(incPDG);
  const Floating E = p.E();
  const Floating Tls = E - m;
  PDG_t outPDG1=-1;
  PDG_t outPDG2=-1;
  const Floating Einit=E+md;
  T3LorentzVector<Floating> outP1 = T3LorentzVector<Floating>(0.0,0.0,0.0,0.0);
  T3LorentzVector<Floating> outP2 = T3LorentzVector<Floating>(0.0,0.0,0.0,0.0);
  for(int i=0; i<aMaterialTable.GetNumberOfIsotopes(matID); ++i)
  {
    const PDG_t isotope=aMaterialTable.GetIsotopes(matID, i);
    if(isotope == deuteronPDG)
    {
      const Floating Tcm = Tls / 2;
      const Floating rat = std::sqrt(Eg / Tcm);
      const Floating g = std::exp(-0.5*rat);
      const Floating nhe31 = 1.0 + std::pow(2.52e-2/Tcm, 1.5);
      const Floating nhe32 = 1.0 + std::pow(4.7e-3/Tcm, 4.5);
      const Floating csnhe3 = 0.21 / nhe31 / nhe32 / (1.+0.23*Tcm) * g * barn;
      const Floating pt1 = 1.0 + std::pow(2.0e-2/Tcm, 1.7);
      const Floating pt2 = 1.0 + std::pow(4.3e-3/Tcm, 4.4);
      const Floating pt3 = 1.0 + std::pow(0.17*Tcm, 0.85);
      const Floating cspt = 0.176/ pt1 / pt2 / pt3 * g * barn;
      const Floating cstot = csnhe3 + cspt;
      const Floating R = RND01(generator_seed);
      if(0.0<=R && R<csnhe3/cstot)       outPDG1 = PDG_t(2112);
      else if(csnhe3/cstot<=R && R<=1.0) outPDG1 = PDG_t(2212);
      
      T3ThreeVector<Floating> InitMomentum = p.Vect();
      const Floating plsinc=InitMomentum.R();
      const Floating Elsinc=sqrt(md2+plsinc*plsinc);
      const Floating s=2*md2+2*Elsinc*md;
      const Floating coef=md/sqrt(s);
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
      const Floating cosThetaCM=GetCosCMFromDDInelasticChannelsApproximationFunctions(outPDG1, generator_seed, Tls);
      tr=cosThetaCM;      
      const Floating phi=2*M_PI*RND01(generator_seed);
      const Floating cosPhi=cos(phi), sinPhi=sin(phi);
      const Floating sinThetaCM=sqrt(1.0-cosThetaCM*cosThetaCM);
      const Floating pcminc=Pcm.R();
      const Floating pcminc2=pcminc*pcminc;
      const Floating pcminc4=pcminc2*pcminc2;
      Floating num=-1.0;
      Floating prodm1=-1.0;
      Floating prodm2=-1.0;
      if(outPDG1==PDG_t(2112))
      { 
        num = ( s - D2Plusn_he3 ) * ( s - D2Minusn_he3 );
        outPDG2=he3PDG;
        prodm1=mn;
        prodm2=mhe3;
      }
      else if(outPDG1==PDG_t(2212))
      { 
        num = ( s - D2Plusp_t ) * ( s - D2Minusp_t );
        outPDG2=tritonPDG;
        prodm1=mp;
        prodm2=mt;
      }
      const Floating pcmfin=sqrt(num)/2/sqrt(s);      
      const Floating ss=pcmfin*sinThetaCM*sinPhi, sc=pcmfin*sinThetaCM*cosPhi;
      T3ThreeVector<Floating> Pcmfin;
      const Floating pcmfinCosTheta=pcmfin*cosThetaCM;
      Pcmfin.SetX(e1.x()*pcmfinCosTheta+e2.x()*ss+e3.x()*sc);
      Pcmfin.SetY(e1.y()*pcmfinCosTheta+e2.y()*ss+e3.y()*sc);
      Pcmfin.SetZ(e1.z()*pcmfinCosTheta+e2.z()*ss+e3.z()*sc);
      T3ThreeVector<Floating> Vcm=InitMomentum/(Elsinc+md);
      const Floating Ecminc=sqrt(md2+pcminc*pcminc);
      const Floating E1cmfin=sqrt(prodm1*prodm1+pcmfin*pcmfin);
      const Floating AbsVcm=Vcm.R();
      const Floating gamma=1.0/std::sqrt(1-AbsVcm*AbsVcm);
      T3ThreeVector<Floating> ParallelProjection=(e1*Pcmfin)*e1;
      T3ThreeVector<Floating> PerpendicularProjection=Pcmfin-ParallelProjection;
      T3ThreeVector<Floating> plsfin1=gamma*(ParallelProjection+Vcm*E1cmfin);
      plsfin1+=PerpendicularProjection;
      T3ThreeVector<Floating> plsfin2=InitMomentum-plsfin1;
      const Floating Absplsfin1=plsfin1.R();
      const Floating Elsfin1=std::sqrt(prodm1*prodm1+Absplsfin1*Absplsfin1);
      const Floating Absplsfin2=plsfin2.R();
      const Floating Elsfin2=std::sqrt(prodm2*prodm2+Absplsfin2*Absplsfin2);
      const Floating Efin=Elsfin1+Elsfin2;
      
      outP1 = T3LorentzVector<Floating>(plsfin1.x(), plsfin1.y(), plsfin1.z(), Elsfin1);
      outP2 = T3LorentzVector<Floating>(plsfin2.x(), plsfin2.y(), plsfin2.z(), Elsfin2);

      T3LorentzVector<Floating> totP=outP1+outP2;

    }
  }

  if(!(outPDG1==protonPDG || outPDG1==neutronPDG) || !(outPDG2==he3PDG || outPDG2==tritonPDG))
  {
    printf("***ERROR: T3Inelasticdd::GetFS(): wrong return PDG: outPDG1=%d outPDG2=%d", outPDG1, outPDG2);
  }
  
  return Four<Floating>(outPDG1, outPDG2, outP1, outP2);
}


}
#endif//T3INELASTICDDIMPL_H
