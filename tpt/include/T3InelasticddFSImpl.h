#pragma once
#ifndef T3INELASTICDDFSIMPL_H
#define T3INELASTICDDFSIMPL_H

#include <random>
#include "T3Defs.h"
#include "T3ThreeVector.h"
#include "T3LorentzVector.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"
#include "T3InelasticddCSImpl.h"

#include "T3Inelasticdd_DB.hh"
#include "T3NSGangular_RW.hh"
#include "T3NSGangular_node.hh"

#include <typeinfo>

namespace t3 {

template <typename Floating = double>
class InelasticddFS {
public:
  InelasticddFS()
  { 
    ParticleTable aParticleTable;
    neutronPDG  = PDG_t(2112);
    protonPDG   = PDG_t(2212);
    deuteronPDG = aParticleTable.makePDGfromZandA(1, 2);
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
    constexpr auto sPDG = 2112;//reaction product=n+3He//this is a neutron from a reaction product.
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

  inline auto GetFS(T3LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID, unsigned int & generator_seed,
                    Floating * csBorderDataFS, Floating & tr,
                    int ind, ParticleTable & aParticleTable, MaterialTable & aMaterialTable,
                    int csBorderSTEP) const;
                    
private:
  PDG_t neutronPDG;
  PDG_t protonPDG;
  PDG_t deuteronPDG;
  PDG_t tritonPDG;
  PDG_t he3PDG;
  const Floating Eg = 0.986 * MeV;
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
  T3Inelasticdd_DB db;//{-1.0, 1.0};
  const Floating Xmin=-1.0;
  const Floating Xmax= 1.0;
};

template <typename Floating>
auto InelasticddFS<Floating>::GetFS(T3LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID,
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
    return Four<Floating>(PDG_t(0), PDG_t(0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0), T3LorentzVector<Floating>(0.0,0.0,0.0,0.0));
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
      const Floating nhe31 = 1.0 + pow(2.52e-2/Tcm, 1.5);
      const Floating nhe32 = 1.0 + pow(4.7e-3/Tcm, 4.5);
      const Floating csnhe3 = 0.21 / nhe31 / nhe32 / (1.+0.23*Tcm) * g * barn;
      const Floating pt1 = 1.0 + pow(2.0e-2/Tcm, 1.7);
      const Floating pt2 = 1.0 + pow(4.3e-3/Tcm, 4.4);
      const Floating pt3 = 1.0 + pow(0.17*Tcm, 0.85);
      const Floating cspt = 0.176/ pt1 / pt2 / pt3 * g * barn;
      const Floating cstot = csnhe3 + cspt;
      const Floating R = RND01(generator_seed);//GenerateSubCanonical<Floating>(engine);      
      if(0.0<=R && R<csnhe3/cstot)       outPDG1 = PDG_t(2112);//neutron.
      else if(csnhe3/cstot<=R && R<=1.0) outPDG1 = PDG_t(2212);//proton.
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
      const Floating cosThetaCM=db.RandomizeCost(generator_seed, Tls);
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
    }
  }
  return Four<Floating>(outPDG1, outPDG2, outP1, outP2);
}
  
}//namespace t3.
#endif // T3INELASTICDDFSIMPL_H
