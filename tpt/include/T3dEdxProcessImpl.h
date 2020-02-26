#pragma once
#ifndef T3DEDXPROCESSIMPL_H
#define T3DEDXPROCESSIMPL_H

#include <cmath>
#include "T3Defs.h"
#include "T3ParticleTable.h"
#include "T3Utility.hh"

namespace t3 {

using namespace units;  
template <typename Floating>
class T3dEdxProcess
{
public:
  T3dEdxProcess() : protonPDG(PDG_t(2212)), deuteronPDG(t3::ParticleTable().makePDGfromZandA(1,2)),
                    titanPDG(t3::ParticleTable().makePDGfromZandA(22,48)),
                    mh(t3::ParticleTable().GetMass(protonPDG)),
                    md(t3::ParticleTable().GetMass(deuteronPDG)),
                    mtitan(t3::ParticleTable().GetMass(titanPDG))
                    
                    
  {
    std::ifstream input1("tpt/include/Electronic_dEdx_of_protons_in_Hydrogen_PSTAR.dat");
    int i=0;
    double Tlsh=0.0;
    double dEdxh=0.0;
    while(input1>>Tlsh>>dEdxh)
    {
      TlspinH[i]=Tlsh * MeV;
      dEdxpinH[i]=dEdxh * MeV*cm*cm/gr;
      ++i;
    }
    input1.close();
    for(int j=0; j<132; ++j) dEdxpinD[j]=dEdxpinH[j]/2;
    std::ifstream input2("tpt/include/Electronic_dEdx_of_protons_in_Titanium_PSTAR.dat");
    i=0;
    double Tlst=0.0;
    double dEdxt=0.0;
    while(input2>>Tlst>>dEdxt)
    {
      TlspinTi[i]=Tlst * MeV;
      dEdxpinTi[i]=dEdxt * MeV*cm*cm/gr;
      ++i;
    }
    input2.close();
    const double A_D=2.0;
    const double A_Ti=47.867;
    const double A_TiD2=2*A_D+A_Ti;
    const double wD=2.0*A_D/A_TiD2;
    const double wTi=A_Ti/A_TiD2;
    for(int j=0; j<132; ++j)
    {
      TlspinTiD2[j]   = TlspinH[j];
      dEdxpinTiD2[j]  = wD * dEdxpinD[j] + wTi * dEdxpinTi[j];
    }
    for(int j=0; j<132; ++j)
    {
      TlsdinTiD2[j] = TlspinTiD2[j] * md/mh;
      dEdxdinTiD2[j]   = dEdxpinTiD2[j] * (ZD/ZH) * (ZD/ZH);
    }
    for(int j=0; j<132; ++j)
    {
      TlstiinTiD2[j]   = TlspinTiD2[j] * mtitan/mh;
      dEdxtiinTiD2[j] = dEdxpinTiD2[j] * (ZTi/ZH) * (ZTi/ZH);
    }
  }

#ifdef OPENACC  
  #pragma acc routine seq
#endif
  inline Floating
  LinhardFunctionOfdEdxAtSmallTls(Floating y/* Tls/m */, int Z/* Z of the particle/ Z of H = 1 */) const;
  
#ifdef OPENACC  
  #pragma acc routine seq
#endif  
  inline Floating
  GetdEdxFromTable(Floating Tls, PDG_t incPDG, MatID_t matID, ParticleTable & aParticleTable) const;

#ifdef OPENACC  
  #pragma acc routine seq
#endif  
  inline Floating
  GetdEdxFromFunction(Floating y/* Tls/m */, int Z/* Z of the particle/ Z of H = 1 */) const;

#ifdef OPENACC  
  #pragma acc routine seq
#endif  
  inline Floating
  GetdEdxFromFunction(Floating Tls, PDG_t incPDG, MatID_t matID, ParticleTable & aParticleTable) const;
  
private:
  const PDG_t protonPDG;
  const PDG_t deuteronPDG;
  const PDG_t titanPDG;
  const Floating mh;
  const Floating md;
  const Floating mtitan;
//--------------------------------------------------------------//  
  const int ZH=1;
  double TlspinH[132]{0.0};
  double dEdxpinH[132]{0.0};
//-----------------------------//
  double dEdxpinD[132]{0.0};  
//-----------------------------//  
  double TlspinTi[132]{0.0};
  double dEdxpinTi[132]{0.0};
//-----------------------------//  
  double TlspinTiD2[132]{0.0};
  double dEdxpinTiD2[132]{0.0};
//--------------------------------------------------------------//  
  const int ZD=1;
  Floating TlsdinTiD2[132]{0.0};
  Floating dEdxdinTiD2[132]{0.0};
//--------------------------------------------------------------//
  const int ZTi=22;
  Floating TlstiinTiD2[132]{0.0};
  Floating dEdxtiinTiD2[132]{0.0};
//--------------------------------------------------------------//  
};


template <typename Floating>
Floating T3dEdxProcess<Floating>::LinhardFunctionOfdEdxAtSmallTls(Floating y/* Tls/m */, int Z/* Z of the particle/ Z of H = 1 */) const
{
  return ( 1.02/( 0.74e-5*std::pow(y, -0.54))*Z*Z * MeV*cm*cm/gr );
}

template <typename Floating>
Floating T3dEdxProcess<Floating>::GetdEdxFromTable(Floating Tls, PDG_t incPDG, MatID_t matID, ParticleTable & aParticleTable) const
{
  if(incPDG!=protonPDG && incPDG!=deuteronPDG && incPDG!=titanPDG) printf("***ERROR: T3dEdxProcess::GetdEdxFromTable: Wrong PDG: incPDG=%i", incPDG);
  Floating dEdx=0.0;
  if(matID==MatID_t(2))
  {
    const size_t SIZE=132;
    if(incPDG==protonPDG)
    {
      if(Tls>TlspinTiD2[131]) printf("***ERROR: T3dEdxProcess::GetdEdxFromTable(): Out of energy range! Tls=%f > Tlsmax=%f", Tls, TlspinTiD2[131]);
      else if(Tls<TlspinTiD2[0])
      {
        const double y=Tls/mh;
        dEdx=LinhardFunctionOfdEdxAtSmallTls(y, ZH);
      }
      else
      {
        const size_t hct = T3Utility::bin_search<SIZE>(TlspinTiD2, Tls) - TlspinTiD2 - 1;
        dEdx=dEdxpinTiD2[hct]+(dEdxpinTiD2[hct+1]-dEdxpinTiD2[hct]) *
          (Tls-TlspinTiD2[hct])/(TlspinTiD2[hct+1]-TlspinTiD2[hct]);
      }
      
    }
    else if(incPDG==deuteronPDG)
    {
      if(Tls>TlsdinTiD2[131]) printf("***ERROR: T3dEdxProcess::GetdEdxFromTable(): Out of energy range! Tls=%f > Tlsmax=%f", Tls, TlsdinTiD2[131]);
      else if(Tls<TlsdinTiD2[0])
      {
        const double y=Tls/md;
        dEdx=LinhardFunctionOfdEdxAtSmallTls(y, ZD);
      }
      else
      {
        const size_t hct = T3Utility::bin_search<SIZE>(TlsdinTiD2, Tls) - TlsdinTiD2 - 1;
        dEdx=dEdxdinTiD2[hct]+(dEdxdinTiD2[hct+1]-dEdxdinTiD2[hct]) *
          (Tls-TlsdinTiD2[hct])/(TlsdinTiD2[hct+1]-TlsdinTiD2[hct]);
      }
    }
    else if(incPDG==titanPDG)
    {
      if(Tls>TlstiinTiD2[131]) printf("***ERROR: T3dEdxProcess::GetdEdxFromTable(): Out of energy range! Tls=%f > Tlsmax=%f", Tls, TlstiinTiD2[131]);
      else if(Tls<TlstiinTiD2[0])
      {
        const double y=Tls/mtitan;
        dEdx=LinhardFunctionOfdEdxAtSmallTls(y, ZTi);
      }
      else
      {
        const size_t hct = T3Utility::bin_search<SIZE>(TlstiinTiD2, Tls) - TlstiinTiD2 - 1;
        dEdx=dEdxtiinTiD2[hct]+(dEdxtiinTiD2[hct+1]-dEdxtiinTiD2[hct]) *
          (Tls-TlstiinTiD2[hct])/(TlstiinTiD2[hct+1]-TlstiinTiD2[hct]);

      }
    }
  }
  return dEdx;
}

template <typename Floating>
Floating T3dEdxProcess<Floating>::GetdEdxFromFunction(Floating y/* Tls/m */, int Z/* Z of the particle/ Z of H = 1 */) const
{
  const Floating term1=1.02/(std::pow(y, 0.746)+0.74e-5*std::pow(y, -0.54));
  Floating phitx=std::log(1.+y);
  for(int j=1; j<=3; ++j) phitx=std::log(1.+phitx);
  const Floating term2=3.56*std::pow(phitx, 1.65);
  const Floating dEdx=(term1+term2)*Z*Z;
  return dEdx * cm*cm/gr;
}

template <typename Floating>
Floating T3dEdxProcess<Floating>::GetdEdxFromFunction(Floating Tls, PDG_t incPDG, MatID_t matID, ParticleTable & aParticleTable) const
{
  if(incPDG!=protonPDG && incPDG!=deuteronPDG && incPDG!=titanPDG) printf("***ERROR: T3dEdxProcess::GetdEdxFromFunction: Wrong PDG: incPDG=%i", incPDG);  
  const Floating m=aParticleTable.GetMass(incPDG);
  const int Z=aParticleTable.GetZ(incPDG);
  const Floating y=Tls/m;
  Floating dEdx=0.0;
  if(matID==MatID_t(2)) dEdx=GetdEdxFromFunction(y,Z);
  return dEdx;
}

}
#endif//T3DEDXPROCESSIMPL_H
