#pragma once
#ifndef T3ELASTICSTRONGIONIONCSIMPL_H
#define T3ELASTICSTRONGIONIONCSIMPL_H

#include <cmath>
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

//**********************************************************************************//
//This class returns integral cross section of elastic D-D scattering approximation //
//from 10^(-3)*2*pcm^2 to 2*pcm^2 in inner units.                                   //
//This is the integral cross section of ElasticStrongIonIon process.                //
//THIS IS YET FOR D-D ONLY!!!                                                       //
//BECAUSE ONLY FOR D-D WE HAVE THE APPROXIMATION OF ELASTIC SCATTERING DIFFRENTIAL  //
//CROSS SECTION (THE INTERFERNCE OF THE ELECTROMAGNETIC AND NYCLEAR AMPLITUDES OF   //
//SCATTERING). IN THE FUTURE IT IS PLANNED TO ADD HERE D-t and some other particles //
//scattering approximation partial sums.                                            //
//**********************************************************************************//

namespace t3 {

using namespace units;  
template <typename Floating>
class T3ElasticStrongIonIonCS
{
public:
  T3ElasticStrongIonIonCS()
  {
    ParticleTable aParticleTable;
    deuteronPDG=aParticleTable.makePDGfromZandA(1,2);
    titanPDG=aParticleTable.makePDGfromZandA(22,48);
//---------------------------------------------------------------------//
    lnEmin=std::log(Emin);
    lnEmax=std::log(Emax);
    //calculate the logarithmic step of energy axis:
    //100 energy values (bin walls) => 99 bins:
    DeltaLnE=(lnEmax-lnEmin)/(100-1);
    const int tgZ = 1;//target deuteron Z=1
    const int tgA = 2;//target deuteron A=2
    //reaction product D+D.
    //this is PDG of 1 of the deuterons in the reaction product.
    const t3::PDG_t sPDG = deuteronPDG;
    const int incZA = 1002;//1000*Z+A//inc deuteron
    //MT index of reaction in ENDF
    //For elastic scattering MT index of reaction in ENDF = 2. See ENDF6.
    const auto rid = "2";
//load the integral cross sections from T3_DATA/ database:
    
    T3R_DDCS rdcs;
    std::cout<<"STEP #1"<<std::endl;

    rdcs.Load_from_CS();

    
    T3R_RW trw=rdcs.GetT3R_RW();

    //*
    std::cout<<"trw: "<<" size="<<trw.size()<<std::endl;
    std::cout<<trw<<std::endl;
    for(int i=0; i<trw.size(); ++i)
    {
      E[i]=trw.NDI(i)->E();
      lnE[i]=std::log(E[i]);      
      CS[i]=trw.NDI(i)->XS();
    }
    //Check:
    std::cout<<"Check cross sections:"<<std::endl;
    for(int i=0; i<100; ++i) std::cout<<"("<<E[i]<<","<<CS[i]<<")   ";
    std::cout<<std::endl;

    //*/

    //exit(0);
    
//loaded.
  }
  
  inline Floating
  GetCS(Floating Tls, PDG_t incPDG, PDG_t targetPDG,
        ParticleTable & aParticleTable) const;
  
  inline Floating
  GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
        ParticleTable & aParticleTable, MaterialTable & aMaterialTable) const;
  
private:
  PDG_t deuteronPDG;
  PDG_t titanPDG;
//-----------------------------------------------------------//
  //E = Tls !!!
  const Floating Emin=30 * keV;
  const Floating Emax=250 * MeV;
  Floating lnEmin;
  Floating lnEmax;
  Floating DeltaLnE;
  Floating lnE[100];
  Floating E[100]; //energies.
  Floating CS[100];//cross sections.
};
  
template <typename Floating>
Floating T3ElasticStrongIonIonCS<Floating>::GetCS(Floating Tls, PDG_t incPDG, PDG_t targetPDG,
                                                  ParticleTable & aParticleTable) const
{
//1. Check if incPDG is deuteronPDG. This class is for getting integral
//cross section approximation for D-D elastic scattering in the energy range
//Tls = 30 * keV - 250 * MeV and for 10^(-3)*2*pcm^2 <= |t| <= 2*pcm^2.
//If incPDG != deuteronPDG, then:
//2. Check if Tls is in the necessary energy range for this proocess:
  const bool cond=(Tls >= 30 * keV) && (Tls <= 250 * MeV);
  if(incPDG != deuteronPDG || !cond)
  {
    std::cout<<"*** Out of energy range or not a deuteron!!!"<<std::endl;
    return 0.0;
  }
//3. find the number of energy bin:
  const Floating lnTls=std::log(Tls);
  const int bin=(lnTls-lnEmin)/DeltaLnE;
//4. interpolate:
  Floating cs=0.0;
  if(bin >= 0 && bin <= 99)
    cs = CS[bin] + (CS[bin+1]-CS[bin]) * (Tls-E[bin]) / (E[bin+1]-E[bin]);
  else
    std::cout<<"*** Out of energy range!!!"<<std::endl;
  return cs;
}

template <typename Floating>
Floating T3ElasticStrongIonIonCS<Floating>::GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
                                                  ParticleTable & aParticleTable,
                                                  MaterialTable & aMaterialTable) const
{
//1. Check if incPDG is deuteronPDG. This class is for getting integral
//cross section approximation for D-D elastic scattering in the energy range
//Tls = 30 * keV - 250 * MeV and for 10^(-3)*2*pcm^2 <= |t| <= 2*pcm^2.
//If incPDG != deuteronPDG, then:
//2. Check if Tls is in the necessary energy range for this proocess:
  const bool cond=(Tls >= 30 * keV) && (Tls <= 250 * MeV); 
  if(incPDG != deuteronPDG || !cond) return 0.0;
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

}
#endif//T3ELASTICSTRONGIONIONCSIMPL_H
