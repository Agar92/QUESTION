#pragma once
#ifndef T3ELASTICEMIONIONCSIMPL_H
#define T3ELASTICEMIONIONCSIMPL_H

#include <cmath>
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

//***************************************************************************************//
//This class returns integral cross section of ion-ion elastic scattering in inner units.//
//1. For equal particles it returns the integral of Rutherford differential cs for       //
//equal particles from |t|min to |t|max-|t|min.                                          //
//2. For different particles it returns the integral of Rutherford differential cross    //
//section from |t|min to |t|max.                                                         //
//3. For D-D:                                                                            //
//1) if(Tls<30 keV)   => Rutherford scattering from |t|min to |t|max-|t|min.             //
//2) if(Tls>= 30 keV) => Rutherford scattering from |t|min to 10^(-3) * 2*pcm^2 =        //
// = 10^(-3) * |t|m.                                                                     //
//***************************************************************************************//

namespace t3 {

using namespace units;  
template <typename Floating>
class T3ElasticEMIonIonCS
{
public:
  T3ElasticEMIonIonCS()
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
  
private:
  PDG_t deuteronPDG;
  PDG_t titanPDG;
  const Floating Eg = 0.986 * MeV;
  const Floating alpha=1.0/137;
  const Floating hc=200.0 * MeV * fm;
  //the border for the D-D approximation partial sums
  //(the partial sums are at kb-:- 1.0)
  const Floating kb=1.0e-3;
  //displacement energies for D and Ti48:
  const Floating Edisplace_deuteron = 10 * eV;//displacement energy for deuteron is 10 eV.
  const Floating Edisplace_titan    = 25 * eV;//displacement energy for Ti48 is 25 eV.
};

template <typename Floating>
Floating T3ElasticEMIonIonCS<Floating>::GetCS(Floating Tls, PDG_t incPDG, PDG_t targetPDG,
                                              ParticleTable & aParticleTable) const
{
  //if the particle is not a nucleus (neutron, proton, electron, positron), return 0:
  if(!aParticleTable.IsNucleus(incPDG)) return 0.0;
  
  //masses:
  const Floating m=aParticleTable.GetMass(incPDG);   // m
  const Floating M=aParticleTable.GetMass(targetPDG);// M
  //charges:
  const int z=aParticleTable.GetZ(incPDG);           // z
  const int Z=aParticleTable.GetZ(targetPDG);        // Z
//Calculate pcm2 and |t|max:
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
//Calculate |t|min:
  Floating tmin=0.0;
  if(targetPDG == deuteronPDG)   tmin=2*Edisplace_deuteron*M;
  else if(targetPDG == titanPDG) tmin=2*Edisplace_titan*M;
//Calculate integral Rutherford (elastic) cross section:
  Floating cs=0.0;
//For equal particles Rutherford scattering:
//For equal particles we integrate Rutherford cross section from
//|t|min to (|t|max - |t|min).
  if(incPDG == targetPDG)
  {
    if(incPDG == deuteronPDG)
    {
      if(Tls < 30 * keV)
      {//pure Rutherford from |t|min to |t|max-|t|min
        const Floating term1=1.0/tmin - 1.0/(tmax-tmin);
        //beta (beta=pcm/Ecm) of the reduced particle in CM:
        const Floating beta_r_cm = 1.0/std::sqrt(1.0+mr*mr/pcm2);
        const Floating term21 = 2.0/3.0/tmax*beta_r_cm/alpha;
        const Floating term22 = alpha/beta_r_cm * std::log( tmax/tmin-1.0 );
        const Floating term2  = term21*std::sin(term22);
        cs = 2*COEF*(term1+term2)*hc*hc;
      }
      else//Tls>=30*keV
      {//pure Rutherford from |t|min to 10^(-3)*|t|m=10^(-3)*2*pcm^2
        const Floating term1=1.0/tmin - 1.0/kb/tm;
        const Floating term2=1.0/(tmax-kb*tm)-1.0/(tmax-tmin);
        //beta (beta=pcm/Ecm) of the reduced particle in CM:
        const Floating beta_r_cm = 1.0/std::sqrt(1.0+mr*mr/pcm2);
        const Floating term3c = 2.0/3.0/tmax*beta_r_cm/alpha;
        const Floating term31 = alpha/beta_r_cm * std::log( kb*tm/(tmax-kb*tm) );
        const Floating term32 = alpha/beta_r_cm * std::log( tmin/(tmax-tmin) );
        const Floating term3  = term3c*(std::sin(term31-std::sin(term32)));
        cs = COEF*(term1+term2+term3)*hc*hc;
      }
    }
    else
    {
      //calculate term1:
      const Floating term1 = 1.0/tmin - 1.0/(tmax-tmin);
      //beta (beta=pcm/Ecm) of the reduced particle in CM:
      const Floating beta_r_cm = 1.0/std::sqrt(1.0+mr*mr/pcm2);
      //calculate spin:
      const Floating spin=aParticleTable.GetSpin(incPDG);
      //calculate sign of the interference term:
      //% - only for integers, fmod(x, y) - for double.
      const bool cond=(std::fmod(spin,1.0))<1.0e-7;
      //check if the spin is integer (0, 1, 2 ...).
      //else spin is 1/2, 3/2, 5/2 ... and the sign is -1.
      const int SIGN=cond?1:-1;
      //calculate interference term2:
      const Floating coef21 = SIGN * 2.0/(2*spin+1)/tmax*beta_r_cm/alpha;
      const Floating term21 = alpha/beta_r_cm * std::log( tmax/tmin-1.0 );
      const Floating term2  = coef21*std::sin(term21);
      cs = 2*COEF*(term1+term2)*hc*hc;
    }
  }
//For different particles Rutherford scattering:
//For different particles we integrate d(sigma)/d|t| from |t|min to |t|max.
  else
  {
    const Floating term = 1.0/tmin - 1.0/tmax;
    cs = COEF*term*hc*hc;
  }
  return cs;
}

template <typename Floating>
Floating T3ElasticEMIonIonCS<Floating>::GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
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

}
#endif//T3ELASTICEMIONIONCSIMPL_H
