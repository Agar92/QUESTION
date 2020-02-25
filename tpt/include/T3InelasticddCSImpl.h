#pragma once
#ifndef T3INELASTICDDCSIMPL_H
#define T3INELASTICDDCSIMPL_H

#include <cmath>
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

namespace t3 {

using namespace units;  
template <typename Floating>
class InelasticddCS {
public:
  InelasticddCS()
  {
    ParticleTable aParticleTable;
    deuteronPDG=aParticleTable.makePDGfromZandA(1,2);
  }

  inline Floating
  GetDDInelasticIntegralCS(Floating Tls) const;

  inline Floating
  GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
        ParticleTable & aParticleTable, MaterialTable & aMaterialTable) const;

private:
  PDG_t deuteronPDG;
  const Floating Eg = 0.986 * MeV;
};

template <typename Floating>
Floating InelasticddCS<Floating>::GetDDInelasticIntegralCS(Floating Tls) const
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
  const Floating cs = csnhe3 + cspt;
  return cs;
}

template <typename Floating>
Floating InelasticddCS<Floating>::GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
                                        ParticleTable & aParticleTable,
                                        MaterialTable & aMaterialTable) const
{
  if(incPDG != deuteronPDG) return 0.;
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

}
#endif//T3INELASTICDDCSIMPL_H
