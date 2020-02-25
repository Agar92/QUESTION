#pragma once
#ifndef T3MATERIALTABLE_H
#define T3MATERIALTABLE_H

#include <utility>
#include "T3Defs.h"

namespace t3 {
//WE HAVE 5 MATERIALS NOW:
constexpr int NumberOfMaterials=5;

class MaterialTable {
public:
  MaterialTable();
  auto const GetIsotopes(MatID_t matID, int i) const {return fMaterials[matID][i];}
  auto GetFractions(MatID_t matID, int i) const {return fFractions[matID][i];}
  auto GetConcentrations(MatID_t matID) const {return fConcentrations[matID];}
  auto GetDensity(MatID_t matID) const {return fDensities[matID];}
  auto GetNumberOfIsotopes(MatID_t matID) const {return NumberOfIsotopes[matID];}
  auto GetNMat() const {return NumberOfMaterials;}
  int GetMaxNumberOfIsotopes() const {return MaxNumberOfIsotopes;}
  int GetcsBorderSTEP() const {return csBorderSTEP;}
private:
  int MaxNumberOfIsotopes;
  int csBorderSTEP;
  int NumberOfIsotopes[NumberOfMaterials];
  PDG_t fMaterials[NumberOfMaterials][2];//maximum 2 isotopes in TiD2 and U mix.
  double fFractions[NumberOfMaterials][2];
  double fConcentrations[NumberOfMaterials];
  double fDensities[NumberOfMaterials];
};
} // namespace t3

#endif // T3MATERIALTABLE_H
