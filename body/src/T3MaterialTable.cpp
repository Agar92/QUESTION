#include "T3Defs.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"
#include <iostream>

namespace t3 {
  MaterialTable::MaterialTable() {
    MaxNumberOfIsotopes=0;
    //number of isotopes in the material:
    NumberOfIsotopes[0]=1;
    NumberOfIsotopes[1]=1;
    NumberOfIsotopes[2]=2;
    NumberOfIsotopes[3]=1;
    NumberOfIsotopes[4]=2;
    //calculate max number of isotopes for T3MultipleScatteringFSImpl.h line 48.
    for(int i=0; i<NumberOfMaterials; ++i)
    {
      if(NumberOfIsotopes[i]>MaxNumberOfIsotopes) MaxNumberOfIsotopes=NumberOfIsotopes[i];
    }
    csBorderSTEP=MaxNumberOfIsotopes+1;
    //elements in materials:
    fMaterials[0][0]=PDG_t{1000010020};//deuterium
    fMaterials[1][0]=PDG_t{1000220480};//titanum
    fMaterials[2][0]=PDG_t{1000010020};//TiD2
    fMaterials[2][1]=PDG_t{1000220480};//TiD2
    fMaterials[3][0]=PDG_t{1000040090};//Be
    fMaterials[4][0]=PDG_t{1000922350};//U
    fMaterials[4][1]=PDG_t{1000922380};//U
    //fractions of elements in materials:
    fFractions[0][0]=1.0;
    fFractions[1][0]=1.0;
    fFractions[2][0]=2.0/3.0;
    fFractions[2][1]=1.0-2.0/3.0;
    fFractions[3][0]=1.0;
    fFractions[4][0]=0.0072;
    fFractions[4][1]=1.0-0.0072;
    //densities of materials:
    fDensities[0]=0.000169 * gcm3;
    fDensities[1]=4.11 * gcm3;
    fDensities[2]=3.91 * gcm3;
    fDensities[3]=1.85 * gcm3;
    fDensities[4]=19.1 * gcm3;
    //concentrations of materials:
    ParticleTable aParticleTable = ParticleTable();
    for(int i=0; i<NumberOfMaterials; ++i)
    {
      auto Atot=0.0;
      for(int j=0; j<NumberOfIsotopes[i]; ++j) Atot += fFractions[i][j] * aParticleTable.GetMass(fMaterials[i][j]);
      fConcentrations[i] = fDensities[i] / Atot;
    }
}
  
} // namespace t3
