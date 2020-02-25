#include "T3ParticleTable.h"

namespace t3 {
ParticleTable::ParticleTable() {
  fMasses[PDG_t(22)]   = 0.0;         //photon
  fMasses[PDG_t(11)]   = electronMass;//electron
  fMasses[PDG_t(-11)]  = electronMass;//positron
  fMasses[PDG_t(2212)] = protonMass;  //proton
  fMasses[PDG_t(2112)] = neutronMass; //neutron
}
} // namespace t3
