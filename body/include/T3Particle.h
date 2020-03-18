#pragma once
#ifndef T3PARTICLE_H
#define T3PARTICLE_H

#include <cmath>
#include "T3Globals.h"
#include "T3NSGangular_node.hh"
#include "T3ParticleTable.h"
#include "T3LorentzVector.h"

namespace t3 {

template <typename Floating>
struct Particle {

Particle()
      : Particle<Floating>(t3::T3LorentzVector<Floating>(0.0,0.0,0.0,0.0),
                           t3::T3LorentzVector<Floating>(1.0/sqrt(3.0),1.0/sqrt(2.0),1.0/sqrt(6.0),G),
                           1.0*MeV, 0.0, t3::PDG_t(2112), 1.0, 0, 1, 1, -1.0/*tr*/) {}


  Particle(T3LorentzVector<Floating> _r, T3LorentzVector<Floating> _p, Floating _m, Floating _de, t3::PDG_t _pdg,
           Floating _wt, unsigned int _rs, int _ir, int _id, Floating _tr);
T3LorentzVector<Floating> r;
T3LorentzVector<Floating> p;
Floating m;
Floating de;
Floating wt;
unsigned int rs;
PDG_t pdg;
int ir;
int id;
void SetEtot(Floating Etot, ParticleTable & aParticleTable)
{
  //if(Etot<m) std::cout<<"***Error: Etot="<<Etot<<" < m="<<m<<std::endl;
  const auto plsnew = std::sqrt(Etot*Etot-m*m);
  SetPxPyPzE(plsnew*vx(), plsnew*vy(), plsnew*vz(), Etot);
}

Floating tr;
void SetPxPyPzE(Floating px, Floating py, Floating pz, Floating E){p.SetPxPyPzE(px, py, pz, E);}
void SetRxRyRzT(Floating rx, Floating ry, Floating rz, Floating t){r.SetPxPyPzE(rx, ry, rz, t);}
Floating mass() const {return m;}
Floating GetdE() const { return de; }
Floating GetEtot() const { return p.E(); }
Floating R() const { return std::sqrt(r.x()*r.x()+r.y()*r.y()+r.z()*r.z());}
Floating P() const { return std::sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());}
Floating vx() const { return p.x() / P(); }
Floating vy() const { return p.y() / P(); }
Floating vz() const { return p.z() / P(); }
int ix() const { return static_cast<int>(std::floor(r.x() / ag)); }
int jy() const { return static_cast<int>(std::floor(r.y() / ag)); }
int kz() const
{
  //std::cout<<"z="<<r.z()<<" ag="<<ag<<std::endl;
  return static_cast<int>(std::floor(r.z() / ag));
}
auto GenerateCanonical() {
   return RND01(rs);
}
};

template <typename Floating>
Particle<Floating>::Particle(T3LorentzVector<Floating> _r, T3LorentzVector<Floating> _p, Floating _m, Floating _de,
                             t3::PDG_t _pdg, Floating _wt, unsigned int _rs, int _ir, int _id, Floating _tr)
  : r{_r}, p{_p}, m(_m), de(_de), pdg(_pdg), wt(_wt), rs(_rs), ir(_ir), id(_id), tr(_tr){}


} // namespace t3

#endif // T3PARTICLE_H
