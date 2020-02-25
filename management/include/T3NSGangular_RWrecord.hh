#pragma once

#ifndef T3NSGANGULAR_RWRECORD_HH
#define T3NSGANGULAR_RWRECORD_HH

#include <vector>
#include "T3NSGangular_RWnode.hh"

class T3NSGangular_RWrecord: public std::vector<T3NSGangular_RWnode>
{
public:
  T3int Get_levN() const {return levN;}
  void Set_levN(T3int val) {levN = val;}
  std::ofstream& save_binary( std::ofstream& out_stream) const;
  std::ifstream& load_binary( std::ifstream& in_stream);
  T3NSGangular_RWnode interpolate(T3double E) const;
private:
  T3int levN;       /// excited level number as in T2NuclearLevels
};

std::ostream& operator<<(std::ostream& os, const T3NSGangular_RWrecord& inst);

#endif
