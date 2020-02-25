#pragma once

#ifndef T3NSGANGULAR_RW_HH
#define T3NSGANGULAR_RW_HH

#include "T3NSGangular_RWrecord.hh"
#include <array>

class T3NSGangular_RW: public std::vector<T3NSGangular_RWrecord>
{
public:
  T3int Get_secPDG() const {return secPDG;}
  void Set_secPDG(T3int val) {secPDG = val;}
  T3int Get_rID() const {return rID;}
  void Set_rID(T3int val) {rID = val;}
  T3int Get_Z() const {return Z;}
  void Set_Z(T3int val) {Z = val;}
  T3int Get_A() const {return A;}
  void Set_A(T3int val) {A = val;}
  void save_binary(T3int tgZ, T3int tgA, const T3String& rid, T3int pdg = 2112,
                   T3int incZA = 1) const;
  T3bool load_binary(T3int tgZ, T3int tgA, const T3String& rid, T3int pdg = 2112,
                     T3int incZA = 1);

  void save_binary( const T3String& fname ) const;      // Save Binary Tables to the file
  T3bool load_binary( const T3String& fname );            // Read binary table from file
  
private:
  T3int Z;
  T3int A;
  T3int secPDG;     /// secondary particle PDGcode
  T3int rID; /// reaction identification
  std::ofstream& save_binary( std::ofstream& out_stream) const;
  std::ifstream& load_binary( std::ifstream& in_stream);
  /////void save_binary( const T3String& fname ) const;      // Save Binary Tables to the file
  /////T3bool load_binary( const T3String& fname );            // Read binary table from file
  T3String default_file( T3int tgZ, T3int tgA,
                         const T3String& rid, T3int pdg,
                         T3int incZA = 1) const;     // Name definition
//std::vector<size_t> index; // map level numbers to record numbers in the vector  
};

std::ostream& operator<<(std::ostream& os, const T3NSGangular_RW& inst);

#endif
