#pragma once
#ifndef T3R_RW_HH
#define T3R_RW_HH

#include <string>
#include <ostream>
#include "T3R_node.hh"
#include "T3TabulatedCS.hh"

class T3R_RW
{
public:
  T3R_RW(const T3TabulatedCS& tcs = T3TabulatedCS(), T3int _tgZ = 0, T3int _tgA = 0);
  //copy constructor:
  T3R_RW(const T3R_RW & right);
  //operator=:
  T3R_RW & operator=(const T3R_RW & right);
  ~T3R_RW();

  T3String default_file_Elastic_approx_e_cs( T3int targZ, T3int targA, T3int pZ, T3int pA) const;

  void save_binary_Elastic_approx_e_cs(T3int targZ, T3int targA, T3int pZ, T3int pA) const;

  T3bool load_binary_Elastic_approx_e_cs(T3int targZ, T3int targA, T3int pZ, T3int pA);
  
//---------------------------------------------------------------------//    

  
  void save_binary(T3int targZ, T3int targA, T3int pZ, T3int pA,
                   const T3String& suffix="unknown") const;
//   void save_binary(G4int targZ, G4int targA, G4int pZ, G4int pA, G4int sZ, G4int sA,
//                    const G4String& suffix="unknown") const;
  void save_binary( T3int targZ, T3int targA, const T3String& sID,
                    const T3String& suffix = "unknown") const;
  void save_binary( const T3String& fname ) const; // @@ Save binary table to file (?)
  T3bool load_binary(T3int targZ, T3int targA,
                     T3int pZ, T3int pA, const T3String& suffix="any");
//   G4bool load_binary(G4int targZ, G4int targA, G4int pZ, G4int pA, G4int sZ, G4int sA,
//                      const G4String& suffix="any");
  T3bool load_binary(T3int targZ, T3int targA, const T3String& rID,
                     const T3String& suffix = "any");

  T3String default_file( T3int targZ, T3int targA, T3int pZ, T3int pA,
                         const T3String& suffix) const;
  T3String default_file( T3int targZ, T3int targA, T3int pZ, T3int pA,
                         T3int sZ, T3int sA, const T3String& suffix) const;
  T3String default_file( T3int targZ, T3int targA, const T3String& sID,
                         const T3String& suffix) const;
  void push_back(T3R_node* N) {ND.push_back(N);}    // add in the end a new R_Node
  size_t size() const       {return ND.size();}     // Get size of the NHG-nodes vector
  const T3R_node* NDI(T3int i) const {return ND.at(i);} // Get i-th NHG-node Pointer
  T3double GetZ() const  {return tgZ;}              // Extract the A-parameter for CS=A+B/p
  T3double GetA() const  {return tgA;}              // Extract the B-parameter for CS=A+B/p
  void SetZ(T3double Z)  {tgZ = Z;}                 // Put the A-parameter for CS=A+B/p
  void SetA(T3double A)  {tgA = A;}                 // Put the B-parameter for CS=A+B/p

private:
  T3bool load_binary( const T3String& fname );      // Read binary table from file

  // Body
  std::vector<T3R_node*> ND;                      // Vector of pointers to (n,h'g) nodes
  T3double tgZ;     // FIXME currently not saved in files (it is in the name of file)
  T3double tgA;     // maybe add to the end for compability ? @@
  T3String T3data;  // path to the data directory
};

std::ostream& operator<<(std::ostream& os, const T3R_RW& inst);

#endif // T3R_RW_HH
