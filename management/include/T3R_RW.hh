#pragma once

//
// ****************************************************************
// * TPT License and Disclaimer                                   *
// *                                                              *
// * The TPT Software  is  copyright  of the Copyright Holders of *
// * the TPT CFAR-VNIIA group. It is provided under the terms and *
// * conditions of  the Software License  (ROSPATENT 2014611928). *
// *                                                              *
// * This code is  NOT  an open code. It is distributed in a form *
// * of compiled libraries. The code implementation is the result *
// * of the  scientific and technical work of the CFAR-VNIIA  TPT *
// * Scientific Group. By using or distributing  the TPT software *
// * (or any work based on the software) you agree to acknowledge *
// * its use  in  resulting scientific publications, and indicate *
// * your acceptance of all terms of the TPT Software License.    *
// ****************************************************************
//
//      ---------------- T2R_RW --------------------
//   Created by Mikhail Kosov, June 2014 from T2NSG_RW.
//  Class header for writing/reading of (h,anything) ENDF6 DB one isotope file
// ---------------------------------------------------------------------------
// Short description: Container for an isotope (n,h'g) Q-levels & E-nodes
// ----------------------------------------------------------------------
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

//---------------------------------------------------------------------//  
//!!!  
//This function saves D-D elastic scattering approximation
//integral cross sections. Yet we chose 1-cos(theta_cm)=1.0e-3
//to be the left border of (E,partial sums) Ox axis. To the
//left from it, there is no difference between red
//approximation and green Rutherford curve.
//The right border is 1-cos(theta_cm)=1.0.
//We sum partial sums from 1-cos(theta_cm)=1.0 to
//1-cos(theta_cm)=1.0e-3 in a logarithmic scale.
//The cross sections we write using this method are
//The partial sum at 1-cos(theta_cm)=1.0e-3 values,
//!!!So, they are not full elastic D-D integral cross sections.!!!  
//!!!
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

  T3bool load_binary( const T3String& fname );      // Read binary table from file
  
private:
  ///T3bool load_binary( const T3String& fname );      // Read binary table from file

  // Body
  std::vector<T3R_node*> ND;                      // Vector of pointers to (n,h'g) nodes
  T3double tgZ;     // FIXME currently not saved in files (it is in the name of file)
  T3double tgA;     // maybe add to the end for compability ? @@
  T3String T3data;  // path to the data directory
};

std::ostream& operator<<(std::ostream& os, const T3R_RW& inst);

#endif // T3R_RW_HH
