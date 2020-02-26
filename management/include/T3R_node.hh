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
//      ------------------ T2R_node header -------------------
//    Created by Mikhail Kossov, Sept 2014 (following T2NSG_node)
//  class header for energy nodes of the (h,anything) reactions in CHIPS.
// --------------------------------------------------------------------------
// Short description: The T2R_node is a part of the CHIPS World. It is
// characterized by the Energy, (h, anything) XS for charged particles inTPT.
// --------------------------------------------------------------------------

#ifndef T3R_node_H
#define T3R_node_H

#include <vector>
#include <fstream>
#include <iostream>
#include "T3Globals.h"
#include "T3Types.h"

class T3R_node
{
public:
  T3R_node(T3double E = 0., T3double XS = 0.);
  void save_binary( std::ofstream* out) const;
  void load_binary( std::ifstream* in); // @@ tmp for old format
  T3bool operator==(const T3R_node rval); // used once
  // Selectors
  T3double E() const {return e;}            // Get the projectile neutron energy
  T3double XS() const {return xs;}          // Get theSummedUp (h,any) reactionCrossSection
  // Modifiers
  void SetE( T3double En )  {e = En;}       // Set the projectile neutron energy
  void SetXS( T3double CS ) {xs = CS;}      // Set the (n,n'g) reaction cross-section
private:
  //body
  T3double e;                // neutron energy
  T3double xs;               // Summed cross-section of all levels + continuum
};

std::ostream& operator<<(std::ostream& os, const T3R_node& inst);

#endif // T3R_node_H
