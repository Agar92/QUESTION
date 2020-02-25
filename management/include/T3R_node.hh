#pragma once

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
