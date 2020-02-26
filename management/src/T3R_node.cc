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
// *** By opening this file you break the TPT License Agreement ***
//
//      ---------------- T2R_node class ----------------
//  Created by Mikhail Kossov, Sept 2012 (following T2NSG_node)
//  Class header for energy nodes of the (h, anything) reactions in TPT.
// ----------------------------------------------------------------------------------------
//       1         2         3         4         5         6         7         8         9
//34567890123456789012345678901234567890123456789012345678901234567890123456789012345678901
// ----------------------------------------------------------------------------------------
// Short description: The T2R_node is a part of the CHIPS World. It is
// characterized by the Energy, (h, anything) XS for charged particles in TPT.
// ---------------------------------------------------------------------------

//#define debug
//#define pdebug

#include "T3R_node.hh"
#include <cstdio>
#include <fstream>
#include <iostream>

T3R_node::T3R_node(T3double Ein, T3double CS) // Ein & XS = 0
{
  e  = Ein;
  xs = CS;
}

void T3R_node::save_binary(std::ofstream* out_stream ) const // New format!
{
  if( !out_stream->good() )
  {
    T3cout<<"-Warning-T3R_node::save_binary:*Bad stream*"<<T3endl;
    return;
  }
  out_stream->write((const char*) &e, sizeof(T3double) );
  out_stream->write((const char*) &xs, sizeof(T3double) );
}

void T3R_node::load_binary(std::ifstream* in_stream ) // @@ old format
{
  if( !in_stream->good() )
  {
    T3cout<<"-Warning-T3R_node::load_binary: *Bad stream*"<<T3endl;
    return;
  }
  in_stream->read((char*) & e, sizeof(T3double) );
  in_stream->read((char*) & xs, sizeof(T3double) );
#ifdef debug
  T3cout<<"T3R_node::load_binary: E="<< e <<", XS="<< xs << T3endl;
#endif
}
std::ostream& operator<<(std::ostream& os, const T3R_node& inst)
{
  os << "E = " << inst.E() << ", xs = " << inst.XS() << T3endl;
  return os;
}


