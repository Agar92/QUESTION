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


