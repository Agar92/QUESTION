#include "T3NSGangular_RW.hh"
#include <map>

#include "T3ParticleTable.h"
#include "unistd.h"

//#define debug

std::ofstream& T3NSGangular_RW::save_binary(std::ofstream& out_stream) const
{
  if( out_stream.good() )
  {
    const unsigned int vector_size = size();
    out_stream.write((const char*) &Z, sizeof(T3int) );
    out_stream.write((const char*) &A, sizeof(T3int) );
    out_stream.write((const char*) &secPDG, sizeof(T3int) );
    out_stream.write((const char*) &rID, sizeof(T3int) );
    out_stream.write((const char*) &vector_size, sizeof(unsigned int) );
    for(T3NSGangular_RW::const_iterator it = begin(); it != end(); ++it)
    {
      it->save_binary(out_stream);
    }
  }
  else
  {
    T3cout<<"-Warning-T3NSGangular_RW::save_binary:*Bad stream*"<<T3endl;
  }
  return out_stream;
}

std::ifstream& T3NSGangular_RW::load_binary(std::ifstream& in_stream )
{
  if( in_stream.good() )
  {
    unsigned int vector_size = 0;
    in_stream.read((char*) &Z, sizeof(T3int) );
    in_stream.read((char*) &A, sizeof(T3int) );
    in_stream.read((char*) &secPDG, sizeof(T3int) );
    in_stream.read((char*) &rID, sizeof(T3int) );
    in_stream.read((char*) &vector_size, sizeof(unsigned int) );
    resize(vector_size);
    for(T3NSGangular_RW::iterator it = begin(); it != end(); ++it)
    {
      it->load_binary(in_stream);
    }
  }
  else
  {
    T3cout<<"-Warning-T3NSGangular_RW::load_binary: *Bad stream*"<<T3endl;
  }

#ifdef debug
  T3cout << *this << T3endl;
#endif
  return in_stream;
}

void T3NSGangular_RW::save_binary( const T3String& fname ) const
{
  std::ofstream out_stream;
  out_stream.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  if(out_stream.good())
  {
    save_binary(out_stream);
  }
  else
  {
#ifdef debug
    T3cout << "-Warning-T3NSGangular_RW::failed to open " << fname << T3endl;
#endif
  }
  out_stream.close();
}

T3bool T3NSGangular_RW::load_binary( const T3String& fname )
{
  T3bool result = false;
  std::ifstream in_stream;
  in_stream.open(fname.c_str(), std::ifstream::in | std::ifstream::binary);
  ///\\\///in_stream.open(fname.c_str(), std::ios::in | std::ios::binary);
  if(in_stream.good())
  {
    load_binary(in_stream);
    result = size();
  }
  else
  {
#ifdef debug
    T3cout << "-Warning-T3NSGangular_RW::failed to open " << fname << T3endl;
#endif
  }
  in_stream.close();

  /*
  std::cout<<"!!!!!!Check filename PS:"<<std::endl;
  std::cout<<"fname="<<fname<<" size="<<fname.size()
           <<"!!!!!!!"<<std::endl;

  std::cout<<"HERE WE CHECK!"<<std::endl;
  sleep(3);
  */
  
  return result;
}

T3String T3NSGangular_RW::default_file(T3int tgZ, T3int tgA,
                                       const T3String& rid, T3int pdg, T3int incZA) const
{
  std::map<T3int, T3String> pnames{{1, "N"},
                                   {1001, "P"},
                                   {1002, "D"},
                                   {1003, "T"},
                                   {2003, "H"},
                                   {2004, "A"}};
//  pnames[2112] = "n";
  char buffer[1023];
  T3String T3data;                // Check that the NNG data exist
  if (getenv("T3_DATA"))
  {
    T3data =  getenv("T3_DATA");
  }
  else
  {
    T3data = T3String("./data");
    T3cout << "-Warning-T3NSGangular_RW::default_file:Cant't get $T3_DATA. Use ./data"
           << T3endl;
  }
  ///\\\///const T3String pname = pnames[pdg];

//Andrey added:  
  t3::ParticleTable aParticleTable;
  T3int secondaryZ=aParticleTable.GetZ(pdg);
  T3int secondaryA=aParticleTable.GetA(pdg);
  const T3int secondaryZA=1000*secondaryZ+secondaryA;
  const T3String pname = pnames[secondaryZA];
//End of Andrey added.
  
  T3String incFolder = (incZA == 1) ? "n" : ("inc" + pnames[incZA]);
  sprintf(buffer, "%s/angular/%s/%s/%s/T3%sSGangular_10%.3d%.3d0.bin", T3data.c_str(),
          pname.c_str(), incFolder.c_str(), rid.c_str(), pnames[incZA].c_str(), tgZ, tgA);
#ifdef debug
  T3cout << "T3NSGangular_RW::default_file: " << buffer << T3endl;
#endif

  return T3String(buffer);
}

T3bool T3NSGangular_RW::load_binary(T3int tgZ, T3int tgA, const T3String& rid, T3int pdg,
                                    T3int incZA)
{
  return load_binary(default_file( tgZ, tgA, rid, pdg, incZA ));
}

void T3NSGangular_RW::save_binary( T3int tgZ, T3int tgA, const T3String& rid, T3int pdg,
                                   T3int incZA) const
{
  save_binary(default_file( tgZ, tgA, rid, pdg, incZA ));
}

std::ostream& operator<<(std::ostream& os, const T3NSGangular_RW& inst)
{
  os << "Z=" << inst.Get_Z() << " A=" << inst.Get_A() << " secPDG = " << inst.Get_secPDG()
     << "rID = " << inst.Get_rID() << ", records:";
  for(size_t ind =0; ind < inst.size(); ++ind) os << "\n[" << inst.at(ind) << "]";
  return os;
}
