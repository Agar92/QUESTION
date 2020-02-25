#include "T3R_RW.hh"
#include <map>

#include "unistd.h"

T3R_RW::T3R_RW(const T3TabulatedCS& tcs, T3int _tgZ, T3int _tgA):
  tgZ(_tgZ), tgA(_tgA)
{
  // TODO getenv once at static initialization ----------------------------------------
  const char* const t3data = getenv("T3_DATA");
  if (t3data)
  {
#ifdef debug
    T3cout<<"T3_DATA="<< t3data <<T3endl;
#endif
    T3data =  T3String(t3data);
  }
  else
  {
    T3data = T3String("./data");
    T3cout<<"-Warning-T3R_RW::default_file: Cant't get $T3_DATA. Using ./data" <<T3endl;
  }
#ifdef debug
  T3cout<<"T3_DATA =  "<< T3data << ", str="<< T3data.c_str() << T3endl;
#endif
  //------------------------------------------------------------------------------------
  for(size_t ind = 0; ind < tcs.Get_size(); ++ind)
  {
    push_back(new T3R_node(tcs.Get_E(ind), tcs.Get_xs(ind)));
  }
  std::cout<<"PRINT: "<<std::endl;
  std::cout<<*this<<std::endl;
}

//copy constructor:
T3R_RW::T3R_RW(const T3R_RW & right)
{
  const T3int iD = ND.size();
  //std::cout<<"iD="<<iD<<std::endl;
  //sleep(3);
  if(iD) for(T3int jd=0; jd<iD; ++jd) if(ND[jd]) delete ND[jd];
  ND.resize(0);
  const T3int iDnew=right.size();
  for(int i=0; i<iDnew; ++i)
  {
    const T3double Ei=right.NDI(i)->E();
    const T3double CSi=right.NDI(i)->XS();
    push_back(new T3R_node(Ei, CSi));
  }
  tgZ=right.tgZ;
  tgA=right.tgA;
  T3data=right.T3data;
}

//operator=:
T3R_RW & T3R_RW::operator=(const T3R_RW & right)
{
  const T3int iD = ND.size();
  std::cout<<"iD="<<iD<<std::endl;
  //sleep(3);
  if(iD) for(T3int jd=0; jd<iD; ++jd) if(ND[jd]) delete ND[jd];
  ND.resize(0);
  const T3int iDnew=right.size();
  for(int i=0; i<iDnew; ++i)
  {
    const T3double Ei=right.NDI(i)->E();
    const T3double CSi=right.NDI(i)->XS();
    push_back(new T3R_node(Ei, CSi));
  }
  tgZ=right.tgZ;
  tgA=right.tgA;
  T3data=right.T3data;
  return *this;
}

T3R_RW::~T3R_RW()
{
  T3int iD = ND.size();
  if(iD) for(T3int jd=0; jd<iD; ++jd) if(ND[jd]) delete ND[jd];
}

void T3R_RW::save_binary(const T3String& fname) const // Write R DB to file
{
  std::ofstream out_stream;
  out_stream.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

  std::cout<<"!!!fname="<<fname<<std::endl;
  
  if(!out_stream.good())
  {
#ifdef debug
    T3cout<<"-Warning-T3R_RW::save_binary:*Bad Stream*"<<T3endl;
#endif
    return;
  }
  unsigned int vector_size = ND.size();
#ifdef debug
  T3cout<<"T3R_RW::save_binary: ND="<<vector_size<<T3endl;
#endif
  out_stream.write(  (const char*) &( vector_size ), sizeof( unsigned int ) );
  for( unsigned int i = 0; i < vector_size; ++i)
  {
#ifdef debug
    T3cout<<"T3R_RW::save_binary: node # "<< i <<", E="<< ND[i]->E() <<", XS="
          << ND[i]->XS() << T3endl;
#endif
    ND[i]->save_binary( &out_stream );
  }

  std::cout<<"VSIZE="<<vector_size<<std::endl;
  
}

T3bool T3R_RW::load_binary(const T3String& fname) // Read NNG DB from file
{
  std::ifstream in_stream;
  std::cout<<"$#$fname="<<fname<<std::endl;
  std::cout<<"$#$GGGGGGGG"<<std::endl;
  ///\\\///in_stream.open( fname.c_str(), std::ios::in | std::ios::binary/*| std::ios::trunc*/);
  in_stream.open( fname.c_str(), std::ifstream::in | std::ifstream::binary);
  std::cout<<"$#$HHHHHHHH"<<std::endl;
  if(!in_stream.good())
  {
#ifdef debug
    T3cout<<"-Warning-T3R_RW::load_binary:*Bad Stream*"<<T3endl;
#endif
    return false;
  }  
  unsigned int vector_size;
  // Now read Nodes
  in_stream.read( (char*) &vector_size, sizeof( unsigned int ) );  
#ifdef debug
  T3cout<<"T3R_RW::load_binary: ND="<< vector_size <<", fname=" << fname << T3endl;
#endif
  T3int iD = ND.size();
  if(iD) for(T3int jd=0; jd<iD; ++jd) delete ND[jd];  
  ND.clear();
  T3R_node* newND = 0;
  for( unsigned int i = 0; i < vector_size; ++i)
  {
    newND = new T3R_node();
#ifdef debug
    T3cout<<"T3R_RW::load_binary: Loading node # "<< i <<T3endl;
#endif
    newND->load_binary( &in_stream ); // @@ temporary old format is used
#ifdef debug
    T3cout<<"T3R_RW::load_binary: node # "<< i <<", E="<< newND->E() <<", XS="
          << newND->XS() << T3endl;
    if(newND->XS() < 0.)
    {
      T3cout << "T3R_RW::load_binary: Warning: xs("<<newND->E()<<") =" << newND->XS()
             << " < 0" << T3endl;
    }
#endif
    ND.push_back(newND);
  }
  in_stream.close();
  return true;
}

T3String T3R_RW::default_file( T3int targZ, T3int targA,
                                 T3int pZ, T3int pA, const T3String& suffix ) const
{
  std::map<T3int, T3String> names;
//   names[1000*0 + 0] = "G";
//   names[1000*0 + 1] = "N";
  names[1000*1 + 1] = "P";
  names[1000*1 + 2] = "D";
  names[1000*1 + 3] = "T";
  names[1000*2 + 3] = "H";
  names[1000*2 + 4] = "A";
  T3String particle_name = names[1000 * pZ + pA];
  return default_file(targZ, targA, particle_name + "any", suffix);
}


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
T3String T3R_RW::default_file_Elastic_approx_e_cs(T3int targZ, T3int targA, T3int pZ, T3int pA) const
{
  std::cout<<"LLL"<<std::endl;
  std::map<T3int, T3String> names;
//   names[1000*0 + 0] = "G";
//   names[1000*0 + 1] = "N";
  names[1000*1 + 1] = "P";
  names[1000*1 + 2] = "D";
  names[1000*1 + 3] = "T";
  names[1000*2 + 3] = "H";
  names[1000*2 + 4] = "A";
  char buffer[1023];
  T3String FolderName = "inc" + names[1000 * pZ + pA];
  T3String base="approx";
  T3String elastic="elastic";

  std::cout<<"10"<<std::endl;
  
  sprintf(buffer, "%s/xs/%s/elastic/T3%s_10%.3d%.3d0_%s.bin", T3data.c_str(),
          FolderName.c_str(), FolderName.c_str(), targZ, targA,
          base.c_str());

  std::cout<<"11"<<std::endl;

  //!!!
  //Not more that 69 symbols in filename, otherwise in_stream.open
  //crashes with seg fault on line 120!!!
  //!!!
  T3String filename = T3String(buffer);
  std::cout<<"#@%Check filename:"<<std::endl;
  std::cout<<"#@%filename="<<filename<<" size="<<filename.size()<<std::endl;

  std::cout<<"12"<<std::endl;
  
  return filename;  
}

void T3R_RW::save_binary_Elastic_approx_e_cs(T3int targZ, T3int targA, T3int pZ, T3int pA) const
{
  save_binary( default_file_Elastic_approx_e_cs(targZ, targA, pZ, pA) );
}

T3bool T3R_RW::load_binary_Elastic_approx_e_cs(T3int targZ, T3int targA, T3int pZ, T3int pA)
{
  return load_binary( default_file_Elastic_approx_e_cs(targZ, targA, pZ, pA) );
}
//---------------------------------------------------------------------//  


void T3R_RW::save_binary( T3int targZ, T3int targA,
                          T3int pZ, T3int pA, const T3String& suffix ) const
{
  save_binary( default_file( targZ, targA, pZ, pA, suffix ) );
}

void T3R_RW::save_binary( T3int targZ, T3int targA,
                          const T3String& sID, const T3String& suffix ) const
{
  save_binary( default_file( targZ, targA, sID, suffix ) );
}

T3bool T3R_RW::load_binary(T3int targZ, T3int targA,
                           T3int pZ, T3int pA, const T3String& suffix)
{
  return load_binary( default_file( targZ, targA, pZ, pA, suffix ) );
}

T3bool T3R_RW::load_binary(T3int targZ, T3int targA,
                             const T3String& rID, const T3String& suffix)
{
  return load_binary( default_file( targZ, targA, rID, suffix ) );
}

T3String T3R_RW::default_file( T3int targZ, T3int targA,
                                 const T3String& sID, const T3String& suffix ) const
{
  const T3int nbases = 4;
  T3String bases[nbases] = {"endf", "tendl", "custom", "unknown"};
  T3bool CSExist = false;
  char buffer[1023];
  T3String filename;
  T3String base;// = bases[0];
  if (suffix == "any")
  {
    for( T3int index = 0; index < nbases; ++index)
    {
      base = bases[index];
      sprintf(buffer, "%s/xs/%s/T3%s_10%.3d%.3d0_%s.bin", T3data.c_str(),
              sID.c_str(), sID.c_str(), targZ, targA, base.c_str());
      filename = T3String(buffer);
      std::ifstream try_file(filename);
      if ( !try_file.good() )
      {
#ifdef debug
        T3cout<<"*WARNING*T3R_RW::CCS:absentCSFile="<<filename<<T3endl;
#endif
        CSExist = false;
        try_file.close();
      }
      else
      {
        CSExist = true;
        break;
      }
    }
    if (CSExist) return filename;
    else
    {
      // @@ TODO: find the file with any suffix
      sprintf(buffer, "%s/xs/%s/T3%s_10%.3d%.3d0_%s.bin", T3data.c_str(),
              sID.c_str(), sID.c_str(), targZ, targA, "<db>");
      filename = T3String(buffer);
#ifdef debug
      T3cout<<"*WARNING*T3R_RW::CCS:(SYNC) no CSFile, return "<<filename<<T3endl;
#endif
      return filename; // TODO return special name or emit warning
    }
  }
  else
  {
    base = suffix;
    sprintf(buffer, "%s/xs/%s/T3%s_10%.3d%.3d0_%s.bin", T3data.c_str(),
            sID.c_str(), sID.c_str(), targZ, targA, base.c_str());
    filename = T3String(buffer);
    std::ifstream try_file(filename);
#ifdef debug
    if(!try_file.good())T3cout<<"T3R_RW::CCS:absentNecessaryCSFile="<<filename<<T3endl;
#endif
    try_file.close();
    return filename;
  }
}

// void T2R_RW::save_binary(G4int targZ, G4int targA, G4int pZ, G4int pA, G4int sZ, G4int sA,
//                          const G4String& suffix) const
// {
//   save_binary( default_file( targZ, targA, pZ, pA, sZ, sA, suffix ) );
// }

// G4bool T2R_RW::load_binary(G4int targZ, G4int targA, G4int pZ, G4int pA, G4int sZ,
//                            G4int sA, const G4String& suffix)
// {
//   return load_binary( default_file( targZ, targA, pZ, pA, sZ, sA, suffix ) );
// }

std::ostream& operator<<(std::ostream& os, const T3R_RW& inst)
{
  for(size_t ind = 0; ind < inst.size(); ++ ind)
  {
    os << "node #" << ind << ' ' << *(inst.NDI(ind)) << T3endl;
  }
  return os;
}
