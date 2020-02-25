#pragma once
#ifndef T3ELASTICSTRONGIONION_DB_HH
#define T3ELASTICSTRONGIONION_DB_HH

#include <iostream>
#include "T3NSGangular_RWnode.hh"

class T3ElasticStrongIonIon_DB
{
public:
  T3ElasticStrongIonIon_DB(T3double xmin=-1.0, T3double xmax=1.0);
  T3ElasticStrongIonIon_DB(const T3ElasticStrongIonIon_DB & nodedb);
  T3ElasticStrongIonIon_DB & operator=(const T3ElasticStrongIonIon_DB & nodedb);
  const size_t _num_point = 127;
  const size_t _num_nodes = 512;

#ifdef OPENACC  
  #pragma acc routine seq
#endif
  T3double RandomizeCost_node(unsigned int & generator_seed, size_t node_num) const;

#ifdef OPENACC
  #pragma acc routine seq
#endif
  T3double RandomizeCost(unsigned int & generator_seed, T3double E) const;

  void Set_V(size_t node_num, size_t point_num, T3double val) {V[node_num][point_num]=val;}
  void Set_a(size_t node_num, size_t point_num, T3double aval) {a[node_num][point_num]=aval;}
  void Set_b(size_t node_num, size_t point_num, T3double bval) {b[node_num][point_num]=bval;}
  void Set_c(size_t node_num, size_t point_num, T3double cval) {c[node_num][point_num]=cval;}
  void Set_Einc(size_t node_num, T3double Eincval) {Einc[node_num]=Eincval;}
  void SetXmin(T3double _Xmin) {Xmin=_Xmin;}
  void SetXmax(T3double _Xmax) {Xmax=_Xmax;}  
  
  T3double Get_V(size_t node_num, size_t point_num) const {return V[node_num][point_num];}
  T3double Get_a(size_t node_num, size_t point_num) const {return a[node_num][point_num];}
  T3double Get_b(size_t node_num, size_t point_num) const {return b[node_num][point_num];}
  T3double Get_c(size_t node_num, size_t point_num) const {return c[node_num][point_num];}
  T3double Get_Einc(size_t node_num) const {return Einc[node_num];}
  T3double GetXmin() const {return Xmin;}
  T3double GetXmax() const {return Xmax;}
  
private:
  T3double V[512][127];
  T3double a[512][127];
  T3double b[512][127];
  T3double c[512][127];
  T3double Einc[512];
  const T3double CosThetaCMmin=1.0e-3;
  const T3double CosThetaCMmax=1.0;   
  T3double Xmin;
  T3double Xmax;
  T3double deltaLnCosThetaCM;
  T3double lncoscm[129];
};

#endif
