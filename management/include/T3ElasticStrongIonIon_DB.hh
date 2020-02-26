#pragma once
#ifndef T3ELASTICSTRONGIONION_DB_HH
#define T3ELASTICSTRONGIONION_DB_HH

#include <iostream>
#include "T3NSGangular_RWnode.hh"

//The scale of the X axis 1-cos( theta_cm ) is logarithmic.
//This means that the scale of the lnX - ln(1-cos( theta_cm )) is linear.
//And the functions of T2 for quadratic interpolation, which used linear scale,
//may be used.!!!
class T3ElasticStrongIonIon_DB
{
public:
  T3ElasticStrongIonIon_DB(T3double xmin=-1.0, T3double xmax=1.0);
  T3ElasticStrongIonIon_DB(const T3ElasticStrongIonIon_DB & nodedb);
  T3ElasticStrongIonIon_DB & operator=(const T3ElasticStrongIonIon_DB & nodedb);
  //number of points in partial sums:
  const size_t _num_point = 127;
  //number of (energy, partial sums) data sets:
  const size_t _num_nodes = 512;

//returns random ln(1-cos(theta_cm))
#ifdef OPENACC  
  #pragma acc routine seq
#endif
  T3double RandomizeCost_node(unsigned int & generator_seed, size_t node_num) const;

//returns random ln(1-cos(theta_cm))  
#ifdef OPENACC
  #pragma acc routine seq
#endif
  T3double RandomizeCost(unsigned int & generator_seed, T3double E) const;

//!!This function is only for CPU, because C array is allocated in it.!!!
/*
//Linear interpolation. This is for tests. Returns random ln(1-cos(theta_cm)).
#ifdef OPENACC
#pragma acc routine seq
#endif
  T3double RandomizeCost_node1(unsigned int & generator_seed, size_t node_num) const;
*/
  
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
  //it is necessary for using this class for different partial sums
  //with different xmin and xmax.
  //TPT2 had, as i understand, xmin=-1, xmax=1 ( cos(theta_cm) from -1 to 1 ).
  //We with our elastic scattering approximation partial sums have now:
  //ln(1-cos(theta_cm))_min=10^(-3) and
  //ln(1-cos(theta_cm))_max=1.
  //So, to use TPT2 quadratic interpolation, we must have a possibility
  //to set different dY, dCosT, ct and fs.
  //For this purpose, we use the following variables:
  //ln( 1-cos(theta_cm) ) axis. For D-D it is from 1-10^(-3) to 1.0 :
  const T3double CosThetaCMmin=1.0e-3;//left bound of 1-cos(theta_cm), which we chose.
  const T3double CosThetaCMmax=1.0;   //corresponds to theta_cm=pi/2 and cos( theta_cm )=0.0 => 1-cos(theta_cm)=1.0.
  T3double Xmin;
  T3double Xmax;
  T3double deltaLnCosThetaCM;
  T3double lncoscm[129];
};

#endif
