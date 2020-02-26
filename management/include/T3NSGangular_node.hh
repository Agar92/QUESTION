#pragma once

#ifndef T3NSGANGULAR_NODE_HH
#define T3NSGANGULAR_NODE_HH

#include "T3NSGangular_RWnode.hh"

//#include "T3RNG.h"
#include <cstring>



namespace t3
{

#ifdef OPENACC
#pragma acc routine seq
#endif
  unsigned int Rand32(unsigned int & xn);

#ifdef OPENACC
#pragma acc routine seq
#endif
  double rndv(unsigned int xn);

#ifdef OPENACC
#pragma acc routine seq
#endif
  double RND01(unsigned int & xn);
  
}//end of namespace t3.



class T3NSGangular_node
{
public:
  T3NSGangular_node(const T3NSGangular_RWnode& rwnode = T3NSGangular_RWnode(), T3double xmin=-1.0, T3double xmax=1.0);
  T3NSGangular_node(T3NSGangular_node const &node);
  static const size_t _num_point = T3NSGangular_RWnode::_num_point;
  T3NSGangular_node& operator=(T3NSGangular_node const &node);
  T3double RandomizeCost(unsigned int & generator_seed);
  T3double Get_V(size_t point_num) const {return V[point_num];}
  T3double Get_a(size_t point_num) const {return a[point_num];}
  T3double Get_b(size_t point_num) const {return b[point_num];}
  T3double Get_c(size_t point_num) const {return c[point_num];}
  T3double Get_pr() const {return pr;}
  T3double Get_sl() const {return sl;}
  void SetXmin(T3double _Xmin){Xmin=_Xmin;}
  void SetXmax(T3double _Xmax){Xmax=_Xmax;}
private:
  T3double V[_num_point];
  T3double a[_num_point];
  T3double b[_num_point];
  T3double c[_num_point];
  T3double pr;
  T3double sl;
  //it is necessary for using this class for different partial sums
  //with different xmin and xmax.
  //SDI had, as i understand, xmin=-1, xmax=1 ( cos(theta_cm) from -1 to 1 ).
  //We with our elastic scattering approximation partial sums have now:
  //ln(1-cos(theta_cm))_min=10^(-3) and
  //ln(1-cos(theta_cm))_max=1.
  //So, to use TPT2 quadratic interpolation, we must have a possibility
  //to set different dY, dCosT, ct and fs.
  //For this purpose, we use the following variables:
  T3double Xmin;
  T3double Xmax;
};

std::ostream& operator<<(std::ostream& os, const T3NSGangular_node& inst);

#endif
