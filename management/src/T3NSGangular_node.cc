#include "T3NSGangular_node.hh"
#include "T3Utility.hh"
#include "unistd.h"

namespace t3
{

#ifdef OPENACC
#pragma acc routine seq
#endif
unsigned int Rand32(unsigned int & xn)
{
  u_quad_t a=0x5DEECE66D;
  u_quad_t c=0xB;
  return (unsigned int)((a*xn+c) & 0xFFFFFFFF);
}

#ifdef OPENACC
#pragma acc routine seq
#endif
double rndv(unsigned int xn)
{
  return (double) xn / (double) 0x100000000LL;
}

#ifdef OPENACC
#pragma acc routine seq
#endif
double RND01(unsigned int & xn)
{
  xn=Rand32(xn);
  return (double)(xn) / (double) 0x100000000LL;
}
  
}//end of namespace t3.

using namespace t3;

T3NSGangular_node::T3NSGangular_node(const T3NSGangular_RWnode& rwnode,
                                     T3double xmin, T3double xmax):
                                     pr(rwnode.Get_pr()), sl(rwnode.Get_sl()),
                                     Xmin(xmin), Xmax(xmax)
{
  std::vector<T3LorentzVector<FloatingType>> vlv = rwnode.CalcXYZ(Xmin,Xmax);

  //std::cout<<"T3NSGangular_node::T3NSGangular_node():"<<std::endl;
  //std::cout<<"Xmin="<<Xmin<<" Xmax="<<Xmax<<std::endl;
  //sleep(1);
  
  for(size_t ind =0; ind < _num_point; ++ind)
  {
    V[ind] = vlv.at(ind).t();
    a[ind] = vlv.at(ind).x();
    b[ind] = vlv.at(ind).y();
    c[ind] = vlv.at(ind).z();
  }
}

T3double T3NSGangular_node::RandomizeCost(unsigned int & generator_seed)
{
  T3double P = RND01(generator_seed);
  if(P < pr)
  {
    return 1. + sl * log(1. - ( 1. - exp(-2 / sl) ) *
                         RND01(generator_seed) );
  }
  static const T3double dY = 2. / (_num_point + 1);
  P = RND01(generator_seed);
  const size_t hct = T3Utility::bin_search<_num_point>(V, P) - V;
  T3int const ct1 = hct ? hct-1 : hct;
  T3int const ct2 = (hct && hct < _num_point) ? hct : -1;
  T3double const r1 = hct ? V[ct1] : 0.;
  T3double const r2 = (hct < _num_point) ? V[hct] : 1.;
  long double cost = a[ct1] + P*( b[ct1] + P * c[ct1] );
  const long double dd = dY / ( r2 - r1);
  if(ct2 >= 0)
  {
    const T3double dr1 = r1 + r1;
    const T3double dr2 = r2 + r2;
    const long double p12 = b[ct1] + dr2 * c[ct1];
    const long double p21 = b[ct2] + dr1 * c[ct2];
    if(p12 < 0. || p21 < 0.)
    {
      long double F  =  a[ct2] + P*( b[ct2] + P * c[ct2] );
      const long double fs = dY * hct - 1.;
      const long double d1 = P - r1;
      const long double d2 = r2 - P;
      if( p12 < 0. )
      {
        T3double z = - ( b[ct1] + r2 * 2 * c[ct1] ) / dd;
        cost = ( d1 * z * ( fs + ( P - r1 ) * dd ) + d2 * cost ) / ( d1 * z + d2 );
      }
      if( p21 < 0. )
      {
        const long double z = - ( b[ct2] + r1 * 2 * c[ct2] ) / dd;
        F = (d2 * z * ( fs + ( P - r1 ) * dd ) + d1 * F ) / ( d1 + d2 * z );
      }
      cost = ( d1 * F + d2 * cost ) / ( d1 + d2);
    }
    else
    {
      const long double p11 = b[ct1] + dr1 * c[ct1];
      const long double p22 = b[ct2] + dr2 * c[ct2];
      long double d1 = P - r1;
      long double d2 = r2 - P;
      long double dr = r2 - r1;
      if(p11 > 0 && p21 > 0 && p11 > p21)
      {
        dr /= p11 / p21;
        if(d1 < dr)
        {
          d2 = r1 + dr - P;
          cost = ( d1 * ( a[ct2] + P * ( b[ct2] + P * c[ct2] ) ) + d2 * cost ) / dr;
        }
        else cost = a[ct2] + P * ( b[ct2] + P * c[ct2] );
      }
      else if (p12 > 0 && p21 > 0 && p22 > p12)
      {
        dr /= p22 / p12;
        if(d2 < dr)
        {
          d1 = P - r2 + dr;
          cost = ( d1 * ( a[ct2] + P * ( b[ct2] + P * c[ct2] ) ) + d2 * cost ) / dr;
        }
      }
      else cost = ( d1 * ( a[ct2] + P * ( b[ct2] + P * c[ct2] ) ) + d2 * cost ) / dr;
    }
  }
  else if(!hct && b[ct1] < 0.)
  {
    const long double z = - b[ct1] / dd;
    const long double d1 = P;
    const long double d2 = r2 - P;
    cost = (d2 * z * (P * dd - 1.) + d1 * cost ) / ( d1 + d2 * z );
  }
  else if(hct == _num_point && b[ct1] + 2 * c[ct1] < 0.)
  {
    const long double z = - ( b[ct1] + 2 * c[ct1] ) / dd;
    const long double d1 = 1. - P;
    const long double d2 = P - r1;
    cost = ( d2 * z * ( 1. - ( 1. - P ) * dd ) + d1 * cost ) / ( d1 + d2 * z );
  }
  return T3double(cost);
}

T3NSGangular_node& T3NSGangular_node::operator=(T3NSGangular_node const &node)
{
  memcpy(this->V, node.V, _num_point * sizeof(T3double));
  memcpy(this->a, node.a, _num_point * sizeof(T3double));
  memcpy(this->b, node.b, _num_point * sizeof(T3double));
  memcpy(this->c, node.c, _num_point * sizeof(T3double));
  sl = node.Get_sl();
  pr = node.Get_pr();
  return (*this);
}


std::ostream& operator<<(std::ostream& os, const T3NSGangular_node& inst)
{
  os << "Node" << " pr = " << inst.Get_pr() << ", sl = "
     << inst.Get_sl() << ", integrated distribution:";
  for(size_t ind =1; ind <= inst._num_point; ++ind) os << " " << inst.Get_V(ind-1);
  return os;
}
