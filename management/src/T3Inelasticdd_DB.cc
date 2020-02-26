#include "T3Inelasticdd_DB.hh"
#include "T3NSGangular_node.hh"
#include "T3Utility.hh"

/*
unsigned int Rand32(unsigned int & xn)
{
  u_quad_t a=0x5DEECE66D;
  u_quad_t c=0xB;
  return (unsigned int)((a*xn+c) & 0xFFFFFFFF);
}

double rndv(unsigned int xn)
{
  return (double) xn / (double) 0x100000000LL;
}

double RND01(unsigned int & xn)
{
  xn=Rand32(xn);
  return (double)(xn) / (double) 0x100000000LL;
}
*/


T3Inelasticdd_DB::T3Inelasticdd_DB(const T3Inelasticdd_DB & nodedb)
{
  for(int i=0; i<512; ++i)
  {
    memcpy(this->V[i], nodedb.V[i], _num_point * sizeof(T3double));
    memcpy(this->a[i], nodedb.a[i], _num_point * sizeof(T3double));
    memcpy(this->b[i], nodedb.b[i], _num_point * sizeof(T3double));
    memcpy(this->c[i], nodedb.c[i], _num_point * sizeof(T3double));
  }
  memcpy(this->pr, nodedb.pr, 512 * sizeof(T3double));
  memcpy(this->sl, nodedb.sl, 512 * sizeof(T3double));
  memcpy(this->Einc, nodedb.Einc, 512 * sizeof(T3double));
  this->Xmin=nodedb.Xmin;
  this->Xmax=nodedb.Xmax;
}

T3Inelasticdd_DB & T3Inelasticdd_DB::operator=(const T3Inelasticdd_DB & nodedb)
{
  for(int i=0; i<512; ++i)
  {
    memcpy(this->V[i], nodedb.V[i], _num_point * sizeof(T3double));
    memcpy(this->a[i], nodedb.a[i], _num_point * sizeof(T3double));
    memcpy(this->b[i], nodedb.b[i], _num_point * sizeof(T3double));
    memcpy(this->c[i], nodedb.c[i], _num_point * sizeof(T3double));
  }
  memcpy(this->pr, nodedb.pr, 512 * sizeof(T3double));
  memcpy(this->sl, nodedb.sl, 512 * sizeof(T3double));
  memcpy(this->Einc, nodedb.Einc, 512 * sizeof(T3double));
  this->Xmin=nodedb.Xmin;
  this->Xmax=nodedb.Xmax;
  return (*this);
}

//#pragma acc routine seq
T3double T3Inelasticdd_DB::RandomizeCost_node(unsigned int & generator_seed, size_t node_num) const
{

  //std::cout<<"RandomizeCost_node():"<<std::endl;
  T3double P = RND01(generator_seed);//t3::GenerateSubCanonical<FloatingType>(generator);
  if(P < pr[node_num])
  {
    return 1. + sl[node_num] * log(1. - ( 1. - exp(-2 / sl[node_num]) ) *
                                   RND01(generator_seed) );
                           //t3::GenerateSubCanonical<FloatingType>(generator) );
  }
  //delta of Ox axis = xmax - xmin:
  const T3double DeltaX=Xmax-Xmin;
  ///\\\///const T3double dY = 2. / (_num_point + 1);
  const T3double dY = DeltaX / (_num_point + 1);
  P = RND01(generator_seed);//t3::GenerateSubCanonical<FloatingType>(generator);
  const size_t SIZE=127;
  const size_t hct = T3Utility::bin_search<SIZE>(V[node_num], P) - V[node_num];
  T3int const ct1 = hct ? hct-1 : hct;                  
  T3int const ct2 = (hct && hct < _num_point) ? hct : -1;
  T3double const r1 = hct ? V[node_num][ct1] : 0.;    
  T3double const r2 = (hct < _num_point) ? V[node_num][hct] : 1.;
  double cost = a[node_num][ct1] + P*( b[node_num][ct1] + P * c[node_num][ct1] );
  const double dd = dY / ( r2 - r1);
  if(ct2 >= 0)
  {
    const T3double dr1 = r1 + r1;
    const T3double dr2 = r2 + r2;
    const double p12 = b[node_num][ct1] + dr2 * c[node_num][ct1];
    const double p21 = b[node_num][ct2] + dr1 * c[node_num][ct2];
    if(p12 < 0. || p21 < 0.)
    {
      double F  =  a[node_num][ct2] + P*( b[node_num][ct2] + P * c[node_num][ct2] );
      ///\\\///const double fs = dY * hct - 1.;
      const double fs = dY * hct + Xmin;
      const double d1 = P - r1;
      const double d2 = r2 - P;
      if( p12 < 0. )                       
      {
        T3double z = - ( b[node_num][ct1] + r2 * 2 * c[node_num][ct1] ) / dd;
        cost = ( d1 * z * ( fs + ( P - r1 ) * dd ) + d2 * cost ) / ( d1 * z + d2 );
      }
      if( p21 < 0. )
      {
        const double z = - ( b[node_num][ct2] + r1 * 2 * c[node_num][ct2] ) / dd;
        F = (d2 * z * ( fs + ( P - r1 ) * dd ) + d1 * F ) / ( d1 + d2 * z );
      }
      cost = ( d1 * F + d2 * cost ) / ( d1 + d2);
    }
    else
    {
      const double p11 = b[node_num][ct1] + dr1 * c[node_num][ct1];
      const double p22 = b[node_num][ct2] + dr2 * c[node_num][ct2];
      double d1 = P - r1;
      double d2 = r2 - P;
      double dr = r2 - r1;
      if(p11 > 0 && p21 > 0 && p11 > p21)
      {
        dr /= p11 / p21;
        if(d1 < dr)
        {
          d2 = r1 + dr - P;
          cost = ( d1 * ( a[node_num][ct2] + P * ( b[node_num][ct2] + P * c[node_num][ct2] ) ) + d2 * cost ) / dr;
        }
        else cost = a[node_num][ct2] + P * ( b[node_num][ct2] + P * c[node_num][ct2] );
      }
      else if (p12 > 0 && p21 > 0 && p22 > p12)
      {
        dr /= p22 / p12;
        if(d2 < dr)
        {
          d1 = P - r2 + dr;
          cost = ( d1 * ( a[node_num][ct2] + P * ( b[node_num][ct2] + P * c[node_num][ct2] ) ) + d2 * cost ) / dr;
        }
      }
      else cost = ( d1 * ( a[node_num][ct2] + P * ( b[node_num][ct2] + P * c[node_num][ct2] ) ) + d2 * cost ) / dr;
    }
  }
  else if(!hct && b[node_num][ct1] < 0.)
  {
    const double z = - b[node_num][ct1] / dd;
    const double d1 = P;
    const double d2 = r2 - P;
    cost = (d2 * z * (P * dd - 1.) + d1 * cost ) / ( d1 + d2 * z );
  }
  else if(hct == _num_point && b[node_num][ct1] + 2 * c[node_num][ct1] < 0.)
  {
    const double z = - ( b[node_num][ct1] + 2 * c[node_num][ct1] ) / dd;
    const double d1 = 1. - P;
    const double d2 = P - r1;
    cost = ( d2 * z * ( 1. - ( 1. - P ) * dd ) + d1 * cost ) / ( d1 + d2 * z );
  }
  return T3double(cost);
}

//#pragma acc routine seq
T3double T3Inelasticdd_DB::RandomizeCost(unsigned int & generator_seed, T3double E/*Tls*/) const
{

  //std::cout<<"RandomizeCost():"<<std::endl;
  
  const T3double Emin = Einc[0];
  const T3double Emax = Einc[_num_nodes - 1];
  const T3double Ediff = Emax - Emin;
  const size_t nbins = _num_nodes - 1;
  const T3double dE = Ediff / nbins;
  if(E <= Emin || E >= Emax)
  {

    //std::cout<<"E="<<E<<std::endl;
    
    const T3double R = RND01(generator_seed);//t3::GenerateSubCanonical<FloatingType>(generator);
    return -1 + 2*R;
  }
  else
  {
    size_t const binn = static_cast<size_t>((E - Emin)*(_num_nodes-1)/Ediff);
    const T3double Elow  = Einc[binn];
    const T3double W = (E - Elow)/(dE);
    const T3double R = RND01(generator_seed);//t3::GenerateSubCanonical<FloatingType>(generator);
    if(R < W)
      return RandomizeCost_node(generator_seed, binn+1);
    else
      return RandomizeCost_node(generator_seed, binn);
  }
}






