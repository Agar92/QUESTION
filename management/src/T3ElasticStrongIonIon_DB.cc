#include "T3ElasticStrongIonIon_DB.hh"
#include "T3NSGangular_node.hh"
#include "T3Utility.hh"

T3ElasticStrongIonIon_DB::T3ElasticStrongIonIon_DB(T3double xmin, T3double xmax):Xmin(xmin),Xmax(xmax)
{
  Xmin=std::log(CosThetaCMmin);//=ln(1.0e-3)=-6.9.
  Xmax=std::log(CosThetaCMmax);//=ln(1.0)=0.
  deltaLnCosThetaCM=(Xmax-Xmin)/(_num_point+1);
  //we have 129 values of 1-costhetacm: 0 ... 1 and 128 bins.
  //but we fill lncoscm only with 127 values excluding 0 and 1:
  for(int i=0; i<129; ++i)
    lncoscm[i]=Xmin+i*deltaLnCosThetaCM;
}

T3ElasticStrongIonIon_DB::T3ElasticStrongIonIon_DB(const T3ElasticStrongIonIon_DB & nodedb)
{
  for(int i=0; i<512; ++i)
  {
    memcpy(this->V[i], nodedb.V[i], _num_point * sizeof(T3double));
    memcpy(this->a[i], nodedb.a[i], _num_point * sizeof(T3double));
    memcpy(this->b[i], nodedb.b[i], _num_point * sizeof(T3double));
    memcpy(this->c[i], nodedb.c[i], _num_point * sizeof(T3double));
  }
  memcpy(this->Einc, nodedb.Einc, 512 * sizeof(T3double));
  memcpy(this->lncoscm, nodedb.lncoscm, 129 * sizeof(T3double));
  this->deltaLnCosThetaCM=nodedb.deltaLnCosThetaCM;
  this->Xmin=nodedb.Xmin;
  this->Xmax=nodedb.Xmax;
}

T3ElasticStrongIonIon_DB & T3ElasticStrongIonIon_DB::operator=(const T3ElasticStrongIonIon_DB & nodedb)
{
  for(int i=0; i<512; ++i)
  {
    memcpy(this->V[i], nodedb.V[i], _num_point * sizeof(T3double));
    memcpy(this->a[i], nodedb.a[i], _num_point * sizeof(T3double));
    memcpy(this->b[i], nodedb.b[i], _num_point * sizeof(T3double));
    memcpy(this->c[i], nodedb.c[i], _num_point * sizeof(T3double));
  }
  memcpy(this->Einc, nodedb.Einc, 512 * sizeof(T3double));
  memcpy(this->lncoscm, nodedb.lncoscm, 129 * sizeof(T3double));
  this->deltaLnCosThetaCM=nodedb.deltaLnCosThetaCM;
  this->Xmin=nodedb.Xmin;
  this->Xmax=nodedb.Xmax;
  return (*this);
}

//this functions return random cos(theta_cm) by random number P.
//V - partial sums.
//a, b, c - TPT2 quadratic interpolation coefficients.
//Einc - incident particle LS energies (Tls).
//Ox must be linear!!!
//Xmin - min value of linear Ox
//Xmax - max value of linear Ox
//Not used:
//deltaLnCosThetaCM - logarithmic step of Ox, not used.
//lncoscm - ln(1-cos(theta_cm)) array - Ox, not used.

//This function takes the partial sums and ln(1-cos(theta_cm))
//linear Ox and interpolates the value of ln(1-cos(theta_cm))
//by the random number P.
//At the end, before return, it makes 1-exp(ln(1-cos(theta_cm)))
//and returns random cos(theta_cm) from 0 to 0.999.

//#pragma acc routine seq
//returns random ln(1-cos(theta_cm))
T3double T3ElasticStrongIonIon_DB::RandomizeCost_node(unsigned int & generator_seed, size_t node_num) const
{
//????//  
  ///\\\///const T3double dY = 2. / (_num_point + 1);
  //Ox approximation sums approximation range:
  const T3double DeltaX=Xmax-Xmin;
  const T3double dY = DeltaX / (_num_point + 1);
//????//


  /*
  std::cout<<"T3ElasticStrongIonIon_DB::RandomizeCost_node():"<<std::endl;
  std::cout<<"STEP #1:\nXmin="<<Xmin<<" Xmax="<<Xmax<<" DeltaX="<<DeltaX
           <<"  _num_point+1="<<(_num_point+1)<<" dY="<<dY
           <<std::endl;
  */
  
  
  
  T3double P = RND01(generator_seed);
  const size_t SIZE=127;
  const size_t hct = T3Utility::bin_search<SIZE>(V[node_num], P) - V[node_num];
  T3int const ct1 = hct ? hct-1 : hct;                  
  T3int const ct2 = (hct && hct < _num_point) ? hct : -1;
  T3double const r1 = hct ? V[node_num][ct1] : 0.;    
  T3double const r2 = (hct < _num_point) ? V[node_num][hct] : 1.;
//!!!
//HERE cost IS NOT cos(theta_cm)!!!
//HERE cost IS ln(1-cos(theta_cm))!!!
//Left the name from TPT2 not to make mistakes.!!!
//So, cost RETURNED FROM THIS FUNCTION GIVES:
// cos(theta_cm) = 1.0-EXP(cost)
//HERE 1-cos(theta_cm) MUST BE FROM 1.0e-3 TO 1.0!!!
//=> ln(1-cos(theta_cm)) IS FROM ln(1.0e-3) TO 0.
//!!!  
  /*long */double cost = a[node_num][ct1] + P*( b[node_num][ct1] + P * c[node_num][ct1] );
//????//
  const /*long */double dd = dY / ( r2 - r1);
//????//


  /*
  std::cout<<"Check partial sums:"<<std::endl;
  std::cout<<"node_num="<<node_num<<std::endl;
  std::cout<<"E="<<Get_Einc(node_num)<<std::endl;
  std::cout<<"V:"<<std::endl;
  for(int i=0; i<_num_point; ++i) std::cout<<V[node_num][i]<<" ";
  std::cout<<std::endl;
  


  std::cout<<"STEP #2:"<<std::endl;
  std::cout<<"P="<<P<<" SIZE="<<SIZE<<std::endl;
  std::cout<<"node_num="<<node_num<<" E["<<node_num-1<<"]="<<Get_Einc(node_num-1)
           <<" E["<<node_num<<"]="<<Get_Einc(node_num)
           <<" E["<<node_num+1<<"]="<<Get_Einc(node_num+1)
           <<std::endl;
  std::cout<<"V:"<<std::endl;
  for(int i=0; i<_num_point; ++i) std::cout<<V[node_num][i]<<" ";
  std::cout<<std::endl;
  */
  
  /*
  std::cout<<"ln(1-coscm):"<<std::endl;
  for(int i=0; i<129; ++i) std::cout<<lncoscm[i]<<" ";
  std::cout<<std::endl;
  */

  /*
  std::cout<<"P="<<P<<" hct="<<hct<<" ct1="<<ct1<<" ct2="<<ct2
           <<" r1="<<r1<<" r2="<<r2
           <<" 1-coscm[ct1]="<<std::exp(lncoscm[ct1])
           <<" 1-coscm[ct1+1]="<<std::exp(lncoscm[ct1+1])
           <<" coscm[ct1]="<<1.0-std::exp(lncoscm[ct1])
           <<" coscm[ct1+1]="<<1.0-std::exp(lncoscm[ct1+1])
           <<" r1="<<r1<<" r2="<<r2
           <<std::endl;
  */

  /*
  std::cout<<"1-coscm:"<<std::endl;
  for(int i=0; i<129; ++i) std::cout<<std::exp(lncoscm[i])<<" ";
  std::cout<<std::endl;
  */

  //std::cout<<"Xmin="<<Xmin<<" Xmax="<<Xmax<<" DeltaX="<<DeltaX
  //         <<" dY="<<dY<<" ct1="<<ct1<<" ct2="<<ct2<<std::endl;
  

  if(ct2 >= 0)
  {

    //std::cout<<"A"<<std::endl;
    
    const T3double dr1 = r1 + r1;
    const T3double dr2 = r2 + r2;
    const /*long */double p12 = b[node_num][ct1] + dr2 * c[node_num][ct1];
    const /*long */double p21 = b[node_num][ct2] + dr1 * c[node_num][ct2];

    //std::cout<<"p12="<<p12<<" p21="<<p21<<std::endl;
    
    if(p12 < 0. || p21 < 0.)
    {

      //std::cout<<"A1"<<std::endl;
      
      /*long */double F  =  a[node_num][ct2] + P*( b[node_num][ct2] + P * c[node_num][ct2] );
//????//
      ///\\\///const /*long */double fs = dY * hct - 1.;
      const /*long */double fs = dY * hct + Xmin;
//????//
      const /*long */double d1 = P - r1;
      const /*long */double d2 = r2 - P;
      if( p12 < 0. )                       
      {
        T3double z = - ( b[node_num][ct1] + r2 * 2 * c[node_num][ct1] ) / dd;
        cost = ( d1 * z * ( fs + ( P - r1 ) * dd ) + d2 * cost ) / ( d1 * z + d2 );
      }
      if( p21 < 0. )
      {
        const /*long */double z = - ( b[node_num][ct2] + r1 * 2 * c[node_num][ct2] ) / dd;
        F = (d2 * z * ( fs + ( P - r1 ) * dd ) + d1 * F ) / ( d1 + d2 * z );
      }
      cost = ( d1 * F + d2 * cost ) / ( d1 + d2);
    }
    else
    {

      //std::cout<<"A2"<<std::endl;
      
      const /*long */double p11 = b[node_num][ct1] + dr1 * c[node_num][ct1];
      const /*long */double p22 = b[node_num][ct2] + dr2 * c[node_num][ct2];
      /*long */double d1 = P - r1;
      /*long */double d2 = r2 - P;
      /*long */double dr = r2 - r1;
      if(p11 > 0 && p21 > 0 && p11 > p21)
      {

        //std::cout<<"A21"<<std::endl;
        
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

        //std::cout<<"A22"<<std::endl;
        
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

    //std::cout<<"B"<<std::endl;
    
    const /*long */double z = - b[node_num][ct1] / dd;
    const /*long */double d1 = P;
    const /*long */double d2 = r2 - P;
    cost = (d2 * z * (P * dd - 1.) + d1 * cost ) / ( d1 + d2 * z );
  }
  else if(hct == _num_point && b[node_num][ct1] + 2 * c[node_num][ct1] < 0.)
  {

    //std::cout<<"C"<<std::endl;
    
    const /*long */double z = - ( b[node_num][ct1] + 2 * c[node_num][ct1] ) / dd;
    const /*long */double d1 = 1. - P;
    const /*long */double d2 = P - r1;
    cost = ( d2 * z * ( 1. - ( 1. - P ) * dd ) + d1 * cost ) / ( d1 + d2 * z );
  }
//HERE cost IS NOT cos(theta_cm)!!!
//HERE cost IS ln(1-cos(theta_cm))!!!
//Left the name from TPT2 not to make mistakes.!!!
//And we need to get cos(theta_cm) from ln(1-cos(theta_cm)):
// cost=1.0-exp(ln(1-cos(theta_cm)))=1.0-exp(cost).
  
  //cost=1.0-std::exp(cost);
  //std::cout<<"Check cost: cost="<<cost<<std::endl;
  
  return T3double(cost);
}


//!!This function is only for CPU, because C array is allocated in it.!!!
/*
//#pragma acc routine seq
T3double T3ElasticStrongIonIon_DB::RandomizeCost_node1(unsigned int & generator_seed, size_t node_num) const
{
  //std::cout<<"AAA"<<std::endl;
  const T3double DeltaX=Xmax-Xmin;
  const T3double dX=DeltaX/(_num_point+1);
  T3double tempar[129]{0.0};
  tempar[0]=0.0;
  for(int j=0; j<127; ++j) tempar[j+1]=V[node_num][j];
  tempar[128]=1.0;
  
  //std::cout<<"tempar:"<<std::endl;
  //for(int j=0; j<129; ++j) std::cout<<"("<<j<<","<<tempar[j]<<") ";
  //std::cout<<std::endl;

  //std::cout<<"lncoscm:"<<std::endl;
  //for(int j=0; j<129; ++j) std::cout<<"("<<j<<","<<lncoscm[j]<<") ";
  //std::cout<<std::endl;
  
  T3double P = RND01(generator_seed);
  const size_t SIZE=129;
  const size_t hct = T3Utility::bin_search<SIZE>(tempar, P) - tempar - 1;  
  const T3double df=P-tempar[hct];
  const T3double f1=tempar[hct];
  const T3double f2=tempar[hct+1];
  const T3double x1=lncoscm[hct];
  const T3double x2=lncoscm[hct+1];
  T3double cost=x1+df/(f2-f1)*(x2-x1);
  //const T3double bb=std::exp(cost);
  //cost=1.0-bb;

  
  //std::cout<<"P="<<P<<" hct="<<hct
  //         <<" SIZE="<<129<<" f1="<<f1<<" f2="<<f2
  //         <<" x1="<<x1<<" x2="<<x2
  //         <<" df="<<df<<std::endl;  
  
  if(std::isinf(cost))
  {
    std::cout<<"cost="<<cost<<" bb="<<bb<<std::endl;
    std::cout<<"Inf***"<<std::endl;
    exit(0);
  }

  //std::cout<<"node_num="<<node_num<<" Xmin="<<Xmin<<" Xmax="<<Xmax
  //         <<" DeltaX="<<DeltaX<<" dX="<<dX<<" P="<<P
  //         <<" SIZE="<<SIZE<<" hct="<<hct<<std::endl;
  
  //std::cout<<"V:"<<std::endl;
  //for(int i=0; i<129; ++i) std::cout<<tempar[i]<<" ";
  //std::cout<<std::endl;

  //std::cout<<"df="<<df<<" f1="<<f1<<" f2="<<f2<<" x1="<<x1<<" x2="<<x2
  //         <<" cost="<<cost<<std::endl;
  //exit(0);

  //std::cout<<"node_num="<<node_num<<" E="<<Get_Einc(node_num)<<std::endl;
  
  
  //std::cout<<"hct="<<hct<<std::endl;
  return T3double(cost);
}
*/


//#pragma acc routine seq
//returns random ln(1-cos(theta_cm))
T3double T3ElasticStrongIonIon_DB::RandomizeCost(unsigned int & generator_seed, T3double E) const
{
  //std::cout<<"RandomizeCost():"<<std::endl;
  const T3double Emin = Einc[0];
  const T3double Emax = Einc[_num_nodes - 1];
  const T3double Ediff = Emax - Emin;
  const size_t nbins = _num_nodes - 1;
  const T3double dE = Ediff / nbins;
  if(E <= Emin || E >= Emax)
  {
    //std::cout<<"E1"<<std::endl;
    const T3double R = RND01(generator_seed);
    return -1 + 2*R;
  }
  else
  {
    //std::cout<<"E2"<<std::endl;
    size_t const binn = static_cast<size_t>((E - Emin)*(_num_nodes-1)/Ediff);
    const T3double Elow  = Einc[binn];
    const T3double W = (E - Elow)/(dE);
    const T3double R = RND01(generator_seed);
    if(R < W)
      return RandomizeCost_node(generator_seed, binn+1);
    else
      return RandomizeCost_node(generator_seed, binn);
  }
}






