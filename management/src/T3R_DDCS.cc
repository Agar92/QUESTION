#include "T3R_DDCS.hh"
#include <cmath>
#include <map>

#include "unistd.h"

namespace t3 {

T3R_DDCS::T3R_DDCS():lnE{0.0}, lncos{0.0}, Energy{0.0},
                     Fk{0.0}, CS{0.0}, ecs{}, tps{}, trw{}
{
  trw.resize(0);
}

void T3R_DDCS::Fill()
{
//prepare E axis for energies in partial sums db:
  const double lnEmin=std::log(Emin);
  const double lnEmax=std::log(Emax);
  deltalnE=(lnEmax-lnEmin)/(Bin2-1);
  for(int b=0; b<Bin2; ++b) lnE.at(b)=lnEmin+b*deltalnE;
//prepare 1-cos(theta_cm) axis:
  deltalncos=(std::log(1.0)-std::log(kb))/Bin3;
  for(int b=0; b<Bin3+1; ++b) lncos.at(b)=std::log(kb)+b*deltalncos;
//here we divide the 1-cos(theta_cm) intevals by SubBin3 subbins:
  std::array<std::array<double, SubBin3+1>, Bin3> SubBinsCS{0.0};
  const double deltalncosSubBins=deltalncos/SubBin3;
  for(int b=0; b<Bin3; ++b)
    for(int sb=0; sb<SubBin3+1; ++sb) SubBinsCS.at(b).at(sb)=lncos.at(b)+sb*deltalncosSubBins;
  
//iterate through all energies:
  for(int i=0; i<Bin2; ++i)
  {
    const double Tls=std::exp(lnE.at(i));//kinetic energy of the particle E[i] in inner units (MeV).
    //calculate tmax:
    const double pls2 = Tls*(2*md+Tls);
    const double pcm2 = pls2/(2.+2.*std::sqrt(1.+pls2/md2));//rat=md/md=1. rat2=1 for D-D.
    const double tm =   2*pcm2;
    //std::cout<<"Tls("<<i<<")="<<Tls<<" pls2="<<pls2/GeV/GeV<<" pcm2="<<pcm2/GeV/GeV<<" tm="<<tm/GeV/GeV<<std::endl;
//fill partial sums.
    double psumb=0.0;//here we will sum sum_0_Bin(dsigma/d|t|*d|t|).   
    Fk.at(i).at(0)=0.0;

    std::cout<<"Tls="<<Tls<<" pls2="<<pls2<<" pcm2="<<pcm2
             <<" tm="<<tm<<std::endl;
    
    for(int b=1; b<Bin3+1; ++b)
    { 
      std::array<double, SubBin3+1> subbinscs=SubBinsCS.at(Bin3-b);
      for(int sb=0; sb<SubBin3; ++sb)
      {
        const double tl=tm*std::exp(subbinscs.at(SubBin3-sb-1));
        const double tr=tm*std::exp(subbinscs.at(SubBin3-sb));
        double dleft=tl*GetdSigmadt(std::exp(lnE.at(i)), tl);
        double dright=tr*GetdSigmadt(std::exp(lnE.at(i)), tr);
        double dlncos=deltalncosSubBins;//=dx
        const double dsigmai=(dleft+dright)*dlncos/2;
        psumb+=dsigmai;    
      }
      Fk.at(i).at(b)=psumb;
    }
  }

//**************************************************************************//
//fill (E,CS) database:
  for(int i=0; i<Bin2; ++i)
  {
    const double Tls=std::exp(lnE.at(i));//kinetic energy of the particle E[i] in inner units.
    const double pls2 = Tls*(2*md+Tls);
    const double pcm2 = pls2/(2.+2.*std::sqrt(1.+pls2/md2));//rat=md/md=1. rat2=1 for D-D.
    const double tm =   2*pcm2;
  
    CS.at(i)=2*Fk.at(i).at(Bin3);
  }
  std::vector<T3double> argE;
  std::vector<T3double> argCS;
  argE.resize(Bin2);
  argCS.resize(Bin2);
  for(int i=0; i<Bin2; ++i)
  {
    argE.at(i)=std::exp(lnE.at(i));
    argCS.at(i)=CS.at(i);
    //ecs.push_back(new T3R_node(argE.at(i), argCS.at(i)));
  }
  T3TabulatedCS tcs(argE,argCS);
  std::cout<<"tcs: "<<" size="<<tcs.Get_size()<<std::endl;
  //std::cout<<tcs<<std::endl;

  std::cout<<"STEP #1"<<std::endl;
  //sleep(3);
  ecs=T3R_RW(tcs, 1, 2);

  std::cout<<"PRINT ecs: size="<<ecs.size()<<std::endl;
  std::cout<<ecs<<std::endl;
  for(int i=0; i<Bin2; ++i)
  {
    int count=0;
    for(int b=0; b<Bin3+1; ++b)
    {
      Fk.at(i).at(b)/=Fk.at(i).at(Bin3);
      ++count;
    }
  }

  tps.resize(0);
  for(int i=0; i<Bin2; ++i)
  {
    std::array<double, _numpoint> PSFk{0.0};
    for(int j=0; j<_numpoint; ++j) PSFk.at(j)=Fk.at(i).at(j+1);
    T3NSGangular_RWnode newnode =
      T3NSGangular_RWnode(std::exp(lnE.at(i)), PSFk);
    tps.push_back(newnode);
  }
  trw.push_back(tps);
  std::cout<<"End of Fill():"<<std::endl;
  std::cout<<" trw size="<<trw.size()<<" tps size="<<tps.size()<<std::endl;
  //std::cout<<"tps: \n"<<tps<<std::endl;
  std::cout<<"trw: \n"<<trw<<std::endl;
  std::cout<<"END ecs: "<<ecs.size()<<std::endl;
}

//returns the differential cross section from D-D approximation in inner units
double T3R_DDCS::GetdSigmadt(double Tls/*LS kinetic energy of the incident deuteron in inner units*/,
                             double t/*|t| - Mandelstahm variable |t| in inner units*/) const
{
//1. in aproximation Tls is in GeV. calculate tmax:
  const double pls2 = Tls*(2*md+Tls);
  const double pcm2 = pls2/(2.+2.*std::sqrt(1.+pls2/md2));//rat=md/md=1. rat2=1 for D-D.
  const double tmax = 4*pcm2;//tmax=4*pcm2.
//if |t|>|t|max, exit:
  if(t>tmax)
  {
    ////std::cout<<"Returns 0.0"<<" |t|="<<t/GeV/GeV<<" tm="
    ////       <<tmax/2/GeV/GeV<<" tmax="<<tmax/GeV/GeV<<std::endl;
    return 0.0;
  }
  if(t>2*pcm2) t=tmax-t;
  //Ar0 without hc
  const double Ar0=md / 137. * std::sqrt(10 * M_PI / pcm2);
  const double pcm4 = 4 * pcm2;
  const double sin2 = t / pcm4;
  const double cos2 = 1. - sin2;
  const double tg2 = sin2 / cos2;
  //Rutherford amplitude for equal particles d-d:
//---------------------------------------------------------------------------
//when pcm/Ecm ~ alpha:
//---------------------------------------------------------------------------  
  const double alpha=1.0/137.0;
  const double pcm=std::sqrt(pcm2);
  const double Ecm=std::sqrt(md2+pcm2);
  const double CA=std::cos(alpha*Ecm/pcm*std::log(tg2));
//---------------------------------------------------------------------------  
  const double Ar = Ar0 * std::sqrt(1. + tg2 * tg2 + (2. / 3.) * tg2 * CA);
  const double Tcm = Tls / 2;//for D-D Tcm=Tls/2.
  const double s = 4 * md2 * (1. + Tcm / md);//here Tcm in inner units - MeV.
  //here Tcm in in inner units - MeV:
  const double g=std::exp(-0.5*std::sqrt(0.986 * MeV / Tcm));//Gamow factor//0.986 MeV - Gamow energy for D-D.
  const double As = 49000.0 * (1. + std::pow(Tcm/1.173, 3.06)) / (1. + std::pow(Tcm/0.64, 2.0)) /
                    (1. + std::pow(Tcm/5.45, 2.0)) /  (1.+std::pow(Tcm/23.7, 1.775)) * g;
//  As=0.49E5/(1.+(Tcm/0.64)**2.0)*
// *(1.+(Tcm/1.173)**3.06)/(1.+(Tcm/5.45)**2.)/
// /(1.+(Tcm/23.7)**1.775)*g
  const double Atu = 1041.0 * (1. + std::pow(Tcm/0.97, 3.89)) /
    (1. + std::pow(Tcm/0.056, 1.29)) / (1. + std::pow(Tcm/1.574, 2.6)) /
    (1. + std::pow(Tcm/16.75, 1.67)) * g;
//  Atu=1041.0/(1.+(Tcm/0.056)**1.29)*
// *(1.+(Tcm/0.97)**(1.29+2.6))/
// /(1.+(Tcm/1.574)**2.6)/
// /(1.+(Tcm/16.75)**1.67)*g  
  const double phases=std::pow(Tcm/0.93, 2.76) / (1. + std::pow(Tcm/1.63, 2.8));
//  phases=(Tcm/0.93)**2.76/
// /(1.+(Tcm/1.63)**2.8)
  const double phasetu=std::pow(Tcm/2.65, 3.1) / (1. + std::pow(Tcm/3.22, 3.245));
//  phasetu=(Tcm/2.65)**3.1/
// /(1.+(Tcm/3.22)**3.245)
  //squared mass of pi-0 mezon = 0.01822 GeV^2.
  const double mpi02=0.01822 * GeV * GeV;
  const double dtg = mpi02 / pcm4;
  const double tg2t = (sin2 + dtg) / (cos2 + dtg);
  const double termAtu = Atu * std::sqrt(1. + tg2t * tg2t + (2. / 3.) * tg2t * CA);
  const double term11 = Ar / t + As * std::cos(phases) / s;
  const double term12 = termAtu * std::cos(phasetu) / (t + mpi02);
  const double term1 = term11 + term12;
  const double term21 = As * std::sin(phases) / s;
  const double term22 = termAtu * std::sin(phasetu) / (t + mpi02);
  const double term2=term21+term22;
  const double hc = 0.2 * GeV * fm;//=200 MeV*fm.  
  double dsigmadt = hc*hc*(term1*term1+term2*term2);// =/mbarn
  dsigmadt/=10;//because std::sqrt(10 * M_PI / pcm2);.
  return dsigmadt;
}

  
//this function returns the value of the integral of the Rutherford equal particles
//Integral(dsigma/d|t|*d|t|) from |t|low to |t|up
double T3R_DDCS::GetRutherfordDDCSIntegral(double Tls/*LS kin energy of the inc deuteron*/,
                                           double tlow/*lower bound*/, double tup/*upper bound*/) const
{
  if(tlow<0.0) return 0.0;
  //calculate |t|min:
  //const double tmin=2*md*Edisplace_deuteron;
  //calculate |t|max:
  const double pls2 = Tls*(2*md+Tls);
  const double pcm2 = pls2/(2.+2.*std::sqrt(1.+pls2/md2));//rat=md/md=1. rat2=1 for D-D.
  const double tmax = 4*pcm2;
  //reduced mass for d-d:
  const double mr=md/2;
  const double alpha=1.0/137.0;
  const double hc=200.0*MeV*fm ;// 200 MeV * fm.
  const double c1=2*mr*alpha*hc;//2*mr*alpha*z*Z | z=Z=1 for d-d.
  const double c2=c1*c1;
  const double COEF=M_PI/pcm2*c2;
  //auto const term1=1.0/tmin-1.0/tup;
  const double term1=1.0/tlow-1.0/tup;
  //auto const term2=1.0/(tmax-tup)-1.0/(tmax-tmin);
  const double term2=1.0/(tmax-tup)-1.0/(tmax-tlow);
  const double pcm=std::sqrt(pcm2);      //CM momentum
  const double Ecm=std::sqrt(pcm2+mr*mr);//CM energy of the reduced particle
  const double beta_r=pcm/Ecm;           //CM beta of the reduced particle
  const double c3=2.0/3.0*1.0/tmax*beta_r/alpha;
  const double term31=std::sin(alpha/beta_r*std::log(tup/(tmax-tup)));
  //auto const term32=std::sin(alpha/beta_r*std::log(tmin/(tmax-tmin)));
  const double term32=std::sin(alpha/beta_r*std::log(tlow/(tmax-tlow)));
  const double term3=c3*(term31-term32);
  const double integral_cs=COEF*(term1+term2+term3);
  //std::cout<<"Tls="<<Tls/MeV<<" t="<<tup/GeV/GeV<<" tlow="<<tlow/GeV/GeV
  //         <<" cs="<<integral_cs/barn<<std::endl;
  // returns the integrated D-D differential cross section
  // from |t|low to |t|up
  //(from |t|min to k*|t|max).
  return integral_cs;
}

//returns the Rutherford differential cross section for d-d (Rutherford formula for equla particles)
double T3R_DDCS::GetdSigmadt_DD_Rutherford(double Tls/*LS kinetic energy of the incident deuteron in inner units*/,
                                           double t/*|t| - Mandelstahm variable |t| in inner units*/) const
{
  const double pls2 = Tls*(2*md+Tls);
  const double pcm2 = pls2/(2.+2.*std::sqrt(1.+pls2/md2));//rat=md/md=1. rat2=1 for D-D.
  const double tmax = 4*pcm2;//tmax=4*pcm2.
  //Ar0 without hc
  const double Ar0=md / 137. * std::sqrt(M_PI / pcm2);
  //Rutherford amplitude for equal particles d-d:
//---------------------------------------------------------------------------
//when pcm/Ecm ~ alpha:
//---------------------------------------------------------------------------
  const double alpha=1.0/137.0;
  const double pcm=std::sqrt(pcm2);
// ////////////////////////////////////////////////
  //Here is an error: Ecm - CM energy of the inc particle or the reduced particle:
  //1) if inc particle - Ecm=std::sqrt(md2+pcm2)
  //2) if reduced particle - Ecm=std::sqrt(md2/4+pcm2)
// ////////////////////////////////////////////////
  const double mr=md/2;
  const double mr2=mr*mr;
  const double Ecm=std::sqrt(mr2+pcm2);
  const double argln=t/(tmax-t);
  const double CA=std::cos(alpha*Ecm/pcm*std::log(argln));
//---------------------------------------------------------------------------
  const double term1 = 1.0/t/t;
  const double term2=1.0/(tmax-t)/(tmax-t);
  const double term3 = 2.0/3.0/t/(tmax-t)*CA;
  const double hc = 0.2 * GeV * fm;//=200 MeV*fm.
  const double coef=Ar0*Ar0;
  double dsigmadt = hc*hc*coef*(term1+term2+term3);// =/mbarn
  return dsigmadt;
}

void T3R_DDCS::Save_to_CS()
{
  //target - deuteron:
  std::cout<<"Check ecs in Save_to_CS():"<<ecs<<std::endl;
  ecs.save_binary_Elastic_approx_e_cs(1,2,1,2);
}

bool T3R_DDCS::Load_from_CS()
{
  //target - deuteron:
  ecs.load_binary_Elastic_approx_e_cs(1,2,1,2);
  //std::cout<<"Check (E,XS) in Load_from_CS()"<<std::endl;
  //std::cout<<ecs<<std::endl;
  return true;
}

void T3R_DDCS::Save_to_PS(/*target Z and A*/int tgZ, int tgA, const T3String & rid,
                          T3int pdg, T3int incZA)
{
  //target - deuteron:
  std::cout<<"Check trw in Save_to_PS():"<<trw<<std::endl;
  trw.Set_Z(tgZ);
  trw.Set_A(tgA);
  trw.save_binary(1,2,rid,pdg,incZA);
}

bool T3R_DDCS::Load_from_PS(/*target Z and A*/int tgZ, int tgA, const T3String & rid,
                          T3int pdg, T3int incZA)
{
  bool debug_Load_from_PS=false;
  //target - deuteron:
  trw.Set_Z(tgZ);
  trw.Set_A(tgA);
  trw.load_binary(1,2,rid,pdg,incZA);
  return true;
}

} // namespace t3
