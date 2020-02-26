#pragma once
#ifndef T3USERSTEPPINGACTIONININJECT_CHECK_DD_ELASTIC_ANG_DISTRIBUTION_H
#define T3USERSTEPPINGACTIONININJECT_CHECK_DD_ELASTIC_ANG_DISTRIBUTION_H

#include <iostream>
#include "T3Globals.h"
#include "T3ParticleTable.h"
#include "T3R_DDCS.hh"

using namespace t3;

//NUMBER OF BINS IN THE HISTOGRAM:
constexpr int HNbinInject_CHECK_DD_ELASTIC=128;

class T3UserSteppingActionInInject_Check_DD_Elastic_ang_distribution
{
public:
  T3UserSteppingActionInInject_Check_DD_Elastic_ang_distribution():
    deuteronPDG(t3::ParticleTable().makePDGfromZandA(1,2)),
    protonPDG(t3::PDG_t(2212)),neutronPDG(t3::PDG_t(2112)),
    titanPDG(t3::ParticleTable().makePDGfromZandA(22,48))
  {
    
    _IsRegistered=false;
    for(int m=0; m<HNbinInject_CHECK_DD_ELASTIC; ++m)
      HistogramElasticDD[m]=0;

    /*
    std::cout<<"T3UserSteppingActionInInject_Check_DD_Elastic_ang_distribution():"<<std::endl;
    std::cout<<"Xmin="<<Xmin<<" Xmax="<<Xmax<<" deltacoscm="<<deltacoscm<<std::endl;
    std::cout<<"deuteronPDG="<<deuteronPDG<<std::endl;
    //sleep(3);
    */
    
  }
  void Register(){_IsRegistered=true;}
  bool IsRegistered(){return _IsRegistered;}
  void UserAction(Particle<FloatingType> & particlei, double tmin, double tm, PDG_t particlePDG)
  {

    //std::cout<<"UserAction():"<<std::endl;

    //std::cout<<"pdg="<<particlei.pdg<<std::endl;
    
    if(particlei.pdg==particlePDG) FillHistogramThetaElasticDD(particlei, tmin, tm, particlePDG);
    if(particlei.pdg==particlePDG) particlei.ir=0;
    
    ////if(particlei.pdg==neutronPDG || particlei.pdg==protonPDG) FillHistogramThetaElasticDD(particlei);

    //the first user action is to kill all the deuterons after their first
    //inelastic d-d reaction:
    //std::cout<<"pdg="<<particlei.pdg<<std::endl;
    //if(particlei.pdg==neutronPDG || particlei.pdg==protonPDG)
    //  particlei.ir=0;
  }
  void FillHistogramThetaElasticDD(Particle<FloatingType> & particlei, double tmin, double tm, PDG_t particlePDG)
  {
    ///const double costcm=particlei.tr;
    const double costcm=log(1.0-particlei.tr);

    //This is for D-D:
    if(particlePDG==deuteronPDG) Xmax=log(1.0e-3);
    //This is for Ti-Ti:
    else if(particlePDG==titanPDG) Xmax=log(1.0);//=0
    
    const double delta=(Xmax-log(tmin/tm))/HNbinInject_CHECK_DD_ELASTIC;
    ///const int bin=(costcm-Xmin)/deltacoscm;
    const int bin=(costcm-log(tmin/tm))/delta;

    /*
    std::cout<<"bin="<<bin<<" tr="<<particlei.tr<<" costcm="<<costcm
             <<" delta="<<delta<<" xmin="<<log(tmin/tm)<<" xmax="
             <<Xmax<<std::endl;
    */
    
    /*
    if(bin>127 || bin<0)
      std::cout<<"bin="<<bin<<" costcm="<<costcm
               <<" Xmin="<<Xmin<<" tmin/tm="<<1.0-tmin/tm
               <<" tr="<<particlei.tr<<std::endl;
    */               
               
    /*
    std::cout<<"coscm="<<costcm<<" costcm-Xmin="<<costcm-Xmin
             <<" deltacoscm="<<deltacoscm<<" bin="<<bin
             <<" 1-tmin/tm="<<1.0-tmin/tm
             <<" tr="<<particlei.tr<<std::endl;
    */
    
    //usleep(100000);

    if(bin>=0 && bin<HNbinInject_CHECK_DD_ELASTIC)
#ifdef OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif
      ++HistogramElasticDD[bin];
  }

  //this function is for writing Rutherford equal particles scattering
  // dsigma/d|t|. It is necessary for choosing how to randomize interference term.
  double GetEqualParticleRutherfordDifferentialCrossSectionValue(PDG_t particlePDG, double Tls, double t/*|t|*/)
  {
    ParticleTable aParticleTable;
    const double m=aParticleTable.GetMass(particlePDG);
    const double hc=0.2*GeV*fm;
    const double alpha=1.0/137.0;
    const double pls2=Tls*(2*m+Tls);
    const double pcm2=pls2/(2.0+2.0*sqrt(1.0+pls2/m/m));
    const double pcm=sqrt(pcm2);
    const double mr=m/2;
    const int    z=aParticleTable.GetZ(particlePDG);
    const double COEF1=2*mr*alpha*hc*z*z;
    const double COEF2=COEF1*COEF1*M_PI/pcm2;
    const double tmax=4*pcm2;
    const double Ecm=sqrt(mr*mr+pcm2);
    const double beta_r_cm=pcm/Ecm;
    const double term1=1.0/t/t;
    const double term2=1.0/(tmax-t)/(tmax-t);
    const double spin=aParticleTable.GetSpin(particlePDG);
    const bool cond=(fmod(spin,1.0))<1.0e-7;
    const int SIGN=(cond)?1:-1;
    const double term3=SIGN*2.0/(2*spin+1.0)/t/(tmax-t)*cos(alpha/beta_r_cm*log(t/(tmax-t)));
    const double dsdt=COEF2*(term1+term2+term3);
    return dsdt;
  }
  
  void Histogram_theta_UserSteppingActionInInject_CHECK_DD_ELASTIC(double CS, double Tls, PDG_t particlePDG)
  {
#ifdef OPENACC
#pragma acc update host(HistogramElasticDD[0:HNbinInject_CHECK_DD_ELASTIC])
#endif
    ParticleTable bParticleTable;
    /*
    const PDG_t deuteronPDG=bParticleTable.makePDGfromZandA(1,2);
    const double md=bParticleTable.GetMass(deuteronPDG);
    const double pls2=Tls*(2*md+Tls);
    const double rat=md/md;
    const double rat2=rat*rat;
    const double pcm2=pls2/(1.0+rat2+2*sqrt(rat2+pls2/md/md));
    const double Edisplace_deuteron=10*eV;
    const double tmin=2*md*Edisplace_deuteron;
    const double tm=2*pcm2;
    const double tmax=4*pcm2;
    */
    
    const double Eded=1.0e-5*MeV;//=10 eV
    const double Edeti=2.5e-5*MeV;//=25 eV
    const double m=bParticleTable.GetMass(particlePDG);
    double tmin;
    if(particlePDG==deuteronPDG) tmin=2*Eded*m;
    else if(particlePDG==titanPDG) tmin=2*Edeti*m;
    const double pls2=Tls*(2*m+Tls);
    const double pcm2=pls2/(2.0+2*sqrt(1.0+pls2/m/m));
    const double tm=2*pcm2;
    //This is for D-D:
    if(particlePDG==deuteronPDG) Xmax=log(1.0e-3);
    //This is for Ti-Ti:
    else if(particlePDG==titanPDG) Xmax=log(1.0);//=0

    /*
    std::cout<<"CS="<<CS/mbarn<<" Tls="<<Tls<<" MeV  mbarn  tmin="<<tmin/GeV/GeV<<" GeV^2"
             <<" tm="<<tm/GeV/GeV<<" GeV^2  tmax="<<tmax/GeV/GeV<<" GeV^2"
             <<std::endl;
    */

//Fill dsigma/d|t| from modelling:    
    std::ofstream foutne_theta;
    //foutne_theta.open("ne_theta_elast_dd.dat");
    foutne_theta.open("t3elastemionion_dd.dat");
    for(int m=0; m<HNbinInject_CHECK_DD_ELASTIC; ++m)
    {
      foutne_theta<<m<<"   ";
      ///const double ti=tm*exp((Xmax-deltacoscm*(m+0.5)));
      const double delta=(Xmax-log(tmin/tm))/HNbinInject_CHECK_DD_ELASTIC;
      const double ti=tm*exp((log(tmin/tm)+delta*(m+0.5)));
      
      /*
      std::cout<<"Xmin="<<Xmin<<" Xmax="<<Xmax
               <<" Nbin="<<HNbinInject_CHECK_DD_ELASTIC
               <<" deltacoscm="<<deltacoscm
               <<" tm="<<tm/MeV/MeV<<" pls2="<<pls2
               <<" pcm2="<<pcm2<<" ti="<<ti/GeV/GeV
               <<" exp(...)="<<exp((Xmax-deltacoscm*(m+0.5)))
               <<std::endl;
      */

      std::cout<<"HistogramElasticDD["<<m<<"]="<<HistogramElasticDD[m]
               <<" CS="<<CS/mbarn<<" Np="<<Np
               <<" deltacoscm="<<deltacoscm<<" ti="<<ti/GeV/GeV<<std::endl;
      
      ///const double dsigmadt=static_cast<double>(CS*HistogramElasticDD[m])/Np/deltacoscm/ti;
      const double dsigmadt_modelling=static_cast<double>(CS*HistogramElasticDD[m])/Np/delta/ti;

      const double dsigmadt_sharp=
        GetEqualParticleRutherfordDifferentialCrossSectionValue(particlePDG,Tls,ti);


      std::cout<<"Tls="<<Tls<<" tmin="<<tmin/GeV/GeV<<" tm="<<tm/GeV/GeV
               <<" 1.0e-3*tm="<<1.0e-3*tm/GeV/GeV<<std::endl;

      std::cout<<"ti="<<ti/GeV/GeV<<" dsigmadt_modelling="
               <<dsigmadt_modelling/mbarn*GeV*GeV<<std::endl;
      
      foutne_theta<<std::setw(8)<<ti/GeV/GeV<<"   "<<HistogramElasticDD[m]<<"   "
                  <<dsigmadt_modelling/mbarn*GeV*GeV
                  <<"   "<<dsigmadt_sharp/mbarn*GeV*GeV;
      foutne_theta<<std::endl;
    }
    foutne_theta.close();

    std::ofstream fout_rat;
    fout_rat.open("t3elastemionion_dd_rat.dat");
    for(int m=0; m<HNbinInject_CHECK_DD_ELASTIC; ++m)
    {
      fout_rat<<m<<"   ";
      const double delta=(Xmax-log(tmin/tm))/HNbinInject_CHECK_DD_ELASTIC;
      const double ti=tm*exp(log(tmin/tm)+delta*(m+0.5));
      
      const double dsigmadt_modelling=static_cast<double>(CS*HistogramElasticDD[m])/Np/delta/ti;
      const double dsigmadt_from_formula_sharp=
        GetEqualParticleRutherfordDifferentialCrossSectionValue(particlePDG,Tls,ti);
      const double ratio=dsigmadt_modelling/dsigmadt_from_formula_sharp;
      ///fout_rat<<std::setw(8)<<1.0-ti/tm<<"  "<<HistogramElasticDD[m]<<"  "<<ratio<<std::endl;

      fout_rat<<std::setw(8)<<ti/GeV/GeV<<"  "<<HistogramElasticDD[m]
              <<"  "<<ratio<<std::endl;


      std::cout<<"m="<<m<<" ti="<<ti/GeV/GeV
               <<" ln(1-ti/tm)="<<log(1.0-ti/tm)<<" rat="<<ratio<<std::endl;
      
    }
    fout_rat.close();

    /*
//Fill ratios:
    T3R_DDCS trdcs;
    std::ofstream fout_rat;
    fout_rat.open("ne_rat.dat");
    for(int m=0; m<HNbinInject_CHECK_DD_ELASTIC; ++m)
    {
      const double minus_coscmi=exp((Xmax-deltacoscm*(m+0.5)));//=1-cos(theta_cm)
      const double ti=tm*exp((Xmax-deltacoscm*(m+0.5)));
      const double dsigmadt_modelling=
        static_cast<double>(CS*HistogramElasticDD[m])/Np/deltacoscm/ti;
      const double dsigmadt_from_approximation_sharp=trdcs.GetdSigmadt(Tls, ti);
      const double ratio=dsigmadt_modelling/dsigmadt_from_approximation_sharp;
      fout_rat<<std::setw(8)<<minus_coscmi<<"  "<<ratio<<std::endl;
    }
    fout_rat.close();
    */
  }

private:
  ///const double Xmin=log(1.0e-3);//1-cos(theta_cm)_min=1.0e-3.
  ///const double Xmax=log(1.0);   //1-cos(theta_cm)_max=1.0.

  const double Xmin=log(1.0e-8);//1-cos(theta_cm)_min=1.0e-3.

  //This is for D-D:
  double Xmax=log(1.0e-3);  
  
  const double deltacoscm=(Xmax-Xmin)/HNbinInject_CHECK_DD_ELASTIC;
  int HistogramElasticDD[HNbinInject_CHECK_DD_ELASTIC];
  bool _IsRegistered;
  const t3::PDG_t deuteronPDG;
  const t3::PDG_t protonPDG;
  const t3::PDG_t neutronPDG;
  const t3::PDG_t titanPDG;
};


#endif//T3USERSTEPPINGACTIONININJECT_CHECK_DD_ELASTIC_ANG_DISTRIBUTION_H
