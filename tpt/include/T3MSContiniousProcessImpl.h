#pragma once
#ifndef T3MSCONTINIOUSPROCESSIMPL_H
#define T3MSCONTINIOUSPROCESSIMPL_H

#include <cmath>
#include "T3Defs.h"
#include "T3ParticleTable.h"
#include "T3LorentzVector.h"

namespace t3 {

using namespace units;  
template <typename Floating>
class T3MSContiniousProcess
{
public:
  T3MSContiniousProcess()
  {
    ParticleTable aParticleTable;
    protonPDG=PDG_t(2212);
    deuteronPDG=aParticleTable.makePDGfromZandA(1,2);
    titanPDG=aParticleTable.makePDGfromZandA(22,48);//Ti48 PDG.
    mh=aParticleTable.GetMass(protonPDG);
    md=aParticleTable.GetMass(deuteronPDG);
    mtitan=aParticleTable.GetMass(titanPDG);
//--------------------------------------------------------------//
    MaterialTable aMaterialTable;
    ntid2=aMaterialTable.GetConcentrations(MatID_t(2));
//--------------------------------------------------------------//
    const Floating Tls=10.0 * MeV;
    MatID_t tid2=MatID_t(2);
    const int zd=aParticleTable.GetZ(deuteronPDG);
    sigma_D_TiD2_10mev_0_tmin=0.0;
    for(int i=0; i<aMaterialTable.GetNumberOfIsotopes(tid2); ++i)
    {
      PDG_t targetPDG = aMaterialTable.GetIsotopes(tid2, i);
      const Floating M=aParticleTable.GetMass(targetPDG);     //M
      const int      Z=aParticleTable.GetZ(targetPDG);        //Z
      Floating tmin;
      if(targetPDG==deuteronPDG)   tmin=2*M*Edisplace_deuteron;
      else if(targetPDG==titanPDG) tmin=2*M*Edisplace_titan;
      sigma_D_TiD2_10mev_0_tmin += aMaterialTable.GetFractions(tid2, i) *
        CalculateMSSigma_Tls_0_tmin(md, M, zd, Z, tmin, Tls, aParticleTable);
    }
  }

#ifdef OPENACC  
#pragma acc routine seq
#endif 
  Floating CalculateMSSigma_Tls_0_tmin(Floating m/*m of the incident particle*/, Floating M/*m of the target particle*/,
                                       int z/*Z of the incident particle*/, int Z/*Z of the target particle*/,
                                       Floating tmin,
                                       /*LS kinetic energy of the inc particle*/
                                       Floating Tls, ParticleTable & aParticleTable) const;

#ifdef OPENACC  
#pragma acc routine seq
#endif
  Floating CalculateSigma_Tls_tmin_tmax(Floating m, Floating M, int z, int Z, Floating tmin,
                                        /*LS kinetic energy of the inc particle*/
                                        Floating Tls, ParticleTable & aParticleTable) const;  

#ifdef OPENACC  
#pragma acc routine seq
#endif
  Floating MSSharpApproximationFunction(PDG_t incPDG,
                                        Floating Tls,/*LS kinetic energy of the inc particle*/
                                        Floating theta,
                                        MatID_t matID,
                                        ParticleTable & aParticleTable) const;

#ifdef OPENACC  
#pragma acc routine seq
#endif
  Floating GetMSStandardRandomAngle1GaussoideFromPDG(
                                            /*the width of the scattering medium, crossed by the
                                              incident particle in inner units(!!!)*/
                                            Floating x,
                                            PDG_t incPDG, Floating Tls/*=p.E()-m*/,
                                            unsigned int & generator_seed,
                                            MatID_t matID,
                                            ParticleTable & aParticleTable) const ;
  

#ifdef OPENACC  
#pragma acc routine seq
#endif
  Floating GetScaledX(PDG_t incPDG,
                      Floating Tls,/*LS kinetic energy of the inc particle*/
                      MatID_t matID,
                      ParticleTable & aParticleTable,
                      MaterialTable & aMaterialTable) const;
  
  
#ifdef OPENACC  
#pragma acc routine seq
#endif
  Floating GetMSRandomAngleOnMeanFreePath(PDG_t incPDG, unsigned int & generator_seed, Floating Tls, MatID_t matID,
                                          ParticleTable & aParticleTable, MaterialTable & aMaterialTable) const;

#ifdef OPENACC  
#pragma acc routine seq
#endif
  void MakeMSRandomAngleScatteringOnParticle4Momentum(PDG_t incPDG,
                                                      /*incident particle 4-momentum vector*/
                                                      T3LorentzVector<Floating> & p,
                                                      unsigned int & generator_seed,
                                                      MatID_t matID,
                                                      ParticleTable & aParticleTable,
                                                      MaterialTable & aMaterialTable);

private:
#ifdef OPENACC  
#pragma acc routine seq
#endif
  Floating W1(Floating z/*the width of the scattering medium, crossed by the incident particle, in mkm*/) const
  {
    Floating t1=33390.0*(1.0+pow(z/17.1, 0.88));
    Floating t2=pow(z, 0.12);
    Floating t3=1.0+pow(z/0.636, 1.16);
    Floating t4=1.0+pow(z/6.128, 2.05);
    Floating w1=t1/t2/t3/t4;
    return w1;
  }

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating W2(Floating z) const
  {
    Floating t1=22979.0/(1.0+0.27/z);
    Floating t2=(1.0+0.0717/z)*(1.0+pow(z/9.628, 0.032));
    Floating w2=t1*t2;
    return w2;
  }

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating W3(Floating z) const
  {
    Floating t1=pow(z, 0.09);
    Floating t2=1.0+pow(0.8884/z, 1.08);
    Floating w3=21859.0/t1/t2;
    return w3;
  }
//---------------------------------------------------------------------------------------------------------//  
#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating D1(Floating z) const
  { 
    Floating t1=0.59e-5*z;
    Floating t2=1.0+0.029/z;
    Floating d1=t1/t2;
    return d1;
  }

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating D2(Floating z) const
  { 
    Floating t1=0.854e-6*(1.0+pow(z/0.25, 1.108));
    Floating d2=D1(z)+t1;
    return d2;
  }

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating D3(Floating z) const
  { 
    Floating t1=0.8e-5*(1.0+pow(z/1.28, 0.77));
    Floating t2=1.0+pow(z/132.37, 4.64);
    Floating d3=D2(z)+t1*t2;
    return d3;
  }
  PDG_t protonPDG;
  PDG_t deuteronPDG;
  PDG_t titanPDG;
  Floating mh;
  Floating md;
  Floating mtitan;
//--------------------------------------------------------------//
  const Floating alpha=1.0/137;
  const Floating hc=200.0 * MeV * fm;
  const Floating Edisplace_deuteron = 10 * eV;
  const Floating Edisplace_titan    = 25 * eV;
  const Floating CTF=0.5*pow(3.0*M_PI/4.0, 2.0/3.0);
  const Floating me=0.511 * MeV;
  const Floating e=sqrt(alpha);
  const Floating a0=1.0/(me*e*e);
  const Floating Lpower=2.0/3.0;
//--------------------------------------------------------------//
  Floating ntid2;
  Floating sigma_D_TiD2_10mev_0_tmin;
//--------------------------------------------------------------//
};

template <typename Floating>
Floating T3MSContiniousProcess<Floating>::CalculateMSSigma_Tls_0_tmin(Floating m/*m of the incident particle*/, Floating M/*m of the target particle*/,
                                                                      int z/*Z of the incident particle*/, int Z/*Z of the target particle*/,
                                                                      Floating tmin,
                                                                      Floating Tls, ParticleTable & aParticleTable) const
{
  const Floating mr=m*M/(m+M);
  const Floating pls2=Tls*(2*m+Tls);
  const Floating rat=m/M;
  const Floating rat2=rat*rat;
  const Floating pcm2=pls2/(1.0+rat2+2*sqrt(rat2+pls2/M/M));
  const Floating pcm=sqrt(pcm2);                            
  const Floating tmax=4*pcm2;                               
  const Floating coef=mr*alpha*z*Z*hc/2/pcm2;
  const Floating COEF=coef*coef*M_PI/pcm2;
  const Floating aI=CTF*a0/sqrt( pow(z, Lpower) + pow(Z, Lpower) );
  const Floating coef1=1.0/(2*pcm*aI);//h=1
  const Floating AS=coef1*coef1*1.13;
  const Floating cs=COEF*tmin/AS/(tmin/tmax+AS);
  return cs;
}

template <typename Floating> 
Floating T3MSContiniousProcess<Floating>::CalculateSigma_Tls_tmin_tmax(Floating m, Floating M, int z, int Z, Floating tmin,
                                                                       Floating Tls, ParticleTable & aParticleTable) const
{
  const Floating mr=m*M/(m+M);
  const Floating pls2=Tls*(2*m+Tls);
  const Floating rat=m/M;
  const Floating rat2=rat*rat;
  const Floating pcm2=pls2/(1.0+rat2+2*sqrt(rat2+pls2/M/M));
  const Floating tmax=4*pcm2;                               
  const Floating coef=2*mr*alpha*z*Z*hc;
  const Floating COEF=coef*coef*M_PI/pcm2;
  Floating cs;
  if(m==M && z==Z) cs=COEF*(tmax-2*tmin)/(tmax-tmin)/tmin;
  else             cs=COEF*(tmax-tmin)/tmax/tmin;
  return cs;
}
 
template <typename Floating>
Floating T3MSContiniousProcess<Floating>::MSSharpApproximationFunction(PDG_t incPDG,
                                                                       Floating Tls,
                                                                       Floating theta,
                                                                       MatID_t matID,
                                                                       ParticleTable & aParticleTable) const
{
  if(matID!=MatID_t(2))
  {
    printf("***ERROR: T3MSContiniousProcess::MSSharpApproximationFunction() Wrong matID=%d! The material is not TiD2!", matID);
    return 0.0;
  }

  const Floating x10=GetScaledX(incPDG, Tls, matID, aParticleTable) / units::um;//scaled x
  
  Floating d1=D1(x10);
  Floating w1=W1(x10);
  Floating t1=2*w1/d1*std::exp(-theta*theta/d1);
  Floating d2=D2(x10);
  Floating w2=W2(x10);
  Floating t2=2*w2/d2*std::exp(-theta*theta/d2);
  Floating d3=D3(x10);
  Floating w3=W3(x10);
  Floating t3=2*w3/d3*std::exp(-theta*theta/d3);
  const Floating app=(t1+t2+t3)/(w1+w2+w3);
  return app;
}

template <typename Floating>
Floating T3MSContiniousProcess<Floating>::GetMSStandardRandomAngle1GaussoideFromPDG(
                                               Floating x,
                                               PDG_t incPDG, Floating Tls/*LS kinetic energy of the incident particle*/,
                                               unsigned int & generator_seed,
                                               MatID_t matID,
                                               ParticleTable & aParticleTable) const
{ 
  const int z_inc=aParticleTable.GetZ(incPDG);                         
  const Floating m=aParticleTable.GetMass(incPDG);                     
  const Floating plsinc2=Tls*(2*m+Tls);                                
  const Floating Elsinc=m+Tls;                                         
  const Floating X0=4.432 * cm;                                        
  const Floating Theta0=13.6 * MeV * Elsinc / plsinc2 * z_inc *        
    sqrt(x/X0) * (1.0 + 0.038 * log(x/X0));                  
  const Floating R=RND01(generator_seed);                              
  const Floating thetams=Theta0*sqrt(-log(R));                         
  return thetams;
}

template <typename Floating>
Floating T3MSContiniousProcess<Floating>::GetScaledX(PDG_t incPDG,
                                                     Floating Tls,/*LS kinetic energy of the inc particle*/
                                                     MatID_t matID,
                                                     ParticleTable & aParticleTable,
                                                     MaterialTable & aMaterialTable) const
{
  
  const Floating m=aParticleTable.GetMass(incPDG);          //m
  const int      z=aParticleTable.GetZ(incPDG);             //z
  Floating sigma_Tls_0_tmin=0.0;
  Floating sigma_Tls_tmin_tmax=0.0;
  for(int i=0; i<aMaterialTable.GetNumberOfIsotopes(matID); ++i)
  {
    const PDG_t targetPDG = aMaterialTable.GetIsotopes(matID, i);
    const Floating M=aParticleTable.GetMass(targetPDG);     //M
    const int Z=aParticleTable.GetZ(targetPDG);             //Z
    Floating tmin;
    if(targetPDG==deuteronPDG)   tmin=2*M*Edisplace_deuteron;
    else if(targetPDG==titanPDG) tmin=2*M*Edisplace_titan;
    
    sigma_Tls_0_tmin += aMaterialTable.GetFractions(matID, i) *
      CalculateMSSigma_Tls_0_tmin(m, M, z, Z, tmin, Tls, aParticleTable);
    
    sigma_Tls_tmin_tmax += aMaterialTable.GetFractions(matID, i) *
      CalculateSigma_Tls_tmin_tmax(m, M, z, Z, tmin, Tls, aParticleTable);
  }

  const Floating N=sigma_Tls_0_tmin/sigma_Tls_tmin_tmax;
  const Floating x10=N/ntid2/sigma_D_TiD2_10mev_0_tmin;
  return x10;
} 

template <typename Floating>
Floating T3MSContiniousProcess<Floating>::GetMSRandomAngleOnMeanFreePath(PDG_t incPDG, unsigned int & generator_seed, Floating Tls, MatID_t matID,
                                                                         ParticleTable & aParticleTable, MaterialTable & aMaterialTable) const
{
  if(matID!=MatID_t(2))
  {
    printf("***ERROR: T3MSContiniousProcess::GetMSRandomAngle() Wrong matID=%d! The material is not TiD2!", matID);
    return 0.0;
  }
  const Floating x10=GetScaledX(incPDG, Tls, matID, aParticleTable, aMaterialTable) / units::um;//scaled x 
  const Floating w1=W1(x10);
  const Floating w2=W2(x10);
  const Floating w3=W3(x10);
  const Floating Wtot=w1+w2+w3;
  const Floating d1=D1(x10);
  const Floating d2=D2(x10);
  const Floating d3=D3(x10);
  const Floating R1=RND01(generator_seed);
  Floating R2=RND01(generator_seed);
  Floating thetams=0.0;
  if(R2<1.0e-10) R2+=1.0e-10;
  if(R1<=w1/Wtot)           thetams=sqrt(-d1*log(R2));
  else if(R1<=(w1+w2)/Wtot) thetams=sqrt(-d2*log(R2));
  else                      thetams=sqrt(-d3*log(R2));
  return thetams;
}

template <typename Floating>
void T3MSContiniousProcess<Floating>::MakeMSRandomAngleScatteringOnParticle4Momentum(PDG_t incPDG,
                                                                                     /*incident particle 4-momentum vector*/
                                                                                     T3LorentzVector<Floating> & p,
                                                                                     unsigned int & generator_seed,
                                                                                     MatID_t matID,
                                                                                     ParticleTable & aParticleTable,
                                                                                     MaterialTable & aMaterialTable)
{
  const T3ThreeVector<Floating> Pls=p.Vect();
  const Floating pls=Pls.R();
  const Floating m=aParticleTable.GetMass(incPDG);
  const Floating Tls=p.E()-m;
  const Floating xy=Pls.x()*Pls.y(), xz=Pls.x()*Pls.z(), yz=Pls.y()*Pls.z();
  const Floating x2=Pls.x()*Pls.x(), y2=Pls.y()*Pls.y(), z2=Pls.z()*Pls.z();
  T3ThreeVector<Floating> e1, e2, e3;
  if(Pls.x() < Pls.y())
  {
    if(Pls.x() < Pls.z())
    {
      e2={0., Pls.z(), -Pls.y()};
      e3={y2+z2, -xy, -xz};
    }
    else
    {
      e2={Pls.y(), -Pls.x(), 0.};
      e3={-xz, -yz, y2+x2};
    }
  }
  else
  {
    if(Pls.y() < Pls.z())
    {
      e2={Pls.z(), 0., -Pls.x()};
      e3={xy, -x2-z2, yz};
    }
    else
    {
      e2={Pls.y(), -Pls.x(), 0.};
      e3={-xz, -yz, y2+x2};
    }
  }
  e1=Pls;
  e1.Unit();
  e2.Unit();
  e3.Unit();
  const Floating thetams=GetMSRandomAngleOnMeanFreePath(incPDG, generator_seed, Tls, matID, aParticleTable, aMaterialTable);
  const Floating cosThetaMS=cos(thetams);
  const Floating sinThetaMS=sin(thetams);
  const Floating phi=2*M_PI*RND01(generator_seed);
  const Floating cosPhi=cos(phi), sinPhi=sin(phi);
  const Floating ss=pls*sinThetaMS*sinPhi, sc=pls*sinThetaMS*cosPhi;
  const Floating plsCosThetaMS=pls*cosThetaMS;
  p.SetX(e1.x()*plsCosThetaMS+e2.x()*ss+e3.x()*sc);
  p.SetY(e1.y()*plsCosThetaMS+e2.y()*ss+e3.y()*sc);
  p.SetZ(e1.z()*plsCosThetaMS+e2.z()*ss+e3.z()*sc);
}

}
#endif//T3MSCONTINIOUSPROCESSIMPL_H
