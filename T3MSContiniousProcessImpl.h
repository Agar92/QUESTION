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
    titanPDG=aParticleTable.makePDGfromZandA(22,48);
    mh=aParticleTable.GetMass(protonPDG);
    md=aParticleTable.GetMass(deuteronPDG);
    mtitan=aParticleTable.GetMass(titanPDG);
    const Floating Tlsd10MeV =10.0 * MeV;
    const Floating Elsd10MeV =md+Tlsd10MeV;
    const Floating Elsd10MeV2=Elsd10MeV * Elsd10MeV;
    const Floating plsd10mev2=Tlsd10MeV*(2*md+Tlsd10MeV);
    const Floating plsd10mev4=plsd10mev2 * plsd10mev2;
    SCALE_COEFFICIENT_FROM_D =ZD * ZD * Elsd10MeV2 / plsd10mev4;
  }

#ifdef OPENACC  
#pragma acc routine seq
#endif
  Floating GetMSStandardRandomAngle1GaussoideFromPDG(
                                            Floating x,
                                            PDG_t incPDG, Floating Elsfull/*=p.E()*/,
                                            unsigned int & generator_seed,
                                            ParticleTable aParticleTable,
                                            MatID_t matID);
  

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating GetMSRandomAngle(Floating z,
                            PDG_t incPDG, T3LorentzVector<Floating> & p,
                            unsigned int & generator_seed,
                            ParticleTable aParticleTable,
                            MatID_t matID) const;

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating MSSharpApproximationFunction(Floating z,
                                        Floating theta,
                                        PDG_t incPDG, T3LorentzVector<Floating> & p,
                                        ParticleTable aParticleTable,
                                        MatID_t matID) const;

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  void MakeMSRandomAngleScatteringOnParticle4Momentum(PDG_t incPDG, Floating x,
                                                      T3LorentzVector<Floating> & p,
                                                      unsigned int & generator_seed, ParticleTable aParticleTable,
                                                      MatID_t matI);

private:
//member functions:

#ifdef OPENACC  
#pragma acc routine seq
#endif
  Floating W1(Floating z) const
  {
    Floating t1=33390.0*(1.0+std::pow(z/17.1, 0.88));
    Floating t2=std::pow(z, 0.12);
    Floating t3=1.0+std::pow(z/0.636, 1.16);
    Floating t4=1.0+std::pow(z/6.128, 2.05);
    Floating w1=t1/t2/t3/t4;
    return w1;
  }

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating W2(Floating z) const
  {
    Floating t1=22979.0/(1+0.27/z);
    Floating t2=(1.0+0.0717/z)*(1.0+std::pow(z/9.628, 0.032));
    Floating w2=t1*t2;
    return w2;
  }

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating W3(Floating z) const
  {
    Floating t1=std::pow(z, 0.09);
    Floating t2=1.0+std::pow(0.8884/z, 1.08);
    Floating w3=21859.0/t1/t2;
    return w3;
  }

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating D1DApp(Floating z) const
  {
    Floating t1=0.59e-5*z;
    Floating t2=1.0+0.029/z;
    Floating d1=t1/t2;
    return d1;
  }

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating D2DApp(Floating z) const
  {
    Floating t1=0.854e-6*(1.0+std::pow(z/0.25, 1.108));
    Floating d2=D1DApp(z)+t1;
    return d2;
  }

#ifdef OPENACC  
#pragma acc routine seq
#endif  
  Floating D3DApp(Floating z) const
  {
    Floating t1=0.8e-5*(1.0+std::pow(z/1.28, 0.77));
    Floating t2=1.0+std::pow(z/132.37, 4.64);
    Floating d3=D2DApp(z)+t1*t2;
    return d3;
  }
  //member variables:
  PDG_t protonPDG;
  PDG_t deuteronPDG;
  PDG_t titanPDG;
  Floating mh;
  Floating md;
  Floating mtitan;
//--------------------------------------------------------------//
  Floating SCALE_COEFFICIENT_FROM_D;
//--------------------------------------------------------------//  
  const int ZH=1;
  const int ZD=1;
  const int ZTi=22;
//--------------------------------------------------------------//
};

template <typename Floating>
Floating T3MSContiniousProcess<Floating>::GetMSStandardRandomAngle1GaussoideFromPDG(
                                               Floating x,
                                               PDG_t incPDG, Floating Elsinc/*=p.E()*/,
                                               unsigned int & generator_seed,
                                               ParticleTable aParticleTable,
                                               MatID_t matID)
{ 
  const int Zinc=aParticleTable.GetZ(incPDG);
  const Floating m=aParticleTable.GetMass(incPDG);
  const Floating Tls=Elsinc-m;
  const Floating plsinc2=Tls*(2*m+Tls);
  const Floating X0=4.432 * cm;
  const Floating Theta0=13.6 * MeV * Elsinc / plsinc2 * Zinc *
    std::sqrt(x/X0) * (1.0 + 0.038 * std::log(x/X0));
  const Floating R=RND01(generator_seed);
  const Floating thetams=Theta0*std::sqrt(-std::log(R));
  return thetams;
}

template <typename Floating>
Floating T3MSContiniousProcess<Floating>::GetMSRandomAngle(Floating z,
                                                           PDG_t incPDG, T3LorentzVector<Floating> & p,
                                                           unsigned int & generator_seed,
                                                           ParticleTable aParticleTable,
                                                           MatID_t matID) const
{
  if(matID!=MatID_t(2))
  {
    printf("***ERROR: T3MSContiniousProcess::GetMSRandomAngle() Wrong matID=%d!", matID);
    return 0.0;
  }
  const int Zinc=aParticleTable.GetZ(incPDG);
  const Floating Elsinc =p.E();
  const Floating Elsinc2=Elsinc*Elsinc;
  const Floating plsinc2=p.x()*p.x()+p.y()*p.y()+p.z()*p.z();
  const Floating plsinc4=plsinc2*plsinc2;
  const Floating COEF   =Zinc*Zinc*Elsinc2/plsinc4;
  const Floating w1=W1(z);
  const Floating w2=W2(z);
  const Floating w3=W3(z);
  const Floating Wtot=w1+w2+w3;
  const Floating D1scaled=D1DApp(z) / SCALE_COEFFICIENT_FROM_D * COEF;
  const Floating D2scaled=D2DApp(z) / SCALE_COEFFICIENT_FROM_D * COEF;
  const Floating D3scaled=D3DApp(z) / SCALE_COEFFICIENT_FROM_D * COEF;
  const Floating R1=RND01(generator_seed);
  const Floating R2=RND01(generator_seed);
  Floating thetams=0.0;
  if(R1<w1/Wtot)           thetams=std::sqrt(-D1scaled*std::log(R2));
  else if(R1<(w1+w2)/Wtot) thetams=std::sqrt(-D2scaled*std::log(R2));
  else                     thetams=std::sqrt(-D3scaled*std::log(R2));
  return thetams;
}

template <typename Floating>
Floating T3MSContiniousProcess<Floating>::MSSharpApproximationFunction(Floating z,
                                                                       Floating theta,
                                                                       PDG_t incPDG, T3LorentzVector<Floating> & p,
                                                                       ParticleTable aParticleTable,
                                                                       MatID_t matID) const
{
  if(matID!=MatID_t(2))
  {
    printf("***ERROR: T3MSContiniousProcess::MSSharpApproximationFunction() Wrong matID=%d!", matID);
    return 0.0;
  }
  const int Zinc=aParticleTable.GetZ(incPDG);
  const Floating Elsinc =p.E();
  const Floating Elsinc2=Elsinc*Elsinc;
  const Floating plsinc2=p.x()*p.x()+p.y()*p.y()+p.z()*p.z();
  const Floating plsinc4=plsinc2*plsinc2;
  const Floating COEF   =Zinc*Zinc*Elsinc2/plsinc4;
  Floating d1=D1DApp(z) / SCALE_COEFFICIENT_FROM_D * COEF;
  Floating w1=W1(z);
  Floating t1=2*w1/d1*std::exp(-theta*theta/d1);
  Floating d2=D2DApp(z) / SCALE_COEFFICIENT_FROM_D * COEF;
  Floating w2=W2(z);
  Floating t2=2*w2/d2*std::exp(-theta*theta/d2);
  Floating d3=D3DApp(z) / SCALE_COEFFICIENT_FROM_D * COEF;
  Floating w3=W3(z);
  Floating t3=2*w3/d3*std::exp(-theta*theta/d3);
  const Floating app=(t1+t2+t3)/(w1+w2+w3);
  return app;
}

template <typename Floating>
void T3MSContiniousProcess<Floating>::MakeMSRandomAngleScatteringOnParticle4Momentum(PDG_t incPDG, Floating x,
                                                                                     T3LorentzVector<Floating> & p,
                                                                                     unsigned int & generator_seed, ParticleTable aParticleTable,
                                                                                     MatID_t matID)
{
  const T3ThreeVector<Floating> Pls=p.Vect();
  const Floating pls=Pls.R();
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
  const Floating xum=x/units::um;
  const Floating thetams=GetMSRandomAngle(xum, incPDG, p, generator_seed,
                                          aParticleTable, matID);
  const Floating cosThetaMS=std::cos(thetams);
  const Floating sinThetaMS=std::sin(thetams);
  const Floating phi=2*M_PI*RND01(generator_seed);
  const Floating cosPhi=std::cos(phi), sinPhi=std::sin(phi);
  const Floating ss=pls*sinThetaMS*sinPhi, sc=pls*sinThetaMS*cosPhi;
  const Floating plsCosThetaCM=pls*cosThetaMS;
  p.SetX(e1.x()*plsCosThetaCM+e2.x()*ss+e3.x()*sc);
  p.SetY(e1.y()*plsCosThetaCM+e2.y()*ss+e3.y()*sc);
  p.SetZ(e1.z()*plsCosThetaCM+e2.z()*ss+e3.z()*sc);
}
  

}
#endif//T3MSCONTINIOUSPROCESSIMPL_H
