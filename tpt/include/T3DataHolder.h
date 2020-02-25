#pragma once
#ifndef T3DATAHOLDER_H
#define T3DATAHOLDER_H

#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include "unistd.h"
#include <type_traits>

#include "T3Globals.h"
#include "T3RNG.h"
#include "T3Particle.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"
#include "T3ThreeVector.h"
#include "T3Process.h"
#include "T3ProcessImplFromCSandFS.h"
#include "T3Utility.hh"

#include "T3InelasticddImpl.h"
#include "T3ElasticStrongIonIonImpl.h"
#include "T3ElasticEMIonIonImpl.h"

using namespace t3;

#define MIN(a, b) ((a<b)?a:b)


//=======================================================================//
using Inelasticdd_scattering_t =
  t3::Process<t3::T3Inelasticdd<FloatingType>>;

Inelasticdd_scattering_t inelasticddProcess =
      Inelasticdd_scattering_t(typename Inelasticdd_scattering_t::Base_t());

using ElasticEMIonIon_scattering_t =
  t3::Process<t3::T3ElasticEMIonIon<FloatingType>>;
                                                                                  
ElasticEMIonIon_scattering_t ElasticEMIonIonProcess =
      ElasticEMIonIon_scattering_t(typename ElasticEMIonIon_scattering_t::Base_t());

using ElasticStrongIonIon_scattering_t =
  t3::Process<t3::T3ElasticStrongIonIon<FloatingType>>;                                         

ElasticStrongIonIon_scattering_t ElasticStrongIonIonProcess =
      ElasticStrongIonIon_scattering_t(typename ElasticStrongIonIon_scattering_t::Base_t());

//=======================================================================//

Particle<FloatingType> particles[GL] __attribute__((aligned(64)));
int ind01[Nbin][BLt] __attribute__((aligned(64)));
int ind23[Nbin][BLt] __attribute__((aligned(64)));
Particle<FloatingType> arr1[GL] __attribute__((aligned(64)));
Particle<FloatingType> arr2[GL] __attribute__((aligned(64)));
Particle<FloatingType> arr3[GL] __attribute__((aligned(64)));
long int MAX_ELEMENT;
unsigned int SHIFT;
int INJECTED_PARTICLES=0;

unsigned int POSITION3;
unsigned int POSITION2;
unsigned int POSITION1;
unsigned int POSITION0;
unsigned int POSITION23;
int LIFE=0;
unsigned int sizep=sizeof(Particle<FloatingType>);
int push=0;
int over=0;
unsigned int Ntop=0;
double SumDGam=0.;
double GAMMA=0.;
double NoNew=0.0;


decltype(INJECTED_PARTICLES) GetNumOfInjectedParticles(){return INJECTED_PARTICLES;}
decltype(LIFE) GetNumOfAliveParticles(){return LIFE;}
decltype(NoNew) GetNoNew(){return NoNew;}
decltype(SumDGam) GetSumDGam(){return SumDGam;}
decltype(Ntop) GetNtop(){return Ntop;}

//-----------------------------------------------------------------------------//
FloatingType csMultipleScattering[GL];
FloatingType csInelasticdd[GL];
FloatingType csElasticEMIonIon[GL];
FloatingType csElasticStrongIonIon[GL];

t3::PDG_t outPDG1Inelasticdd[GL];
t3::T3LorentzVector<FloatingType> outP1Inelasticdd[GL];
t3::PDG_t outPDG2Inelasticdd[GL];
t3::T3LorentzVector<FloatingType> outP2Inelasticdd[GL];

t3::PDG_t outPDG1ElasticEMIonIon[GL];
t3::T3LorentzVector<FloatingType> outP1ElasticEMIonIon[GL];
t3::PDG_t outPDG1ElasticStrongIonIon[GL];
t3::T3LorentzVector<FloatingType> outP1ElasticStrongIonIon[GL];

t3::PDG_t outPDG2[GL];
t3::T3LorentzVector<FloatingType> outP2[GL];
FloatingType csBorderDataFS[GL];
//-----------------------------------------------------------------------------//

int mini[Nbin];
int count01[Nbin];
int count23[Nbin];
int count0[Nbin];
int count1[Nbin];
int count2[Nbin];
int pos2[Nbin];
int count3[Nbin];
int ii1[Nbin];
int ii3[Nbin];
int ii23[Nbin];

int init[Nbin];
int fin[Nbin];

int pointer1[Nbin];
int pointer2[Nbin];
int pointer3[Nbin];
int dL;
int DL;
int n;
int numbin;


//********************************************************//
long int KILLED_DEUTERON_COUNT=0;
int COUNT_IR3_IN_PROPAGATE=0;
int COUNT_IR3_IN_REACT=0;
//--------------------------------------------------------//
int DEUTERON_COUNT=0;
int PROTON_COUNT=0;
int TRITON_COUNT=0;
int He3_COUNT=0;
long int Ti48_COUNT=0;
int NEUTRON_COUNT=0;
//--------------------------------------------------------//
const int NINT=1000;
int NEVENTS[NINT]{0};
double LF[NINT]{0.0};
const double TlsMin=10 * keV;
const double TlsMax=100 * keV;
const double deltaTls=(TlsMax-TlsMin)/NINT;
double ENINT[NINT]{0.0};
//********************************************************//

#ifdef OPENACC  
  #pragma acc routine seq
#endif
template<typename Floating>
bool UserKillParticle(PDG_t pdg, Floating Tls,
                      int & neutron_count, int & proton_count,
                      int & triton_count, int & he3_count,
                      long int & ti48_count,
                      long int & killed_deuteron_count,
                      const PDG_t & neutron_PDG, const PDG_t & proton_PDG,
                      const PDG_t & triton_PDG, const PDG_t & he3_PDG,
                      const PDG_t & ti48_PDG, const PDG_t & deuteron_PDG)
{
  bool kill=false;
  if(pdg==neutron_PDG)
  {
#ifdef OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif    
    ++neutron_count;
    kill=true;
  }
  else if(pdg==proton_PDG)
  {
#ifdef OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif    
    ++proton_count;
    kill=true;
  }
  else if(pdg==triton_PDG)
  {
#ifdef OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif   
    ++triton_count;
    kill=true;
  }
  else if(pdg==he3_PDG)
  {
#ifdef OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif    
    ++he3_count;
    kill=true;
  }   
  else if(pdg==ti48_PDG)
  {
#ifdef OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif    
    ++ti48_count;
    kill=true;
  }
  else if(pdg==deuteron_PDG)
  {
    if(Tls<Tls_CUT)
    {
#ifdef OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif      
      ++killed_deuteron_count;
      kill=true;
    }
  }
  return kill;
}

#ifdef OPENACC  
  #pragma acc routine seq
#endif
template<typename Floating>
Floating GetEnergyLossWithFluctuations(unsigned int & generator_seed,
                                       T3LorentzVector<Floating> & p,
                                       Floating loss, Floating _2me)
{
  const Floating beta_LS=p.R()/p.E();
  const Floating sigma=sqrt( loss * _2me ) * beta_LS;
  const Floating eloss=
    T3Utility::RandomizeNormalDistribution(generator_seed, 0.0, sigma);
  return eloss;
}

template <typename Floating>
class DataHolder
{
public:
  DataHolder(): NEUTRON_PDG(PDG_t(2112)), PROTON_PDG(PDG_t(2212)),
                TRITON_PDG(aParticleTable.makePDGfromZandA(1,3)),
                DEUTERON_PDG(aParticleTable.makePDGfromZandA(1,2)),
                He3_PDG(aParticleTable.makePDGfromZandA(2,3)),
                Ti48_PDG(aParticleTable.makePDGfromZandA(22,48)),
                //electron mass (electron PDG = 11):
                me(aParticleTable.GetMass(PDG_t(11))),
                _2me(2*me)
  {}
  void Propagate();
  void React();
  void Compress();
  void Inject();  
private:
  ParticleTable aParticleTable;
  MaterialTable aMaterialTable;
  //---PDG's-for-UserKillParticle()-function:---------------//
  const PDG_t NEUTRON_PDG;
  const PDG_t PROTON_PDG;
  const PDG_t TRITON_PDG;
  const PDG_t DEUTERON_PDG;
  const PDG_t He3_PDG;
  const PDG_t Ti48_PDG;
  //-END-of-PDG's-for-UserKillParticle()-function:----------//
  //Electron mass for energy losses fluctuations: 
  const Floating me;
  const Floating _2me;
  
  T3ElasticEMIonIon<FloatingType>     ElasticEMIonIonImpl;
  T3ElasticStrongIonIon<FloatingType> ElasticStrongIonIonImpl;
  
};

template <typename Floating>
void DataHolder<Floating>::Propagate()
{
  const t3::MatID_t material=t3::MatID_t(2u);

  inelasticddProcess.GetCS(particles, material, LIFE,
                           csInelasticdd, aParticleTable, aMaterialTable);

  ElasticEMIonIonProcess.GetCS(particles, material, LIFE,
                               csElasticEMIonIon, aParticleTable, aMaterialTable);

  ElasticStrongIonIonProcess.GetCS(particles, material, LIFE,
                                   csElasticStrongIonIon, aParticleTable, aMaterialTable);
  
  Floating da;
  if(std::is_same<double, Floating>::value)     da=ag * 1.0e-10;
  else if(std::is_same<float, Floating>::value) da=ag * 1.0e-7;
  constexpr Floating lMAX = std::numeric_limits<Floating>::max();
  const Floating ntid2 = aMaterialTable.GetConcentrations(material);
  const Floating rho_tid2 = aMaterialTable.GetDensity(material);
  
#ifdef OPENACC
#pragma acc parallel loop gang vector copy(LIFE) present(particles,csInelasticdd,csElasticEMIonIon,csElasticStrongIonIon) copy(lMAX,ntid2,rho_tid2) reduction(+:NEUTRON_COUNT,PROTON_COUNT,TRITON_COUNT,He3_COUNT,Ti48_COUNT,KILLED_DEUTERON_COUNT,COUNT_IR3_IN_PROPAGATE)
#else
#pragma omp parallel for simd reduction(+:NEVENTS,LF,NEUTRON_COUNT,PROTON_COUNT,TRITON_COUNT,He3_COUNT,Ti48_COUNT,KILLED_DEUTERON_COUNT,COUNT_IR3_IN_PROPAGATE)
#endif
  for(int i=0; i<LIFE; ++i)
  {
    const Floating En = particles[i].GetEtot();
    if(particles[i].ir > 0)
    {
      const Floating csInelasticddi = csInelasticdd[i];
      const Floating csElasticEMIonIoni = csElasticEMIonIon[i];
      const Floating csElasticStrongIonIoni = csElasticStrongIonIon[i];      
      T3ThreeVector<FloatingType> r0(particles[i].r.x(), particles[i].r.y(), particles[i].r.z());
      Floating  l1x =
        (particles[i].vx() == 0.)
            ? lMAX/1.0e3
            : ((particles[i].vx() > 0.)
                     ? ((particles[i].ix() + 1) * ag - particles[i].r.x())
                     : (particles[i].ix() * ag - particles[i].r.x())) /
                        particles[i].vx() +
                    da;
      Floating  l1y =
        (particles[i].vy() == 0.)
            ? lMAX/1.0e3
            : ((particles[i].vy() > 0.)
                     ? ((particles[i].jy() + 1) * ag - particles[i].r.y())
                     : (particles[i].jy() * ag - particles[i].r.y())) /
                        particles[i].vy() +
                     da;
      Floating  l1z =
        (particles[i].vz() == 0.)
            ? lMAX/1.0e3
            : ((particles[i].vz() > 0.)
                     ? ((particles[i].kz() + 1) * ag - particles[i].r.z())
                     : (particles[i].kz() * ag - particles[i].r.z())) /
                        particles[i].vz() +
                     da;

      const Floating m=aParticleTable.GetMass(particles[i].pdg);
      const Floating Tls=particles[i].p.E() - m;
      const Floating dEdxfull = 1000.0 * MeV * cm * cm / gr;
      const Floating lambdar_Inelasticdd =
        (csInelasticddi<1.0e-7*mbarn)?lMAX/1.0e3:1.0/ntid2/csInelasticddi;
      const Floating lambdar_ElasticEMIonIon =
        (csElasticEMIonIoni<1.0e-7*mbarn)?lMAX/1.0e3:1.0/ntid2/csElasticEMIonIoni;
      const Floating lambdar_ElasticStrongIonIon =
        (csElasticStrongIonIoni<1.0e-7*mbarn)?lMAX/1.0e3:1.0/ntid2/csElasticStrongIonIoni;
      Floating l1=MIN(MIN(l1x, l1y), MIN(l1y, l1z));
      Floating R1 = RND01(particles[i].rs);
      if(R1<1.0e-10) R1+=1.0e-10;
      Floating R2 = RND01(particles[i].rs);
      if(R2<1.0e-10) R2+=1.0e-10;
      Floating R3 = RND01(particles[i].rs);
      if(R3<1.0e-10) R3+=1.0e-10;
      const Floating l_ElasticEMIonIon = fabs(lambdar_ElasticEMIonIon * log(R1));
      const Floating l_ElasticStrongIonIon = fabs(lambdar_ElasticStrongIonIon * log(R2));
      const Floating l3 = fabs(lambdar_Inelasticdd * log(R3));
      const Floating l2=MIN(l_ElasticEMIonIon, l_ElasticStrongIonIon);
      int irc=-1;
      Floating l=-1.0;
      particles[i].rindex=-1;
      if(l1 < l2 && l1 < l3)
      {
        irc=1;
        l=l1;
      }
      else
      {
        if(l2<l3)
        {
          irc=2;
          l=l2;
          particles[i].rindex=
            (l_ElasticEMIonIon<l_ElasticStrongIonIon)?1:2;
        }
        else
        {
          irc=3;
          l=l3;
        }
      }
      Floating dl = fabs(l);
      const int BinF=(Tls-dEdxfull*dl/2-TlsMin)/deltaTls;
      if(BinF>=0 && BinF<NINT)
      {
#ifdef OPENACC
      #pragma acc atomic update
#else
      #pragma omp atomic update
#endif        
        ++NEVENTS[BinF];
#ifdef OPENACC
      #pragma acc atomic update
#else
      #pragma omp atomic update
#endif        
        LF[BinF]+=dl;
      }
              
      particles[i].r.SetPxPyPzE(particles[i].r.x()+particles[i].vx() * (dl + da),
                                particles[i].r.y()+particles[i].vy() * (dl + da),
                                particles[i].r.z()+particles[i].vz() * (dl + da), 0.0);//t=0.0
      
      bool const out = (particles[i].ix() >= cuba || particles[i].jy() >= cuba ||
                        particles[i].kz() >= cuba || particles[i].ix() < -cuba ||
                        particles[i].jy() < -cuba || particles[i].kz() < -cuba);

      Floating loss = 0.;
      if(aParticleTable.IsNucleus(particles[i].pdg))
      {
        loss = dEdxfull * dl;
        const Floating deltaE=GetEnergyLossWithFluctuations(particles[i].rs,particles[i].p,loss,_2me);
        loss+=deltaE;
        if (loss >= particles[i].GetEtot() - aParticleTable.GetMass(particles[i].pdg))
        {
          particles[i].de += particles[i].GetEtot() - particles[i].m;
          particles[i].SetEtot(aParticleTable.GetMass(particles[i].pdg), aParticleTable);       
          particles[i].ir = 0;
        }
        else
        {
          const Floating old_energy = particles[i].GetEtot();
          const Floating new_energy = old_energy - loss;
          particles[i].SetEtot(new_energy, aParticleTable);          
          particles[i].de += loss;
          particles[i].ir = irc;
        }
      }
      if(out)
      {
        particles[i].ir = 0;
      }
    }
  }
  
#ifdef OPENACC
#pragma acc parallel loop gang vector present(particles,aParticleTable) reduction(+:KILLED_DEUTERON_COUNT,DEUTERON_COUNT,PROTON_COUNT,TRITON_COUNT,He3_COUNT,Ti48_COUNT)
#else
#pragma omp parallel for simd reduction(+:KILLED_DEUTERON_COUNT,DEUTERON_COUNT,PROTON_COUNT,TRITON_COUNT,He3_COUNT,Ti48_COUNT)
#endif
  for(int i=0; i<LIFE; ++i)
  {
    if(particles[i].pdg==DEUTERON_PDG)
    {
      const Floating tls=particles[i].GetEtot() - aParticleTable.GetMass(particles[i].pdg);
      if(UserKillParticle(particles[i].pdg,tls,
                          NEUTRON_COUNT,PROTON_COUNT,TRITON_COUNT,
                          He3_COUNT,Ti48_COUNT,KILLED_DEUTERON_COUNT,
                          NEUTRON_PDG,PROTON_PDG,TRITON_PDG,
                          He3_PDG,Ti48_PDG,DEUTERON_PDG))
        particles[i].ir=0;
    }    
  }
}

template <typename Floating>
void DataHolder<Floating>::React()
{
  const int CSBORDERSTEP=aMaterialTable.GetcsBorderSTEP();
  const auto material=t3::MatID_t(2);
  
  inelasticddProcess.GetFS(particles, material, POSITION3,
                           outPDG1Inelasticdd, outP1Inelasticdd, outPDG2Inelasticdd, outP2Inelasticdd,
                           csBorderDataFS, aParticleTable, aMaterialTable,
                           CSBORDERSTEP);
  
#ifdef OPENACC
#pragma acc parallel loop gang vector present(particles,aParticleTable)
#else
#pragma omp parallel for simd reduction(+:COUNT_IR3_IN_REACT)
#endif
  for (unsigned int i = 0; i < POSITION3; ++i)
  {
    particles[i].p = outP1Inelasticdd[i];
    particles[i].pdg = outPDG1Inelasticdd[i];
    particles[i].m=aParticleTable.GetMass(particles[i].pdg);
    particles[i+SHIFT].p=outP2Inelasticdd[i];
    particles[i+SHIFT].pdg=outPDG2Inelasticdd[i];
    particles[i+SHIFT].m=aParticleTable.GetMass(particles[i+SHIFT].pdg);

    const Floating Tls1sec=particles[i].p.E() - particles[i].m;
    if(UserKillParticle(particles[i].pdg, Tls1sec,
                        NEUTRON_COUNT,PROTON_COUNT,TRITON_COUNT,
                        He3_COUNT,Ti48_COUNT,KILLED_DEUTERON_COUNT,
                        NEUTRON_PDG,PROTON_PDG,TRITON_PDG,
                        He3_PDG,Ti48_PDG,DEUTERON_PDG))
          particles[i].ir=0;
    else  particles[i].ir=3;

    const Floating Tls2sec=particles[i+SHIFT].p.E() - particles[i+SHIFT].m;
    if(UserKillParticle(particles[i+SHIFT].pdg, Tls2sec,
                        NEUTRON_COUNT,PROTON_COUNT,TRITON_COUNT,
                        He3_COUNT,Ti48_COUNT,KILLED_DEUTERON_COUNT,
                        NEUTRON_PDG,PROTON_PDG,TRITON_PDG,
                        He3_PDG,Ti48_PDG,DEUTERON_PDG))
          particles[i+SHIFT].ir=0;
    else  particles[i+SHIFT].ir=3;

    ++COUNT_IR3_IN_REACT;
    
  }

  
#ifdef OPENACC
#pragma acc parallel loop gang vector present(particles,ElasticEMIonIonImpl,ElasticStrongIonIonImpl,aParticleTable)
#else
#pragma omp parallel for simd
#endif
  for (unsigned int i = POSITION3; i < POSITION23; ++i)
  {   
    if(particles[i].rindex==1)
    {
      Four<FloatingType> four1=ElasticEMIonIonImpl.GetFS(particles[i].p, particles[i].pdg,
                               material, particles[i].rs, csBorderDataFS,
                               particles[i].tr, i,
                               aParticleTable, aMaterialTable,
                               CSBORDERSTEP);
      particles[i].p = four1.P1;
      particles[i].pdg = four1.pdg1;
      particles[i].m=aParticleTable.GetMass(particles[i].pdg);    
      
      particles[i+SHIFT].p=four1.P2;
      particles[i+SHIFT].pdg=four1.pdg2;
      particles[i+SHIFT].m=aParticleTable.GetMass(particles[i+SHIFT].pdg);
    }
    else if(particles[i].rindex==2)
    {
      Four<FloatingType> four2=ElasticStrongIonIonImpl.GetFS(particles[i].p, particles[i].pdg,
                               material, particles[i].rs, csBorderDataFS,
                               particles[i].tr, i,
                               aParticleTable, aMaterialTable,
                               CSBORDERSTEP); 
      particles[i].p = four2.P1;
      particles[i].pdg = four2.pdg1;
      particles[i].m=aParticleTable.GetMass(particles[i].pdg);
      
      particles[i+SHIFT].p=four2.P2;
      particles[i+SHIFT].pdg=four2.pdg2;
      particles[i+SHIFT].m=aParticleTable.GetMass(particles[i+SHIFT].pdg);
    }
    
    const Floating Tls1sec=particles[i].p.E() - particles[i].m;
    if(UserKillParticle(particles[i].pdg, Tls1sec,
                        NEUTRON_COUNT,PROTON_COUNT,TRITON_COUNT,
                        He3_COUNT,Ti48_COUNT,KILLED_DEUTERON_COUNT,
                        NEUTRON_PDG,PROTON_PDG,TRITON_PDG,
                        He3_PDG,Ti48_PDG,DEUTERON_PDG))
          particles[i].ir=0;
    else  particles[i].ir=2;
    
    const Floating Tls2sec=particles[i+SHIFT].p.E() - particles[i+SHIFT].m;
    if(UserKillParticle(particles[i+SHIFT].pdg, Tls2sec,
                        NEUTRON_COUNT,PROTON_COUNT,TRITON_COUNT,
                        He3_COUNT,Ti48_COUNT,KILLED_DEUTERON_COUNT,
                        NEUTRON_PDG,PROTON_PDG,TRITON_PDG,
                        He3_PDG,Ti48_PDG,DEUTERON_PDG))
          particles[i+SHIFT].ir=0;
    else  particles[i+SHIFT].ir=2;
    
  }
}

template <typename Floating>
void DataHolder<Floating>::Compress(){}

template <typename Floating> void DataHolder<Floating>::Inject()
{   
  static int push = 0;
  if (LIFE > Ntop) Ntop = LIFE;
  if (push < Np)
  {
    if(LIFE < K)
    {
      int NUMBER_OF_PARTICLES_TO_PUSH=INJ;
      if(push+INJ>Np) NUMBER_OF_PARTICLES_TO_PUSH=Np-push;
      
#ifdef OPENACC
#pragma acc parallel loop gang vector present(particles,aParticleTable)
#else
#pragma omp parallel for        
#endif
      for(int i=0; i<NUMBER_OF_PARTICLES_TO_PUSH; ++i)
#ifdef PROBLEM
        InitParticle(LIFE+i, MAX_ELEMENT+i, aParticleTable);
#else
      {
        t3::PDG_t initPDG = aParticleTable.makePDGfromZandA(1, 2);
        auto const m = aParticleTable.GetMass(initPDG);
        auto const InitEls= m+TLS;
        const auto pls=sqrt(InitEls*InitEls-m*m);
        particles[LIFE+i] = Particle<Floating>(
        t3::T3LorentzVector<Floating>(InitParticlex0*ag,
                                      InitParticley0*ag,
                                      -cuba*ag,0.0),
        t3::T3LorentzVector<Floating>(0.0,0.0,pls,InitEls),
        m, 0.0, initPDG, 1., MAX_ELEMENT+i, 1, MAX_ELEMENT+i, -1.0, -1);
      }
#endif
      push+=NUMBER_OF_PARTICLES_TO_PUSH;
      LIFE+=NUMBER_OF_PARTICLES_TO_PUSH;
      MAX_ELEMENT+=NUMBER_OF_PARTICLES_TO_PUSH;
      INJECTED_PARTICLES+=NUMBER_OF_PARTICLES_TO_PUSH;
    }
  }
 }

#endif//T3DATAHOLDER_H
