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
#include "T3MultipleScatteringCSImpl.h"
#include "T3MultipleScatteringFSImpl.h"

#include "T3UserSteppingActionInInject_Check_DD_Elastic_ang_distribution.h"

#include "T3InelasticddImpl.h"
#include "T3ElasticStrongIonIonImpl.h"
#include "T3ElasticEMIonIonImpl.h"

#include "T3dEdxProcessImpl.h"
#include "T3MSContiniousProcessImpl.h"

using namespace t3;

#define MIN(a, b) ((a<b)?a:b)

//#define DEBUG

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
//Arrays fo storing cross sections calculated in Propagate.
FloatingType csMultipleScattering[GL];
FloatingType csInelasticdd[GL];
FloatingType csElasticEMIonIon[GL];
FloatingType csElasticStrongIonIon[GL];
//End of arrays fo storing cross sections calculated in Propagate.

//for incident particles after interaction:
t3::PDG_t outPDG1Inelasticdd[GL];
t3::T3LorentzVector<FloatingType> outP1Inelasticdd[GL];
//for secondary particles:
t3::PDG_t outPDG2Inelasticdd[GL];
t3::T3LorentzVector<FloatingType> outP2Inelasticdd[GL];

t3::PDG_t outPDG1ElasticEMIonIon[GL];
t3::T3LorentzVector<FloatingType> outP1ElasticEMIonIon[GL];
t3::PDG_t outPDG1ElasticStrongIonIon[GL];
t3::T3LorentzVector<FloatingType> outP1ElasticStrongIonIon[GL];


t3::PDG_t outPDG2[GL];
t3::T3LorentzVector<FloatingType> outP2[GL];
//array for T3MultipleScatteringFSImpl.h for storing csBorder arrays for each particle.
FloatingType csBorderDataFS[GL];//this is for GetFS() in Reactor().
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
//Number of killed deuterons:
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
//For finding fluence:
const int NINT=1000;
//events in energy interval:
int NEVENTS[NINT]{0};
//sum of all events dl's:
double LF[NINT]{0.0};//here fluence shoul be.
//Tls scale of the histogram:
const double TlsMin=10 * keV;
const double TlsMax=100 * keV;
const double deltaTls=(TlsMax-TlsMin)/NINT;
//Tls axis of the histogram:
double ENINT[NINT]{0.0};
//********************************************************//

int IR1=0;
int IR2=0;

//decides whether to kill or not to kill the particle:
#ifdef OPENACC  
  #pragma acc routine seq
#endif
template<typename Floating>
bool UserKillParticle(PDG_t pdg/*particle PDG code*/, Floating Tls,/*LS kinetic energy if the particle*/
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
                                       Floating loss/*=dE/dx*l*/, Floating _2me/*=2*me*/)
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
  const Floating me;
  const Floating _2me;
  T3dEdxProcess<FloatingType> dedxprocess;
  T3MSContiniousProcess<FloatingType> mscontiniousprocess;
  T3ElasticEMIonIon<FloatingType>     ElasticEMIonIonImpl;
  T3ElasticStrongIonIon<FloatingType> ElasticStrongIonIonImpl;
};

template <typename Floating>
void DataHolder<Floating>::Propagate()
{
//CALCULATE CROSS SECTIONS:
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
#pragma acc parallel loop gang vector copy(LIFE) present(particles,csInelasticdd,csElasticEMIonIon,csElasticStrongIonIon,dedxprocess) copy(lMAX,ntid2,rho_tid2) reduction(+:IR1,IR2,NEUTRON_COUNT,PROTON_COUNT,TRITON_COUNT,He3_COUNT,Ti48_COUNT,KILLED_DEUTERON_COUNT,COUNT_IR3_IN_PROPAGATE)
#else
#pragma omp parallel for simd reduction(+:IR1,IR2,/*NEVENTS,LF,*/NEUTRON_COUNT,PROTON_COUNT,TRITON_COUNT,He3_COUNT,Ti48_COUNT,KILLED_DEUTERON_COUNT,COUNT_IR3_IN_PROPAGATE)
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
      const Floating Tls=En-m;
      const Floating dEdx = dedxprocess.GetdEdxFromTable(Tls, particles[i].pdg,
                                                         material,aParticleTable);

#ifdef DEBUG
      if(i%10==0) printf("i=%d pdg=%d\n", i, particles[i].pdg);
#endif      
      
      
      
      const Floating dEdxfull = dEdx*rho_tid2;
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
          particles[i].rindex=(l_ElasticEMIonIon<l_ElasticStrongIonIon)?1:2;
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
      if(irc==3) ++COUNT_IR3_IN_PROPAGATE;
      if(irc==1) ++IR1;
      if(irc==2) ++IR2;
      
      mscontiniousprocess.MakeMSRandomAngleScatteringOnParticle4Momentum(
                                           particles[i].pdg,
                                           particles[i].p, particles[i].rs,
                                           material, aParticleTable, aMaterialTable);

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
    }//End if ir>0    
  }//End of i-particle loop
  
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

  std::cout<<"NEUTRON_COUNT="<<NEUTRON_COUNT<<" He3_COUNT="<<He3_COUNT
           <<" PROTON_COUNT="<<PROTON_COUNT<<" TRITON_COUNT="<<TRITON_COUNT
           <<" COUNT_IR3_IN_PROPAGATE="<<COUNT_IR3_IN_PROPAGATE
           <<" COUNT_IR3_IN_REACT="<<COUNT_IR3_IN_REACT<<std::endl;

  
}//End of Propagator

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
}//End of Reactor

template <typename Floating>
void DataHolder<Floating>::Compress()
{
  double sdg=SumDGam;
#ifdef OPENACC
#pragma acc parallel loop gang vector reduction(+:sdg) present(particles)
#else
#pragma omp parallel for simd reduction(+:sdg)
#endif
  for(int i=0; i<LIFE; ++i) sdg += particles[i].de;
  SumDGam=sdg;
  dL=LIFE/Nbin;
  DL=dL+1;  
  unsigned int const n = Nbin - LIFE % Nbin;
    
  POSITION0=POSITION1=POSITION2=POSITION3=POSITION23=0;
#ifdef OPENACC
#pragma acc parallel loop gang vector copy(GL1,dL,DL,n,count01,count23,count0,count1,count2,count3,init,fin)
#else
#pragma omp parallel for simd
#pragma distribute_point
#endif
    for(int b=0; b<Nbin; ++b)    
    {
      count01[b]=GL1;
      count23[b]=0;
      count0[b]=GL1;
      count1[b]=0;
      count2[b]=GL1;
      count3[b]=0;
      if(b<n)
      {
        init[b]=b*dL;
        fin[b]=(b+1)*dL;
      }
      else if(b==n)
      {
        init[b]=n*dL;
        fin[b]=n*dL+DL;
      }
      else if(b>n)
      {
        init[b]=n*dL+DL*(b-n);
        fin[b]=n*dL+DL*(b-n+1);
      }
    }

#ifdef OPENACC
#pragma acc parallel loop gang vector copy(count01,count23,init,fin) present(particles,ind23)
#endif
  {
#ifndef OPENACC
#pragma omp parallel for
#endif
    for(int b=0; b<Nbin; ++b)
    {
#ifdef OPENACC
      //#pragma acc loop vector reduction(+:count01,count23)
#else
      //#pragma omp simd reduction(+:count01,count23)
#endif
      for(int i=init[b]; i<fin[b]; ++i)
      {
        if(particles[i].ir<2) ind23[b][count01[b]--]=i;
        else                  ind23[b][count23[b]++]=i;
      }
    }
  }

#ifdef OPENACC
#pragma acc parallel loop gang copy(count0,count1,count2,count3,count01,count23,mini,ii23,init,fin) present(particles,ind01,ind23)
#else
#pragma omp parallel for
#pragma distribute_point
#endif
    for(int b=0; b<Nbin; ++b)    
    {
      ii23[b]=count23[b]-1;
      mini[b]=GL1-count01[b];
      if(count23[b]<mini[b]) mini[b]=count23[b];
      int js=0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+:js)
#endif
      for(int j=0; j<mini[b]; ++j)
        if (ind23[b][ii23[b] - j] > ind23[b][GL1 - j]) ++js;
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
      for(int j=0; j<js; ++j) std::swap(particles[ind23[b][ii23[b]-j]],particles[ind23[b][GL1-j]]);
      
      for(int i=init[b]; i<fin[b]; ++i)
      {
        if     (particles[i].ir==0) ind01[b][count0[b]--]=i;
        else if(particles[i].ir==1) ind01[b][count1[b]++]=i;
        else if(particles[i].ir==2) ind23[b][count2[b]--]=i;
        else                        ind23[b][count3[b]++]=i;
      }
    }
#ifdef OPENACC
#pragma acc parallel loop gang copy(count0,count1,count2,count3,mini,ii1,ii3,GL1) present(particles,ind01,ind23)
#else
#pragma omp parallel for
#pragma distribute_point
#endif
    for(int b=0; b<Nbin; ++b)    
    {
      ii1[b]=count1[b]-1;
      mini[b]=GL1-count0[b];
      if(count1[b]<mini[b]) mini[b]=count1[b];
      int js=0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+:js)
#endif
      for(int j=0; j<mini[b]; ++j)
        if (ind01[b][ii1[b] - j] > ind01[b][GL1 - j]) ++js;
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
      for(int j=0; j<js; ++j) std::swap(particles[ind01[b][ii1[b]-j]],particles[ind01[b][GL1-j]]);
      ii3[b]=count3[b]-1;
      mini[b]=GL1-count2[b];
      if(count3[b]<mini[b]) mini[b]=count3[b];
      js=0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+:js)
#endif
      for(int j=0; j<mini[b]; ++j)
        if (ind23[b][ii3[b] - j] > ind23[b][GL1 - j]) ++js;
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
      for(int j=0; j<js; ++j) std::swap(particles[ind23[b][ii3[b]-j]],particles[ind23[b][GL1-j]]);
    }
  
#ifdef OPENACC
#pragma acc parallel loop gang vector reduction(+:POSITION0,POSITION1,POSITION2,POSITION3,POSITION23)
#else
#pragma omp parallel for simd reduction(+:POSITION0,POSITION1,POSITION2,POSITION3,POSITION23)
#endif
  for(int b=0; b<Nbin; ++b)
  {
    count0[b]=GL1-count0[b];
    count2[b]=GL1-count2[b];
    POSITION0+=count0[b];
    POSITION1+=count1[b];
    POSITION2+=count2[b];
    POSITION3+=count3[b];
    POSITION23+=count23[b];
  }

  auto prevLIFE = LIFE;
  SHIFT = LIFE - POSITION0;

  if(HISTOGRAM) LIFE = LIFE - POSITION0;
  else          LIFE = LIFE + POSITION23 - POSITION0;

  pointer1[0]=pointer2[0]=pointer3[0]=0;
#ifdef OPENACC
#pragma acc serial loop copy(pointer1,pointer2,pointer3)
#endif
  for(int b=0; b<Nbin-1; ++b)
  {
    pointer1[b+1]=pointer1[b]+count1[b];
    pointer2[b+1]=pointer2[b]+count2[b];
    pointer3[b+1]=pointer3[b]+count3[b];
  }

  for(int b=0; b<Nbin; ++b)
  {
#ifdef OPENACC
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
    {
      cudaMemcpy(&arr1[pointer1[b]],&particles[init[b]+count3[b]+count2[b]],count1[b]*sizep,cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr2[pointer2[b]],&particles[init[b]+count3[b]],count2[b]*sizep,cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr3[pointer3[b]],&particles[init[b]],count3[b]*sizep,cudaMemcpyDeviceToDevice);
    }
#else
    memcpy(&arr1[pointer1[b]],&particles[init[b]+count3[b]+count2[b]],count1[b]*sizep);
    memcpy(&arr2[pointer2[b]],&particles[init[b]+count3[b]],count2[b]*sizep);
    memcpy(&arr3[pointer3[b]],&particles[init[b]],count3[b]*sizep);
#endif
  }
  
#ifdef OPENACC
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
  {
    cudaMemcpy(&particles[0],&arr3[0],POSITION3*sizep,cudaMemcpyDeviceToDevice);
    cudaMemcpy(&particles[POSITION3],&arr2[0],POSITION2*sizep,cudaMemcpyDeviceToDevice);
    cudaMemcpy(&particles[POSITION23],&arr1[0],POSITION1*sizep,cudaMemcpyDeviceToDevice);
  }
#else
  memcpy(&particles[0],&arr3[0],POSITION3*sizep);
  memcpy(&particles[POSITION3],&arr2[0],POSITION2*sizep);
  memcpy(&particles[POSITION23],&arr1[0],POSITION1*sizep);
#endif

  if(!HISTOGRAM)
  {
#ifdef OPENACC
#pragma acc parallel loop gang vector copy(SHIFT,MAX_ELEMENT) present(particles)
#else
#pragma omp parallel for simd
#endif
    for(int i=0; i<POSITION3; ++i)
    {
      int i2=i+i;
      particles[i+SHIFT]=particles[i];
      particles[i].rs=static_cast<unsigned int>(MAX_ELEMENT + i2);
      particles[i+SHIFT].rs=static_cast<unsigned int>(MAX_ELEMENT + i2 + 1);
      particles[i].id=MAX_ELEMENT+i2;
      particles[i].ir=3;
      particles[i+SHIFT].id=MAX_ELEMENT+i2+1;
      particles[i+SHIFT].ir=-1;
      particles[i+SHIFT].de=0.0;
      particles[i+SHIFT].wt=1.0;
      particles[i].de=0.0;
      particles[i].wt=1.0;
    }    
    MAX_ELEMENT+=POSITION3+POSITION3;
#ifdef OPENACC
#pragma acc parallel loop gang vector copy(POSITION3,POSITION23,SHIFT,MAX_ELEMENT) present(particles)
#else
#pragma omp parallel for simd
#endif
    for(int i=POSITION3; i<POSITION23; ++i)
    {
      particles[i+SHIFT] = particles[i];
      particles[i+SHIFT].rs=static_cast<unsigned int>(MAX_ELEMENT + i);

      particles[i].ir=2;
      particles[i+SHIFT].id=MAX_ELEMENT+i;
      particles[i+SHIFT].ir=-1;
    }
    MAX_ELEMENT+=POSITION23-POSITION3;
  }//end of (!HISTOGRAM)
}//End of Compressor

template <typename Floating> void DataHolder<Floating>::Inject()
{

  double sg = 0.;
  double wt = 0.;
  double sgCompensation = 0.;
  double wtCompensation = 0.;
#ifdef OPENACC
#pragma acc parallel loop gang vector reduction(+ : sg, wt)
#else
#pragma omp parallel for simd reduction(+ : sg, wt)
#endif
  for (unsigned int j = 0; j < LIFE; ++j)
  {
    auto sgi = particles[j].GetEtot() * static_cast<double>(particles[j].wt)-sgCompensation;
    auto sgtemp = sg + sgi;
    sgCompensation = (sgtemp - sg) - sgi;
    sg = sgtemp;
    auto wti = static_cast<double>(particles[j].wt) - wtCompensation;
    auto wttemp = wt + wti;
    wtCompensation = (wttemp - wt) - wti;
    wt = wttemp;
  }
  GAMMA = sg;
  NoNew = wt;

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

      std::cout<<"##################################################"<<std::endl;
      std::cout<<"push="<<push<<" INJ="<<INJ<<" N="<<N<<std::endl;
      std::cout<<"##################################################"<<std::endl;      

    }
  }
 }//End of Injector

#endif//T3DATAHOLDER_H
