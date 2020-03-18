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
#include "T3Particle.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"
#include "T3ThreeVector.h"
#include "T3Process.h"


#include "T3UserSteppingActionInInject.h"
#include "T3InelasticddImpl.h"
#include "T3AllocateData.h"


using namespace data;
using namespace t3;

#define MIN(a, b) ((a<b)?a:b)

//=======================================================================//
//For T3InelastiddImpl.h:
using Inelasticdd_scattering_t =
  t3::Process<t3::Inelasticdd<FloatingType>>;

Inelasticdd_scattering_t inelasticddProcess =
      Inelasticdd_scattering_t(typename Inelasticdd_scattering_t::Base_t());
//=======================================================================//

//-----------------------------------------------------------------------------//

template <typename Floating>
class DataHolder
{
public:
  DataHolder(){}
  void Propagate();
  void React();
  void Compress();
  void Inject();
  void RegisterUserSteppingActionInInject()
  {
    actionInInject.Register();
  }
  void Histogram_theta_UserSteppingActionInInject()
  {
    actionInInject.Histogram_theta_UserSteppingActionInInject();
  }
private:
  ParticleTable aParticleTable;
  MaterialTable aMaterialTable;
  UserSteppingActionInInject    actionInInject;
};

template <typename Floating>
void DataHolder<Floating>::Propagate()
{  
  const MatID_t material=t3::MatID_t(2u);

  inelasticddProcess.GetCS(material, LIFE, csInelasticdd, aParticleTable, aMaterialTable);

  Floating da;
  if(std::is_same<double, Floating>::value)     da=ag * 1.0e-10;
  else if(std::is_same<float, Floating>::value) da=ag * 1.0e-7;
  
  constexpr Floating rho0 = 1.0;
  constexpr auto lMAX = std::numeric_limits<Floating>::max();
  Floating const ntid2 = aMaterialTable.GetConcentrations(t3::MatID_t(2u));
  Floating const rho_tid2 = aMaterialTable.GetDensity(t3::MatID_t(2u));
  Floating DE=0.0;
  
#ifdef OPENACC
#pragma acc parallel loop gang vector copy(LIFE) present(particles,csMultipleScattering) copy(lMAX,ntid2,rho_tid2)
#else
#pragma omp parallel for simd
#endif
  for(int i=0; i<LIFE; ++i)
  {
    const Floating En = particles[i].GetEtot();    
    if(particles[i].ir > 0)
    {
      const Floating csInelasticddi = csInelasticdd[i];      
      T3ThreeVector<FloatingType> r0(particles[i].r.x(), particles[i].r.y(), particles[i].r.z());
      //l1 l2 l3
      Floating  l1x =
        (particles[i].vx() == 0.)
            ? lMAX
            : ((particles[i].vx() > 0.)
                     ? ((particles[i].ix() + 1) * ag - particles[i].r.x())
                     : (particles[i].ix() * ag - particles[i].r.x())) /
                        particles[i].vx() +
                    da;
      Floating  l1y =
        (particles[i].vy() == 0.)
            ? lMAX
            : ((particles[i].vy() > 0.)
                     ? ((particles[i].jy() + 1) * ag - particles[i].r.y())
                     : (particles[i].jy() * ag - particles[i].r.y())) /
                        particles[i].vy() +
                     da;
      Floating  l1z =
        (particles[i].vz() == 0.)
            ? lMAX
            : ((particles[i].vz() > 0.)
                     ? ((particles[i].kz() + 1) * ag - particles[i].r.z())
                     : (particles[i].kz() * ag - particles[i].r.z())) /
                        particles[i].vz() +
                     da;

      Floating const dEdx = 2.0 * MeV * cm * cm / gr;
      Floating const dEdxfull = dEdx*rho_tid2;
      Floating const range = (particles[i].GetEtot() - particles[i].m)/dEdxfull;
      Floating l0 = range;
      const Floating CSInelasticddi = std::numeric_limits<FloatingType>::max();
      Floating const lambdar = 1.0/ntid2/CSInelasticddi; 
      Floating R = particles[i].GenerateCanonical();
      if(R<1.0e-10) R+=1.0e-10;

      Floating l2 = fabs(lambdar * log(R));
      Floating l1=MIN(MIN(l1x, l1y), MIN(l1y, l1z));
      int irc = (l0 < l2 && l0 < l1) ? 0 : ((l2 < l1) ? 2 : 1);
      Floating l = MIN(MIN(l1, l2), MIN(l2, l0));
      Floating dl = fabs(l);

      //std::cout<<"irc="<<irc<<std::endl;
      
      particles[i].r.SetPxPyPzE(particles[i].r.x()+particles[i].vx() * (dl + da),
                                particles[i].r.y()+particles[i].vy() * (dl + da),
                                particles[i].r.z()+particles[i].vz() * (dl + da), 0.0);//t=0.0

      bool const out = (particles[i].ix() >= cuba || particles[i].jy() >= cuba ||
                        particles[i].kz() >= cuba || particles[i].ix() < -cuba ||
                        particles[i].jy() < -cuba || particles[i].kz() < -cuba);
      
      Floating loss = 0.;
      if (aParticleTable.IsNucleus(particles[i].pdg))
      {
        loss = dEdxfull * dl;
        if (loss >= particles[i].GetEtot() - particles[i].m)
        {
          particles[i].de += particles[i].GetEtot() - particles[i].m;
          particles[i].SetEtot(aParticleTable.GetMass(particles[i].pdg), aParticleTable);       
          particles[i].ir = 0;
        }
        else
        {
          Floating const old_energy = particles[i].GetEtot();
          Floating const new_energy = old_energy - loss;
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
}//End of Propagator

template <typename Floating>
void DataHolder<Floating>::React()
{
  const int CSBORDERSTEP=aMaterialTable.GetcsBorderSTEP();
  const auto material=t3::MatID_t(2);

  std::cout<<"P23="<<POSITION23<<" SHIFT="<<SHIFT<<std::endl;
  
  inelasticddProcess.GetFS(material, POSITION23, SHIFT,
                           aParticleTable, aMaterialTable);
    
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
  for(int i=0; i<LIFE; ++i)
  {
    sdg += particles[i].de;
  }
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
      //#pragma acc loop vector reduction(+:count01,count23) ///???????????? vector//
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
      particles[i+SHIFT] = particles[i];
      particles[i].rs=static_cast<unsigned int>(MAX_ELEMENT + i2);
      particles[i+SHIFT].rs=static_cast<unsigned int>(MAX_ELEMENT + i2 + 1);

      particles[i].id=MAX_ELEMENT+i2;
      particles[i].ir=-1;
      particles[i+SHIFT].id=MAX_ELEMENT+i2+1;
      particles[i+SHIFT].ir=-1;

      if (particles[i].pdg == aParticleTable.makePDGfromZandA(1,2))
      {
        particles[i].pdg = aParticleTable.makePDGfromZandA(1,2);
        particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(1,2);
      }
      else if (particles[i].pdg == aParticleTable.makePDGfromZandA(22,48))
      {
        particles[i].pdg = aParticleTable.makePDGfromZandA(22,48);
        particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(22,48);
      }
      else if (particles[i].pdg == t3::PDG_t(22))
      {
        particles[i].pdg = t3::PDG_t(11);
        particles[i+SHIFT].pdg = t3::PDG_t(-11);
      } else if (particles[i].pdg == t3::PDG_t(11))
      {
        particles[i].pdg = t3::PDG_t(22);
        particles[i+SHIFT].pdg = t3::PDG_t(11);
      } else if (particles[i].pdg == t3::PDG_t(-11))
      {
        particles[i].pdg = t3::PDG_t(22);
        particles[i+SHIFT].pdg = t3::PDG_t(22);
      } else if (particles[i].pdg == t3::PDG_t(2112))
      {
        particles[i].pdg = t3::PDG_t(2112);
        particles[i+SHIFT].pdg = t3::PDG_t(2112);
      }
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

      particles[i].ir=-1;
      particles[i+SHIFT].id=MAX_ELEMENT+i;
      particles[i+SHIFT].ir=-1;

      if (particles[i].pdg == aParticleTable.makePDGfromZandA(1,2))
        particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(1,2);
      else if (particles[i].pdg == aParticleTable.makePDGfromZandA(22,48))
        particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(22,48);
      else if (particles[i].pdg == t3::PDG_t(22))
        particles[i+SHIFT].pdg = t3::PDG_t(11);
      else
        particles[i+SHIFT].pdg = t3::PDG_t(22);
    }
    MAX_ELEMENT+=POSITION23-POSITION3;
  }//end of (!HISTOGRAM)
}//End of Compressor

template <typename Floating> void DataHolder<Floating>::Inject() {  
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

//----------------------------------------------------------------------------//
  if(actionInInject.IsRegistered())
  {
#ifdef OPENACC
#pragma acc parallel loop gang vector present(particles)
#else
#pragma omp parallel for simd
#endif
    for(int i=0; i<LIFE; ++i)
    {
      actionInInject.UserAction(particles[i]);
    }
  }
//----------------------------------------------------------------------------//

  static int push = 0;
  if (LIFE > Ntop) Ntop = LIFE;
  if (push < N)
  {
    if (LIFE < K)
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
        /*
        particles[LIFE+i] = Particle<Floating>(
        t3::T3LorentzVector<Floating>(InitParticlex0*ag,
                                    InitParticley0*ag,
                                    -cuba*ag,0.0),
        t3::T3LorentzVector<Floating>(0.0,0.0,pls,InitEls),
        m, 0.0, initPDG, 1., MAX_ELEMENT+i, 1, MAX_ELEMENT+i-1, -1.0);
        */
        particles[LIFE+i].r=t3::T3LorentzVector<Floating>(InitParticlex0*ag,InitParticley0*ag,-cuba*ag,0.0);
        particles[LIFE+i].p=t3::T3LorentzVector<Floating>(0.0,0.0,pls,InitEls);
        particles[LIFE+i].m=m;
        particles[LIFE+i].de=0.0;
        particles[LIFE+i].pdg=initPDG;
        particles[LIFE+i].wt=1.0;
        particles[LIFE+i].rs=MAX_ELEMENT+i;
        particles[LIFE+i].ir=1;
        particles[LIFE+i].id=MAX_ELEMENT+i-1;
        particles[LIFE+i].id=-1.0;
      }
#endif
      push+=NUMBER_OF_PARTICLES_TO_PUSH;
      LIFE+=NUMBER_OF_PARTICLES_TO_PUSH;
      MAX_ELEMENT+=NUMBER_OF_PARTICLES_TO_PUSH;
      INJECTED_PARTICLES+=NUMBER_OF_PARTICLES_TO_PUSH;


#ifdef DEBSEED
      int COUNT_EQUAL_SEEDS=0;
      int INDEXI=-1;
      int INDEXJ=-1;
#ifdef OPENACC
#pragma acc parallel loop gang vector reduction(+:COUNT_EQUAL_SEEDS) present(particles)
//#pragma acc serial loop reduction(+:COUNT_EQUAL_SEEDS) present(particles)      
#else
#pragma omp parallel for simd reduction(+:COUNT_EQUAL_SEEDS)       
#endif      
      for(int i=LIFE-NUMBER_OF_PARTICLES_TO_PUSH; i<LIFE; ++i)
      {
        for(int j=i+1; j<LIFE; ++j)
        {
          if(particles[i].rs==particles[j].rs)
          {
            ++COUNT_EQUAL_SEEDS;
            break;
          }
        }
      }
  #ifdef OPENACC
  #pragma acc update host(particles[LIFE-NUMBER_OF_PARTICLES_TO_PUSH:LIFE])
  #endif      
      std::cout<<"There are "<<COUNT_EQUAL_SEEDS<<" equal seeds"
               <<" INDEXI="<<INDEXI<<" INDEXJ="<<INDEXJ
               <<" "<<particles[INDEXI].rs<<" "
               <<particles[INDEXJ].rs
               <<std::endl;
      /*
      std::cout<<"Print seeds: LIFE="<<LIFE<<" NUMBER_OF_PARTICLES_TO_PUSH="<<NUMBER_OF_PARTICLES_TO_PUSH<<std::endl;
      for(int i=LIFE-NUMBER_OF_PARTICLES_TO_PUSH; i<LIFE; ++i) std::cout<<particles[i].rs<<" ";
      std::cout<<std::endl;
      sleep(1);
      */
#endif          
      
    }
  }
 }//End of Injector

#endif//T3DATAHOLDER_H
