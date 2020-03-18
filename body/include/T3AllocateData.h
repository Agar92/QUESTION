#pragma once
#ifndef T3ALLOCATEDATA_H
#define T3ALLOCATEDATA_H

namespace data {

//=======================================================================//

Particle<double> particles[GL] __attribute__((aligned(64)));

#ifdef OPENACC
#pragma acc declare create(particles[0:GL])
#endif

//=======================================================================//

double csBorderDataFS[GL]{0.0};//this is for GetFS() in Reactor().

/*
#ifdef OPENACC
#pragma acc declare create(csBorderDataFS[0:GL])
#endif
*/

//=======================================================================//

int ind01[Nbin][BLt] __attribute__((aligned(64)));
int ind23[Nbin][BLt] __attribute__((aligned(64)));
Particle<FloatingType> arr1[GL] __attribute__((aligned(64)));
Particle<FloatingType> arr2[GL] __attribute__((aligned(64)));
Particle<FloatingType> arr3[GL] __attribute__((aligned(64)));

//=======================================================================//
//For storing cross sections:
//=======================================================================//
FloatingType csMultipleScattering[GL];
FloatingType csInelasticdd[GL];
FloatingType csElasticEMIonIon[GL];
FloatingType csElasticStrongIonIon[GL];
//=======================================================================//

long int MAX_ELEMENT;
int SHIFT;
int INJECTED_PARTICLES=0;

int POSITION3;
int POSITION2;
int POSITION1;
int POSITION0;
int POSITION23;
int LIFE=0;
unsigned int sizep=sizeof(Particle<FloatingType>);
int push=0;
int over=0;
unsigned int Ntop=0;
double SumDGam=0.;
double GAMMA=0.;
double NoNew=0.0;

//=======================================================================//

decltype(INJECTED_PARTICLES) GetNumOfInjectedParticles(){return INJECTED_PARTICLES;}
decltype(LIFE) GetNumOfAliveParticles(){return LIFE;}
decltype(NoNew) GetNoNew(){return NoNew;}
decltype(SumDGam) GetSumDGam(){return SumDGam;}
decltype(Ntop) GetNtop(){return Ntop;}

//=======================================================================//

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

//=======================================================================//

}//namespace data.  

#endif//T3ALLOCATEDATA_H
