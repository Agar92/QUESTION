#include <chrono>
#include <iostream>
#include <unistd.h>
#include <cmath>

#include <accelmath.h>
#include <openacc.h>
#include <cuda.h>
#include <cuda_runtime.h>

//cmake . -DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++
//-DCMAKE_CXX_FLAGS="-acc -Minfo=all -mcmodel=medium
//-ta=tesla:cc30 -Mcuda=cuda10.1"

const unsigned int GL=2500000;
const int N=100000;
const unsigned int Nbin=1;
const int BLt=GL/Nbin;

struct Particle
{
  float        x;
  unsigned int rs;
  int          ir;
  int       count;
};

Particle particles[GL] __attribute__((aligned(64)));
int ind01[Nbin][BLt] __attribute__((aligned(64)));
int ind23[Nbin][BLt] __attribute__((aligned(64)));
Particle arr1[GL] __attribute__((aligned(64)));
Particle arr2[GL] __attribute__((aligned(64)));
Particle arr3[GL] __attribute__((aligned(64)));

long int MAX_ELEMENT;
unsigned int POSITION3;
unsigned int POSITION2;
unsigned int POSITION1;
unsigned int POSITION0;
unsigned int POSITION23;
int LIFE=0;

unsigned int sizep=sizeof(Particle);

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
int GL1=BLt-1;
int dL;
int DL;
int n;
int numbin;

const double UP_BORDER = 1.0;
const int SIZE=UP_BORDER/0.1*5.0;
int HIST[SIZE]{0};

#ifdef OPENACC
#pragma acc routine seq
#endif
unsigned int Rand32(unsigned int & xn)
{
  u_quad_t a=0x5DEECE66D;
  u_quad_t c=0xB;
  return (unsigned int)((a*xn+c) & 0xFFFFFFFF);
}

#ifdef OPENACC
#pragma acc routine seq
#endif
double RND01(unsigned int & xn)
{
  xn=Rand32(xn);
  return (double)(xn) / (double) 0x100000000LL;
}

void propagator()
{
#pragma acc parallel loop gang vector copy(LIFE) present(particles)
  for(int i=0; i<LIFE; ++i)
  {
    const double R1=RND01(particles[i].rs);
    if(R1>2.0/3.0)      particles[i].ir=3;
    else if(R1>1.0/3.0) particles[i].ir=2;
    else                particles[i].ir=1;
    const double R2=RND01(particles[i].rs);
    particles[i].x += R2/10.0;
    ++particles[i].count;
    if(particles[i].x > UP_BORDER)
    {
      particles[i].ir=0;
      ++HIST[particles[i].count];
    }
  }
}

void make_particle()
{
#pragma acc parallel num_gangs(1) vector_length(1) present(particles[0:GL]) copy(LIFE,MAX_ELEMENT)
  {
    particles[LIFE].rs=MAX_ELEMENT+27;
    particles[LIFE].x=0.0f;
    particles[LIFE].ir=-1;
    particles[LIFE].count=0;
#pragma acc atomic update
    ++LIFE;
#pragma acc atomic update
    ++MAX_ELEMENT;
  }
}

void compressor()
{
  dL=LIFE/Nbin;
  DL=dL+1;
  if(LIFE%Nbin==0) n=Nbin;
  else     n=Nbin*DL-LIFE;
  POSITION0=POSITION1=POSITION2=POSITION3=POSITION23=0;
  #pragma acc parallel loop copy(GL1,dL,DL,n,init,fin)
  {
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
  }

#pragma acc parallel loop copy(count01,count23,init,fin) present(particles,ind23)
  {
    for(int b=0; b<Nbin; ++b)    
    {
      for(int i=init[b]; i<fin[b]; ++i)
      {
        if(particles[i].ir<2) ind23[b][count01[b]--]=i;
        else                  ind23[b][count23[b]++]=i;
      }
    }
  }
  
#pragma acc parallel loop copy(count0,count1,count2,count3,count01,count23,mini,ii23,init,fin) present(particles,ind01,ind23)
  for(int b=0; b<Nbin; ++b)    
  {
    ii23[b]=count23[b]-1;
    mini[b]=GL1-count01[b];
    if(count23[b]<mini[b]) mini[b]=count23[b];
    int js=0;
#pragma acc loop vector reduction(+:js)
    for(int j=0; j<mini[b]; ++j) js+=int(ind23[b][ii23[b]-j]>ind23[b][GL1-j]);
#pragma acc loop vector
    for(int j=0; j<js; ++j) std::swap(particles[ind23[b][ii23[b]-j]],particles[ind23[b][GL1-j]]);    
    for(int i=init[b]; i<fin[b]; ++i)
    {
      if     (particles[i].ir==0) ind01[b][count0[b]--]=i;
      else if(particles[i].ir==1) ind01[b][count1[b]++]=i;
      else if(particles[i].ir==2) ind23[b][count2[b]--]=i;
      else                        ind23[b][count3[b]++]=i;
    }
  }
  
#pragma acc parallel loop copy(count0,count1,count2,count3,mini,ii1,ii3,GL1) present(particles)
  for(int b=0; b<Nbin; ++b)    
  {
    ii1[b]=count1[b]-1;
    mini[b]=GL1-count0[b];
    if(count1[b]<mini[b]) mini[b]=count1[b];
    int js=0;
#pragma acc loop vector reduction(+:js)
    for(int j=0; j<mini[b]; ++j) js+=int(ind01[b][ii1[b]-j]>ind01[b][GL1-j]);
#pragma acc loop vector      
    for(int j=0; j<js; ++j) std::swap(particles[ind01[b][ii1[b]-j]],particles[ind01[b][GL1-j]]);  
    ii3[b]=count3[b]-1;
    mini[b]=GL1-count2[b];
    if(count3[b]<mini[b]) mini[b]=count3[b];
    js=0;
#pragma acc loop vector reduction(+:js)
    for(int j=0; j<mini[b]; ++j) js+=int(ind23[b][ii3[b]-j]>ind23[b][GL1-j]);
#pragma acc loop vector
    for(int j=0; j<js; ++j) std::swap(particles[ind23[b][ii3[b]-j]],particles[ind23[b][GL1-j]]);
  }
  
#pragma acc parallel loop reduction(+:POSITION0,POSITION1,POSITION2,POSITION3,POSITION23)
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
  LIFE-=POSITION0;
  pointer1[0]=pointer2[0]=pointer3[0]=0;
#pragma acc parallel num_gangs(1) vector_length(1) copy(pointer1[0:Nbin],pointer2[0:Nbin],pointer3[0:Nbin])
  for(int b=0; b<Nbin-1; ++b)
  {
    pointer1[b+1]=pointer1[b]+count1[b];
    pointer2[b+1]=pointer2[b]+count2[b];
    pointer3[b+1]=pointer3[b]+count3[b];
  }
  
  for(int b=0; b<Nbin; ++b)
  {
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
    {
      cudaMemcpy(&arr1[pointer1[b]],&particles[init[b]+count3[b]+count2[b]],count1[b]*sizep,cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr2[pointer2[b]],&particles[init[b]+count3[b]],count2[b]*sizep,cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr3[pointer3[b]],&particles[init[b]],count3[b]*sizep,cudaMemcpyDeviceToDevice);
    }
  }
  
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
  {
    cudaMemcpy(&particles[0],&arr3[0],POSITION3*sizep,cudaMemcpyDeviceToDevice);
    cudaMemcpy(&particles[POSITION3],&arr2[0],POSITION2*sizep,cudaMemcpyDeviceToDevice);
    cudaMemcpy(&particles[POSITION23],&arr1[0],POSITION1*sizep,cudaMemcpyDeviceToDevice);
  }
}

int main(int argc, char **argv)
{
	auto begin=std::chrono::steady_clock::now();
  LIFE=0;
  MAX_ELEMENT=0;
  int step=0;
#pragma acc data create(particles,ind01,ind23,arr1,arr2,arr3,HIST) copy(HIST[0:SIZE])
  {
    for(int i=0; i<N; ++i) make_particle();
    while(LIFE>0)
    {
      propagator();
      compressor();
      ++step;
    }
  }
  std::cout<<"Check the histogram:"<<std::endl;
  for(int i=0; i<SIZE; ++i) std::cout<<HIST[i]<<"  ";
  std::cout<<std::endl;
	auto end=std::chrono::steady_clock::now();
	auto elapsed_ms=std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
	std::cout<<"time="<<elapsed_ms.count()<<" ms"<<std::endl;
}
