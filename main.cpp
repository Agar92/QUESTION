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

constexpr int NUMBER_OF_GPUS=1; 

const unsigned int GL=2500000;
const int N=10000;
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

long int MAX_ELEMENT[NUMBER_OF_GPUS];
unsigned int POSITION3[NUMBER_OF_GPUS];
unsigned int POSITION2[NUMBER_OF_GPUS];
unsigned int POSITION1[NUMBER_OF_GPUS];
unsigned int POSITION0[NUMBER_OF_GPUS];
unsigned int POSITION23[NUMBER_OF_GPUS];
int LIFE[NUMBER_OF_GPUS]{0};

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
int dL[NUMBER_OF_GPUS];
int DL[NUMBER_OF_GPUS];
int n[NUMBER_OF_GPUS];
int numbin[NUMBER_OF_GPUS];

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
  for(int i=0; i<NUMBER_OF_GPUS; ++i)
  {
    acc_set_device_num(i,acc_device_nvidia);
#pragma acc parallel loop gang vector copy(LIFE[i]) present(particles)
    for(int i=0; i<LIFE[i]; ++i)
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
}

void make_particle(int N)
{
  const int PORTION=N / NUMBER_OF_GPUS;
  int PUSHED=0; 
  int PUSH=PORTION;
  for(int gpu=0; gpu<NUMBER_OF_GPUS; ++gpu)
  {
#pragma acc set device_type(acc_device_nvidia) device_num(gpu)    
    if(gpu==NUMBER_OF_GPUS-1) PUSH=N-PUSHED;
    PUSHED+=PUSH;
    std::cout<<"PORTION="<<PORTION<<" PUSHED="<<PUSHED<<std::endl;
    std::cout<<"STEP #5"<<std::endl;
    std::cout<<"Current device #"<<acc_get_device_num(acc_device_nvidia)<<std::endl;
    for(int j=0; j<PUSH; ++j)
    {
#pragma acc parallel num_gangs(1) vector_length(1) present(particles[0:GL]) copy(LIFE[gpu],MAX_ELEMENT[gpu])
      {
        particles[LIFE[gpu]].rs=MAX_ELEMENT[gpu]+27;
        particles[LIFE[gpu]].x=0.0f;
        particles[LIFE[gpu]].ir=-1;
        particles[LIFE[gpu]].count=0;
#pragma acc atomic update
        ++LIFE[gpu];
#pragma acc atomic update
        ++MAX_ELEMENT[gpu];
      }
    }
    std::cout<<"STEP #6"<<std::endl;
  }
}

void compressor()
{
  for(int gpu=0; gpu<NUMBER_OF_GPUS; ++gpu)
  {
    dL[gpu]=LIFE[gpu]/Nbin;
    DL[gpu]=dL[gpu]+1;
    if(LIFE[gpu]%Nbin==0) n[gpu]=Nbin;
    else                n[gpu]=Nbin*DL[gpu]-LIFE[gpu];
    POSITION0[gpu]=POSITION1[gpu]=POSITION2[gpu]=POSITION3[gpu]=POSITION23[gpu]=0;
    acc_set_device_num(gpu,acc_device_nvidia); 
    #pragma acc parallel loop copy(GL1,dL[gpu],DL[gpu],n[gpu],init,fin)
    {
      for(int b=0; b<Nbin; ++b)    
      {
        count01[b]=GL1;
        count23[b]=0;
        count0[b]=GL1;
        count1[b]=0;
        count2[b]=GL1;
        count3[b]=0;
        if(b<n[gpu])
        {
          init[b]=b*dL[gpu];
          fin[b]=(b+1)*dL[gpu];
        }
        else if(b==n[gpu])
        {
          init[b]=n[gpu]*dL[gpu];
          fin[b]=n[gpu]*dL[gpu]+DL[gpu];
        }
        else if(b>n[gpu])
        {
          init[b]=n[gpu]*dL[gpu]+DL[gpu]*(b-n[gpu]);
          fin[b]=n[gpu]*dL[gpu]+DL[gpu]*(b-n[gpu]+1);
        }
      }
    }
    acc_set_device_num(gpu,acc_device_nvidia); 
#pragma acc parallel loop copy(count01,count23,init,fin) present(particles,ind23)
  {
    for(int b=0; b<Nbin; ++b)    
    {
      for(int i=init[b]; i<fin[b]; ++i)
      {
        if(particles[gpu].ir<2) ind23[b][count01[b]--]=i;
        else                  ind23[b][count23[b]++]=i;
      }
    }
  }
  acc_set_device_num(gpu,acc_device_nvidia); 
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
  acc_set_device_num(gpu,acc_device_nvidia); 
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

  int P0=0, P1=0, P2=0, P3=0, P23=0;
  acc_set_device_num(gpu,acc_device_nvidia); 
#pragma acc parallel loop reduction(+:P0,P1,P2,P3,P23)  //reduction(+:POSITION0,POSITION1,POSITION2,POSITION3,POSITION23)
  for(int b=0; b<Nbin; ++b)
  {
    count0[b]=GL1-count0[b];
    count2[b]=GL1-count2[b];
    P0+=count0[b];
    P1+=count1[b];
    P2+=count2[b];
    P3+=count3[b];
    P23+=count23[b];
  }
  POSITION0[gpu]+=P0;
  POSITION1[gpu]+=P1;
  POSITION2[gpu]+=P2;
  POSITION3[gpu]+=P3;
  POSITION23[gpu]+=P23;
  LIFE[gpu]-=POSITION0[gpu];
  pointer1[0]=pointer2[0]=pointer3[0]=0;

  acc_set_device_num(gpu,acc_device_nvidia); 
#pragma acc parallel num_gangs(1) vector_length(1) copy(pointer1[0:Nbin],pointer2[0:Nbin],pointer3[0:Nbin])
  for(int b=0; b<Nbin-1; ++b)
  {
    pointer1[b+1]=pointer1[b]+count1[b];
    pointer2[b+1]=pointer2[b]+count2[b];
    pointer3[b+1]=pointer3[b]+count3[b];
  }

  acc_set_device_num(gpu,acc_device_nvidia); 
  for(int b=0; b<Nbin; ++b)
  {
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
    {
      cudaMemcpy(&arr1[pointer1[b]],&particles[init[b]+count3[b]+count2[b]],count1[b]*sizep,cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr2[pointer2[b]],&particles[init[b]+count3[b]],count2[b]*sizep,cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr3[pointer3[b]],&particles[init[b]],count3[b]*sizep,cudaMemcpyDeviceToDevice);
    }
  }
  
  acc_set_device_num(gpu,acc_device_nvidia); 
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
  {
    cudaMemcpy(&particles[0],&arr3[0],POSITION3[gpu]*sizep,cudaMemcpyDeviceToDevice);
    cudaMemcpy(&particles[POSITION3[gpu]],&arr2[0],POSITION2[gpu]*sizep,cudaMemcpyDeviceToDevice);
    cudaMemcpy(&particles[POSITION23[gpu]],&arr1[0],POSITION1[gpu]*sizep,cudaMemcpyDeviceToDevice);
  }


  }//END OF NUMBER_OF_GPUS loop.

  
}

int main(int argc, char **argv)
{
  std::cout<<"Num of GPUs="<<acc_get_num_devices(acc_device_nvidia)<<std::endl;
	auto begin=std::chrono::steady_clock::now();
  for(int gpu=0; gpu<NUMBER_OF_GPUS; ++gpu)
  {
    LIFE[gpu]=0;
    MAX_ELEMENT[gpu]=0;
  }
  std::cout<<"STEP #1"<<std::endl;
  int step=0;
  /*
  for(int i=0; i<NUMBER_OF_GPUS; ++i) acc_set_device_num(i,acc_device_nvidia);
  */
  std::cout<<"STEP #2"<<std::endl;
#pragma acc data create(particles,ind01,ind23,arr1,arr2,arr3,HIST) copy(HIST[0:SIZE]) 
  {
    std::cout<<"STEP #3"<<std::endl;
    make_particle(N);
    std::cout<<"STEP #4"<<std::endl;
    int life=0;
    for(int gpu=0; gpu<NUMBER_OF_GPUS; ++gpu) life+=LIFE[gpu];
    while(life>0)
    {
      propagator();
      compressor();
      ++step;
      life=0;
      for(int gpu=0; gpu<NUMBER_OF_GPUS; ++gpu) life+=LIFE[gpu];
    }
  }
  std::cout<<"Check the histogram:"<<std::endl;
  for(int i=0; i<SIZE; ++i) std::cout<<HIST[i]<<"  ";
  std::cout<<std::endl;
	auto end=std::chrono::steady_clock::now();
	auto elapsed_ms=std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
	std::cout<<"time="<<elapsed_ms.count()<<" ms"<<std::endl;
}
