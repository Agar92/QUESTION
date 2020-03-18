#pragma once

//
// ****************************************************************
// * TPT License and Disclaimer                                   *
// *                                                              *
// * The TPT Software  is  copyright  of the Copyright Holders of *
// * the TPT CFAR-VNIIA group. It is provided under the terms and *
// * conditions of  the Software License  (ROSPATENT 2014611928). *
// *                                                              *
// * This code is  NOT  an open code. It is distributed in a form *
// * of compiled libraries. The code implementation is the result *
// * of the  scientific and technical work of the CFAR-VNIIA  TPT *
// * Scientific Group. By using or distributing  the TPT software *
// * (or any work based on the software) you agree to acknowledge *
// * its use  in  resulting scientific publications, and indicate *
// * your acceptance of all terms of the TPT Software License.    *
// ****************************************************************
//

#ifndef T3R_DDCS_HH
#define T3R_DDCS_HH

#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include <array>
#include <vector>
#include "T3Types.h"
#include "T3R_RW.hh"
#include "T3TabulatedFunction.hh"
#include "T3TabulatedCS.hh"
#include "T3NSGangular_RWnode.hh"
#include "T3NSGangular_RWrecord.hh"
#include "T3NSGangular_RW.hh"

namespace t3 {
using namespace units;

//--------------------------------------------------------------
//  In this file there is a d-d elastic scattering approximation
//  from /users/ALPHA_ELAST/FINAL_CORRECT_dd .
//--------------------------------------------------------------

//*************************************************************//
//DESCRIPTION:                                                 //
//This file is only for filling and Reading the database of    //
//energy/partial sums of D-D elastic scattering in the energy  //
//range 30 keV - 250 MeV                                       //
//in /home/70-gaa/T3data/ on i7.                               //
//*************************************************************//
  
class T3R_DDCS
{
public:
  //constructor:
  T3R_DDCS();
  //destructor:
  ~T3R_DDCS(){std::cout<<"~T3R_DDCS(): "<<std::endl;}
  //these functions fills the (E, partial sums) table in T3_DATA:
  void Fill();
  
  //returns the value of D-D differential cross-section (in inner units - MeV)
  //from  our D-D elastic scattering approximation.
  double GetdSigmadt(double Tls/*LS kinetic energy of the incident deuteron in inner units (MeV)*/,
                     double t/*|t| - Mandelstahm variable in inner units*/) const;

  //returns the Rutherford D-D differential cross section.
  double GetdSigmadt_DD_Rutherford(double Tls/*LS kinetic energy of the incident deuteron in inner units (MeV)*/,
                                   double t/*|t| - Mandelstahm variable in inner units*/) const;

  //returns thye integral Rutherford cross section in inner units.
  double GetRutherfordDDCSIntegral(double Tls/*LS kinetic energy of the incident deuteron*/,
                                   double tlow/*lower bound*/, double tup/*upper bound*/) const;
  
private:
  //number of energies in (E, partial sums) database:
  static constexpr int Bin2=100;
  //129 values (walls of bins)=>128 bins.
  //number of bins in partial sums of (E, partial sums):
  static constexpr int Bin3=128;
  //129 values (walls of bins) => without 0 and 1 there will be 127 values.
  //number of values in partial sums:
  static constexpr int _numpoint=127;
  //number of subbins, which we use for more accurate
  //integration of dsigma/d|t| over the Bin3 bins.
  static constexpr int SubBin3=4;
public:
  static int GetNumberOfBins2() {return Bin2;}
  static int GetNumberOfBins3() {return Bin3;}
  T3R_RW & GetT3R_RW()
  {
    std::cout<<"BEFORE seg fault"<<std::endl;
    return ecs;
  }
  T3NSGangular_RW GetT3NSGangular_RW(){return trw;}
  //write to a file (binary)
//WRITE TO FILE
//1. save the file (E,CS):
  void Save_to_CS();//CS -cross section
//1. read the file (E,CS):
  bool Load_from_CS();
//2. save the file (E, partial sums):
  void Save_to_PS(/*target Z and A*/int tgZ, int tgA, const T3String & rid,
                          T3int pdg, T3int incZA);//PS - partial sums
//2. read the file (E, partial sums):
  bool Load_from_PS(/*target Z and A*/int tgZ, int tgA, const T3String & rid,
                          T3int pdg, T3int incZA);//PS - partial sums
  //may be called after Fill() was called:
  std::array<double, Bin2> GetE() const {return Energy;}
  //may be called after Fill() was called:
  double GetCS(double Tls/*LS kinetic energy of the inc deuteron in inner units (MeV)*/) const;

private:
  const double Edisplace_deuteron = 10 * eV;//displacement energy for deuteron is 10 eV.
  const double md=1.8756*GeV;
  const double md2=md*md;
  double tmin=2*md*Edisplace_deuteron;
//------------------------------------------------------------------------------//
//1. We found that the difference of elastic D-D scattering differential cross section
//from Rutherford differential cross section at the scattering angle in CM
//theta_cm=90 grad becomes 1% at Tls=30 keV
//(Tls - the LS kineticenergy of the incident particle).
//That is why we chose Tls_min=30 keV.
//In /users/ALPHA_ELAST/CHECK_DD_APPROXIMATION_IN_TPROC/
// /SEE_ALL_DIFF_CS/Make01_1_minus_costhetacm_axis/
// /FIND_AT_WHICH_ENERGY_THE_DIFFERENCE_BETWEEN_RUTHERFORD_AND_APPROXIMATION_BEGINS
//
//At Tls<30 keV there is only Rutherford scattering.
//2. We decided that the lower bound of 1-cos(theta_cm) is 1.0e-3.
//In /users/ALPHA_ELAST/WORK_ON_FILLING_ELASTIC_DD_APPROXIMATION_PARTIAL_SUMS_DATABASE_IN_TPROC/
// /FIND_Point_of_difference_e_dependency
//------------------------------------------------------------------------------//
  const double kb=1.0e-3;//lower bound of 1-cos(theta_cm) axis for all energies.
  double Emin=30.0 * keV;//minimum LS kinetic energy (for elastic scattering, if less, than Rutherford scattering) of the inc deuteron.
  double Emax=250.0 * MeV;//maximum LS kinetic energy (for elastic scattering) of the inc deuteron.
  double deltalncos;//the step of logarithmic 1-cos(theta_cm) axis
  double deltalnE;//the step of logarithmic lnE axis for energies in partial sums
  std::array<double, Bin2> lnE;//make logarithmic E axis for partial sums db.
  std::array<double, Bin3+1> lncos;//make logarithmic 1-cos(theta_cm) for partial sums db.
  std::array<double, Bin2> Energy;
  std::array<std::array<double, Bin3+1>, Bin2> Fk;//for normalized partial sums
//-----------------------------------------------------------------------------//
  //for integral cross sections:
  std::array<double, Bin2> CS;
//-----------------------------------------------------------------------------//
  T3R_RW ecs;
  T3NSGangular_RWrecord tps;
  T3NSGangular_RW trw;
};

} // namespace t3

#endif//T3R_DDCS_HH
