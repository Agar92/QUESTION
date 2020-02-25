#pragma once
#ifndef T3DEFS_H
#define T3DEFS_H

#include <cstdint>

namespace t3 {

// TODO calculate with needed precision at compile time
constexpr auto PI = 3.14159265358979323846;
constexpr auto TWOPI = PI * 2;

namespace units {
//energy is measured in MeV's:
constexpr double MeV = 1.0;         // 1 MeV = 1.0e6 eV;
constexpr double keV = 1.0e-3 * MeV;// 1 keV = 1.0e3 eV;
constexpr double eV  = 1.0e-6 * MeV;// 1 eV
constexpr double GeV = 1.0e+3 * MeV; // 1 GeV = 1.0e9 eV.
constexpr double J   = 6.24150934e+18 * eV; // 1 Joule = 6.24115093e18 eV.
constexpr double erg = 1.0e-7 * J; // 1 GeV = 1.0e9 eV.
//mass is measured in grams:
constexpr double gr = 5.60958884e+32 * eV; //1 g   = 5.609585e26  MeV.
constexpr double ug = 1.0e-6 * gr; // 1 ug = 1.0e-6 * gr; microgram.
constexpr double mg = 1.0e-3 * gr; // 1 mg = 1.0e-3 * gr; milligram.
constexpr double kg = 1.0e+3 * gr;  // 1 kg= 1.0e3 g; kilogram.
//distance is measured in mm's:
constexpr double mm = 1.0;          //1 millimeter = 1.0e-3 m.
constexpr double cm = 10.0 * mm;    //1 cantimeter = 1.0e-2 m.
constexpr double m  = 1.0e+3 * mm;  //1 meter      = 1.0 m.
constexpr double um = 1.0e-3 * mm;  //1 micron     = 1.0e-6 m.
constexpr double nm = 1.0e-6 * mm;  //1 nanometer  = 1.0e-9 m.
constexpr double fm = 1.0e-12 * mm; //1 femtometer = 1.0e-15 m.
constexpr double in = 0.0254 * m;   //1 inch = 0.0254 m.
constexpr double A  = 0.1 * nm;     //1 Angstrem  = 0.1 nm.
//time is measured in nanoseconds:
constexpr double ns = 1.0;          //1 nanosecond  = 1.0e-9 s.
constexpr double us = 1.0e+3 * ns;  //1 microsecond = 1.0e-6 s.
constexpr double ms = 1.0e+6 * ns;  //1 millisecond = 1.0e-3 s.
constexpr double s  = 1.0e+9 * ns;  //1 second      = 1.0 s.
//velocity is measured in units v/c (dimensionless).
//surface:
constexpr double barn  = 1.0e-28 * m * m;//1 barn = 1.0e-28 * m^2.
constexpr double mbarn = 1.0e-3 * barn; //1 millibarn = 1.0e-3 * barn.
//magnetic field:
constexpr double T = 1.0;           //1 Tesla = 1.0.
constexpr double Gauss = 1.0e-4 * T;//1 Gauss = 1.0e-4 * T.
//force:
constexpr double Newton  = 1.0;         //1 Newton = 1.0.
constexpr double dyne = 1.0e-5 * Newton;//1 dyne = 1.0e-5 * N.
//temperature:
constexpr double Kelvin = 1.0;//grad Kelvin.
constexpr double ograd = 273.15 * Kelvin;//0 grad Celsius = 273.15K.
constexpr double kT = 1.0 / 38.68173135 * eV;//at T=300K kT=(38.68173135)^-1 * eV.
//pressure:
constexpr double Pa  = 1.0;//1 Pascal.
constexpr double atm = 101325.0 * Pa;//1 atmosphere = 101325 Pascal.
constexpr double Torr  = 1.0 * atm / 760;//1 Torr = 1 atmosphere / 760.
//---------------------------//
//coefficient to properly write density (use in MaterialTable.h):
//[density]=gr/cm^3
constexpr double  cm3 = cm * cm * cm;
constexpr double gcm3 = gr / cm3;
//r - mm=1.0. cm=10.0 * mm. m=1000.0*mm. nm, fm
//E - MeV eV=1.0e-6 * eV.
//t - nanoseconds. 1s=10^9ns.
}

} // namespace t3

#endif // T3DEFS_H
