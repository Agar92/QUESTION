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
//
// Created by Dmitry Savin 2014
// The last revision:
// ---------------------------------------------------------------------------
// Short description:
// ---------------------------------------------------------------------------

#ifndef T3TABULATEDCS_HH
#define T3TABULATEDCS_HH

#include "T3TabulatedFunction.hh"

class T3TabulatedCS
{
public:
  T3TabulatedCS(std::vector<T3double> arg = std::vector<T3double>(),
                std::vector<T3double> val = std::vector<T3double>()): _data(arg,val){}
  T3double Get_E(size_t n) const {return _data.Get_arg(n);}
  T3double Get_xs(size_t n) const {return _data.Get_val(n);}
  std::vector<T3double> Get_E() const {return _data.Get_arg();}
  std::vector<T3double> Get_xs() const {return _data.Get_val();}
  const T3TabulatedFunction<T3double, T3double>& Get_data() const {return _data;}
  void Set_E(size_t n, T3double E) {_data.Set_arg(n, E);}
  void Set_xs(size_t n, T3double xs) {_data.Set_val(n, xs);}
  size_t Get_size() const {return _data.Get_size();}
  T3double Calc_xs(T3double arg, T3double by_def = 0.) const
    {return _data.Calc_val(arg, by_def);}
  void Clear() {_data.Clear();}
  void Erase(size_t ind_min, size_t ind_after_max) {_data.Erase(ind_min, ind_after_max);}
  void Resize(size_t size){_data.Resize(size);}
  T3bool Is_sorted() const {return _data.Is_sorted();} // Monotone Energy check
  T3double Get_min_E() const {return Get_size() ? _data.Get_min_arg() : 0;}
  T3double Get_max_E() const {return Get_size() ? _data.Get_max_arg() : 0;}
  void Add_point(T3double arg, T3double val, T3double delta = 1e-11)
    {return _data.Add_point(arg, val, delta);}
  void Push_back(T3double arg, T3double val) {_data.Push_back(arg, val);}
  std::vector<T3double> Get_equidistant(T3double Emin, T3double Emax, size_t nsteps)
    const {return _data.Get_equidistant(Emin, Emax, nsteps);}
  void Save_to(std::ostream& os) const {_data.Save_to(os);}
  void Load_from(std::istream& is) {_data.Load_from(is);}
  inline T3TabulatedCS& operator=(const T3TabulatedCS& tab);
  inline const T3TabulatedCS operator*(const T3double mul) const;
  inline const T3TabulatedCS operator/(const T3double mul) const;
  inline const T3TabulatedCS operator+(const T3TabulatedCS& tab) const;
  inline const T3TabulatedCS operator-(const T3TabulatedCS& tab) const;
  inline T3bool operator==(const T3TabulatedCS& rval) const;
private:
  T3TabulatedFunction<T3double, T3double> _data;
};

inline std::ostream& operator<<(std::ostream& os, const T3TabulatedCS& inst)
{
  os << "(E, xs): " << inst.Get_data();
  return os;
}

inline T3TabulatedCS& T3TabulatedCS::operator=(const T3TabulatedCS& tab)
{
  _data = tab._data;
  return *this;
}

inline const T3TabulatedCS T3TabulatedCS::operator+(const T3TabulatedCS& tab) const
{
  T3TabulatedCS result;
  result._data = _data + tab._data;
  return result;
}

inline const T3TabulatedCS T3TabulatedCS::operator-(const T3TabulatedCS& tab) const
{
  return (*this) + (tab * (-1));
}

inline const T3TabulatedCS T3TabulatedCS::operator*(const T3double mul) const
{
  T3TabulatedCS result;
  result._data = _data * mul;
  return result;
}

inline const T3TabulatedCS T3TabulatedCS::operator/(const T3double den) const
{
  return (*this) * (1./den);
}

inline T3bool T3TabulatedCS::operator==(const T3TabulatedCS& rval) const
{
  T3bool temp = true;
  for (size_t ind = 0; ind < Get_size(); ++ind)
  {
    temp = temp && (fabs (Get_E(ind) - rval.Get_E(ind)) <= fabs(Get_E(ind)) * 1e-7 )
    && (fabs (Get_xs(ind) - rval.Get_xs(ind)) <= fabs(Get_xs(ind)) * 1e-7 );
  }
  return temp;
}

#endif
