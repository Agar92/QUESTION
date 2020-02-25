#include "T3TabulatedFunction.hh"
#include <algorithm>

// #define DEBUG

template<class T1, class T2> size_t T3TabulatedFunction<T1, T2>::Get_size() const
{
#ifdef debug
  if( _arg.size() != _val.size() )
  {
    T3cout << "T2IsotopeFission::Get_size Error: numbers of arguments and values differ: "
           << "#arg = " << _arg.size() << ", #val = " << _val.size() << T3endl;
    rexit(-1);
  }
#endif
  return _arg.size();
}

template<class T1, class T2> T3bool T3TabulatedFunction<T1, T2>::Is_sorted() const
{
  T3bool result = true;
  for(size_t index = 0; index + 1 < _arg.size(); ++ index)
    result = result && Get_arg(index) <= Get_arg(index+1);
  return result;
}
template<class T1, class T2>
void T3TabulatedFunction<T1, T2>:: Erase(size_t ind_min, size_t ind_after_max)
{
  if(ind_min >= ind_after_max)
  {
    T3cout << "Warning: wrong index order: " << ind_min << "<=" << ind_after_max << T3endl;
  }
  else if(ind_min >= Get_size())
  {
    T3cout << "Warning: min_index out of range: " << ind_min << "=>" << Get_size()
           << T3endl;
  }
  else if(ind_after_max > Get_size())
  {
    T3cout << "Warning: min_index out of range: " << ind_after_max << ">" << Get_size()
           << "Erasing to the end" << T3endl;
           _arg.erase(_arg.begin() + ind_min, _arg.end());
           _val.erase(_val.begin() + ind_min, _val.end());
  }
  else
  {
    _arg.erase(_arg.begin() + ind_min, _arg.begin() + ind_after_max);
    _val.erase(_val.begin() + ind_min, _val.begin() + ind_after_max);
  }
}

template<class T1, class T2> void T3TabulatedFunction<T1, T2>::Add_point(T1 arg, T2 val,
                                                                         T1 delta)
{
  if (!Get_size()) Push_back(arg, val);
  else if( arg > Get_max_arg()) Push_back(arg, val);
  else
  {
#ifdef debug
    if (arg < Get_max_arg())
      T3cout << "T3TabulatedFunction::Add_point: Warning: non-ascending argument order, "
             << " max_arg = " << Get_max_arg() << ", new arg = " << arg << T3endl;
#endif
    size_t index = 0;
    while (Get_arg(index) < arg) index += 1; // no pre-increment operator ??
    if( Get_arg(index) == arg )
    {
#ifdef debug
      T3cout << "Double-value point arg = " << arg << T3endl;
#endif
      if(delta) Add_point(arg + delta, val, delta);
      else
      {
        _arg.insert( _arg.begin() + index, arg);
        _val.insert(     _val.begin() + index, val);
      }
    }// FIXED sometimes added twice
    else
    {
      _arg.insert( _arg.begin() + index, arg);
      _val.insert(     _val.begin() + index, val);
    }
  }
}

template<class T1, class T2> T2 T3TabulatedFunction<T1, T2>::Calc_val(T1 arg,
                                                                      T2 by_def) const
{
//   G4cout << Get_min_arg() << ' ' << arg << ' ' << Get_max_arg() << G4endl;
  if (!Get_size()) return by_def;  // return the default value if empty
  else if( Get_min_arg() > arg || arg > Get_max_arg() ) return by_def; // or out of range
  else
  {
    const size_t ind_ins = std::upper_bound(_arg.begin(), _arg.end(), arg) - _arg.begin();
    const size_t ind_upp = Get_size() == ind_ins ? ind_ins - 1 : ind_ins;
    const size_t ind_low = ind_upp ? ind_upp - 1 : ind_upp;
//     G4cout << ind_upp << ' ' << ind_low << ' ' << Get_size() << G4endl;
    const T1 arg_low = Get_arg(ind_low);
    const T1 arg_upp = Get_arg(ind_upp);
    const T1 arg_diff = arg_upp - arg_low;
    const T3double val_low = Get_val(ind_low);
    const T3double val_upp = Get_val(ind_upp);
#ifdef debug
    T3cout << "T2TabulatedFunction::Calc_xs(" << arg << "): E:xs(ind_low=" << ind_low
    << ") = " << arg_low << " : " << val_low << ", E:xs(ind_upp="
    << ind_upp << ") = " << arg_upp<< " : " << val_low << T3endl;
#endif
    if ( arg_diff > 1e-8) // FIXME some resonable T1 value, currently for float only
    {
//       ret_val = ( val_upp*(arg - arg_low)
//                 + val_low*(arg_upp - arg)*(-1))
//                 *(1/(arg_diff ));

      return val_low + (val_upp - val_low)*(arg - arg_low)/(arg_diff);
    }
    else
    {
      return (val_low + val_low) / 2;
    }
  }
}

template<class T1, class T2> const T3TabulatedFunction<T1, T2>
T3TabulatedFunction<T1, T2>::operator*(const T1& t1) const
{
  T3TabulatedFunction<T1, T2> result;
  result.Resize(Get_size());
  for(size_t ind = 0; ind < Get_size(); ++ind)
  {
    result.Set_arg(ind, Get_arg(ind));
    result.Set_val(ind, Get_val(ind)*t1);
  }
  return result;
}

template<class T1, class T2> const T3TabulatedFunction<T1, T2>
T3TabulatedFunction<T1, T2>::operator+(const T3TabulatedFunction<T1, T2>& t) const
{
  T3TabulatedFunction<T1, T2> result;
  size_t ind1 =0;
  size_t ind2 = 0;
  while( ind1 < Get_size() && ind2 < t.Get_size())
  {
    if ( Get_arg(ind1) < t.Get_arg(ind2) ){
      result.Push_back(Get_arg(ind1), Get_val(ind1) + t.Calc_val( Get_arg(ind1) ) );
      ++ind1;
    }
    else if ( Get_arg(ind1) > t.Get_arg(ind2) ){
      result.Push_back(t.Get_arg(ind2), t.Get_val(ind2) + Calc_val( t.Get_arg(ind2) ) );
      ++ind2;
    }
    else
    {
      result.Push_back(Get_arg(ind1), Get_val(ind1) + t.Get_val(ind2) );
      ++ind1;
      ++ind2;
    }
  }
  while( ind1 < Get_size() )
  {
    result.Push_back(Get_arg(ind1), Get_val(ind1));
    ++ind1;
  }
  while( ind2 < t.Get_size() )
  {
    result.Push_back(t.Get_arg(ind2), t.Get_val(ind2));
    ++ind2;
  }
  return result;
}

template<class T1, class T2> std::vector<T2> const
T3TabulatedFunction<T1, T2>::Get_equidistant(T1 arg_min, T1 arg_max, size_t nsteps) const
{
  std::vector<T2> result;
  result.reserve(nsteps);
  for(T1 arg = arg_min; arg <= arg_max; arg += (arg_max-arg_min)/nsteps)
  {
    result.push_back(Calc_val(arg));
  }
  return result;
}

template<class T1, class T2>
void T3TabulatedFunction<T1, T2>::Save_to(std::ostream& os) const
{
  unsigned int vector_size = Get_size();
#ifdef DEBUG
  T3cout << "T2TabulatedFunction<T1, T2>::Save_to: " << *this;
  T3cout << "T2TF::Save_to: size " << this->Get_size() << " " << vector_size
  << "sofT1="<<sizeof(T1) << " sofT2="<<sizeof(T2)<< T3endl;
#endif
  os.write( (const char*) &( vector_size ), sizeof( unsigned int ) );
  os.write( (const char*) &( _arg[0] ), vector_size*sizeof(T1));
  os.write( (const char*) &( _val[0] ), vector_size*sizeof(T2));
}

template<class T1, class T2>
void T3TabulatedFunction<T1, T2>::Load_from(std::istream& is)
{
  unsigned int vector_size;
  is.read( (char*) &vector_size, sizeof( unsigned int ) );
  Resize(vector_size);
  is.read((char*) &_arg[0], vector_size * sizeof( T1 ) );
  is.read((char*) &_val[0], vector_size * sizeof( T2 ) );
#ifdef DEBUG
  T3cout << "T2TabulatedFunction<T1, T2>::Load_from: " << *this;
  T3cout << "T2TF::L_f: size " << this->Get_size() << " " << vector_size
  << "sofT1="<<sizeof(T1) << " sofT2="<<sizeof(T2)<< T3endl;
#endif
}

template<class T1, class T2> T3TabulatedFunction<T1, T2>&
T3TabulatedFunction<T1, T2>::operator=(const T3TabulatedFunction<T1, T2>& tab)
{
  _arg = tab._arg;
  _val = tab._val;
  return *this;
}

template<class T1, class T2> T3bool
T3TabulatedFunction<T1, T2>::operator==(const T3TabulatedFunction<T1, T2>& tab) const
{
  return _arg == tab._arg && _val == tab._val;
}

template class T3TabulatedFunction<T3double, T3double>;
// template class T2TabulatedFunction<G4double, T2TabulatedFunction<G4double, G4double> >;
