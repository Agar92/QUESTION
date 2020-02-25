#include "T3NSGangular_RWnode.hh"

namespace
{
  T3LorentzVector<FloatingType> CalcXYZ(T3double xmin, T3double xmax, T3double p1, T3double p2, T3double p3, size_t i2);
}

T3NSGangular_RWnode::T3NSGangular_RWnode(T3double e,
                                         const std::array<T3double, _num_point> src,
                                         T3double newpr, T3double newsl) : E(e),
                                         sl(newsl), pr(newpr), V(src)
{}

T3NSGangular_RWnode::T3NSGangular_RWnode(T3NSGangular_RWnode const &node): E(node.E),
                                         sl(node.sl), pr(node.pr), V(node.V){}

T3NSGangular_RWnode& T3NSGangular_RWnode::operator=(T3NSGangular_RWnode const &node)
{
  E = node.E;
  sl = node.sl;
  pr = node.pr;
  V = node.V;
  return *this;
}

T3double T3NSGangular_RWnode::Get_V(size_t point_num) const
{
  if (point_num > _num_point + 1)
  {
    return 0;
  }
  else if (0 == point_num)              return 0;
  else if (_num_point + 1 == point_num) return 1;
  else                                  return V.at(point_num - 1);
}

void T3NSGangular_RWnode::Set_V(size_t point_num, T3double val)
{
  if ( _num_point < point_num || 0 == point_num)
  {
    return;
  }
  else V.at(point_num - 1) = val;
}

void T3NSGangular_RWnode::Set_V(const std::array<T3double, _num_point> src)
{
  V = src;
}

std::ofstream& T3NSGangular_RWnode::save_binary(std::ofstream& out_stream) const
{
  if( out_stream.good() )
  {
    unsigned int vector_size = _num_point;
    out_stream.write((const char*) &E, sizeof(T3double) );
    out_stream.write((const char*) V.data(), vector_size * sizeof(T3double));
    out_stream.write((const char*) &pr, sizeof(T3double) );
    out_stream.write((const char*) &sl, sizeof(T3double) );
  }
  else{}
  return out_stream;
}

std::ifstream& T3NSGangular_RWnode::load_binary(std::ifstream& in_stream )
{
  if( in_stream.good() )
  {
    unsigned int vector_size = _num_point;
    in_stream.read((char*) & E, sizeof(T3double) );
    in_stream.read((char*) V.data(), vector_size * sizeof( T3double) );
    in_stream.read((char*) & pr, sizeof(T3double) );
    in_stream.read((char*) & sl, sizeof(T3double) );
  }
  return in_stream;
}

T3LorentzVector<FloatingType> T3NSGangular_RWnode::CalcXYZ(size_t point_num, T3double xmin, T3double xmax) const
{
  if(0 < point_num && point_num <= _num_point)
    return ::CalcXYZ(xmin,xmax,Get_V(point_num - 1),Get_V(point_num),Get_V(point_num + 1),point_num);
  else
  {
    return T3LorentzVector<FloatingType>();
  }
}

std::vector<T3LorentzVector<FloatingType>> T3NSGangular_RWnode::CalcXYZ(T3double xmin, T3double xmax) const
{
  std::vector<T3LorentzVector<FloatingType>> result;
  result.reserve(_num_point);
  for(size_t ind = 1; ind <= _num_point; ++ind )
  {
    result.push_back(CalcXYZ(ind,xmin,xmax));
  }
  return result;
}


std::array<T3double, T3NSGangular_RWnode::_num_point> T3NSGangular_RWnode::isotropic()
{
  std::array<T3double, _num_point> array;
  static const T3double nbins = _num_point + 1;
  for (size_t ind = 0; ind < _num_point; ++ ind) array.at(ind) = T3double(ind+1)/nbins;
  return array;
}

namespace
{
  T3LorentzVector<FloatingType> CalcXYZ(T3double xmin, T3double xmax, T3double p1, T3double p2, T3double p3,  size_t i2)
  {
    static const size_t   _num_point = T3NSGangular_RWnode::_num_point;
    const T3double DeltaX=xmax-xmin;
    /*static */const T3double dCosT= DeltaX / ( _num_point + 1 );
    static const T3double tCosT= 2 * dCosT;
    static const T3double dop= 1.E-6;
    static const T3double eql= 1.E-7;
    const long double ct = i2 * dCosT + xmin;
    T3LorentzVector<FloatingType> res(0., 0., 0., p2);
    if(std::fabs(p2-p1) < eql && std::fabs(p3-p2) < eql )
    {
      p2=p1+eql;
      p3=p2+eql;
    }
    else if(std::fabs(p2-p1) < eql)
    {
      p2=p1+eql;
    }
    else if(std::fabs(p3-p2) < eql)
    {
      p3=p2+eql;
    }
    const long double r12 = p2-p1;
    const long double r23 = p3-p2;
    const long double r13 = p3-p1;
    
    const long double s12 = p1 + p2;
    const long double s23 = p2 + p3;
    const long double s13 = p1 + p3;
    const long double d12 = p1*p2*r12;
    const long double d23 = p2*p3*r23;
    const long double d13 = p1*p3*r13;
    const long double D   = d12 + d23 - d13;
    if(std::fabs(D) < 1.E-13 * d13 || D == 0.)
    {
      long double a = ct;
      long double b = tCosT;
      if(r13 > dop)
      {
        a -= dCosT*s13/r13;
        b /= r13;
      }
      else if(r23 > r12 && r23 > dop)
      {
        a -= dCosT*s23/r23;
        b /= r23;
      }
      else if (r12 > dop)
      {
        a -= dCosT*s12/r12;
        b /= r12;
      }
      else if (r13 > 0.)
      {
        b /= r13;
        a -= b*(s13+p2)/3;
      }
      res.SetX(T3double(a));
      res.SetY(T3double(b));
      return res;
    }
    const long double r = dCosT / D;
    const long double a = ct + r * (d12-d23);
    const long double b = r * (r23*s23 - r12*s12);
    const long double c = r * (r12 - r23);
    res.SetX(T3double(a));
    res.SetY(T3double(b));
    res.SetZ(T3double(c));
    return res;
  }
}



std::ostream& operator<<(std::ostream& os, const T3NSGangular_RWnode& inst)
{
  os << "Node e = " << inst.Get_E() << ", pr = " << inst.Get_pr() << ", sl = "
     << inst.Get_sl() << ", integrated distribution:";
  for(size_t ind =1; ind <= inst._num_point; ++ind) os << " " << inst.Get_V(ind);
  return os;
}

