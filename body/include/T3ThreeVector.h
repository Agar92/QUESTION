#pragma once

#ifndef T3THREEVECTOR_H
#define T3THREEVECTOR_H

#include <cmath>

namespace t3
{

//************************************************************//
//These forward declarations are necessary for defining friend//
//functions in template class ThreeVector<T>:                 //
//************************************************************//
template <typename T>
class T3ThreeVector;

//scalar product:
template <typename T>
T operator*(const T3ThreeVector<T> & lhsvector, const T3ThreeVector<T> & rhsvector);

template <typename T>
T3ThreeVector<T> operator*(const T a, T3ThreeVector<T> const & rhsside);

template <typename T>
T3ThreeVector<T> operator*(T3ThreeVector<T> const & lhsside, const T a);

template <typename T>
T3ThreeVector<T> cross(T3ThreeVector<T> const & lhsvector, T3ThreeVector<T> const & rhsvector);  
//***************************************************************//

template <typename T>
class T3ThreeVector
{
private:
  T fx, fy, fz;
public:
//contsructors:
//\\//T3ThreeVector()=delete;
  T3ThreeVector():fx(0.0), fy(0.0), fz(0.0){}
  T3ThreeVector(T x, T y, T z):fx(x),fy(y),fz(z){}
  T3ThreeVector(T3ThreeVector const & vector);
//return components and vector lenghth:
//getters:
  T x() const {return fx;}
  T y() const {return fy;}
  T z() const {return fz;}  
  T R() const {return sqrt(fx*fx+fy*fy+fz*fz);}
//normalizes 3-vector:
  void Unit()
  {
    const T vector_length=sqrt(fx*fx+fy*fy+fz*fz);
    this->fx/=vector_length, this->fy/=vector_length, this->fz/=vector_length;
  }
//setters:
  void SetX(T X){fx=X;}
  void SetY(T Y){fy=Y;}
  void SetZ(T Z){fz=Z;}
//operator overloading:  
  T3ThreeVector & operator=(T3ThreeVector const & rhsvector);
  T3ThreeVector operator+(T3ThreeVector const & rhsvector);
  T3ThreeVector operator-(T3ThreeVector const & rhsvector);
  //unary operator, returns a negative copy of *this:
  T3ThreeVector operator-();
  T3ThreeVector & operator+=(const T3ThreeVector & rhsvector);
  T3ThreeVector & operator-=(const T3ThreeVector & rhsvector);
  //scalar product:
  friend T operator*<>(const T3ThreeVector & lhsvector, const T3ThreeVector & rhsvector);
  friend T3ThreeVector operator*<>(const T a, T3ThreeVector const & rhsside);
  friend T3ThreeVector operator*<>(T3ThreeVector const & lhsside, const T a);
  friend T3ThreeVector cross<>(T3ThreeVector const & lhsvector, T3ThreeVector const & rhsvector);  
  T3ThreeVector & operator*=(const T a);
  T3ThreeVector operator/(const T a);
  T3ThreeVector & operator/=(const T a);
};

template <typename T>
T3ThreeVector<T>::T3ThreeVector(T3ThreeVector const & vector):fx(vector.fx), fy(vector.fy), fz(vector.fz){}

template <typename T>
T3ThreeVector<T> & T3ThreeVector<T>::operator=(const T3ThreeVector & rhsvector)
{
  if(this==&rhsvector)
  {
    return *this;
  }
  this->fx=rhsvector.fx, this->fy=rhsvector.fy, this->fz=rhsvector.fz;
  return *this;
}

template <typename T>
T3ThreeVector<T> T3ThreeVector<T>::operator+(T3ThreeVector const & rhsvector)
{
  T3ThreeVector<T> result(this->fx+rhsvector.fx,this->fy+rhsvector.fy,this->fz+rhsvector.fz);
  return result;
}

template <typename T>
T3ThreeVector<T> T3ThreeVector<T>::operator-(T3ThreeVector const & rhsvector)
{
  T3ThreeVector<T> result(this->fx-rhsvector.fx,this->fy-rhsvector.fy,this->fz-rhsvector.fz);
  return result;
}

//unary operator, returns a negative copy of *this:
template <typename T>
T3ThreeVector<T> T3ThreeVector<T>::operator-()
{
  T3ThreeVector<T> result(-this->fx,-this->fy,-this->fz);
  return result;
}

template <typename T>
T3ThreeVector<T> & T3ThreeVector<T>::operator+=(T3ThreeVector const & rhsvector)
{
  this->fx+=rhsvector.fx, this->fy+=rhsvector.fy, this->fz+=rhsvector.fz;
  return *this;
}

template <typename T>
T3ThreeVector<T> & T3ThreeVector<T>::operator-=(T3ThreeVector const & rhsvector)
{
  this->fx-=rhsvector.fx, this->fy-=rhsvector.fy, this->fz-=rhsvector.fz;
  return *this;
}

//scalar product:
template <typename T>
T operator*(T3ThreeVector<T> const & lhsvector, const T3ThreeVector<T> & rhsvector)
{
  T result=lhsvector.fx*rhsvector.fx + lhsvector.fy*rhsvector.fy + lhsvector.fz*rhsvector.fz;
  return result;
}

template <typename T>
T3ThreeVector<T> operator*(const T a, T3ThreeVector<T> const & rhsside)
{
  T3ThreeVector<T> result(rhsside.fx*a, rhsside.fy*a, rhsside.fz*a);
  return result;
}

template <typename T>
T3ThreeVector<T> operator*(T3ThreeVector<T> const & lhsside, const T a)
{
  T3ThreeVector<T> result(lhsside.fx*a, lhsside.fy*a, lhsside.fz*a);
  return result;
}

template <typename T>
T3ThreeVector<T> cross(T3ThreeVector<T> const & lhsvector, T3ThreeVector<T> const & rhsvector)
{
  T3ThreeVector<T> result;
  result.fx=lhsvector.fy*rhsvector.fz-lhsvector.fz*rhsvector.fy;
  result.fy=lhsvector.fz*rhsvector.fx-lhsvector.fx*rhsvector.fz;
  result.fz=lhsvector.fx*rhsvector.fy-lhsvector.fy*rhsvector.fx;
  return result;
}

template <typename T>
T3ThreeVector<T> & T3ThreeVector<T>::operator*=(const T a)
{
  this->fx*=a, this->fy*=a, this->fz*=a;
  return *this;
}

template <typename T>
T3ThreeVector<T> T3ThreeVector<T>::operator/(const T a)
{
  T3ThreeVector<T> result(this->fx/a, this->fy/a, this->fz/a);
  return result;
}

template <typename T>
T3ThreeVector<T> & T3ThreeVector<T>::operator/=(const T a)
{
  this->fx/=a, this->fy/=a, this->fz/=a;
  return *this;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const T3ThreeVector<T> & vector)
{
  os<<vector.x()<<" "<<vector.y()<<" "<<vector.z()<<" ";
  return os;
}
  
}
#endif  //T3THREEVECTOR
