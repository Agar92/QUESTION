#pragma once
#ifndef T3RNG_H
#define T3RNG_H

#include <cmath>
#include <cstdint>
#include <iostream>
//------------------------//
//Random number generator.//
//------------------------//

namespace t3 {

class RNDGenerator
{
private:
  unsigned int fseed;
public:
  void seed(unsigned int seed){fseed=seed;}
  using result_type=unsigned int;
  result_type GetSeed(){return fseed;}
  RNDGenerator():fseed(1){}
  RNDGenerator(unsigned int seed):fseed(seed){}
  result_type min() const {return 0u;}
  result_type max() const {return 0xFFFFFFFF;}
  result_type Rand32(unsigned int xn)
  {
    u_quad_t a=0x5DEECE66D;
    u_quad_t c=0xB;
    return (unsigned int)((a*xn+c) & 0xFFFFFFFF);
  }
  result_type operator()()
  {
    fseed=Rand32(fseed);
    return fseed;
  }
};
  
//From T3Util.h:
template <typename Floating, typename RandomEngine>
auto GenerateSubCanonical(RandomEngine &engine) {
  auto const randomNumber = engine();
  auto maxRandomNumber = engine.max();
  auto maxRandomNumberAsFloating =
      static_cast<Floating>(maxRandomNumber);
  auto const result =
      static_cast<Floating>(randomNumber) / 0xFFFFFFFF;
  return result;
}

template <typename Floating>
bool isnan(Floating a){return (!(a>1.0) && !(a<2.0));}

} // namespace t3

#endif // T3RNG_H
