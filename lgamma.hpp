#ifndef LGAMMA_TEMPLATE_H
#define LGAMMA_TEMPLATE_H

namespace l2func {
  template<class F> F LowerGamma(int n, F z) {
    if(n <= 0)
      throw "n must be positive";

    if(n == 1) {
      return 1.0-exp(-z);
    } else {
      F surface = -F(pow(z, n-1))*exp(-z);
      F prev = LowerGamma(n-1, z);
    
      return surface + F(n-1)*prev;
    }
  }
}
#endif
