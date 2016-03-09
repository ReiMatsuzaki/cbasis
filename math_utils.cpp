#include "math_utils.hpp"
#include<complex>

typedef std::complex<double> CD;

namespace l2func {

  int Factorial(int n) {

    if(n < 0) {
      throw std::invalid_argument("n must be zero or positive in Factorial.");
    }

    int acc = 1;
    for(int i = n; i > 0; i--)
      acc *= i;

    return acc;
  }
  double DFactorial(int n) {
    return Factorial(n) * 1.0;
  }
  int DoubleFactorial(int n) {

    if(n < 0) {
      throw std::invalid_argument("n must be zero or positive in DoubleFactorial.");
    }


    int acc = 1;
    for(int i = n; i > 0; i-=2) {
      acc *= i;
    }
    
    return acc;
  }
  double DDoubleFactorial(int n) {
    return DoubleFactorial(n) * 1.0;
  }

  template<class F> F ConjugateIfPossible(F x) { return x; }
  template<> double ConjugateIfPossible(double x) { return x; }
  template<> CD ConjugateIfPossible<CD>(CD x) { return std::conj(x); }
}
