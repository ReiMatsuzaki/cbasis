#include <string>
#include "fact.hpp"

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
}
