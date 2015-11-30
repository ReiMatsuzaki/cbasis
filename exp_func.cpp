//#include <boost/static_assert.hpp>
//#include <boost/type_traits.hpp>
#include "exp_func.hpp"
#include <complex.h>
#include "math_utils.hpp"

namespace l2func {

  // ==== Methods =====
  template<class F, int m>
  F ExpFunc<F,m>::at(F x) const {
    return this->c() * pow(x, this->n()) * exp(-this->z() * pow(x, m));
  }

  template<class F, int m>
  void ExpFunc<F,m>::SetComplexConjugate() {
    ExpFunc<F,m> a(ConjugateIfPossible<F>(c()),
		   n(),
		   ConjugateIfPossible<F>(z()));
    *this = a;

  }

  template<class F, int m>
  void ExpFunc<F,m>::SetDerivParam() {
    
    ExpFunc<F,m> a(-this->c(), this->n() + m, this->z());
    *this = a;
    
  }

  template<class F, int m>
  void ExpFunc<F,m>::SetScalarProd(F c) {
    this->set_c(this->c() * c);
  }

  template<class F, int m>
  void ExpFunc<F,m>::SetRmProd(int n) {
    this->set_n(this->n() + n);
  }

  template<class F, int m>
  typename ExpFunc<F,m>::FuncDerivOne ExpFunc<F,m>::DerivParamOne() {
    FuncDerivOne s(*this);
    s.SetDerivParam();
    return s;
  }
  template<class F, int m>
  typename ExpFunc<F,m>::FuncDerivTwo ExpFunc<F,m>::DerivParamTwo() {
    FuncDerivOne s(*this);
    s.SetDerivParam();
    s.SetDerivParam();
    return s;
  }

  // ==== Externals ====
  template<class F, int m>
  std::ostream& operator << (std::ostream& os, const ExpFunc<F,m>& a) {
  F c = a.c();

    if(m == 1) 
      os << (char *)"STO(" ;
    else
      os << (char *)"GTO(";
    char* comma = (char*)", ";
    os << c;
    os << comma << a.n() << comma << a.z() << (char*)")";
    return os;
  } 

  // ==== Explicit Decralation ====
  template class ExpFunc<double, 1>;
  template class ExpFunc<std::complex<double>, 1>;
  template class ExpFunc<double, 2>;
  template class ExpFunc<std::complex<double>, 2>;

  template std::ostream& operator <<<std::complex<double>,1> 
  (std::ostream& os, const CSTO& a);
  template std::ostream& operator <<<double,1> 
  (std::ostream& os, const RSTO& a);
  template std::ostream& operator <<<std::complex<double>,2> 
  (std::ostream& os, const CGTO& a);
  template std::ostream& operator <<<double,2>
  (std::ostream& os, const RGTO& a);


}
