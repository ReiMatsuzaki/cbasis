//#include <boost/static_assert.hpp>
//#include <boost/type_traits.hpp>
#include "exp_func.hpp"
#include <complex.h>
#include "math_utils.hpp"

namespace l2func {

  // ==== Constructors ====
  template<class F, int m>
  ExpFunc<F,m>::ExpFunc(): c_(0), 
			   n_(0), 
			   z_(0) {}

  template<class F, int m>
  ExpFunc<F,m>::ExpFunc(F _c, int _n, F _z): c_(_c), 
					     n_(_n), 
					     z_(_z) {}
  
  template<class F, int m>
  ExpFunc<F,m>::ExpFunc(const ExpFunc<F,m>& o): c_(o.c_), 
						n_(o.n_), 
						z_(o.z_) {}
    

  // ==== Methods =====
  template<class F, int m>
  F ExpFunc<F,m>::at(F x) const {
    return c_ * pow(x, n_) * exp(-z_ * pow(x, m));
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
