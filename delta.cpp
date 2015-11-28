#include "delta.hpp"

namespace l2func {

  // ==== Constructors ====
  template<class F>
  DiracDelta<F>::DiracDelta() :
    r0_(0.0), h_(0.001) {}
  template<class F>
  DiracDelta<F>::DiracDelta(double _r0) :
    r0_(_r0), h_(0.001) {}
  template<class F>
  DiracDelta<F>::DiracDelta(const DiracDelta& o) :
    r0_(o.r0()), h_(0.001) {}

  // ==== Method ====
  template<class F> F DiracDelta<F>::at(F x) const {
    
    if( std::abs(x-this->r0()) < this->h_) 
      return F(1);
    else
      return F(0);

  }

  // ==== Externals ====
  template<class F>
  std::ostream& operator << (std::ostream& os, const DiracDelta<F>& a) {
    os << "Delta(r0=" << a.r0() << ")";
    return os;
  }

  // ==== Explicit Decralation ====
  template class DiracDelta<double>;
  template class DiracDelta<std::complex<double> >;
}
