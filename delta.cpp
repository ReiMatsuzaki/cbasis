#include "delta.hpp"

namespace l2func {

  // ==== Big 3 ====
  template<class F>
  DiracDelta<F>::DiracDelta() :
    r0_(0.0), h_(0.001) {}
  template<class F>
  DiracDelta<F>::DiracDelta(double _r0) :
    r0_(_r0), h_(0.001) {}
  template<class F>
  DiracDelta<F>::DiracDelta(const DiracDelta& o) :
    r0_(o.r0()), h_(0.001) {}
  
  template<class F> F DiracDelta<F>::at(F x) const {
    
    if( abs(x-this->r0()) < this->h_) 
      return F(1);
    else
      return F(0);

  }
  
  template<class F> F CIP(const DiracDelta<F>& a, const CSTO& b) {

    return b.at(a.r0());

  }
  template<class F> F CIP(const DiracDelta<F>& a, const CGTO& b) {

    return b.at(a.r0());

  }
  template<class F> F CIP(const CSTO& b, const DiracDelta<F>& a) {

    return CIP(a, b);

  }
  template<class F> F CIP(const CGTO& b, const DiracDelta<F>& a) {

    return CIP(a, b);
    
  }

  // ==== Explicit Decralation ====
  template class DiracDelta<CD>;
  template CD CIP(const DiracDelta<CD>& a, const CSTO& b);
  template CD CIP(const CSTO& b, const DiracDelta<CD>& a);
}
