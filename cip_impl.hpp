#ifndef CIP_IMPL_TEMPLATE_H
#define CIP_IMPL_TEMPLATE_H

#include "delta.hpp"


namespace l2func {

  // ==== CIP for Delta ====
  template<class F, class A, class B> 
  F _CIP(const A& a, const B& b, delta_tag, func_tag) {

    return b.at(a.r0());

  }
  template<class F, class A, class B> 
  F _CIP(const A& a, const B& b, func_tag, delta_tag) {
    return _CIP(b, a, delta_tag(), func_tag());
  }

  // ==== normalization ====
  template<class A>
  typename A::Field CNorm(const A& a) {
    return sqrt(CIP(a, a));
  }
  template<class A>
  void CNormalize(A *a) {
    typedef typename A::Field F;
    F cc = CNorm(*a);
    a->SetScalarProd(F(1)/cc);
  }
}
#endif
