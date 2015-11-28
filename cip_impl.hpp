#ifndef CIP_IMPL_TEMPLATE_H
#define CIP_IMPL_TEMPLATE_H

#include "cip.hpp"
#include "delta.hpp"
#include "exp_func.hpp"

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

  //  template double CIP(const DiracDelta<double>&, const RSTO&);

  template<class F, class A, class TagA>
  typename DiracDelta<F>::Field CIP(const A& a, DiracDelta<F>, TagA, delta_tag);
  template<class F, class A, class TagA>
  typename DiracDelta<F>::Field CIP(DiracDelta<F>, const A& a, delta_tag,  TagA);

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
