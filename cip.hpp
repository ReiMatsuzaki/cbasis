#ifndef CIP_TEMPLATE_H
#define CIP_TEMPLATE_H

namespace l2func {

  // ==== Inner product  ====
  template<class FuncA, class FuncB>
  typename FuncA::Field CIP(const FuncA& a, const FuncB& b);

  // ==== Other Function ====
  template<class FuncA> typename FuncA::Field CNorm(const FuncA& a);
  template<class FuncA> void CNormalize(FuncA *a);
  
  // ==== matrix element of operator ====
  template<class FuncA, class OpT, class FuncB>
  typename FuncA::Field CIP(const FuncA& a, const OpT& o, const FuncB& b);

}
#endif
