#ifndef CIP_TEMPLATE_H
#define CIP_TEMPLATE_H

namespace l2func {

  // ==== Utilities ====
  //  int int_pow(int base, unsigned int expo);
  //  template<class F> F STO_Int(F z, int n);
  //  template<class F> F GTO_Int(F z, int n);

  // ==== Inner product of two function ====
  template<class A, class B>
  typename A::Field CIP(const A& a, const B& b);

  // ==== Other Function ====
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

  // ==== Inner product of operator ====
/*
  template<class F, class A, class O, class B>
  F _CIP_op(const A& a, const O& o, const B& b, 
	    sto_tag,    rm_tag,     sto_tag) {

    return STO_Int(a.z() + b.z(), a.n() + o.m() + b.n()) * a.c() * b.c();
  
  }
  template<class F, class A, class O, class B>
  F _CIP_op(const A& a, const O& o, const B& b, 
	    gto_tag,    rm_tag,     gto_tag) {

    return STO_Int(a.z() + b.z(), a.n() + o.m() + b.n()) * a.c() * b.c();
    
  }
  template<class F, class A, class O, class B>
  F _CIP_op(const A& a, const O& o, const B& b, 
	    sto_tag,    d2_tag,     sto_tag) {
    F c = a.c() * b.c();
    F z = a.z() + b.z();
    int n = a.n() + b.n();
  
    F acc(0);
    
    if(b.n() > 1)
      acc += F(b.n()*b.n()-b.n()) * STO_Int(z, n-2);
    
    if(b.n() > 0)
      acc += -F(b.n())*b.z() * STO_Int(z, n-1);
    
    acc += b.z() * b.z() * STO_Int(z, n);
    
    acc *= c;
    
    return acc;
  
}



  template<class F, class A, class O, class B> 
  F _CIP_op(const A& a,  const O& o, const B& b,
	    linfunc_tag, op_tag,     func_tag) {

  F acc(0);
  for(typename A::const_iterator it = a.begin();
      a != a.end(); 
      a++) {

    acc += it->first * CIP(a->second, o, b);

  }

  return acc;
}
  template<class F, class A, class O, class B> 
  F _CIP_op(const A& a,  const O& o, const B& b,
	    func_tag,     op_tag,     linfunc_tag) {
  return CIP(b, o, a);
}
*/

/*
  template<class F, class A, class B> struct Summer{
  
  // Field
  const A& a_;
  const B& b_;

  // Constructor
  Summer(const A& _a, const B& _b): a_(_a), b_(_b) {}
    
  // function called in boost::fusion::fold
  template<class OpT>
  F operator()(F x0, const std::pair<F, OpT>& co)  const {
    return x0 + co.first * CIP(a_, co.second, b_);
  }
};
  template<class F,  class A, class O, class B>
  F _CIP_op(const A& a, const O& o, const B& b,
	  func_tag, linop_tag, func_tag) {
  Summer<F, A, B> summer(a, b);
  F res = boost::fusion::fold(o.CoefOp(), F(0), summer);
  return res;
  }

  template<class A, class O, class B>
  typename A::Field CIP(const A& a, const O& o, const B& b) {

  is_l2func<A>(); is_l2func<B>(); is_op<O>();
  typedef typename A::Field FA;
  typedef typename B::Field FB;
  static_assert(std::is_same<FA, FB>::value, 
		"Field of A and B must be same");
  return _CIP_op<FA, A, O, B>(a, o, b, 
			      typename func_traits<A>::func_tag(),
			      typename op_traits<O>::op_tag(),
			      typename func_traits<B>::func_tag());

}
  */
}
#endif
