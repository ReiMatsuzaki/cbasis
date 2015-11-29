#ifndef CIP_IMPL_TEMPLATE_H
#define CIP_IMPL_TEMPLATE_H

#include "cip.hpp"

#include "func.hpp"
#include "exp_func.hpp"
#include "cut_exp.hpp"
#include "delta.hpp"
#include "lin_func.hpp"

#include "op.hpp"

//#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace l2func {

  typedef std::complex<double> CD;

  // ==== Calculation ====
  template<class F> F STO_Int(F z, int n);
  template<class F> F GTO_Int(F z, int n);
  template<class F> F STO_GTO_Int(F as, F ag, int n);
  template<class F> F ExpExp_Int(F zA, F zB, int n, sto_tag, sto_tag) {
    return STO_Int(zA+zB, n);
  }
  template<class F> F ExpExp_Int(F zA, F zB, int n, sto_tag, gto_tag) {
    return STO_GTO_Int(zA, zB, n);
  }
  template<class F> F ExpExp_Int(F zA, F zB, int n, gto_tag, sto_tag) {
    return STO_GTO_Int(zB, zA, n);
  }
  template<class F> F ExpExp_Int(F zA, F zB, int n, gto_tag, gto_tag) {
    return GTO_Int(zB+zA, n);
  }

  template<class F> F CutSTO_Int(F z, int n, double r0);


  // ==== Inner Product ====
  // ---- STO/GTO ----
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, exp_func_tag, exp_func_tag) {

    return a.c() * b.c() * 
      ExpExp_Int(a.z(), b.z(), a.n() + b.n(),
		 typename func_traits<A>::func_tag(), 
		 typename func_traits<B>::func_tag());

  }

  // ---- CutSTO/CutGTO ----
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, cut_sto_tag, cut_sto_tag) {
    double r0 = a.r0() < b.r0() ? a.r0() : b.r0();
    return a.c() * b.c() * CutSTO_Int(a.z()+b.z(), a.n()+b.n(), r0);
  }
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, cut_sto_tag, sto_tag) {
    double r0 = a.r0();
    return a.c() * b.c() * CutSTO_Int(a.z()+b.z(), a.n()+b.n(), r0);
  }
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, sto_tag, cut_sto_tag) {
    return _CIP<F, B, A>(b, a, cut_sto_tag(), sto_tag());
  }

  // ----CIP for Delta ----
  template<class F, class A, class B> 
  F _CIP(const A& a, const B& b, delta_tag, func_tag) {

    return b.at(a.r0());

  }
  template<class F, class A, class B> 
  F _CIP(const A& a, const B& b, func_tag, delta_tag) {
    return _CIP<F, B, A>(b, a, delta_tag(), func_tag());
  }

  // ---- CIP for LinFunc ----
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, linfunc_tag, func_tag, 
	 typename boost::disable_if<
	    boost::is_same<
	    typename func_traits<B>::func_tag, linfunc_tag> >::type* =0
	 ) {
    
    /**  Memo:
	 boost::disable_if is for ambiguous error when A and B are both LinFunc
    */

    F acc(0);
    for(typename A::const_iterator it = a.begin(); it != a.end(); ++it) {

      acc += it->first * CIP(it->second, b);

    }
    return acc;

  }
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, func_tag, linfunc_tag) {

    F acc(0);
    for(typename B::const_iterator it = b.begin(); it != b.end(); ++it) {

      acc += it->first * CIP(it->second, a);

    }
    return acc;

  }

  // ---- static interface ----
  template<class A, class B>
  typename A::Field CIP(const A& a, const B& b) {
    
    is_l2func<A>(); is_l2func<B>();
    typedef typename A::Field FA;
    typedef typename B::Field FB;
    //    BOOST_STATIC_ASSERT(boost::has_logical_not<boost::is_same<FA, FB> >::value, "Field of A and B must be same");
			
			
    return _CIP<FA, A, B>(a, b,
			  typename func_traits<A>::func_tag(), 
			  typename func_traits<B>::func_tag());

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


  // ==== Matrix Element of Operator ====
  // ---- Rm ----
  template<class ExpFuncA, class ExpFuncB>
  typename ExpFuncA::Field RmMat(const ExpFuncA& a, const ExpFuncB& b) {

    

  }


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
	    sto_tag,    rm_tag,     gto_tag) {

    return STO_GTO_Int(a.z(), b.z(), a.n() + b.n()) * a.c() * b.c();
    
  }
  template<class F, class A, class O, class B>
  F _CIP_op(const A& a, const O& o, const B& b, 
	    gto_tag,    rm_tag,     sto_tag) {

    return _CIP(b, o, a, sto_tag(), rm_tag(), gto_tag());
    
  }

  template<class F, class A, class O, class B>
  F _CIP_op(const A& a, const O& o, const B& b, 
	    cut_sto_tag, rm_tag,    cut_sto_tag) {
    double r0 = a.r0() < b.r0() ? a.r0() : b.r0();
    return a.c() * b.c() * CutSTO_Int(a.z()+b.z(), a.n()+b.n()+o.m(), r0);
  }
  template<class F, class A, class O, class B>
  F _CIP_op(const A& a, const O& o, const B& b, 
	    cut_sto_tag, rm_tag,    sto_tag) {
    double r0 = a.r0();
    return a.c() * b.c() * CutSTO_Int(a.z()+b.z(), a.n()+b.n()+o.m(), r0);
  }
  template<class F, class A, class O, class B>
  F _CIP_op(const A& a, const O& o, const B& b, 
	    sto_tag,    rm_tag,     cut_sto_tag) {
    return _CIP_op(b, o, a, cut_sto_tag(), rm_tag(), sto_tag());
  }

  // ---- D1 ----
  template<class ExpFuncA, class ExpFuncB>
  typename ExpFuncA::Field D1Mat(const ExpFuncA& a, const ExpFuncB& b) {
    
    typedef typename ExpFuncA::Field F;
    int n = b.n();
    int m = b.exp_power;

    return F(n) * CIP(a, OpRm(-1), b) + F(-b.z()*m) * CIP(a, OpRm(m-1), b);

  }

  template<class F, class A, class O, class B>
  F _CIP_op(const A& a,   const O&, const B& b,
	    exp_func_tag, d1_tag,   exp_func_tag) {

    return D1Mat(a, b);

  }

  template<class F, class A, class O, class B>
  F _CIP_op(const A& a, const O&, const B& b,
	    cut_exp_tag, d1_tag, cut_exp_tag) {

    return D1Mat(a, b);
    
}

  // ---- D2 ----
  template<class F, class A, class O, class B>
  F _CIP_op(const A& a, const O&, const B& b, 
	    sto_tag,    d2_tag,   sto_tag) {
    F c = a.c() * b.c();
    F z = a.z() + b.z();
    int n = a.n() + b.n();
  
    F acc(0);
    
    if(b.n() > 1)
      acc += F(b.n()*b.n()-b.n()) * STO_Int(z, n-2);
    
    if(b.n() > 0)
      acc += -F(2*b.n())*b.z() * STO_Int(z, n-1);
    
    acc += b.z() * b.z() * STO_Int(z, n);
    
    acc *= c;
    
    return acc;
  
  }

  // ---- linfunc ----
  template<class F, class A, class O, class B> 
  F _CIP_op(const A& a,  const O& o, const B& b,
	    linfunc_tag, op_tag,     func_tag,
	    typename boost::disable_if<
	    boost::is_same<
	    typename func_traits<B>::func_tag, linfunc_tag> >::type* =0) {

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
  return _CIP_op<FA, A, O, B>(a, o, b, 
			      typename func_traits<A>::func_tag(),
			      typename op_traits<O>::op_tag(),
			      typename func_traits<B>::func_tag());

}

}
#endif
