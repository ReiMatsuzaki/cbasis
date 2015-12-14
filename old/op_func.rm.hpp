#ifndef OP_FUNC_TEMPLATE_H
#define OP_FUNC_TEMPLATE_H

/**
   Operand operator OpT to function FuncA.
 */

#include "op.hpp"

#include "exp_func.hpp"
#include "cut_exp.hpp"
#include "lin_func.hpp"

namespace l2func {
  
  // ==== r^m operator ====
  template<class FuncRes, class OpT, class FuncA>
  FuncRes _OP(const OpT& op, const FuncA& a, 
	      rm_tag, func_tag,
	      typename boost::disable_if< is_compound<FuncA> >::type* = 0) {
    typedef typename FuncA::Field F;
    LinFunc<FuncA> res;
    FuncA aa(a); aa.SetRmProd(op.m());
    res.Add(F(1), aa);
    return res;
  }

  // ==== d/dr operator ====
  template<class FuncA>
  LinFunc<FuncA> OperateD1_ExpFunc(const FuncA& a) {
    // c(r^n)e^(-zr^m)  -> { n r^(-1) - z*m r^(m-1) } u
    LinFunc<FuncA> res;
    
    typedef typename FuncA::Field F;

    int n = a.n();
    int m = a.exp_power;
    
    FuncA a1(a); a1.SetRmProd(-1);
    res.Add(F(n), a1);

    FuncA a2(a); 
    if(m != 1 )
      a2.SetRmProd(m-1);
    res.Add(-a.z()*F(m), a2);

    return res;
  }

  template<class FuncRes, class OpT, class FuncA>
  LinFunc<FuncA> _OP(const OpT&, const FuncA& a,
		     d1_tag, exp_func_tag) {
    return OperateD1_ExpFunc(a);
  }

  template<class FuncRes, class OpT, class FuncA>
  LinFunc<FuncA> _OP(const OpT&, const FuncA& a,
		     d1_tag, cut_exp_tag) {
    return OperateD1_ExpFunc(a);
  }


  // ==== d2/dr2 operator ====
  template<class FuncA>
  LinFunc<FuncA> OperateD2_ExpFunc(const FuncA& a) {
    // r^n e^(-zr^m) -> n r^(n-1) -zm r^(n+m-1)
    //               -> n(n-1)r^(n-2) -zm(n +n+m-1)r^(n+m-2) +zzmm r^(n+2m-2)

    LinFunc<FuncA> res;

    typedef typename FuncA::Field F;

    int n = a.n();
    int m = a.exp_power;
    F z = a.z();
    
    if( n != 1) {
      FuncA a1(a); a1.SetRmProd(-2);
      res.Add(n*n-n, a1);
    }

    FuncA a2(a); a2.SetRmProd(m-2);
    res.Add(-m*z*(2*n+m-1), a2);

    FuncA a3(a); a3.SetRmProd(2*m-2);
    res.Add(z*z*m*m, a3);

    return res;
  }

  template<class FuncRes, class OpT, class FuncA>
  FuncRes _OP(const OpT&, const FuncA& a,
	      d2_tag, exp_func_tag) {
    return OperateD2_ExpFunc(a);
  }

  template<class FuncRes, class OpT, class FuncA>
  LinFunc<FuncA> _OP(const OpT&, const FuncA& a,
		     d2_tag, cut_exp_tag) {
    return OperateD2_ExpFunc(a);
  }

  // ==== Scalar Prod of Operator ====
  template<class FuncRes, class OpT, class FuncA>
  FuncRes _OP(const OpT& op,   const FuncA& a,
	      scalar_prod_tag, func_tag) {
    
    typedef typename FuncA::Field F;

    LinFunc<FuncA> res = OP(op.op, a);
    res.SetScalarProd(op.c);

    return res;
  }

  // ==== Add ====
  template<class FuncRes, class OpT, class FuncA>
  FuncRes _OP(const OpT& op, const FuncA& f,
	      op_add_tag,    func_tag) {

    LinFunc<FuncA> Af = OP(op.opA, f);
    LinFunc<FuncA> Bf = OP(op.opB, f);

    Af.Add(Bf);

    return Af;

  }

  // ==== LinFunc ====
  template<class FuncRes, class OpT, class FuncA>
  FuncRes _OP(const OpT& o, const FuncA& a,
	      op_tag,       linfunc_tag, 
	      typename boost::disable_if< is_compound<OpT> >::type* = 0) {

    typedef typename FuncA::Field F;
    LinFunc<FuncA> res;
    typedef typename FuncA::const_iterator IT;
    for(IT it = a.begin(); 
	it != a.end(); ++it) {
      F c = it->first;
      FuncA op_i = OP(o, it->second);
      res.Add(c, op_i);
    }

    return res;

  }


  // ==== Interface ====
  template<class OpT, class FuncA>
  LinFunc<FuncA> OP(const OpT& op, const FuncA& a) {
    
    is_l2func<FuncA>(); is_op<OpT>();
    typedef typename FuncA::Field F;
    return _OP<LinFunc<FuncA>, OpT, FuncA>(op, a, 
					   typename op_traits<OpT>::op_tag(),
					   typename func_traits<FuncA>::func_tag());

  }

}

#endif
