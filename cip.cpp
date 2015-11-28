//#include <boost/static_assert.hpp>
// #include <boost/type_traits.hpp>
#include "cip.hpp"
#include "lgamma.hpp"
#include "erfc.hpp"
#include "fact.hpp"

#include "func.hpp"
#include "exp_func.hpp"
#include "cut_exp.hpp"

using namespace erfc_mori;
using namespace std;

namespace l2func {

  // ==== Calculation ====
  const CD operator*(const CD& a, int b) {
    return CD(a.real() * b, a.imag() * b); }
  const CD operator*(int b, const CD& a) {
    return CD(a.real() * b, a.imag() * b); }  

  int int_pow(int base, unsigned int expo) {

  int acc = 1;
  for(unsigned int i = 0; i < expo; i++) 
    acc *= base;
  return(acc);
    
}
  template<class F> F STO_Int(F z, int n) {
    return pow(z, -n-1.0) * (1.0 * fact::Factorial(n));
  }
  template<class F> F GTO_Int(F z, int n) {
    
    F res;
    if(n % 2 == 0) {
      
    int nn = n/2;
    res = fact::DoubleFactorial(2*nn-1) * sqrt(M_PI) /
      (F(int_pow(2, nn+1)) * pow(sqrt(z), 2*nn+1));
    
    } else {
      
      int nn = (n-1)/2;
      res = F(fact::Factorial(nn)) / (F(2) * pow(z, nn+1));
    
    }
    return res;
  } 
  template<class F> F sto_gto_int_0(F as, F ag) {
    F erfcVal, expVal, sqrtPi,pi,res;
    ErfcCalcData data;
    Erfc(as/(2*sqrt(ag)),erfcVal,data);
    expVal=exp(as*as/(4*ag));
    pi=M_PI;
    sqrtPi=sqrt(pi);
    res = (erfcVal*expVal*sqrtPi)/(2*sqrt(ag));
    return (res);
  }
  template<class F> F sto_gto_int_1(F as, F ag) {
    F exp2erfc, sqrtPi, sqrt_ag, pi, res;
    ErfcCalcData data;
    sqrt_ag = sqrt(ag);
    pi = M_PI;
    sqrtPi = sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfc, data);

    res = F(1)/(2*ag) - (as * exp2erfc * sqrtPi) /
      (4*ag*sqrt_ag);  
  
    return(res);
  }
  template<class F> F sto_gto_int_2(F as, F ag) {
    F erfcVal, expVal, sqrtPi,pi,res;
    ErfcCalcData data;
    Erfc(as/(2*sqrt(ag)),erfcVal,data);
    expVal=exp(as*as/(4*ag));
    pi=M_PI;
    sqrtPi=sqrt(pi);
    res = (-2*sqrt(ag)*as + (2*ag + pow(as,2))*erfcVal*expVal*sqrtPi)/(8*pow(sqrt(ag),5));

    return (res);
  }
  template<class F> F sto_gto_int_3(F as, F ag) {
    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfcVal, data);
    
    res = (4*ag + pow(as,2))/(8*pow(ag,3)) 
      -(as*(6*ag + pow(as,2)) * exp2erfcVal * sqrtPi) / 
      (16 * sqrt_ag * ag * ag * ag);

    return (res);
  }
  template<class F> F sto_gto_int_4(F as, F ag) {
    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfcVal, data);
    
    res = (-2*sqrt_ag*as*(10*ag + as*as)
	   + (12*ag*ag + 12*ag*as*as + as*as*as*as)*exp2erfcVal*sqrtPi)/(32*ag*ag*ag*ag*sqrt_ag);
    
    return (res);
  }
  template<class F> F sto_gto_int_5(F as, F ag) {
    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    Exp2Erfc(as/(2*sqrt(ag)), exp2erfcVal, data);
    pi=M_PI;
    sqrtPi=sqrt(pi);
    sqrt_ag = sqrt(ag);
    F x = sqrt_ag;
 
    res = (2*sqrt_ag*(2*ag + as*as)*(16*ag + as*as) -
	   as*(60*ag*ag + 20*ag*as*as + as*as*as*as)*
	   exp2erfcVal*sqrtPi)/
      (64*x*ag*ag*ag*ag*ag);
    return (res);
  }
  template<class F> F sto_gto_int_6(F as, F ag) {

    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    F a = as/(2*sqrt_ag);
    Exp2Erfc(a, exp2erfcVal, data);

    if(! data.convergence) {
      string msg("is not convergence in sto_gto_int_6");
      throw msg;
    }
    
    res = (-2*sqrt_ag*as*(6*ag + as*as)*(22*ag + as*as) +
	   (120*ag*ag*ag + 180*ag*ag*as*as + 30*ag*pow(as,4) + 
	  pow(as,6))*exp2erfcVal*sqrtPi)/
      (128*sqrt_ag*ag*ag*ag*ag*ag*ag);

    return (res);
  }
  template<class F> F sto_gto_int_7(F as, F ag) {
    
    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfcVal, data);
  
    res = (2*sqrt(ag)*(384*pow(ag,3) + 348*pow(ag,2)*pow(as,2) + 40*ag*pow(as,4) + pow(as,6)) - as*(840*pow(ag,3) + 420*pow(ag,2)*pow(as,2) + 42*ag*pow(as,4) + pow(as,6))*exp2erfcVal*sqrtPi)/
    (256*sqrt_ag * ag* ag* ag* ag* ag* ag* ag);

    return (res);
  }
  template<class F> F sto_gto_int_8(F as, F ag) {
    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfcVal, data);
    
    res = (-2*sqrt(ag)*as*(2232*pow(ag,3) + 740*pow(ag,2)*pow(as,2) + 54*ag*pow(as,4) + pow(as,6)) + (1680*pow(ag,4) + 3360*pow(ag,3)*pow(as,2) + 840*pow(ag,2)*pow(as,4) + 56*ag*pow(as,6) + pow(as,8))*exp2erfcVal*sqrtPi)/(512*sqrt_ag*pow(ag,8));

    return (res);
  }
  template<class F> F sto_gto_int_9(F as, F ag) {

    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfcVal, data);

    res = (2*sqrt_ag*(6144*pow(ag,4) + 7800*pow(ag,3)*pow(as,2) + 1380*pow(ag,2)*pow(as,4) + 70*ag*pow(as,6) + pow(as,8)) - as*(15120*pow(ag,4) + 10080*pow(ag,3)*pow(as,2) + 1512*pow(ag,2)*pow(as,4) + 72*ag*pow(as,6) + pow(as,8))*exp2erfcVal*sqrtPi)/(1024*sqrt_ag*pow(ag,9));


    return (res);
  }
  template<class F> F STO_GTO_Int(F as, F ag, int n) {
    F res;
    switch(n) {
    case 0:
      res = sto_gto_int_0(as, ag);
      break;
    case 1:
      res = sto_gto_int_1(as, ag);
      break;
    case 2:
      res = sto_gto_int_2(as, ag);
      break;
    case 3:
      res = sto_gto_int_3(as, ag);
      break;
    case 4:
      res = sto_gto_int_4(as, ag);
      break;
    case 5:
      res = sto_gto_int_5(as, ag);
      break;
    case 6:
      res = sto_gto_int_6(as, ag);
      break;
    case 7:
      res = sto_gto_int_7(as, ag);
      break;
    case 8:
      res = sto_gto_int_8(as, ag);
      break;
    case 9:
      res = sto_gto_int_9(as, ag);
      break;      
    default:
      string msg;
      msg = "this is not supported in STO_GTO_int";
      throw msg;
    }

    return res;
  }
  template<class F> F CutSTO_Int(F z, int n, double r0) {
    return LowerGamma<F>(n+1, z*r0) / pow(z, n+1);
  }

  // ==== Inner Product ====

  // ---- STO/GTO ----
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, sto_tag, sto_tag) {
    return a.c() * b.c() * STO_Int(a.z()+b.z(), a.n()+b.n());
  }
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, gto_tag, gto_tag) {
    return a.c() * b.c() * GTO_Int(a.z()+b.z(), a.n()+b.n());
  }
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, sto_tag, gto_tag) {
    return a.c() * b.c() * STO_GTO_Int(a.z(), b.z(), a.n()+b.n());
  }
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, gto_tag, sto_tag) {
    return CIP(b, a);
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


/*
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, linfunc_tag, func_tag) {
    
    
    F acc(0);
    for(typename A::const_iterator it = a.begin();
	it != a.end();
	it++) {
    
      acc += it->first * CIP(a, it->second, b);
      
    }
    
    return acc;
    
  }
  template<class F, class A, class B>
  F _CIP(const A& a, const B& b, func_tag, linfunc_tag) {
    return CIP(b, a);
  }
*/


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

  // ---- explicit instance ----

  template double CIP(const RSTO&, const RSTO&);
  template std::complex<double> CIP(const CSTO&, const CSTO&);
  template double CIP(const RGTO&, const RGTO&);
  template std::complex<double> CIP(const CGTO&, const CGTO&);
  template double CIP(const RSTO&, const RGTO&);
  template std::complex<double> CIP(const CSTO&, const CGTO&);
  template double CIP(const RGTO&, const RSTO&);
  template std::complex<double> CIP(const CGTO&, const CSTO&);

  template double CIP(const CutRSTO&, const CutRSTO&);
  template std::complex<double> CIP(const CutCSTO&, const CutCSTO&);
  template double CIP(const CutRSTO&, const RSTO&);
  template std::complex<double> CIP(const CutCSTO&, const CSTO&);
  template double CIP(const RSTO&, const CutRSTO&);
  template std::complex<double> CIP(const CSTO&, const CutCSTO&);

}
