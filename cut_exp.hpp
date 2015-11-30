#ifndef CUT_EXP_TEMPLATE_H
#define CUT_EXP_TEMPLATE_H

#include "exp_func.hpp"

/**
   Cutted STO/GTO
*/

namespace l2func {

  template<class F, int m>
  class CutExpFunc {
    
  public:
    // ---- type ----
    typedef F Field;
    enum EExpPower {exp_power=m};
    typedef CutExpFunc<F,m> FuncDerivOne;
    typedef CutExpFunc<F,m> FuncDerivTwo;

  private:
    // ---- Member Field ----
    ExpFunc<F, m> func;
    double r0_;          // R-matrix radius

  public:
    // ---- Constructors ----
    CutExpFunc();
    CutExpFunc(F _c, int _n, F _z, double r0);
    CutExpFunc(const CutExpFunc<F, m>& o);
    template<class U> CutExpFunc(const CutExpFunc<U,m>& o):
      func(o.func), r0_(o.r0_)  {}

    // ---- Accessor ----
    F c() const { return func.c(); }
    int n() const { return func.n(); }
    F z() const { return func.z(); }
    double r0() const { return this->r0_; }
    
    void set_c(F c) { func.set_c(c); }
    void set_n(int n) { func.set_n(n);}
    void set_z(F z) { func.set_z(z);}
    void set_r0(double r0) { this->r0_ = r0;}

    const ExpFunc<F,m>& get_func() const { return func; }
    
    // ---- Method ----
    F at(F x) const;
    void SetComplexConjugate();
    void SetDerivParam();
    void SetScalarProd(F c);
    void SetRmProd(int n);

    FuncDerivOne DerivParamOne();
    FuncDerivTwo DerivParamTwo();
  };


  // ==== External ====
  template<class F, int m>
  std::ostream& operator << (std::ostream& os, const ExpFunc<F,m>& a);

  // ==== STO or GTO ====
  struct cut_exp_tag :public func_tag {};
  template<class F, int m> struct is_fundamental<CutExpFunc<F, m> > {};

  // ==== STO ====
  struct cut_sto_tag :public cut_exp_tag {};
  template<class F> struct func_traits<CutExpFunc<F, 1> > {
    typedef cut_sto_tag func_tag;
  };
  template<class F> struct is_l2func<CutExpFunc<F, 1> > {};

  typedef CutExpFunc<double, 1> CutRSTO;
  typedef CutExpFunc<std::complex<double>, 1> CutCSTO;

  // ==== GTO ====
  struct cut_gto_tag :public cut_exp_tag {};
  template<class F> struct func_traits<CutExpFunc<F, 2> > {
    typedef cut_gto_tag func_tag;
  };
  template<class F> struct is_l2func<CutExpFunc<F, 2> > {};
  
  typedef CutExpFunc<double, 2> CutRGTO;
  typedef CutExpFunc<std::complex<double>, 2> CutCGTO;
}
#endif
