#ifndef EXP_FUNC_TEMPLATE_H
#define EXP_FUNC_TEMPLATE_H

#include <iostream>
#include <complex>
#include "linspace.hpp"
#include "func.hpp"

namespace l2func {

  // ==== STO or GTO ====
  template<class F, int m>
  class ExpFunc {
  public:
    // ---- type ----
    typedef F Field;
    enum EExpPower { exp_power=m };
    typedef ExpFunc<F,m> FuncDerivOne;
    typedef ExpFunc<F,m> FuncDerivTwo;

  private:    
    // ---- Member Field ----
    F c_;
    int n_;
    F z_;

  public:
    // ---- Constructors ----
    ExpFunc();
    ExpFunc(F c, int n, F z);
    ExpFunc(const ExpFunc<F, m>& o);
    template<class F2>
    ExpFunc(const ExpFunc<F2,m>& o): c_(o.c()), n_(o.n()), z_(o.z()) {}
    ~ExpFunc();

    // ---- Accessor ----
    F c() const { return this->c_; }
    int n() const {return this->n_; }
    F z() const { return this->z_; }
    
    void set_c(F c) { this->c_ = c; }
    void set_n(int n) { this->n_ = n;}
    void set_z(F z) { this->z_ = z;}
    
    // ---- Method ----
    F at(F x) const;
    std::string str() const;
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
  struct exp_func_tag :public func_tag {};
  template<class F, int m> 
  struct is_fundamental<ExpFunc<F, m> > : public boost::true_type {};
  template<class F, int m> 
  struct is_compound<ExpFunc<F, m> > : public boost::false_type {};

  // ==== STO ====
  struct sto_tag :public exp_func_tag {};
  template<class F> struct func_traits<ExpFunc<F, 1> > {
    typedef sto_tag func_tag;
  };
  template<class F> struct is_l2func<ExpFunc<F, 1> > {};

  typedef ExpFunc<double, 1> RSTO;
  typedef ExpFunc<std::complex<double>, 1> CSTO;

  // ==== GTO ====
  struct gto_tag :public exp_func_tag {};
  template<class F> struct func_traits<ExpFunc<F, 2> > {
    typedef gto_tag func_tag;
  };
  template<class F> struct is_l2func<ExpFunc<F, 2> > {};
  
  typedef ExpFunc<double, 2> RGTO;
  typedef ExpFunc<std::complex<double>, 2> CGTO;
}

#endif
