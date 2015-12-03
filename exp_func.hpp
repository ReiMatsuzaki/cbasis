#ifndef EXP_FUNC_TEMPLATE_H
#define EXP_FUNC_TEMPLATE_H

//#include <type_traits>
//#include <math.h>
#include <iostream>
#include <complex>
#include "linspace.hpp"
#include "func.hpp"
//#include <boost/share_ptr.hpp>

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
    class Impl {
    public:
      F c_;
      int n_;
      F z_;
      Impl() : c_(0), n_(0), z_(0) {}
      Impl(F c, int n, F z) : c_(c), n_(n), z_(z) {}
      Impl(const Impl& o) : c_(o.c_), n_(o.n_), z_(o.z_) {}
    };

    Impl *impl_;
    
  public:
    // ---- Constructors ----
    ExpFunc() : impl_(new Impl()) {};
    ExpFunc(F _c, int _n, F _z): impl_(new Impl(_c, _n, _z)) {}
    ExpFunc(const ExpFunc<F,m>& o): impl_(new Impl(o.c(), o.n(), o.z())) {}
    template<class F2>
    ExpFunc(const ExpFunc<F2,m>& o): impl_(new Impl(o.c(), o.n(), o.z())) {}
    ~ExpFunc() {
      delete impl_;
    }
    ExpFunc<F,m>& operator=(const ExpFunc<F,m>& o) {
      if(this == &o) return *this;
      delete impl_;
      impl_ = new Impl(*o.impl_);
      return *this;
    }

    // ---- Accessor ----
    F c() const { return impl_->c_; }
    int n() const {return impl_->n_; }
    F z() const { return impl_->z_; }
    
    void set_c(F c) { impl_->c_ = c; }
    void set_n(int n) { impl_->n_ = n;}
    void set_z(F z) { impl_->z_ = z;}
    
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
