#ifndef H_TEMPLATE_H
#define H_TEMPLATE_H

#include <iostream>
#include <complex>
// #include <boost/array.hpp>

#include "math_utils.hpp"
#include "linspace.hpp"
#include "op.hpp"

typedef l2func::array3<double> d3;
typedef l2func::array3<int>    i3;
typedef l2func::array3<std::complex<double> > c3;

namespace l2func {

  template<class F, class FC>
  class CartGTO :public Func<F, l2func::array3<FC> > {
  public:
    // ---- type ----
    typedef F Field;
    typedef FC FieldCoord;
    typedef l2func::array3<FC> FC3;
    
  private:
    // ---- Member Field ----
    F c_;
    i3 nml_;
    FC3 xyz_;
    F z_;

  public:
    // ---- Constructors ----
    CartGTO(F c, const i3& nml, const FC3& xyz, F z) : c_(c), nml_(nml), xyz_(xyz), z_(z) {}
    template<class F2, class FC2>
    CartGTO(const CartGTO<F2, FC2>& o): c_(o.c()), nml_(o.nml()), xyz_(o.xyz()), z_(o.z()) {}
    
    // ---- Accessors ----
    F c() const { return this->c_; }
    i3 nml() const { return this->nml_; }
    d3 xyz() const { return this->xyz_; }
    F z() const { return this->z_; }

    void set_c(F c) {this->c_ = c; } 
    void set_nml(i3 nml) {this->nml_ = nml; } 
    void set_xyz(d3 xyz) {this->xyz_ = xyz; } 
    void set_z(F z) {this->z_ = z; } 

    // ---- Method ----
    F at(FC3 xyz) const;
    std::string str() const;
    
    void SetComplexConjugate();
    void SetScalarProd(F c);
  };

  // ==== External ====
  template<class F, class FC>
  std::ostream& operator << (std::ostream& os, const CartGTO<F, FC>& a);

  typedef CartGTO<double, double> CartRGTO;
  typedef CartGTO<std::complex<double>, double> CartCGTO;

  // ==== Operator ====
  class OpKE {};
  template<class F>
  class OpNA {
  private:
    F z_;
    array3<F> xyz_;
  public:
    OpNA(F z, const array3<F>& xyz): z_(z), xyz_(xyz) {}
    F z() const { return z_;}
    array3<F> xyz() const { return xyz_;}
  };

  // ==== Inner product ====

  template<class F, class FC>
  F CIP(const CartGTO<F, FC>& a, const CartGTO<F, FC>& b) {
    return 1.10;
  }
  template<class F, class FC>
  F CIP(const CartGTO<F, FC>& a, const OpKE, const CartGTO<F, FC>& b) {
    return 1.10;
  }
  template<class F, class FC>
  F CIP(const CartGTO<F, FC>& a, const OpNA<FC>& v, const CartGTO<F, FC>& b) {
    return 1.10;
  }

}

#endif
