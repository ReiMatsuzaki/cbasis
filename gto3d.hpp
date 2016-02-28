#ifndef H_TEMPLATE_H
#define H_TEMPLATE_H

#include <iostream>
#include <complex>
#include <sstream>
// #include <boost/exception/exception.hpp>
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>

#include "math_utils.hpp"
#include "linspace.hpp"
#include "op.hpp"
#include "cints.hpp"

typedef l2func::array3<double> d3;
typedef l2func::array3<int>    i3;
typedef l2func::array3<std::complex<double> > c3;

namespace l2func {

  // ==== Cartesian GTO ====
  template<class F, class FC>
  class CartGTO :public Func<F, FC> {
  public:
    // ---- type ----
    typedef F Field;
    typedef FC FieldCoord;
    //typedef l2func::array3<FC> FC3;
    
  private:
    // ---- Member Field ----
    F c_;
    i3 nml_;
    FC xyz_;
    F zeta_;

  public:
    // ---- Constructors ----
    CartGTO(F c, const i3& nml, const FC& xyz, F z) : c_(c), nml_(nml), xyz_(xyz), zeta_(z) {}
    template<class F2, class FC2>
    CartGTO(const CartGTO<F2, FC2>& o): c_(o.c()), nml_(o.nml()), xyz_(o.xyz()), zeta_(o.z()) {}
    
    // ---- Accessors ----
    F c() const { return this->c_; }
    const i3& nml() const { return this->nml_; }
    int n() const {return this->nml_[0];}
    int m() const {return this->nml_[1];}
    int l() const {return this->nml_[2];}
    const FC& xyz() const { return this->xyz_; }
    FC x() const {return this->xyz_[0];}
    FC y() const {return this->xyz_[1];}
    FC z() const {return this->xyz_[2];}
    F zeta() const { return this->zeta_; }

    void set_c(F c) {this->c_ = c; } 
    void set_nml(i3 nml) {this->nml_ = nml; } 
    void set_xyz(FC xyz) {this->xyz_ = xyz; } 
    void set_z(F z) {this->z_ = z; } 

    // ---- Method ----
    F at(FC xyz) const;
    std::string str() const;
    
    void SetComplexConjugate() {
      this->c_ = conjug(this->c_);
      this->z_ = conjug(this->z_);
    }
    void SetScalarProd(F c) {
      this->c_ *= c;
    }
    void SetNormalize() {
      F norm2 = CIP(*this, *this);
      this->SetScalarProd(F(1)/sqrt(norm2));
    }
  };

  // ==== Spherical GTO ====
  /*
  template<class F, class FC>
  class SpheGTO :public Func<F, l2func::array3<FC> > {
  public:
    // ---- type ----
    typedef F Field;
    typedef FC FieldCoord;
    typedef l2func::array3<FC> FC3;

  private:
    // ---- Member Field ----
    F c_;
  };
  */
  // ---- Exception class ----
  class ExceptionBadYlm : public std::exception, public boost::exception {
  public:
    int L_, M_;
    std::string msg_;
    ExceptionBadYlm(int L, int M): std::exception(), L_(L), M_(M) {
      stringstream ss;
      ss << "Unphysical (L, M) pair. (L, M) = (" << L_ << ", " << M_ << ")";
      msg_ = ss.str();
    }
    ~ExceptionBadYlm() throw() {}
    virtual const char* what() const throw() {
      return msg_.c_str();
    }

  };
  //  typedef boost::error_info<struct tag_errno_code, int> error_code;

  template<class Field, class Coord>
  void SetSphericalGTO(int L, int M, Coord xyz, Field zeta,
		       LinFunc<CartGTO<Field, Coord> >* target) {
    
    if(L == 0 && M == 0) {
      CartGTO<Field, Coord> gto(Field(1), i3(0, 0, 0), xyz, zeta);
      gto.SetNormalize();
      target->Add(Field(1), gto);
    } else if(L == 1 && M == 0) {
      CartGTO<Field, Coord> g3(Field(1), i3(0, 0, 1), xyz, zeta);
      g3.SetNormalize(); target->Add(1, g3);
    } else if(L == 1 && M == 1) {
      CartGTO<Field, Coord> g1(1, i3(1, 0, 0), xyz, zeta);
      g1.SetNormalize(); target->Add(1, g1);
    } else if(L == 1 && M == -1) {
      CartGTO<Field, Coord> g1(Field(1), i3(0, 1, 0), xyz, zeta);
      g1.SetNormalize(); target->Add(1, g1);
    } else if(L == 2 && M == 1) {
      CartGTO<Field, Coord> g(Field(1), i3(1, 0, 1), xyz, zeta); 
      g.SetNormalize(); target->Add(Field(1), g);
    } else if(L == 2 && M == -1) {
      CartGTO<Field, Coord> g(Field(1), i3(0, 1, 1), xyz, zeta); 
      g.SetNormalize(); target->Add(Field(1), g);
    } else if(L == 2 && M == 2) {
      CartGTO<Field, Coord> g1(Field(1), i3(2, 0, 0), xyz, zeta); 
      g1.SetNormalize(); target->Add(Field(1), g1);
      CartGTO<Field, Coord> g2(Field(1), i3(0, 2, 0), xyz, zeta); 
      g2.SetNormalize(); target->Add(-Field(1), g2);	
    } else if(L == 2 && M == -2) {
      CartGTO<Field, Coord> g1(Field(1), i3(1, 1, 0), xyz, zeta); 
      g1.SetNormalize(); target->Add(1, g1);
    } else if(L == 2 && M == 0) {
      CartGTO<Field, Coord> g1(Field(1), i3(2, 0, 0), xyz, zeta); 
      g1.SetNormalize(); target->Add(-Field(1), g1);
      CartGTO<Field, Coord> g2(Field(1), i3(0, 2, 0), xyz, zeta); 
      g2.SetNormalize(); target->Add(-Field(1), g2);
      CartGTO<Field, Coord> g3(Field(1), i3(0, 0, 2), xyz, zeta); 
      g3.SetNormalize(); target->Add(Field(2), g3);
    } else {
      BOOST_THROW_EXCEPTION(ExceptionBadYlm(L, M));
    }
    target->SetNormalize();
  }
  
  // ==== External ====
  template<class F, class FC>
  std::ostream& operator << (std::ostream& os, const CartGTO<F, FC>& a);

  typedef CartGTO<double, double> CartRGTO;
  typedef CartGTO<std::complex<double>, double> CartCGTO;

  // ==== Operator ====
  template<class Field, class Coord>
  class OpKE : public Op<Field, Coord> {};
  template<class Field, class Coord>
  class OpNA :public Op<Field, Coord>  {
  private:
    Field q_;
    Coord xyz_;
  public:
    OpNA(Field q, const Coord& xyz): q_(q), xyz_(xyz) {}
    Field q() const { return q_;}
    const Coord& xyz() const { return xyz_;}
  };
  template<class Field, class Coord>
  class OpXyz :public Op<Field, Coord>  {
  private:
    i3 nml_;
  public:
    OpXyz(int n, int m, int l) : nml_(n, m, l) {}
    const i3& nml() const { return nml_; }
    int n() const { return nml_[0]; }
    int m() const { return nml_[1]; }
    int l() const { return nml_[2]; }
  };

  // ==== Inner product ====
  template<class F, class FC>
  F CIP_impl_prim(const CartGTO<F, FC>& a, const One&, const CartGTO<F, FC>& b) {
    return a.c() * b.c() * overlap(a.zeta(),
				   a.nml()[0], a.nml()[1], a.nml()[2],
				   a.xyz()[0], a.xyz()[1], a.xyz()[2],
				   b.zeta(),
				   b.nml()[0], b.nml()[1], b.nml()[2],
				   b.xyz()[0], b.xyz()[1], b.xyz()[2]);
  }
  template<class F, class FC>
  F CIP_impl_prim(const CartGTO<F, FC>& a, const OpKE<F, FC>&, const CartGTO<F, FC>& b) {
    return a.c() * b.c() * kinetic(a.zeta(),
				   a.nml()[0], a.nml()[1], a.nml()[2],
				   a.xyz()[0], a.xyz()[1], a.xyz()[2],
				   b.zeta(),
				   b.nml()[0], b.nml()[1], b.nml()[2],
				   b.xyz()[0], b.xyz()[1], b.xyz()[2]);
  }
  template<class F, class FC>
  F CIP_impl_prim(const CartGTO<F, FC>& a, const OpNA<F, FC>& v, const CartGTO<F, FC>& b) {

    return v.q() *  nuclear_attraction(a.xyz()[0], a.xyz()[1], a.xyz()[2],
				       a.c(),
				       a.nml()[0], a.nml()[1], a.nml()[2],
				       a.zeta(),
				       b.xyz()[0], b.xyz()[1], b.xyz()[2],
				       b.c(),
				       b.nml()[0], b.nml()[1], b.nml()[2],
				       b.zeta(),
				       v.xyz()[0], v.xyz()[1], v.xyz()[2]);
  }
  template<class F, class FC>
  F CIP_impl_prim(const CartGTO<F, FC>& a, const OpXyz<F, FC> & op, const CartGTO<F, FC>& b) {
    return a.c() * b.c() * overlap(a.zeta(),
				   a.nml()[0] + op.nml()[0],
				   a.nml()[1] + op.nml()[1],
				   a.nml()[2] + op.nml()[2],
				   a.xyz()[0], a.xyz()[1], a.xyz()[2],
				   b.zeta(),
				   b.nml()[0], b.nml()[1], b.nml()[2],
				   b.xyz()[0], b.xyz()[1], b.xyz()[2]);    
  }
}

#endif
