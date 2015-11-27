#ifndef CUT_PRIME_TEMPLATE_H
#define CUT_PRIME_TEMPLATE_H

#include "delta.hpp"
#include "prim.hpp"
#include "delta.hpp"

namespace l2func {

  // ============= Inner Product ==============
  template<class F> F LowerGamma(int n, F z);
  template<class F> F CutSTO_Int(F z, int n, double r0);
  //  template<class F> F CutGTO_Int(F z, int n, double r0);

  // ============= Cutted STO or GTO =============
  // m==1 => STO, m==2=>GTO,
  // value is 0 on (r0, oo) region
  // this function is used in R-matrix theory
  template<class F, int m>
  class CutExpBasis {

  public:
    // ------- typedef -------
    typedef F Field;
    enum EExpPower {exp_power=m};

  private:
    // ------- Field Member -----
    ExpBasis<F, m>  basis;
    double r0_;              // R-matrix radius

  public:
    // ------- Big3 -------
    CutExpBasis();
    CutExpBasis(F _c, int _n, F _z, double r0);
    CutExpBasis(int _n, F _z, double r0, ENormalized);

    template<class U> CutExpBasis(const CutExpBasis<U,m>& o):
      r0_(o.r0_), basis(o.basis) {}

    // ------- Accessor ------
    F c() const { return this->basis.c();}
    int n() const { return basis.n();}
    F z() const { return basis.z();}
    double r0() const { return this->r0_; }
    void set_z(F z) { basis.set_z(z);}
    F& ref_z() { return basis.ref_z();}
    F at(F x) const;

    // ------- Operation ------
    void SetComplexConjugate(const CutExpBasis<F,m>& o);
    CutExpBasis<F, m> ComplexConjugate() const;    
  };

  // ========== type =========
  typedef CutExpBasis<CD, 1> CutCSTO;
  template<> struct IsPrimitive<CutCSTO> {};

  template<class F>
  F CIP(const CutExpBasis<F, 1>& a, const ExpBasis<F, 1>& b);
  template<class F>
  F CIP(const CutExpBasis<F, 1>& a, const CutExpBasis<F, 1>& b);
  template<class F>
  F CIP(const ExpBasis<F, 1>& a, const CutExpBasis<F, 1>& b);

  template<class F>
  F CIP(const CutExpBasis<F, 1>& a, const DiracDelta<F>& b);
  template<class F>
  F CIP(const CutExpBasis<F, 2>& a, const DiracDelta<F>& b);
  template<class F>
  F CIP(const DiracDelta<F>& b, const CutExpBasis<F, 1>& a);
  template<class F>
  F CIP(const DiracDelta<F>& b, const CutExpBasis<F, 2>& a);

  // =========== External Function ===========
  template<class F>
  F CIP(const CutExpBasis<F, 1>&, const CutExpBasis<F, 1>&);
  template<class F>
  F CIP(const CutExpBasis<F, 1>&, const ExpBasis<F, 1>&);
  template<class F>
  F CIP(const ExpBasis<F, 1>&, const CutExpBasis<F, 1>&);


  template<class F, int m>
  ostream& operator << (ostream& os, const ExpBasis<F,m>& a);

  template<class Prim>
  Prim OperateCst(typename Prim::Field c, const Prim& f);

  template<int m>
  CutExpBasis<CD, m> ComplexConj(const CutExpBasis<CD, m>& f) {
    CutExpBasis<CD, m> res(conj(f.c()), f.n(), conj(f.z()), f.r0);
    return res;
  }


}
#endif
