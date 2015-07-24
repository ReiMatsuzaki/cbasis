#ifndef PRIM_HPP
#define PRIM_HPP

#include <complex>

namespace {
  using std::ostream;
  typedef std::complex<double> CD;
}

namespace l2func {

  // ==================== static ======================
  enum ENormalized {
    Normalized
  };
  template<typename T> struct raise_error;
  
  // ============= inner product ==========
  const CD operator*(const CD& a, int b);
  const CD operator*(int b, const CD& a);
  int int_pow(int base, unsigned int expo);
  template<class F> F STO_Int(F z, int n);
  template<class F> F GTO_Int(F z, int n);
       
  // ==================== STO or GTO =================
  // represent STO or GTO.
  // m==1 => STO, m==2 => GTO
  // value of m can be called exp_power;
  template<class F, int m>
  class ExpBasis {
    
  public:
    // ---------- typedef ---------------------------
    typedef F Field;
    enum EExpPower { exp_power=m };
    
    // ---------- Field Member ----------------------
  private:
    F c_;    // coefficient
    int n_;  // principle number
    F z_;    // orbital exponent
    
  public:
    // ----------- Constructors ---------------------
    ExpBasis();
    ExpBasis(F _c, int _n, F _z);
    ExpBasis(int _n, F _z, ENormalized);
    template<class U>
    ExpBasis(const ExpBasis<U, m>& o):
      c_(o.c()), n_(o.n()), z_(o.z()) {}      

    // ----------- Accessors ------------------------
    F c() const { return c_; }
    int n() const { return n_; }
    F z() const { return z_; }
    void set_z(F z) { z_ = z; }    
    F at(F x) const;
    // ------------ Operation -----------------------
    void SetComplexConjugate(const ExpBasis<F,m>& a);
    ExpBasis<F, m> ComplexConjugate() const;

  };

  // =========== typedef =========================
  typedef ExpBasis<double, 1> RSTO;
  typedef ExpBasis<double, 2> RGTO;
  typedef ExpBasis<CD, 1> CSTO;
  typedef ExpBasis<CD, 2> CGTO;
  
  template<class F>
  F CIP(const ExpBasis<F, 1>& a, const ExpBasis<F, 1>& b);
  template<class F>
  F CIP(const ExpBasis<F, 1>& a, const ExpBasis<F, 2>& b);
  template<class F>  
  F CIP(const ExpBasis<F, 2>& a, const ExpBasis<F, 1>& b);
  template<class F>
  F CIP(const ExpBasis<F, 2>& a, const ExpBasis<F, 2>& b);

  // =========== raise error if not primitive =========
  template<class Prim> struct IsPrimitive;
  template<> struct IsPrimitive<RSTO> {};
  template<> struct IsPrimitive<CSTO> {};
  template<> struct IsPrimitive<RGTO> {};
  template<> struct IsPrimitive<CGTO> {};

  // ============= operation =================
  template<class F, int m>
  ostream& operator << (ostream& os, const ExpBasis<F,m>& a);

  template<class Prim>
  Prim OperateRm( int m, const Prim& f);
  template<class Prim>
  Prim OperateCst(typename Prim::Field c, const Prim& f);

  template<int m>
  ExpBasis<CD, m> ComplexConj(const ExpBasis<CD, m>& f) {
    ExpBasis<CD, m> res(conj(f.c()), f.n(), conj(f.z()));
    return res;
  }

}

#endif
