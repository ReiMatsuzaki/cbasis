#include "cut_prim.hpp"

namespace l2func {

  // ============= Inner Product ==============  
  template<class F> F LowerGamma(int n, F z) {
    if(n <= 0)
      throw "n must be positive";

    if(n == 1) {
      return 1.0-exp(-z);
    } else {
      F surface = -pow(z, n-1)*exp(-z);
      F prev = LowerGamma(n-1, z);
    
      return surface + (n-1)*prev;
    }
  }
  template<class F> F CutSTO_Int(F z, int n, double r0) {
    return LowerGamma<F>(n+1, z*r0) / pow(z, n+1);
  }

  template<class F, int m>
  CutExpBasis<F,m>::CutExpBasis(): 
    basis(ExpBasis<F,m>()), r0_(10.0) {}
  template<class F, int m>
  CutExpBasis<F,m>::CutExpBasis(F _c, int _n, F _z, double _r0): 
    basis(_c, _n, _z), r0_(_r0) {}
  template<class F, int m>
  CutExpBasis<F,m>::CutExpBasis(int _n, F _z, double _r0, ENormalized): 
    basis(1, _n, _z), r0_(_r0) {
    
    F c2;
    if(m == 1) 
      c2 = CutSTO_Int(2*_z, 2*_n, _r0);
    else
      throw std::runtime_error("not implemented yet in CutExpBasis<F,m>::CutExpBasis");
    ExpBasis<F,m> u(F(1) / sqrt(c2), _n, _z);
    this->basis = u;
  }

  template<class F, int m>
  F CutExpBasis<F,m>::at(F x) const {
    if(abs(x) < this->r0_)
      return this->basis.at(x);
    else
      return F(0);
  }
  template<class F, int m>
  void CutExpBasis<F,m>::SetComplexConjugate(const CutExpBasis<F,m>& o) {
    this->basis.SetComplexConjugate(o.basis);
  }
  template<class F, int m>
  CutExpBasis<F,m> CutExpBasis<F,m>::ComplexConjugate() const {
    CutExpBasis<F,m> res(*this);    // copu this object
    res.SetComplexConjugate(*this); // set complex conjugate values
    return res;    
  }

/*
  template<class F, int m>
  std::ostream& operator << (std::ostream& os, 
			     const CutExpBasis<F,m>& a) {
    //    int m = Prim::exp_power;
    //    typename Prim::Field c = a.c();
    //typename Prim::Field z = a.z();
    F c = a.c();

    if(m == 1) 
      os << (char *)"STO(" ;
    else
      os << (char *)"GTO(";
    char* comma = (char*)", ";
    os << c;
    os << comma << a.n() << comma << a.z() << (char*)")";
    
    os << (char *)"Cut(at " << a.r0 << ")";
    os << a.basis;
    return os;
  }
  template std::ostream& operator <<<CD,1> (std::ostream& os, const CutCSTO a);
*/
  template<>
  CutCSTO OperateRm<CutCSTO>(int m, const CutCSTO& o) {
    double r0 = o.r0();
    CutCSTO res(o.c(), o.n()+m, o.z(), r0);
    return res;
  }
  template<>
  CutCSTO OperateCst(CD c, const CutCSTO& o) {
    double r0 = o.r0();
    CutCSTO res(o.c()*c, o.n(), o.z(), r0);
    return res;
  }

  template<class F>
  F CIP(const CutExpBasis<F,1>& a, const CutExpBasis<F,1>& b) {
    double r0 = a.r0() < b.r0() ? a.r0() : b.r0();
    return CutSTO_Int<F>(a.z()+b.z(), a.n()+b.n(), r0) * a.c() * b.c();
  }
  template<class F>
  F CIP(const CutExpBasis<F,1>& a, const ExpBasis<F,1>& b) {
    double r0 = a.r0();
    return CutSTO_Int<F>(a.z()+b.z(), a.n()+b.n(), r0) * a.c() * b.c();
  }
  template<class F>
  F CIP(const ExpBasis<F,1>& a, const CutExpBasis<F,1>& b) {
    return CIP(b, a);
  }

  template<class F>
  F CIP(const CutExpBasis<F,1>& a, const DiracDelta<F>& b) {
    return a.at(b.r0());
  }

  template<class F>
  F CIP(const DiracDelta<F>& a, const CutExpBasis<F,1>& b) {
    return CIP(b, a);
  }

  // ========= Explicit Decralation ===========
  template class CutExpBasis<CD, 1>;
  template double LowerGamma<double>(int n, double z);
  template CD LowerGamma<CD>(int n, CD z);
  template CD CutSTO_Int<CD>(CD z, int n, double r0);
  template double CutSTO_Int<double>(double z, int n, double r0);
  // template CD CutSTO_Int<CD>(CD z, int n, double r0);

  template CD CIP(const CutCSTO&, const CutCSTO&);
  template CD CIP(const CutCSTO&, const CSTO&);
  template CD CIP(const CSTO&, const CutCSTO&);
  template CD CIP(const CutCSTO& a, const DiracDelta<CD>& b);
  //  template CD CIP(const CutCGTO& a, const DiracDelta<CD>& b);
  template CD CIP(const DiracDelta<CD>& b, const CutCSTO& a);
  //  template CD CIP(const DiracDelta<CD>& b, const CutCGTO& a);

}
