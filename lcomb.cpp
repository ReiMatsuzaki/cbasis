#include "lcomb.hpp"

namespace l2func {

  template<class Prim>
  LinearComb<Prim>::LinearComb() { IsPrimitive<Prim>();}
  template<class Prim>
  LinearComb<Prim>::LinearComb(int n):
    cf_list_(n) {}
  template<class Prim>
  LinearComb<Prim>::LinearComb(const Prim& prim):
    cf_list_(1) {
    cf_list_[0] = make_pair(1.0, prim);
  }

  template<class Prim>
  LinearComb<Prim> LinearComb<Prim>::Clone() const {
    LinearComb<Prim> res(*this);
    return res;
  }

  template<class Prim>
  typename Prim::Field LinearComb<Prim>::at(Field x) const {

    Field acc(0);

    typedef typename LinearComb<Prim>::const_iterator IT;
    for(IT it = this->begin(), end_it = this->end(); it != end_it; ++it) 
      acc += it->first * it->second.at(x);

    return acc;    
  }

  template<class Prim>
  void LinearComb<Prim>::AddCoefPrim
  (Field c, const Prim& o) {
    cf_list_.push_back(make_pair(c, o));
  }
  template<class Prim>
  void LinearComb<Prim>::AddPair(pair<Field, Prim> cf) {
    cf_list_.push_back(cf);
  }
  template<class Prim>
  void LinearComb<Prim>::AddOther
  (const LinearComb<Prim>& o) {

    for(const_iterator it = o.begin(),
	  end = o.end(); it != end; ++it) {
      Field c = it->first;
      Prim  f = it->second;
      cf_list_.push_back(make_pair(c, f));
    }

  }

  template<class Prim>
  void LinearComb<Prim>::Add(Field c, const Prim& o) {
    this->AddCoefPrim(c, o);
  }
  template<class Prim>
  void LinearComb<Prim>::Add(pair<Field, Prim> cf) {
    this->AddPair(cf);
  }
  template<class Prim>
  void LinearComb<Prim>::Add(const LinearComb<Prim>& o) {
    this->AddOther(o);
  }
 
  template<class Prim>
  void LinearComb<Prim>::operator += (pair<Field, Prim> cf) {
    this->AddPair(cf);
  }
  template<class Prim>
  void LinearComb<Prim>::operator += (const Prim& f) {
    this->AddCoefPrim(1, f);
  }
  template<class Prim>
  void LinearComb<Prim>::operator += (const LinearComb<Prim>& o) {
    this->AddOther(o);
  }

  template<class Prim>
  void LinearComb<Prim>::ScalarProduct(Field c) {
    for(typename VFP::iterator it = cf_list_.begin(),
	  end = cf_list_.end(); it != end; ++it) {
      it->first *=c;
    }

  }

  // ------ explicit instance -------
  template class LinearComb<RSTO>;
  template class LinearComb<CSTO>;
  template class LinearComb<RGTO>;
  template class LinearComb<CGTO>;


  template<class Prim>
  LinearComb<Prim> OperateDDr(const Prim& f) {
    
    typedef typename Prim::Field F;
    F c = f.c();
    int n = f.n();
    F z = f.z();
    int m = Prim::exp_power;
    
    Prim a(c * F(n),     n-1, z);
    Prim b(c * (-z*F(m)), n+m-1, z);

    LinearComb<Prim> res;
    res += F(1) * a;
    res += F(1) * b;
    return(res);
  }
  template<class Prim>
  LinearComb<Prim> OperateDDr2 (const Prim& f) {
    
    // f  = r^n exp(-zr^m)
    // df = (nr^{n-1} -mzr^{n+m-1}) exp(-zr^m)
    // ddf= (n(n-1)r^{n-2} -mz(n+m-1)r^{n+m-2}
    //     +-zmnr^{n+m-2} +mmzzr^{n+2m-2} ) exp(-zr^m)
    
    typedef typename Prim::Field F;
    F   c = f.c();
    int n = f.n();
    int m = Prim::exp_power;
    F   z = f.z();
    
    LinearComb<Prim> res;
    res += c * Prim(n*(n-1),        n-2,     z);
    res += c * Prim(-m*z*(2*n+m-1), n+m-2,   z);
    res += c * Prim(m*m*z*z,        n+2*m-2, z);
    return res;
  }
  template LinearComb<RSTO> OperateDDr(const RSTO&);
  template LinearComb<CSTO> OperateDDr(const CSTO&);
  template LinearComb<RGTO> OperateDDr(const RGTO&);
  template LinearComb<CGTO> OperateDDr(const CGTO&);
  template LinearComb<RSTO> OperateDDr2(const RSTO&);
  template LinearComb<CSTO> OperateDDr2(const CSTO&);
  template LinearComb<RGTO> OperateDDr2(const RGTO&);
  template LinearComb<CGTO> OperateDDr2(const CGTO&);

  // --------- derivative basis -------------
  // assuming a is normalized basis set.
  // above formula is from Sotsuron.
  // STO  
  // a(r) = N r^n Exp[-zr]
  //    N = 1/sqrt((2z)^{-2n-1} (2n)!) = z^{n+1/2} / sqrt(2n!)
  //   dN = (n+1/2)z^{n-1/2} / sqrt(2n!)
  //  ddN = (n-1/2)(n+1/2)z^{n-3/2} / sqrt(2n!)
  // dN/N = (n+1/2)/z
  //ddN/N = (n-1/2)(n+1/2)/(z^2)
  // GTO
  // a(r) = N r^n Exp[-zr^2]
  // dN/N = (1/4 + n/2) / z
  //ddN/N = (-3/4+n/2)(1/4+n/2)/(z^2)
  template<class F>
  LinearComb<ExpBasis<F,1> > D1Normalized(const ExpBasis<F, 1>& a) {
    LinearComb<ExpBasis<F, 1> > res;
    F   c = a.c(); // normalization term
    int n = a.n(); // principle number
    F   z = a.z(); // orbital expoent

    F   cp = (n + 0.5) / z * c;
    
    res += F(1) * ExpBasis<F, 1>(cp, n,   z);
    res += F(1) * ExpBasis<F, 1>(-c, n+1, z);

    return res;
  }
  template<class F>
  LinearComb<ExpBasis<F,2> > D1Normalized(const ExpBasis<F, 2>& a) {
    LinearComb<ExpBasis<F, 2> > res;
    F   c = a.c(); // normalization term
    int n = a.n(); // principle number
    F   z = a.z(); // orbital expoent

    F   cp = (n * 0.5 + 0.25) / z * c;
    
    res += 1.0 * ExpBasis<F, 2>(cp, n,   z);
    res += 1.0 * ExpBasis<F, 2>(-c, n+2, z);

    return res;
  }
  template<class F>
  LinearComb<ExpBasis<F,1> > D2Normalized(const ExpBasis<F, 1>& a) {
    LinearComb<ExpBasis<F, 1> > res;
    F   c = a.c(); // normalization term
    int n = a.n(); // principle number
    F   z = a.z(); // orbital expoent

    F   cp = (n + 0.5) / z * c;
    F  cpp = F(4 * n * n - 1) / (F(4) * z * z) * c;
    
    res += F(1) * ExpBasis<F, 1>(cpp,        n,   z);
    res += F(1) * ExpBasis<F, 1>(-F(2) * cp, n+1, z);
    res += F(1) * ExpBasis<F, 1>(c,          n+2, z);

    return res;
  }
  template<class F>
  LinearComb<ExpBasis<F,2> > D2Normalized(const ExpBasis<F, 2>& a) {
    LinearComb<ExpBasis<F, 2> > res;
    F   c = a.c(); // normalization term
    int n = a.n(); // principle number
    F   z = a.z(); // orbital expoent

    F   cp = (n * 0.5 + 0.25) / z * c;
    F  cpp = (-3.0/4.0 + n / 2.0) * (0.25 + n * 0.5) / (z*z) * c;    
    
    res += 1.0 * ExpBasis<F, 2>(cpp,       n,   z);
    res += 1.0 * ExpBasis<F, 2>(-F(2)*cp , n+2, z);
    res += 1.0 * ExpBasis<F, 2>(c,         n+4, z);

    return res;
  }    
  template LinearComb<RSTO> D1Normalized(const RSTO&);
  template LinearComb<CSTO> D1Normalized(const CSTO&);
  template LinearComb<RGTO> D1Normalized(const RGTO&);
  template LinearComb<CGTO> D1Normalized(const CGTO&);
  template LinearComb<RSTO> D2Normalized(const RSTO&);
  template LinearComb<CSTO> D2Normalized(const CSTO&);
  template LinearComb<RGTO> D2Normalized(const RGTO&);
  template LinearComb<CGTO> D2Normalized(const CGTO&);

}
