#include <stdexcept>
#include <iostream>
#include "../math/erfc.hpp"
#include "../math/int_exp.hpp"
#include "../utils/fact.hpp"
#include "../utils/macros.hpp"
#include "r1basis.hpp"
#include "r1_lc.hpp"

using namespace std;
using namespace Eigen;
using namespace erfc_mori;

namespace cbasis {

  // ==== Basic calculation ====    
  template<> dcomplex EXPInt<1,1>(int n, dcomplex a, dcomplex b) {
    return STOInt_Rplus(n, a+b);
  }  
  template<> dcomplex EXPInt<2,2>(int n, dcomplex a, dcomplex b) {
    return GTOInt_Rplus(n, a+b);
  }  
  template<> dcomplex EXPInt<1,2>(int n, dcomplex a, dcomplex b) {
    return STO_GTOInt_Rplus(n, a, b);
  }  
  template<> dcomplex EXPInt<2,1>(int n, dcomplex a, dcomplex b) {
    return STO_GTOInt_Rplus(n, b, a);
  }
  
  
  template<> dcomplex EXPIntD2<1, 1>(int ni, dcomplex zi, int nj, dcomplex zj) {

    //
    // dr2 r^n exp(-zr) = { zz -2nzr^{-1} + (nn-n)r^{-2}} r^n exp[-zr]
    //

    dcomplex acc(0);
    acc += zj*zj      * EXPInt<1,1>(ni+nj,   zi, zj);
    acc += -2.0*nj*zj * EXPInt<1,1>(ni+nj-1, zi, zj);
    if(nj > 1)
      acc += dcomplex(nj*nj-nj) * EXPInt<1,1>(ni+nj-2, zi, zj);

    return acc;

  }  
  template<> dcomplex EXPIntD2<2, 1>(int ni, dcomplex zi, int nj, dcomplex zj) {

    //
    // dr2 r^n exp(-zr) = { zz -2nzr^{-1} + (nn-n)r^{-2}} r^n exp[-zr]
    //

    dcomplex acc(0);
    acc += zj*zj      * EXPInt<2,1>(ni+nj,   zi, zj);
    acc += -2.0*nj*zj * EXPInt<2,1>(ni+nj-1, zi, zj);
    if(nj > 1)
      acc += dcomplex(nj*nj-nj) * EXPInt<2,1>(ni+nj-2, zi, zj);

    return acc;

  }
  template<> dcomplex EXPIntD2<1, 2>(int ni, dcomplex zi, int nj, dcomplex zj) {
    
    //
    // dr2 r^n exp(-zr^2) = { 4zz r^{2} -2z(2n+1) + (nn-n)r^{-2}} r^n exp[-zr^2]
    //
    
    dcomplex acc(0);
    acc += ( +4.0*zj*zj * EXPInt<1,2>(ni+nj+2, zi, zj)
	     -2.0*(2*nj+1)*zj* EXPInt<1,2>(ni+nj, zi, zj));
    if(nj > 1)
      acc += dcomplex(nj*nj-nj) * EXPInt<1,2>(ni+nj-2, zi, zj);
    return acc;
    
  } 
  template<> dcomplex EXPIntD2<2, 2>(int ni, dcomplex zi, int nj, dcomplex zj) {
    
    //
    // dr2 r^n exp(-zr^2) = { 4zz r^{2} -2z(2n+1) + (nn-n)r^{-2}} r^n exp[-zr^2]
    //
    
    dcomplex acc(0);
    acc += ( +4.0*zj*zj * EXPInt<2,2>(ni+nj+2, zi, zj)
	     -2.0*(2*nj+1)*zj* EXPInt<2,2>(ni+nj, zi, zj));
    if(nj > 1)
      acc += dcomplex(nj*nj-nj) * EXPInt<2,2>(ni+nj-2, zi, zj);
    return acc;
    
  }

  template<int m1, int m2>
  dcomplex EXPIntLC(typename _EXPs<m1>::LC_EXPs a,
		    int m,
		    typename _EXPs<m2>::LC_EXPs b) {
    dcomplex acc(0);
    for(int i = 0; i < a->size(); i++)
      for(int j = 0; j < b->size(); j++) {
	dcomplex c(a->c(i) * b->c(j));
	int      n(a->n(i) + b->n(j));
	acc +=  c * EXPInt<m1, m2>(n+m, a->z(i), b->z(j));
      }
    return acc;    
  }
  template dcomplex EXPIntLC<1,1>(_EXPs<1>::LC_EXPs a, int m, _EXPs<1>::LC_EXPs b);
  template dcomplex EXPIntLC<1,2>(_EXPs<1>::LC_EXPs a, int m, _EXPs<2>::LC_EXPs b);
  template dcomplex EXPIntLC<2,1>(_EXPs<2>::LC_EXPs a, int m, _EXPs<1>::LC_EXPs b);
  template dcomplex EXPIntLC<2,2>(_EXPs<2>::LC_EXPs a, int m, _EXPs<2>::LC_EXPs b);
      
  // ==== member field ====
  template<int m>
  _EXPs<m>::_EXPs() {
  }
  template<int m>
  VectorXcd _EXPs<m>::AtR(const VectorXcd& rs, const VectorXcd& cs) const {

    if(cs.size() != this->size()) {
      string msg;
      SUB_LOCATION(msg);
      msg += " : size mismatch.";
      throw runtime_error(msg);
    }

    if(!this->HasCoefAll()) {
      string msg;
      SUB_LOCATION(msg);
      msg = "\n" + msg + " : coef is not set";
    }

    int num(this->size());
    VectorXcd ys = VectorXcd::Zero(rs.size());

    for(int i = 0; i < num; i++) {
      VectorXcd ys0 = this->basis(i)->AtR(rs);
      ys += cs[i] * ys0;
    }

    return ys;

  }
  template<int m>
  VectorXcd _EXPs<m>::DAtR(const VectorXcd& rs, const VectorXcd& cs) const {

    if(cs.size() != this->size()) {
      string msg;
      SUB_LOCATION(msg);
      msg += " : size mismatch.";
      throw runtime_error(msg);
    }

    if(!this->HasCoefAll()) {
      string msg;
      SUB_LOCATION(msg);
      msg = "\n" + msg + " : coef is not set";
    }    

    int num(this->size());
    VectorXcd ys = VectorXcd::Zero(rs.size());

    for(int i = 0; i < num; i++) {
      VectorXcd ys0 = this->basis(i)->DAtR(rs);
      ys += ys0;
    }

    return ys;

  }
  template<int m>
  dcomplex _EXPs<m>::AtR_One(dcomplex r, const VectorXcd& cs) const {

    VectorXcd rs(1);
    rs[0] = r;
    return this->AtR(rs, cs)[0];

  }
  template<int m>
  string _EXPs<m>::str() const {
    ostringstream oss; 
    if(m == 1) {
      oss << "==== STOs ====" << endl;
      oss << "STO basis set" << endl;
    } else if (m == 2) {
      oss << "==== GTOs ====" << endl;
      oss << "GTO basis set" << endl;
    }

    for(int i = 0; i < this->size(); i++) {
      oss << i << " : " << this->basis(i)->str() << endl;
    }
    
    oss << "==============" << endl;
    
    return oss.str();
  }
  template<int m>
  bool _EXPs<m>::IsPrim(int i) const {
    if(i < 0 || this->size() <= i ) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + " : index out of range";
      throw runtime_error(msg);
    }
    return (this->basis(i)->size() == 1);
  }
  template<int m>
  bool _EXPs<m>::IsPrimAll() const {
    bool acc = true;
    typedef vector<int>::const_iterator It;
    for(int i = 0; i < this->size(); i++) {
      acc = acc && (this->IsPrim(i));
    }
    return acc;
  }
  template<int m>
  bool _EXPs<m>::IsNormal(int i) const {
    if(i < 0 || this->size() <= i ) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + " : index out of range";
      throw runtime_error(msg);
    }
    return (coef_type_[i] == COEF_NORMAL);
  }
  template<int m>
  bool _EXPs<m>::IsNormalAll() const {
    bool acc = true;
    typedef vector<int>::const_iterator It;
    for(int i = 0; i < this->size(); i++) {
      acc = acc && (this->IsNormal(i));
    }
    return acc;
  }
  template<int m>
  bool _EXPs<m>::HasCoef(int i) const {
    if(i < 0 || this->size() <= i ) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + " : index out of range";
      throw runtime_error(msg);
    }
    return (coef_type_[i] == COEF_NORMAL ||
	    coef_type_[i] == COEF_NOT_NORMAL);
  }  
  template<int m> 
  bool _EXPs<m>::HasCoefAll() const {
    bool acc = true;
    typedef vector<int>::const_iterator It;
    for(int i = 0; i < this->size(); i++) {
      acc = acc && (this->HasCoef(i));
    }
    return acc;
  }

  template<int m>
  _EXPs<m>* _EXPs<m>::AddPrim(int n, dcomplex z) {

    LC_EXPs g = Create_LC_EXPs<m>();
    g->Add(1.0, n, z);
    this->basis_.push_back(g);
    this->coef_type_.push_back(COEF_NO);
    return this;
    
  }
  template<int m>
  _EXPs<m>* _EXPs<m>::AddPrims(int n, Eigen::VectorXcd zs) {

    for(int i = 0; i < zs.size(); i++) 
      this->AddPrim(n, zs[i]);    
    return this;
    
  }
  template<int m>
  _EXPs<m>* _EXPs<m>::AddLC(LC_EXPs lc) {
    
    this->basis_.push_back(lc);
    this->coef_type_.push_back(COEF_NO);
    return this;
    
  }
  template<int m>
  _EXPs<m>* _EXPs<m>::AddNotNormalPrim(dcomplex c, int n, dcomplex z) {

    LC_EXPs g = Create_LC_EXPs<m>();
    g->Add(c, n, z);
    this->basis_.push_back(g);
    this->coef_type_.push_back(COEF_NOT_NORMAL);
    return this;
    
  }
  template<int m>
  _EXPs<m>* _EXPs<m>::AddNotNormalLC(LC_EXPs lc) {

    this->basis_.push_back(lc);
    this->coef_type_.push_back(COEF_NOT_NORMAL);
    return this;
    
  }      
  template<int m>
  _EXPs<m>* _EXPs<m>::SetUp() {

    if(this->HasCoefAll())
      return this;

    int num(this->size());
    
    for(int i = 0; i < num; i++) {

      if(this->coef_type_[i] == COEF_NO) {
	LC_EXPs bi = this->basis(i);
	dcomplex nterm(1.0/sqrt(EXPIntLC<m, m>(bi, 0, bi)));
	for(int ii = 0; ii < bi->size(); ii++) {
	  bi->c(ii) *= nterm;
	}
	this->coef_type_[i] = COEF_NORMAL;
      }
    }

    return this;
    
  }

  template<int m>
  _EXPs<m>* _EXPs<m>::ReplaceLC(int i, LC_EXPs lc) {

    if(i < 0 || i >= this->size()) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + " : index out of range";
      throw runtime_error(msg);
    }

    this->basis_[i] = lc->Clone();
    this->coef_type_[i] = COEF_NO;
    return this;
    
  }
  
  // ---- Create ----
  template<int m>
  typename _EXPs<m>::EXPs _EXPs<m>::self() {
    EXPs ptr = boost::enable_shared_from_this<_EXPs<m> >::shared_from_this();
    return ptr;
  }
  template<int m> void _EXPs<m>::Clone_Symbolic(typename _EXPs<m>::EXPs *o) const {

    EXPs ptr = Create_EXPs<m>();
    for(int i = 0; i < this->size(); i++) {
      LC_EXPs lc;
      this->basis(i)->Clone_Symbolic(&lc);
      ptr->AddLC(lc);
    }
    *o = ptr;
    
  }
  template<int m> void _EXPs<m>::Clone(int reuse, typename _EXPs<m>::EXPs *o) const {
    
    if(reuse == INITIAL) {
      this->Clone_Symbolic(o);
    } else if(reuse != REUSE) {
      THROW_ERROR("reuse <- {INITIAL, REUSE}");
    }

    EXPs ptr = *o;
    
    int num(this->size());
    for(int i = 0; i < num; i++) {
      this->basis(i)->Clone_Numeric(ptr->basis(i));
      ptr->coef_type_[i] = this->coef_type_[i];
    }
    
  }  
  template<> void _EXPs<1>::DerivOneZeta(int reuse, EXPs *other) const {
    
    if(!this->IsPrimAll()) {
      THROW_ERROR("all basis must be primitive");
    }
    if(!this->IsNormalAll()) {
      THROW_ERROR("all basis is not normalized");
    }

    if(reuse == INITIAL) {
      *other = Create_EXPs<1>(this->size(), 2, COEF_NOT_NORMAL);
    } else if (reuse != REUSE) {
      THROW_ERROR("reuse <- {INITIAL, REUSE} ");
    } 

    for(int i = 0; i < this->size(); i++) {
      LC_EXPs a = this->basis(i);
      LC_EXPs b = (*other)->basis(i);
      dcomplex N(a->c(0));
      dcomplex z(a->z(0));
      int      n(a->n(0));
      dcomplex N1((n+0.5)/z);
      b->c(0) = N1*N;
      b->n(0) = n;
      b->z(0) = z;
      b->c(1) = -N;
      b->n(1) = n+1;
      b->z(1) = z;      
    }
  }   
  template<> void _EXPs<1>::DerivTwoZeta(int reuse, EXPs *other) const {
    
    if(!this->IsPrimAll()) {
      THROW_ERROR("all basis must be primitive");
    }
    if(!this->IsNormalAll()) {
      THROW_ERROR("all basis is not normalized");
    }

    if(reuse == INITIAL) {
      *other = Create_EXPs<1>(this->size(), 3, COEF_NOT_NORMAL);
    } else if (reuse != REUSE) {
      THROW_ERROR("reuse <- {INITIAL, REUSE} ");
    } 

    for(int i = 0; i < this->size(); i++) {
      LC_EXPs a = this->basis(i);
      LC_EXPs b = (*other)->basis(i);
      dcomplex N(a->c(0));
      dcomplex z(a->z(0));
      int      n(a->n(0));
      dcomplex N1((n+0.5)/z);
      dcomplex N2((n*n-1.0/4.0)/(z*z));
      b->c(0) = N2*N;
      b->n(0) = n;
      b->z(0) = z;
      b->c(1) = -2.0*N1*N;
      b->n(1) = n+1;
      b->z(1) = z;
      b->c(2) = N;
      b->n(2) = n+2;
      b->z(2) = z;      
    }
  }  
  template<> void _EXPs<2>::DerivOneZeta(int reuse, EXPs *other) const {
    if(!this->IsPrimAll()) {
      THROW_ERROR("all basis must be primitive");
    }
    if(!this->IsNormalAll()) {
      THROW_ERROR("all basis is not normalized");
    }

    if(reuse == INITIAL) {
      EXPs ptr = Create_EXPs<2>(this->size(), 2, COEF_NOT_NORMAL);
      *other = ptr;
    } else if (reuse != REUSE) {
      THROW_ERROR("reuse <- {INITIAL, REUSE} ");
    }
    
    for(int i = 0; i < this->size(); i++) {
      LC_EXPs a = this->basis(i);
      LC_EXPs b = (*other)->basis(i);
      dcomplex N(a->c(0));
      dcomplex z(a->z(0));
      int      n(a->n(0));
      dcomplex N1((2*n+1.0)/(4.0*z));
      b->c(0) = N1*N;
      b->n(0) = n;
      b->z(0) = z;
      b->c(1) = -N;
      b->n(1) = n+2;
      b->z(1) = z;
    }    
  }    
  template<> void _EXPs<2>::DerivTwoZeta(int reuse, EXPs *other) const {
    
    if(!this->IsPrimAll()) {
      THROW_ERROR("all basis must be primitive");
    }
    if(!this->IsNormalAll()) {
      THROW_ERROR("all basis is not normalized");
    }

    if(reuse == INITIAL) {
      *other = Create_EXPs<2>(this->size(), 3, COEF_NOT_NORMAL);
    } else if (reuse != REUSE) {
      THROW_ERROR("reuse <- {INITIAL, REUSE} ");
    } 

    for(int i = 0; i < this->size(); i++) {
      LC_EXPs a = this->basis(i);
      LC_EXPs b = (*other)->basis(i);
      dcomplex N(a->c(0));
      dcomplex z(a->z(0));
      int      n(a->n(0));
      dcomplex N1((n+0.5)/z);
      dcomplex N2((n*n-1.0/4.0)/(z*z));
      b->c(0) = N2*N;
      b->n(0) = n;
      b->z(0) = z;
      b->c(1) = -2.0*N1*N;
      b->n(1) = n+1;
      b->z(1) = z;
      b->c(2) = N;
      b->n(2) = n+2;
      b->z(2) = z;      
    }
    
  }

  template<int m> typename _EXPs<m>::EXPs _EXPs<m>::Clone() const {
    
    EXPs ptr = Create_EXPs<m>();

    int num(this->size());
    for(int i = 0; i < num; i++) {
      LC_EXPs lc = this->basis_[i]->Clone();
      ptr->basis_.push_back(lc);
      ptr->coef_type_.push_back(this->coef_type_[i]);
    }
    ptr->SetUp();
    return ptr;
  }  
  template<int m> typename _EXPs<m>::EXPs _EXPs<m>::Conj() const {

    EXPs ptr = Create_EXPs<m>();

    int num(this->size());
    for(int i = 0; i < num; i++) {
      LC_EXPs lc = this->basis_[i]->Conj();
      ptr->basis_.push_back(lc);
      ptr->coef_type_.push_back(this->coef_type_[i]);
    }
    ptr->SetUp();
    return ptr;
    
  }  
  
  // -- vec/matrix
  template<int m>
  void _EXPs<m>::CalcVec(LC_STOs stos, int scall, VectorXcd *vec) {
    cbasis::CalcVec<m,1>(this->self(), stos, scall, vec);
  }
  template<int m>
  void _EXPs<m>::CalcVec(LC_GTOs gtos, int scall, VectorXcd *vec) {
    cbasis::CalcVec<m,2>(this->self(), gtos, scall, vec);
  }
  template<int m>
  void _EXPs<m>::CalcRmMat(int M, int scall, MatrixXcd *mat) {
    cbasis::CalcRmMat<m, m>(this->self(), M, this->self(), scall, mat);
  }
  template<int m>
  void _EXPs<m>::CalcD2Mat(int scall, MatrixXcd *mat) {
    cbasis::CalcD2Mat<m, m>(this->self(), this->self(), scall, mat);
  }

  // -- for compatible --
  template<int m>
  MatrixXcd _EXPs<m>::CalcRmMat(int M) const {
    throw runtime_error("do not use");
  }
  template<int m>
  MatrixXcd _EXPs<m>::CalcD2Mat() const {
    throw runtime_error("do not use");
  }
  template<int m>
  VectorXcd _EXPs<m>::CalcVecSTO(LC_STOs o) const {
    throw runtime_error("do not use");
  }
  template<int m>
  VectorXcd _EXPs<m>::CalcVecGTO(LC_GTOs gtos) const {
    throw runtime_error("do not use");
  }
  
  
  // ---- realize ----
  template class _EXPs<1>;
  template class _EXPs<2>;

  // ==== external etc =====
  template<int m> typename _EXPs<m>::EXPs Create_EXPs() {
    
    typedef typename _EXPs<m>::EXPs EXPs;
    EXPs ptr(new _EXPs<m>);
    return ptr;

  }
  template<int m> typename _EXPs<m>::EXPs Create_EXPs(int nbasis, int nprim, int coef_type) {
    typename _EXPs<m>::EXPs ptr = Create_EXPs<m>();
    for(int i = 0; i < nbasis; i++) {
      typename _EXPs<m>::LC_EXPs lc = Create_LC_EXPs<m>(nprim);
      if(coef_type == COEF_NORMAL) {
	ptr->AddLC(lc);
      } else if(coef_type == COEF_NOT_NORMAL) {
	ptr->AddNotNormalLC(lc);
      } else {
	THROW_ERROR("coef_type <- {COEF_NOT_NORMAL, COEF_NORMAL}");
      }
    }
    return ptr;
  }
  STOs Create_STOs() {

    return Create_EXPs<1>();
    
  }
  GTOs Create_GTOs() {

    return Create_EXPs<2>();
    
  }
  
  // ==== Vector ====
  template<int m1>
  void CalcVec_Symbolic(typename _EXPs<m1>::EXPs a, VectorXcd *m) {
    int n(a->size());
    *m = VectorXcd::Zero(n);
  }
  template<int m1, int m2>
  void CalcVec_Numeric(typename _EXPs<m1>::EXPs a, typename _EXPs<m2>::LC_EXPs b,
		       VectorXcd& vec) {
    if(vec.size() != a->size()) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + " : size mismatch";
      throw runtime_error(msg);
    }

    if(!a->HasCoefAll()) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "EXP does not have coef. Call SetUp or use no_normal method only";
      throw runtime_error(msg);
    }

    int n = vec.size();
    for(int i = 0; i < n; i++) {
      typename _EXPs<m1>::LC_EXPs bi = a->basis(i);
      dcomplex ele = EXPIntLC<m1, m2>(bi, 0, b);      
      vec.coeffRef(i) = ele;
    }
  }
  template<int m1, int m2>
  void CalcVec(typename _EXPs<m1>::EXPs a, typename _EXPs<m2>::LC_EXPs b,
	       int reuse, VectorXcd *vec) {
    if(reuse == INITIAL) {
      CalcVec_Symbolic<m1>(a, vec);
    } else if(reuse != REUSE) {
      THROW_ERROR("scall <- {INITIAL, REUSE}");
    }
    cbasis::CalcVec_Numeric<m1, m2>(a, b, *vec);
  }
  template void CalcVec<1,1>(_EXPs<1>::EXPs, _EXPs<1>::LC_EXPs, int, VectorXcd*);
  template void CalcVec<1,2>(_EXPs<1>::EXPs, _EXPs<2>::LC_EXPs, int, VectorXcd*);
  template void CalcVec<2,1>(_EXPs<2>::EXPs, _EXPs<1>::LC_EXPs, int, VectorXcd*);
  template void CalcVec<2,2>(_EXPs<2>::EXPs, _EXPs<2>::LC_EXPs, int, VectorXcd*);

  
  // ==== Matrix====
  template<int m1, int m2>  
  void CalcMat_Symbolic(typename _EXPs<m1>::EXPs a,
			typename _EXPs<m2>::EXPs b, Eigen::MatrixXcd *mat) {
    *mat = MatrixXcd::Zero(a->size(), b->size());
  } 
  template<int m1, int m2>
  void CalcRmMat_Numeric(typename _EXPs<m1>::EXPs a, int M,
			 typename _EXPs<m2>::EXPs b,
			 MatrixXcd& mat) {
    int na(a->size());
    int nb(b->size());
    if(mat.rows() != na || mat.cols() != nb) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + " : size mismatch";
      throw runtime_error(msg);
    }

    if(!a->HasCoefAll() || !b->HasCoefAll()) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "EXP does not have coef.";
      throw runtime_error(msg);
    }
    
    for(int i = 0; i < na; i++) {
      typename _EXPs<m1>::LC_EXPs bi = a->basis(i);
      for(int j = 0; j < nb; j++) {
	typename _EXPs<m2>::LC_EXPs bj = b->basis(j);
	mat.coeffRef(i, j) = EXPIntLC<m1, m2>(bi, M, bj);
      }
    }
    
  }  
  template<int m1, int m2>
  void CalcD2Mat_Numeric(typename _EXPs<m1>::EXPs a, typename _EXPs<m2>::EXPs b,
			 MatrixXcd& mat) {
    int na(a->size());
    int nb(b->size());
    if(mat.rows() != na || mat.cols() != nb) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + " : size mismatch";
      throw runtime_error(msg);
    }

    if(!a->HasCoefAll() || !b->HasCoefAll()) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "EXP does not have coef.";
      throw runtime_error(msg);
    }

    for(int i = 0; i < na; i++) {
      typename _EXPs<m1>::LC_EXPs bi = a->basis(i);
      for(int j = 0; j < nb; j++) {
	typename _EXPs<m2>::LC_EXPs bj = b->basis(j);
	
	dcomplex acc(0);
	for(int ii = 0; ii < bi->size(); ii++)
	  for(int jj = 0; jj < bj->size(); jj++) {
	    dcomplex c(bi->c(ii) * bj->c(jj));
	    int      ni(bi->n(ii));
	    int      nj(bj->n(jj));
	    dcomplex zi(bi->z(ii));
	    dcomplex zj(bj->z(jj));
	    acc += c * EXPIntD2<m1, m2>(ni, zi, nj, zj);
	  }
	mat.coeffRef(i, j) = acc;
      }
    }
  }
  template<int m1, int m2>
  void CalcRmMat(typename _EXPs<m1>::EXPs a, int M,
		 typename _EXPs<m2>::EXPs b, int reuse,
		 Eigen::MatrixXcd *mat) {
    if(reuse == INITIAL) {
      CalcMat_Symbolic<m1, m2>(a, b, mat);
    } else if(reuse != REUSE) {
      THROW_ERROR("scall <- {INITIAL, REUSE}");
    }
    CalcRmMat_Numeric<m1, m2>(a, M, b, *mat);
  }
  template<int m1, int m2>
  void CalcD2Mat(typename _EXPs<m1>::EXPs a, typename _EXPs<m2>::EXPs b, 
		 int scall, Eigen::MatrixXcd *mat) {
    if(scall == INITIAL) {
      CalcMat_Symbolic<m1, m2>(a, b, mat);
    }
    if(scall != INITIAL && scall != REUSE) {
      THROW_ERROR("scall <- {INITIAL, REUSE}");
    }
    CalcD2Mat_Numeric<m1, m2>(a, b, *mat);
  }
  template void CalcRmMat<1,1>(_EXPs<1>::EXPs, int, _EXPs<1>::EXPs, int, MatrixXcd*);
  template void CalcRmMat<1,2>(_EXPs<1>::EXPs, int, _EXPs<2>::EXPs, int, MatrixXcd*);
  template void CalcRmMat<2,1>(_EXPs<2>::EXPs, int, _EXPs<1>::EXPs, int, MatrixXcd*);
  template void CalcRmMat<2,2>(_EXPs<2>::EXPs, int, _EXPs<2>::EXPs, int, MatrixXcd*);
  template void CalcD2Mat<1,1>(_EXPs<1>::EXPs a, _EXPs<1>::EXPs b, int, MatrixXcd*);
  template void CalcD2Mat<1,2>(_EXPs<1>::EXPs a, _EXPs<2>::EXPs b, int, MatrixXcd*);
  template void CalcD2Mat<2,1>(_EXPs<2>::EXPs a, _EXPs<1>::EXPs b, int, MatrixXcd*);
  template void CalcD2Mat<2,2>(_EXPs<2>::EXPs a, _EXPs<2>::EXPs b, int, MatrixXcd*);

  
  // ==== for compatible ====
  template<int m1, int m2>
  VectorXcd CalcVec(typename _EXPs<m1>::EXPs a,
		    typename _EXPs<m2>::LC_EXPs b) {
    
    if(!a->HasCoefAll()) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "EXP does not have coef. Call SetUp or use no_normal method only";
      throw runtime_error(msg);
    }
    
    int numa(a->size());
    VectorXcd vec(numa);
    
    for(int i = 0; i < numa; i++) {
      typename _EXPs<m1>::LC_EXPs bi = a->basis(i);
      vec(i) = EXPIntLC<m1, m2>(bi, 0, b);
    }
    return vec;
  }
  template VectorXcd CalcVec<1,1>(_EXPs<1>::EXPs a, _EXPs<1>::LC_EXPs b);
  template VectorXcd CalcVec<1,2>(_EXPs<1>::EXPs a, _EXPs<2>::LC_EXPs b);
  template VectorXcd CalcVec<2,1>(_EXPs<2>::EXPs a, _EXPs<1>::LC_EXPs b);
  template VectorXcd CalcVec<2,2>(_EXPs<2>::EXPs a, _EXPs<2>::LC_EXPs b);

  template<int m1, int m2>
  MatrixXcd CalcRmMat(typename _EXPs<m1>::EXPs a,
		      int M,
		      typename _EXPs<m2>::EXPs b) {
    
    if(!a->HasCoefAll() || !b->HasCoefAll()) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "EXP does not have coef. Call SetUp or use no_normal method only";
      throw runtime_error(msg);
    }
    
    int numa(a->size());
    int numb(b->size());
    MatrixXcd mat(numa, numb);
    
    for(int i = 0; i < numa; i++)
      for(int j = 0; j < numb; j++) {

	typename _EXPs<m1>::LC_EXPs bi = a->basis(i);
	typename _EXPs<m2>::LC_EXPs bj = b->basis(j);

	dcomplex acc(0);
	for(int ii = 0; ii < bi->size(); ii++) {
	  for(int jj = 0; jj < bj->size(); jj++) {
	    dcomplex c(bi->c(ii) * bj->c(jj));
	    int      n(bi->n(ii) + bj->n(jj));
	    dcomplex zi(bi->z(ii));
	    dcomplex zj(bj->z(jj));
	    acc +=  c * EXPInt<m1, m2>(n+M, zi, zj);
	  }
	}
	mat(i, j) = acc;
      }

    return mat;
  }
  template MatrixXcd CalcRmMat<1,1>(_EXPs<1>::EXPs a, int M, _EXPs<1>::EXPs b);
  template MatrixXcd CalcRmMat<1,2>(_EXPs<1>::EXPs a, int M, _EXPs<2>::EXPs b);
  template MatrixXcd CalcRmMat<2,1>(_EXPs<2>::EXPs a, int M, _EXPs<1>::EXPs b);
  template MatrixXcd CalcRmMat<2,2>(_EXPs<2>::EXPs a, int M, _EXPs<2>::EXPs b);
  
  template<int m1, int m2>
  Eigen::MatrixXcd CalcD2Mat(typename _EXPs<m1>::EXPs a,
			     typename _EXPs<m2>::EXPs b) {

    if(!a->HasCoefAll() || !b->HasCoefAll()) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "No coefficient.";
      throw runtime_error(msg);
    }

    int numa(a->size());
    int numb(b->size());
    MatrixXcd mat(numa, numb);
    
    for(int i = 0; i < numa; i++) {
      for(int j = 0; j < numb; j++) {
	
	typename _EXPs<m1>::LC_EXPs bi = a->basis(i);
	typename _EXPs<m2>::LC_EXPs bj = b->basis(j);
	
	dcomplex acc(0);
	for(int ii = 0; ii < bi->size(); ii++)
	  for(int jj = 0; jj < bj->size(); jj++) {
	    dcomplex c(bi->c(ii) * bj->c(jj));
	    int      ni(bi->n(ii));
	    int      nj(bj->n(jj));
	    dcomplex zi(bi->z(ii));
	    dcomplex zj(bj->z(jj));
	    acc += c * EXPIntD2<m1, m2>(ni, zi, nj, zj);
	  }
	mat(i, j) = acc;
      }
    }
    
    return mat;    

  }
  template MatrixXcd CalcD2Mat<1,1>(_EXPs<1>::EXPs a, _EXPs<1>::EXPs b); 
  template MatrixXcd CalcD2Mat<1,2>(_EXPs<1>::EXPs a, _EXPs<2>::EXPs b); 
  template MatrixXcd CalcD2Mat<2,1>(_EXPs<2>::EXPs a, _EXPs<1>::EXPs b); 
  template MatrixXcd CalcD2Mat<2,2>(_EXPs<2>::EXPs a, _EXPs<2>::EXPs b);
  
}
