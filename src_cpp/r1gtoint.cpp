#include <iostream>
#include <exception>
#include <algorithm>
#include "r1gtoint.hpp"
#include "macros.hpp"
#include "exp_int.hpp"
#include "mult_array.hpp"

using namespace Eigen;
using namespace std;

namespace l2func {  

  typedef MultArray<dcomplex, 2> A2;

  // ==== External ====
  dcomplex NPrimeGTO(dcomplex nterm, int n, dcomplex z) {
    /*
      nterm : normalization term
      n : principle number
      z : orbital exponents
     */
    return (0.5+n)/(2.0*z) * nterm;
    
  }
  dcomplex NDoublePrimeGTO(dcomplex nterm, int n, dcomplex z) {
    /*
      nterm : normalization term
      n : principle number
      z : orbital exponents
     */
    return dcomplex((-3+2*n)*(1+2*n)) / (16.0*z*z) * nterm;
    
  }
  dcomplex TMat(int ni, int nj, dcomplex zj, int L, dcomplex* gs) {
    dcomplex cumsum(0);
    cumsum += zj * dcomplex(2*nj+1) * gs[ni+nj];
    cumsum += -0.5 * dcomplex(nj*nj-nj-L*L-L)*gs[ni+nj-2];
    cumsum += -2.0*zj*zj*gs[ni+nj+2];
    return cumsum;
  }
  
  // ==== Data ====
  const MatrixXcd& MatVecMap::mat(string lbl) const {

    map<string, MatrixXcd>::const_iterator it = mat_.find(lbl);
    if(it == this->mat_.end()) {
      string msg; SUB_LOCATION(msg);
      msg += ": failed to find matrix. lbl=" + lbl;
      throw runtime_error(msg);
    }
    return it->second;
    
  }
  MatrixXcd& MatVecMap::mat(string lbl) {

    //    return 
    //      const_cast<MatrixXcd&>(static_cast<const MatVecMap>(*this).mat(lbl));

    map<string, MatrixXcd>::const_iterator it = mat_.find(lbl);
    if(it == this->mat_.end()) {
      string msg; SUB_LOCATION(msg);
      msg += ": failed to find matrix. lbl=" + lbl;
      throw runtime_error(msg);
    }
    return mat_[lbl];

  }
  void MatVecMap::InitMatIfNecessary(string lbl, int ni, int nj) {
    if(!this->exist_mat(lbl)) {
      this->mat_[lbl] = MatrixXcd::Zero(ni, nj);
    } else {
      MatrixXcd& m = this->mat_[lbl];
      if(m.rows() != ni || m.cols() != nj) {
	m = MatrixXcd::Zero(ni, nj);
      }
    }
  }
  const VectorXcd& MatVecMap::vec(string lbl) const {

    map<string, VectorXcd>::const_iterator it = vec_.find(lbl);
    if(it == this->vec_.end()) {
      string msg; SUB_LOCATION(msg);
      msg += ": failed to find matrix. lbl=" + lbl;
      throw runtime_error(msg);
    }
    return it->second;
    
  }
  VectorXcd& MatVecMap::vec(string lbl) {
    
    //return 
    //      const_cast<VectorXcd&>(static_cast<const MatVecMap>(*this).vec(lbl));
    map<string, VectorXcd>::iterator it = vec_.find(lbl);
    if(it == this->vec_.end()) {
      string msg; SUB_LOCATION(msg);
      msg += ": failed to find matrix. lbl=" + lbl;
      throw runtime_error(msg);
    }
    return it->second;
  }
  void MatVecMap::InitVecIfNecessary(string lbl, int n) {
    if(!this->exist_vec(lbl)) {
      this->vec_[lbl] = VectorXcd::Zero(n);
    }else {
      VectorXcd& v = this->vec_[lbl];
      if(v.size() != n) {
	v = VectorXcd::Zero(n);
      }
    }
  }

  // ==== R1GTO/R1STO ====
  R1GTO::R1GTO(dcomplex _c, int _n, dcomplex _z): c(_c), n(_n), z(_z) {
    if(n < 1) {
      string msg; SUB_LOCATION(msg);
      msg += "n must be positive.";
      throw runtime_error(msg);
    }
  }
  R1STO::R1STO(dcomplex _c, int _n, dcomplex _z): c(_c), n(_n), z(_z) {
    /*
    if(n < 0) {
      string msg; SUB_LOCATION(msg);
      msg += "n must be positive or 0.";
      throw runtime_error(msg);
    }
    */
  }

  void CalcGTOInt(int maxn, dcomplex a, dcomplex* res) {

    if(maxn < 0) {
      std::string msg; SUB_LOCATION(msg);
      msg += "maxn must be bigger than 0";
      throw std::runtime_error(msg);
    }

    res[0] = sqrt(M_PI)/(2.0*sqrt(a));
    
    if(maxn == 0)
      return;

    res[1] = 0.5 / a;

    if(maxn == 1)
      return;

    for(int n = 2; n <= maxn; n++)
      res[n] = dcomplex(n-1)/(2.0*a)*res[n-2];

  }
  ostream& operator<<(ostream& out, const R1GTO& basis) {
    out << "R1GTO("
	<< basis.c << ", "
	<< basis.n << ", "
	<< basis.z << ")";
    return out;
  }
  ostream& operator<<(ostream& out, const R1STO& basis) {
    out << "R1STO("
	<< basis.c << ", "
	<< basis.n << ", "
	<< basis.z << ")";
    return out;
  }

  //R1GTOs::R1GTOs(int _L):
  //    coef_set_q_(false),
  //    coef_type_("Nothing"),
  //    L_(_L) { }
  R1GTOs::R1GTOs(): setup_q_(false), coef_type_("Nothing") {}

  // ---- Utils ----
  int  R1GTOs::max_n() const {

    int maxn(0);
    for(cItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      for(cItPrim it_p = it->prim.begin(), end_p = it->prim.end();
	  it_p != end_p; ++it_p) {
	if(maxn < it_p->n)
	  maxn = it_p->n;
      }
    }

    return maxn;    
  }
  bool R1GTOs::IsSameStructure(const R1GTOs& o) const {
    if(this->conts_.size() != o.conts_.size())
      return false;
    
    for(cItCont it = this->conts_.begin(), it_end = this->conts_.end(),
	  jt = o.conts_.begin(); it != it_end; ++it, ++jt) {
      if(it->size_prim() != jt->size_prim())
	return false;
      if(it->size_basis() != jt->size_basis())
	return false;
    }

    return true;
      
  }

  // ---- Getter ----
  int R1GTOs::calc_size_basis() const {
    int cumsum(0);
    for(vector<Contraction>::const_iterator it = conts_.begin(), end = conts_.end();
	it != end; ++it) {
      cumsum += it->coef.rows();
    }
    return cumsum;
  }
  int R1GTOs::calc_size_prim() const {
    int cumsum(0);
    for(vector<Contraction>::const_iterator it = conts_.begin(), end = conts_.end();
	it != end; ++it) {
      cumsum += it->coef.cols();
    }
    return cumsum;
  }
  const R1GTOs::Prim& R1GTOs::prim(int i) const {

    typedef vector<Contraction>::const_iterator ContIt;
    int cumsum(0);

    for(ContIt it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      if(cumsum <= i && i < cumsum + it->size_prim() )
	return it->prim[i-cumsum];
      cumsum += it->size_prim();
    }

    string msg; SUB_LOCATION(msg);
    ostringstream oss; 
    oss << msg << endl << ": Failed to find index i." << endl;
    oss << "i: " << i << endl;
    throw runtime_error(oss.str());
    
  }
  R1GTOs::Prim& R1GTOs::prim(int i) {
    typedef vector<Contraction>::iterator ContIt;
    int cumsum(0);

    for(ContIt it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      if(cumsum <= i && i < cumsum + it->size_prim() )
	return it->prim[i-cumsum];
      cumsum += it->size_prim();
    }

    string msg; SUB_LOCATION(msg);
    ostringstream oss; 
    oss << msg << endl << ": Failed to find index i." << endl;
    oss << "i: " << i << endl;
    throw runtime_error(oss.str());
  }

  // ---- Setter ----
  void R1GTOs::Add(int _n, dcomplex _zeta) {
    this->setup_q_ = false;
    this->coef_type_ = "Nothing";
    Contraction cont;
    cont.prim.push_back(Prim(_n, _zeta));
    cont.coef = MatrixXcd::Ones(1, 1);    
    conts_.push_back(cont);    
  }
  void R1GTOs::Add(int n, const VectorXcd& zs) {
    this->setup_q_ = false;
    this->coef_type_ = "Nothing";

    for(int i = 0; i < zs.size(); i++) {
      this->Add(n, zs[i]);
    }
    
  }
  void R1GTOs::Add(int n, const VectorXcd& zs, const MatrixXcd& coef) {

    if(zs.size() != coef.cols()) {
      string msg; SUB_LOCATION(msg);
      throw runtime_error(msg);
    }

    this->setup_q_ = false;
    this->coef_type_ = "Nothing";

    Contraction cont;
    for(int i = 0; i < zs.size(); i++) 
      cont.prim.push_back(Prim(n, zs(i)));
    cont.coef = coef;
    cont.offset = this->calc_size_basis();
    this->conts_.push_back(cont);
    
  }
  void R1GTOs::Set(const Eigen::VectorXcd& zs) {
    
    if(zs.size() != this->size_prim()) {
      string msg; SUB_LOCATION(msg);
      msg += ": size mismatch.";
      ostringstream oss;
      oss << msg << endl;
      oss << "size_prim: " << this->size_prim() << endl;
      oss << "zs.size:   " << zs.size() << endl;
      throw runtime_error(oss.str());
    }

    int i(0);
    for(vector<Contraction>::iterator it = conts_.begin(), end = conts_.end();
	it != end; ++it) {
      for(int ip = 0; ip < it->size_prim(); ip++, i++) {
	it->prim[ip].z = zs[i];
      }
    }

    this->setup_q_ = false;
    this->coef_type_ = "Nothing";

  }
  void R1GTOs::Set(int n, const Eigen::VectorXcd& zs) {
    
    if(zs.size() != this->size_prim()) {
      string msg; SUB_LOCATION(msg);
      msg += ": size mismatch.";
      throw runtime_error(msg);
    }

    int i(0);
    for(vector<Contraction>::iterator it = conts_.begin(), end = conts_.end();
	it != end; ++it) {
      for(int ip = 0; ip < it->size_prim(); ip++, i++) {
	it->prim[ip].n = n;
	it->prim[ip].z = zs[i];
      }
    }

    this->setup_q_ = false;
    this->coef_type_ = "Nothing";

  }

  // ---- basis convert ----  
  void R1GTOs::SetConj(const R1GTOs& o) {
    if(!o.setup_q()) {
      string msg; SUB_LOCATION(msg);
      msg += ": Call SetUp first";
      throw runtime_error(msg);
    }

    if(!this->IsSameStructure(o)) {

      *this = o; // copy

    }

    ItCont it = this->conts_.begin();
    ItCont end = this->conts_.end();
    cItCont jt =  o.conts_.begin();
    for(; it != end; ++it, ++jt) {
      ItPrim iprim = it->prim.begin();
      ItPrim prim0 = it->prim.end();
      cItPrim  jprim = jt->prim.begin();
      for(; iprim != prim0; ++iprim, ++jprim) {
	iprim->z = conj(iprim->z);
      }
      it->coef = it->coef.conjugate();
    }
    
    
    this->coef_type_  = o.coef_type_;

  }
  void R1GTOs::SetOneDeriv(const R1GTOs& o) {

    if(o.coef_type() != "Normalized" || !o.setup_q()) {
      string msg; SUB_LOCATION(msg);
      msg += ": o must be normalized.\n";
      msg += "coef_type: " + o.coef_type() + "\n";
      msg += "setup_q: " + string(o.setup_q()? "Yes" : "No" )+ "\n";
      throw runtime_error(msg);
    }

    for(cItCont it = o.conts_.begin(), end = o.conts_.end(); it != end; ++it) {
      if(it->size_prim() != 1 || it->size_basis() != 1) {
	string msg; SUB_LOCATION(msg);
	msg += ": Only support non contracted basis set.";
	throw runtime_error(msg);
      }
    }

    // -- we assume if setup_q is true, then this object has same structure as o.
    if(!this->setup_q()) {
      vector<Contraction> conts;
      this->conts_.swap(conts);
      for(cItCont it = o.conts_.begin(), end = o.conts_.end(); it != end; ++it) {
	VectorXcd zs = VectorXcd::Zero(2);
	MatrixXcd cs = MatrixXcd::Ones(1, 2);
	this->Add(it->prim[0].n, zs, cs);
      }
      this->SetUp();
    }

    ItCont it = this->conts_.begin();
    for(cItCont jt = o.conts_.begin(), end = o.conts_.end(); jt != end; ++jt, ++it) {
      dcomplex c(jt->coef(0, 0));
      int      n(jt->prim[0].n);
      dcomplex z(jt->prim[0].z);
      it->coef(0, 0) = NPrimeGTO(c, n, z);
      it->prim[0].n =  n;
      it->prim[0].z =  z;

      it->coef(0, 1) = -c;
      it->prim[1].n  = n+2;
      it->prim[1].z  = z;
	
    }
    
    this->coef_type_  = "D1Normalized";
  }
  void R1GTOs::SetTwoDeriv(const R1GTOs& o) {
    
    if(o.coef_type() != "Normalized" || !o.setup_q()) {
      string msg; SUB_LOCATION(msg);
      msg += ": o must be normalized.";
      throw runtime_error(msg);
    }

    for(cItCont it = o.conts_.begin(), end = o.conts_.end(); it != end; ++it) {
      if(it->size_prim() != 1 || it->size_basis() != 1) {
	string msg; SUB_LOCATION(msg);
	msg += ": Only support non contracted basis set.";
	throw runtime_error(msg);
      }
    }
    
    // -- we assume if setup_q is true, then this object has same structure as o.
    if(!this->setup_q()) {
      vector<Contraction> conts;
      this->conts_.swap(conts);
      for(cItCont it = o.conts_.begin(), end = o.conts_.end(); it != end; ++it) {
	VectorXcd zs = VectorXcd::Zero(3);
	MatrixXcd cs = MatrixXcd::Ones(1, 3);
	this->Add(it->prim[0].n, zs, cs);
      }
      this->SetUp();
    }

    ItCont it = this->conts_.begin();
    for(cItCont jt = o.conts_.begin(), end = o.conts_.end(); jt != end; ++jt, ++it) {
      dcomplex c(jt->coef(0, 0));
      int      n(jt->prim[0].n);
      dcomplex z(jt->prim[0].z);

      it->coef(0, 0) = NDoublePrimeGTO(c, n, z);
      it->prim[0].n =  n;
      it->prim[0].z =  z;

      it->coef(0, 1) = -2.0*NPrimeGTO(c, n, z);
      it->prim[1].n  = n+2;
      it->prim[1].z  = z;

      it->coef(0, 2) = c;
      it->prim[2].n  = n+4;
      it->prim[2].z  = z;
      
    }    
    
    this->coef_type_  = "D2Normalized";
  }

  // ---- Calculation ----
  void R1GTOs::CreateMat(const R1GTOs& o, Eigen::MatrixXcd& m) const {
    int ni(this->size_basis());
    int nj(o.size_basis());
    if(m.rows() != ni || m.cols() != nj) {
      m = MatrixXcd::Zero(ni, nj);
    }
  }
  void R1GTOs::CreateVec(Eigen::VectorXcd& m) const {
    int n(this->size_basis());
    if(m.size() != n) 
      m = VectorXcd::Zero(n);
  }
  
  void R1GTOs::CalcMatSTV(const R1GTOs& o, int L, MatVecMap& mat_vec,
			  string s_lbl, string t_lbl, string v_lbl) const {

    int num = this->size_basis();
    mat_vec.InitMatIfNecessary(s_lbl, num, num);
    mat_vec.InitMatIfNecessary(t_lbl, num, num);
    mat_vec.InitMatIfNecessary(v_lbl, num, num);

    MatrixXcd& s = mat_vec.mat(s_lbl);
    MatrixXcd& t = mat_vec.mat(t_lbl);
    MatrixXcd& v = mat_vec.mat(v_lbl);

    this->CalcMatSTV(o, L, s, t, v);

  }
  void R1GTOs::CalcMatSTV(const R1GTOs& o, int L,
			  MatrixXcd& S, MatrixXcd& T, MatrixXcd& V) const {
    if(!this->setup_q()) {
      string msg; SUB_LOCATION(msg);
      msg += ": coef is not set.";
      throw runtime_error(msg);
    }

    static MultArray<dcomplex, 2> s_prim_(1000);
    static MultArray<dcomplex, 2> t_prim_(1000);
    static MultArray<dcomplex, 2> v_prim_(1000);

    int maxn = 2*this->max_n() + 2 + 1;
    static vector<dcomplex> buf(maxn);
    if((int)buf.capacity() < maxn)
      buf.reserve(maxn);
    dcomplex* gs = &buf[0];

    this->CreateMat(o, S); this->CreateMat(o, T); this->CreateMat(o, V);

    for(cItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      for(cItCont jt = o.conts_.begin(); jt != o.conts_.end(); ++jt) {

	s_prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	t_prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	v_prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	for(int i = 0; i < it->size_prim(); i++) {
	  for(int j = 0; j < jt->size_prim(); j++) {
	    
	    int      ni(it->prim[i].n);
	    int      nj(jt->prim[j].n);
	    dcomplex zi(it->prim[i].z);
	    dcomplex zj(jt->prim[j].z);
	    CalcGTOInt(ni+nj+2, zi+zj, gs);
	    s_prim_(i, j) = gs[ni+nj];
	    v_prim_(i, j) = -gs[ni+nj-1];
	    t_prim_(i, j) = TMat(ni, nj, zj, L, gs);;
	    
	  }
	}

	int idx(it->offset);
	for(int ib = 0; ib < it->size_basis(); ib++, idx++) {
	  int jdx(jt->offset);
	  for(int jb = 0; jb < jt->size_basis(); jb++, jdx++) {
	    dcomplex s0(0), t0(0), v0(0);
	    for(int i = 0; i < it->size_prim(); i++) {
	      for(int j = 0; j < jt->size_prim(); j++) {
		dcomplex cc(it->coef(ib, i) * jt->coef(jb, j));
		s0 += cc * s_prim_(i, j);
		t0 += cc * t_prim_(i, j);
		v0 += cc * v_prim_(i, j);
	      }
	    }
	    S(idx, jdx) = s0;
	    T(idx, jdx) = t0;
	    V(idx, jdx) = v0;
	  }
	}
      }
    }        

  }
  void R1GTOs::CalcMatSTV(int L, MatVecMap& mat_vec,
			  string s_lbl, string t_lbl, string v_lbl) const {
    this->CalcMatSTV(*this, L, mat_vec, s_lbl, t_lbl, v_lbl);
  }
  void R1GTOs::CalcMatSTV(int L, MatrixXcd& S, 
			  MatrixXcd& T, MatrixXcd& V) const {
    this->CalcMatSTV(*this, L, S, T, V);
  }

  void R1GTOs::CalcMatSTO(const R1GTOs& o, const R1STOs& sto,
			  MatVecMap& res, string label) const {
    int num = size_basis();    
    res.InitMatIfNecessary(label, num, num);
    MatrixXcd& s = res.mat(label);
    this->CalcMatSTO(o, sto, s);
  }
  void R1GTOs::CalcMatSTO(const R1STOs& v, MatVecMap& res, string lbl) const {
    this->CalcMatSTO(*this, v, res, lbl);
  }
  void R1GTOs::CalcMatSTO(const R1GTOs& o, const R1STOs& v, MatrixXcd& s) const {
    if(!this->setup_q()) {
      string msg; SUB_LOCATION(msg);
      msg += ": coef is not set.";
      throw runtime_error(msg);
    }

    static MultArray<dcomplex, 2> prim_(1000);

    this->CreateMat(o, s);

    for(cItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {      
      for(cItCont jt = o.conts_.begin(); jt != o.conts_.end(); ++jt) {

	prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	for(int i = 0; i < it->size_prim(); i++) {
	  for(int j = 0; j < jt->size_prim(); j++) {
	    int      n(it->prim[i].n + jt->prim[j].n);
	    dcomplex z(it->prim[i].z + jt->prim[j].z);
	    dcomplex cumsum(0);
	    for(int io = 0; io < v.size_basis(); io++ ) {
	      dcomplex sto_c(v.basis(io).c);
	      int      sto_n(v.basis(io).n);
	      dcomplex sto_z(v.basis(io).z);
	      cumsum += sto_c * STO_GTO_Int(sto_z, z, sto_n+n);
	    }
	    prim_(i, j) = cumsum;
	  }
	}

	int idx(it->offset);
	for(int ib = 0; ib < it->size_basis(); ib++, idx++) {
	  int jdx(jt->offset);
	  for(int jb = 0; jb < jt->size_basis(); jb++, jdx++) {
	    dcomplex s0(0);
	    for(int i = 0; i < it->size_prim(); i++) {
	      for(int j = 0; j < jt->size_prim(); j++) {
		dcomplex cc(it->coef(ib, i) * jt->coef(jb, j));
		s0 += cc * prim_(i, j);
	      }
	    }
	    s(idx, jdx) = s0;
	  }
	}
      }
    }

  }
  void R1GTOs::CalcMatSTO(const R1STOs& v, MatrixXcd& s) const {
    this->CalcMatSTO(*this, v, s);
  }

  void R1GTOs::CalcVec(const R1STOs& o, MatVecMap& mat_vec, string label) const {
    int num = this->size_basis();
    mat_vec.InitVecIfNecessary(label, num);
    VectorXcd& m = mat_vec.vec(label);
    this->CalcVec(o, m);
  }
  void R1GTOs::CalcVec(const R1STOs& o, Eigen::VectorXcd& m) const {

    static MultArray<dcomplex, 1> m_prim_(20);

    if(!this->setup_q()) {
      string msg; SUB_LOCATION(msg); 
      msg += ": R1GTOs is not set up.";
      throw runtime_error(msg);
    }          
    
    this->CreateVec(m);

    for(cItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {

      int num_prim(it->size_prim());
      m_prim_.SetRange(0, num_prim);
      for(int i = 0; i < num_prim; i++ ) {
	dcomplex cumsum(0);
	for(int io = 0; io < o.size_basis(); io++ ) {
	  dcomplex c(o.basis(io).c);
	  int      n(it->prim[i].n + o.basis(io).n);
	  dcomplex z(it->prim[i].z);
	  cumsum += c * STO_GTO_Int(o.basis(io).z, z, n);
	}
	m_prim_(i) = cumsum;

	int idx(it->offset);
	for(int ib = 0; ib < it->size_basis(); ib++, idx++) {
	  dcomplex mval(0);
	  for(int i = 0; i < it->size_prim(); i++) {
	    dcomplex cc(it->coef(ib, i));
	    mval += cc * m_prim_(i);
	  }
	  m(idx) = mval;
	}
      }
    }
    
  }

  dcomplex HmES(int ni, int nj, dcomplex zj, int L, double energy, dcomplex *gs) {
    return (TMat(ni,   nj, zj, L, gs)
	    -gs[ni+nj-1]
	    -energy * gs[ni+nj]);
  }
  /*
  void R1GTOs::CalcDerivMat(double energy) {

    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      if(it->size_prim() != 1 || it->size_basis() != 1) {
	string msg; SUB_LOCATION(msg);
	msg += ": Only support non contracted basis set.";
	throw runtime_error(msg);
      }
    }

    int maxn = 2*this->max_n() + 6 + 1;
    static vector<dcomplex> buf(maxn);
    if((int)buf.capacity() < maxn)
      buf.reserve(maxn);
    dcomplex* gs = &buf[0];

    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg); 
      msg += ": coef is not set up.";
      throw runtime_error(msg);
    }

    int num = this->size_basis();    
    if(mat_.find("d10_hmes") == mat_.end())
      mat_["d10_hmes"] = MatrixXcd::Zero(num, num);
    if(mat_.find("d20_hmes") == mat_.end())
      mat_["d20_hmes"] = MatrixXcd::Zero(num, num);
    if(mat_.find("d11_hmes") == mat_.end())
      mat_["d11_hmes"] = MatrixXcd::Zero(num, num);

    MatrixXcd& d10 = mat_["d10_hmes"];
    MatrixXcd& d20 = mat_["d20_hmes"];
    MatrixXcd& d11 = mat_["d11_hmes"];

    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it)
      for(ItCont jt = conts_.begin(); jt != end; ++jt) {

	dcomplex ci(it->coef(0, 0));
	int      ni(it->prim[0].n);
	dcomplex zi(it->prim[0].z);
	dcomplex ci1(NPrimeGTO(ci, ni, zi));
	dcomplex ci2(NDoublePrimeGTO(ci, ni, zi));
	dcomplex cj(jt->coef(0, 0));
	int      nj(jt->prim[0].n);
	dcomplex zj(jt->prim[0].z);
	dcomplex cj1(NPrimeGTO(cj, nj, zj));	

	CalcGTOInt(ni+nj+6, zi+zj, gs);
	int idx(it->offset); int jdx(jt->offset);

       
//	     u  = (N  r^n                         ) Exp[-zr^2]
//	     u' = (N' r^n - N r^(n+2)             ) Exp[-zr^2]
//	     u''= (N''r^n -2N'r^(n+2) + N r^(n+4) ) Exp[-zr^2]

	d10(idx, jdx) = (+HmES(ni,   nj, zj, L_, energy, gs) * ci1 * cj
			 -HmES(ni+2, nj, zj, L_, energy, gs) * ci  * cj);
	d20(idx, jdx) = (+    HmES(ni,   nj, zj, L_, energy, gs) * ci2 * cj
			 -2.0*HmES(ni+2, nj, zj, L_, energy, gs) * ci1 * cj
			 +    HmES(ni+4, nj, zj, L_, energy, gs) * ci  * cj);
			 
	d11(idx, jdx) = (+HmES(ni,   nj,   zj, L_, energy, gs) * ci1 * cj1
			 -HmES(ni,   nj+2, zj, L_, energy, gs) * ci1 * cj
			 -HmES(ni+2, nj,   zj, L_, energy, gs) * ci  * cj1
			 +HmES(ni+2, nj+2, zj, L_, energy, gs) * ci  * cj);
      }
  }
*/
  /*
  void R1GTOs::CalcVec(const R1GTOs& o, string label) {
    string msg; SUB_LOCATION(msg);
    msg += ": Not implemented yet.";
    msg += "label: " + label;
    throw runtime_error(msg);
    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg); 
      msg += ": coef is not set up.";
      throw runtime_error(msg);
    }      
    if(o.coef_set_q_) {
      string msg; SUB_LOCATION(msg); 
      msg += ": coef is not set up.";
      throw runtime_error(msg);
    }      
  }
*/
/*
  void R1GTOs::CalcVec(const R1STOs& o, string label) {

    static MultArray<dcomplex, 1> m_prim_(20);

    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg); 
      msg += ": coef is not set up.";
      throw runtime_error(msg);
    }          
    
    int num(this->size_basis());
    if(vec_.find(label) == vec_.end())
      vec_[label] = VectorXcd::Zero(num);
    VectorXcd& m = vec_[label];
    
    typedef vector<Contraction>::iterator ContIt;
    for(ContIt it = conts_.begin(), end = conts_.end(); it != end; ++it) {

      int num_prim(it->size_prim());
      m_prim_.SetRange(0, num_prim);
      for(int i = 0; i < num_prim; i++ ) {
	dcomplex cumsum(0);
	for(int io = 0; io < o.size_basis(); io++ ) {
	  dcomplex c(o.basis(io).c);
	  int      n(it->prim[i].n + o.basis(io).n);
	  dcomplex z(it->prim[i].z);
	  cumsum += c * STO_GTO_Int(o.basis(io).z, z, n);
	}
	m_prim_(i) = cumsum;

	int idx(it->offset);
	for(int ib = 0; ib < it->size_basis(); ib++, idx++) {
	  dcomplex mval(0);
	  for(int i = 0; i < it->size_prim(); i++) {
	    dcomplex cc(it->coef(ib, i));
	    mval += cc * m_prim_(i);
	  }
	  m(idx) = mval;
	}
      }
    }
  }
*/
  /*
  void R1GTOs::CalcDerivVec(const R1STOs& o) { 
    
    
    //u  = (N  r^n                         ) Exp[-zr^2]
    //u' = (N' r^n - N r^(n+2)             ) Exp[-zr^2]
    //u''= (N''r^n -2N'r^(n+2) + N r^(n+4) ) Exp[-zr^2]
    
    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg); 
      msg += ": coef is not set up.";
      throw runtime_error(msg);
    }      

    //    typedef vector<Contraction>::iterator ContIt;
    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      if(it->size_prim() != 1 || it->size_basis() != 1) {
	string msg; SUB_LOCATION(msg);
	msg += ": CalcDerivVec only support non contracted basis set.";
	throw runtime_error(msg);
      }
    }

    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg); 
      msg += ": coef is not set up.";
      throw runtime_error(msg);
    }      

    int num = this->size_basis();    
    if(vec_.find("d1_m") == vec_.end())
      vec_["d1_m"] = VectorXcd::Zero(num);
    if(vec_.find("d2_m") == vec_.end())
      vec_["d2_m"] = VectorXcd::Zero(num);
    VectorXcd& d1_m = vec_["d1_m"];
    VectorXcd& d2_m = vec_["d2_m"];
    
    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      dcomplex cumsum_d1(0);
      dcomplex cumsum_d2(0);
      for(int i = 0; i < o.size_basis(); i++) {
	
	int      ng(it->prim[0].n);
	dcomplex zg(it->prim[0].z);
	dcomplex c0(it->coef(0, 0));
	dcomplex c1(NPrimeGTO(      c0, ng, zg));
	dcomplex c2(NDoublePrimeGTO(c0, ng, zg));
	dcomplex cs(o.basis(i).c);
	int      ns(o.basis(i).n);
	dcomplex zs(o.basis(i).z);
	cumsum_d1 += (+STO_GTO_Int(zs, zg, ns+ng  ) * cs * c1
		      -STO_GTO_Int(zs, zg, ns+ng+2) * cs * c0);
	cumsum_d2 += (+STO_GTO_Int(zs, zg, ns+ng  ) * cs * c2
		      -STO_GTO_Int(zs, zg, ns+ng+2) * cs * c1 * 2.0
		      +STO_GTO_Int(zs, zg, ns+ng+4) * cs * c0);
      }
      int idx(it->offset);
      d1_m(idx) = cumsum_d1;
      d2_m(idx) = cumsum_d2;
      
    }
    
  }
  */
  void R1GTOs::SetUp() {

    int cumsum(0);
    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      it->offset = cumsum;
      cumsum += it->size_basis();
    }
    num_basis_ = this->calc_size_basis();
    num_prim_ = this->calc_size_prim();
    this->setup_q_ = true;
    this->coef_type_ = "NotNormalized";
  }
  void R1GTOs::Normalize() {

    static MultArray<dcomplex, 1> m_prim_(20);

    this->SetUp();

    int maxn = this->max_n() + this->max_n() + 2 + 1;
    static vector<dcomplex> buf(maxn);
    if((int)buf.capacity() < maxn)
      buf.reserve(maxn);
    dcomplex* gs = &buf[0];

    typedef vector<Contraction>::iterator ItCont;
    typedef vector<R1GTO>::iterator ItPrim;

    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {

      for(int ibasis = 0; ibasis < it->size_basis(); ++ibasis) {
	dcomplex norm2(0);
	for(int i = 0; i < it->size_prim(); i++) {
	  for(int j = 0; j < it->size_prim(); j++) {
	    dcomplex c(it->coef(ibasis, i) * it->coef(ibasis, j));
	    int      n(it->prim[i].n + it->prim[j].n);
	    dcomplex z(it->prim[i].z + it->prim[j].z);
	    CalcGTOInt(n, z, gs);
	    norm2 += c*gs[n];
	  }
	}
	dcomplex scale(1.0/sqrt(norm2));
	for(int iprim = 0; iprim < it->size_prim(); iprim++) {
	  it->coef(ibasis, iprim) *= scale;
	}
      }
    }

    this->coef_type_ = "Normalized";
  }
  void R1GTOs::AtR(const VectorXcd& cs, const VectorXcd& rs,
		   VectorXcd* ys) {

    if(!this->setup_q()) {
      string msg; SUB_LOCATION(msg); 
      msg += ": this object is not set up.";
      throw runtime_error(msg);
    }      


    if(cs.size() != this->size_basis()) {
      string msg; SUB_LOCATION(msg);
      msg += ": \n size of cs must be equal to basis size";
      throw runtime_error(msg);
    }

    typedef vector<Contraction>::iterator ContIt;

    VectorXcd res = VectorXcd::Zero(rs.size());
    for(int i = 0; i < rs.size(); i++) {
      dcomplex r(rs[i]);
      dcomplex cumsum(0);

      for(ContIt it = conts_.begin(), end = conts_.end();
	  it != end; ++it) {
	for(int ib = 0; ib < it->size_basis(); ib++) {
	  for(int ip = 0; ip < it->size_prim(); ip++) {
	    dcomplex c(it->coef(ib, ip)); 
	    int      n(it->prim[ip].n); 
	    dcomplex z(it->prim[ip].z); 
	    cumsum += cs[it->offset+ib] * c * pow(r, n) * exp(-z*r*r);
	  }
	}
      }
      res(i) = cumsum;
    }
    ys->swap(res);

  }
  void R1GTOs::AtR(const VectorXcd& rs, VectorXcd* ys) {
    VectorXcd cs = VectorXcd::Ones(this->size_basis());
    this->AtR(cs, rs, ys);
  }  
  VectorXcd* R1GTOs::AtR(const VectorXcd& cs, const VectorXcd& rs) {
    VectorXcd* ys = new VectorXcd();
    this->AtR(cs, rs, ys);
    return ys;
  }
  VectorXcd* R1GTOs::AtR(const VectorXcd& rs) {
    VectorXcd cs = VectorXcd::Ones(this->size_basis());
    return this->AtR(cs, rs);
  }
  void R1GTOs::DerivAtR(const VectorXcd& cs, const VectorXcd& rs,
			VectorXcd* ys) const {
    if(!this->setup_q()) {
      string msg; SUB_LOCATION(msg);
      msg += ": this object is not set up.";
      throw runtime_error(msg);
    }
    if(cs.size() != this->size_basis()) {
      string msg; SUB_LOCATION(msg);
      msg += ": \n size of cs must be equal to basis size";
      throw runtime_error(msg);
    }
    
    VectorXcd res = VectorXcd::Zero(rs.size());
    for(int i = 0; i <rs.size(); i++) {
      dcomplex r(rs[i]);
      dcomplex cumsum(0);
      for(cItCont it = conts_.begin(), end = conts_.end();
	  it != end; ++it) {
	for(int ib = 0; ib < it->size_basis(); ib++) {
	  for(int ip = 0; ip < it->size_prim(); ip++) {
	    dcomplex c(it->coef(ib, ip));
	    int      n(it->prim[ip].n);
	    dcomplex z(it->prim[ip].z);
	    cumsum += cs[it->offset+ib] * c *
	      (dcomplex(n)*pow(r, n-1) -2.0*z*pow(r,n+1)) * exp(-z*r*r);
	  }
	}
      }
      res(i) = cumsum;
    }
    ys->swap(res);
  }
  void R1GTOs::DerivAtR(const VectorXcd& rs, VectorXcd* ys) const {
    VectorXcd cs = VectorXcd::Ones(this->size_basis());
    this->DerivAtR(cs, rs, ys);
  } 
  VectorXcd* R1GTOs::DerivAtR(const VectorXcd& cs, const VectorXcd& rs) const {
    VectorXcd* ys = new VectorXcd();
    this->DerivAtR(cs, rs, ys);
    return ys;
  }
  VectorXcd* R1GTOs::DerivAtR(const VectorXcd& rs) const {
    VectorXcd cs = VectorXcd::Ones(this->size_basis());
    return this->DerivAtR(cs, rs);
  }
  void R1GTOs::Deriv2AtR(const VectorXcd& cs, const VectorXcd& rs, VectorXcd* ys) const {
    if(!this->setup_q()) {
      string msg; SUB_LOCATION(msg);
      msg += ": this object is not set up";
      throw runtime_error(msg);
    }
    if(cs.size() != this->size_basis()) {
      string msg; SUB_LOCATION(msg);
      msg += ": size mismatch";
      throw runtime_error(msg);
    }

    VectorXcd res = VectorXcd::Zero(rs.size());
    for(int i = 0; i < rs.size(); i++) {
      dcomplex r(rs[i]);
      dcomplex cumsum(0);
      for(cItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
	for(int ib = 0; ib < it->size_basis(); ib++) {
	  for(int ip = 0; ip < it->size_prim(); ip++) {
	    dcomplex c(it->coef(ib, ip));
	    int      n(it->prim[ip].n);
	    dcomplex z(it->prim[ip].z);
	    dcomplex tmp(0);
	    tmp += 4.0 * z * z               * pow(r, n+2);
	    tmp += -2.0 * z * dcomplex(2*n+1) * pow(r, n  );
	    if(n > 1)
	      tmp += dcomplex(n*n-n) * pow(r, n-2);
	    tmp *= cs[it->offset+ib] * c * exp(-z*r*r);
	    cumsum += tmp;
	  }
	}
      }
      res(i) = cumsum;
    }
    ys->swap(res);
  }
  void R1GTOs::Deriv2AtR(const VectorXcd& rs, VectorXcd* ys) const {
    VectorXcd cs = VectorXcd::Ones(this->size_basis());
    this->Deriv2AtR(cs, rs, ys);
  }
  VectorXcd* R1GTOs::Deriv2AtR(const VectorXcd& cs, const VectorXcd& rs) const {
    VectorXcd* ys = new VectorXcd();
    this->Deriv2AtR(cs, rs, ys);
    return ys;
  }
  VectorXcd* R1GTOs::Deriv2AtR(const VectorXcd& rs) const {
    VectorXcd cs = VectorXcd::Ones(this->size_basis());
    return this->Deriv2AtR(cs, rs);
  }
  void R1GTOs::swap(R1GTOs& o) {
    
    ::swap(this->setup_q_, o.setup_q_);
    ::swap(this->coef_type_, o.coef_type_);
    this->conts_.swap(o.conts_);
  }

  void R1STOs::Add(dcomplex c, int n, dcomplex zeta) {
    stos_.push_back(R1STO(c, n, zeta));
  }
  void R1STOs::Add(int n, dcomplex zeta) {
    stos_.push_back(R1STO(1, n, zeta));
  }
  void R1STOs::AtR(const Eigen::VectorXcd& rs, Eigen::VectorXcd* ys) const {

    int num = rs.size();

    typedef vector<R1STO>::const_iterator cIt;
    VectorXcd res(num);
    for(int i = 0; i < num; i++) {
      dcomplex r(rs[i]);
      dcomplex cumsum(0);
      for(cIt it = stos_.begin(), end = stos_.end(); it != end; ++it) {
	cumsum += it->c * pow(r, it->n) * exp(-it->z * r);
      }
      res[i] = cumsum;
    }
    ys->swap(res);
  }
  Eigen::VectorXcd* R1STOs::AtR(const Eigen::VectorXcd& rs)       const {
    VectorXcd* ys = new VectorXcd();
    this->AtR(rs, ys);
    return ys;
  }

  ostream& operator<<(ostream& out, const R1GTOs& us) {
    out << "R1GTOs" << endl;
    out << "setup? : "
	<< (us.setup_q() ? "Yes" : "No")
	<< endl;
    out << "coef type : " << us.coef_type_ << endl;
    out << "num(contractions): " << us.conts_.size() << endl;
    
    for(R1GTOs::cItCont it = us.conts_.begin();
	it != us.conts_.end(); ++it) {
      //      if(it->size_prim() == 1 && it->size_basis() == 1) {
      //	cout << it->prim[0].n << ", " << it->prim[0].z << endl;
      //      } else {
      out << "offset:" << it->offset << endl;
      out << "prim:" << endl;
      for(R1GTOs::cItPrim iprim = it->prim.begin();
	  iprim != it->prim.end(); ++iprim) {
	
	out << iprim->n << ", " << iprim->z << endl;
	  
      }
      out << "coef:" << endl;
      out << it->coef << endl;
    }
    //    }

    return out;
  }
  ostream& operator<<(ostream& out, const R1STOs& us) {
    out << "R1STOs" << endl;
    out << "normalized_q : " << (us.normalized_q() ? "Yes" : "No");
    for(int i = 0; i < us.size_basis(); ++i) {
      out << out << us.basis(i) << endl;
    }
    return out;
  }

}
