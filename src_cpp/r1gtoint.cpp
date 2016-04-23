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

  R1GTOs::R1GTOs(int _L):
    coef_set_q_(false),
    coef_type_("Nothing"),
    L_(_L) { }

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
  int R1GTOs::size_basis() const {
    int cumsum(0);
    for(vector<Contraction>::const_iterator it = conts_.begin(), end = conts_.end();
	it != end; ++it) {
      cumsum += it->coef.rows();
    }
    return cumsum;
  }
  int R1GTOs::size_prim() const {
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
    this->coef_set_q_ = false;
    this->coef_type_ = "Nothing";
    Contraction cont;
    cont.prim.push_back(Prim(_n, _zeta));
    cont.coef = MatrixXcd::Ones(1, 1);
    cont.offset = this->size_basis();
    conts_.push_back(cont);    
  }
  void R1GTOs::Add(int n, const VectorXcd& zs) {
    this->coef_set_q_ = false;
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

    this->coef_set_q_ = false;
    this->coef_type_ = "Nothing";

    Contraction cont;
    for(int i = 0; i < zs.size(); i++) 
      cont.prim.push_back(Prim(n, zs(i)));
    cont.coef = coef;
    cont.offset = this->size_basis();
    this->conts_.push_back(cont);
    
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

    this->coef_set_q_ = false;
    this->coef_type_ = "Nothing";

  }

  // ---- basis convert ----  
  void R1GTOs::SetConj(const R1GTOs& o) {
    if(!o.coef_set_q_) {
      string msg; SUB_LOCATION(msg);
      msg += ": coef is not set.";
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
    
    this->coef_set_q_ = o.coef_set_q_;
    this->coef_type_  = o.coef_type_;
    this->L_          = o.L_;

  }
  
  void R1GTOs::Normalize() {

    static MultArray<dcomplex, 1> m_prim_(20);

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

    this->coef_set_q_ = true;
    this->coef_type_ = "Normalized";
  }
  void R1GTOs::CalcMatSTV(const R1GTOs& o, int L, MatVecMap& mat_vec,
			  string s_lbl, string t_lbl, string v_lbl) const {
    if(!this->coef_set_q_) {
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

    int num = size_basis();    
    
    if(mat_vec.mat.find(s_lbl) == mat_vec.mat.end())
      mat_vec.mat[s_lbl] = MatrixXcd::Zero(num, num);
    if(mat_vec.mat.find(t_lbl) == mat_vec.mat.end())
      mat_vec.mat[t_lbl] = MatrixXcd::Zero(num, num);
    if(mat_vec.mat.find(v_lbl) == mat_vec.mat.end())
      mat_vec.mat[v_lbl] = MatrixXcd::Zero(num, num);

    MatrixXcd& s = mat_vec.mat[s_lbl];
    MatrixXcd& t = mat_vec.mat[t_lbl];
    MatrixXcd& v = mat_vec.mat[v_lbl];

    for(cItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      for(cItCont jt = o.conts_.begin(); jt != end; ++jt) {

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
	    s(idx, jdx) = s0;
	    t(idx, jdx) = t0;
	    v(idx, jdx) = v0;
	  }
	}
      }
    }    
  }
  void R1GTOs::CalcMatSTV(int L, MatVecMap& mat_vec,
			  string s_lbl, string t_lbl, string v_lbl) const {
    this->CalcMatSTV(*this, L, mat_vec, s_lbl, t_lbl, v_lbl);
  }
  void R1GTOs::CalcMatSTO(const R1GTOs& o, const R1STOs& sto,
			  MatVecMap& res, string label) const {

    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg);
      msg += ": coef is not set.";
      throw runtime_error(msg);
    }

    static MultArray<dcomplex, 2> prim_(1000);

    int num = size_basis();    

    if(res.mat.find(label) == res.mat.end())
      res.mat[label] = MatrixXcd::Zero(num, num);

    MatrixXcd& s = res.mat[label];

    for(cItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {      
      for(cItCont jt = o.conts_.begin(); jt != end; ++jt) {

	prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	for(int i = 0; i < it->size_prim(); i++) {
	  for(int j = 0; j < jt->size_prim(); j++) {
	    int      n(it->prim[i].n + jt->prim[j].n);
	    dcomplex z(it->prim[i].z + jt->prim[j].z);
	    dcomplex cumsum(0);
	    for(int io = 0; io < sto.size_basis(); io++ ) {
	      dcomplex sto_c(sto.basis(io).c);
	      int      sto_n(sto.basis(io).n);
	      dcomplex sto_z(sto.basis(io).z);
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
  void R1GTOs::CalcMatSTO(const R1STOs& v, MatVecMap& res, string lbl) const {
    this->CalcMatSTO(*this, v, res, lbl);
  }
  void R1GTOs::CalcVec(const R1STOs& o, MatVecMap& mat_vec, string label) const {
    static MultArray<dcomplex, 1> m_prim_(20);

    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg); 
      msg += ": coef is not set up.";
      throw runtime_error(msg);
    }          
    
    int num(this->size_basis());
    if(mat_vec.vec.find(label) == mat_vec.vec.end())
      mat_vec.vec[label] = VectorXcd::Zero(num);
    VectorXcd& m = mat_vec.vec[label];

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

  void R1GTOs::CalcMat() {
    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg); 
      msg += ": coef is not set up.";
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

    int num = size_basis();    

    if(mat_.find("s") == mat_.end())
      mat_["s"] = MatrixXcd::Zero(num, num);
    if(mat_.find("t") == mat_.end())
      mat_["t"] = MatrixXcd::Zero(num, num);
    if(mat_.find("v") == mat_.end())
      mat_["v"] = MatrixXcd::Zero(num, num);

    MatrixXcd& s = mat_["s"];
    MatrixXcd& t = mat_["t"];
    MatrixXcd& v = mat_["v"];

    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      for(ItCont jt = conts_.begin(); jt != end; ++jt) {

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
	    t_prim_(i, j) = TMat(ni, nj, zj, L_, gs);;
	    
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
	    s(idx, jdx) = s0;
	    t(idx, jdx) = t0;
	    v(idx, jdx) = v0;
	  }
	}
      }
    }
  }
  void R1GTOs::CalcMatH() {
    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg); 
      msg += ": coef is not set up.";
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

    int num = size_basis();    

    if(mat_.find("hs") == mat_.end())
      mat_["hs"] = MatrixXcd::Zero(num, num);
    if(mat_.find("ht") == mat_.end())
      mat_["ht"] = MatrixXcd::Zero(num, num);
    if(mat_.find("hv") == mat_.end())
      mat_["hv"] = MatrixXcd::Zero(num, num);

    MatrixXcd& s = mat_["hs"];
    MatrixXcd& t = mat_["ht"];
    MatrixXcd& v = mat_["hv"];

    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      for(ItCont jt = conts_.begin(); jt != end; ++jt) {

	s_prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	t_prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	v_prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	for(int i = 0; i < it->size_prim(); i++) {
	  for(int j = 0; j < jt->size_prim(); j++) {
	    
	    int      ni(it->prim[i].n);
	    int      nj(jt->prim[j].n);
	    dcomplex zi(conj(it->prim[i].z));
	    dcomplex zj(jt->prim[j].z);
	    CalcGTOInt(ni+nj+2, zi+zj, gs);
	    dcomplex sval = gs[ni+nj];
	    dcomplex tval(0.0);
	    tval -= zj * dcomplex(4*nj+2) * gs[ni+nj];
	    if(ni+nj >= 2)
	      tval += dcomplex(nj*nj-nj-L_*L_-L_) * gs[ni+nj-2];
	    tval += 4.0*zj*zj*gs[ni+nj+2];
	    tval *= -0.5;
	    dcomplex vval = -gs[ni+nj-1];
	    s_prim_(i, j) = sval;
	    v_prim_(i, j) = vval;
	    t_prim_(i, j) = tval;
	    
	  }
	}

	int idx(it->offset);
	for(int ib = 0; ib < it->size_basis(); ib++, idx++) {
	  int jdx(jt->offset);
	  for(int jb = 0; jb < jt->size_basis(); jb++, jdx++) {
	    dcomplex s0(0), t0(0), v0(0);
	    for(int i = 0; i < it->size_prim(); i++) {
	      for(int j = 0; j < jt->size_prim(); j++) {
		dcomplex cc(conj(it->coef(ib, i)) * jt->coef(jb, j));
		s0 += cc * s_prim_(i, j);
		t0 += cc * t_prim_(i, j);
		v0 += cc * v_prim_(i, j);
	      }
	    }
	    s(idx, jdx) = s0;
	    t(idx, jdx) = t0;
	    v(idx, jdx) = v0;
	  }
	}
      }
    }
  }
  void R1GTOs::CalcMat(const R1STOs& sto, string label) {
    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg);
      msg += ": coef is not set.";
      throw runtime_error(msg);
    }

    static MultArray<dcomplex, 2> prim_(1000);

    int num = size_basis();    

    if(mat_.find(label) == mat_.end())
      mat_[label] = MatrixXcd::Zero(num, num);

    MatrixXcd& s = mat_[label];

    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {      
      for(ItCont jt = conts_.begin(); jt != end; ++jt) {

	prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	for(int i = 0; i < it->size_prim(); i++) {
	  for(int j = 0; j < jt->size_prim(); j++) {
	    int      n(it->prim[i].n + jt->prim[j].n);
	    dcomplex z(it->prim[i].z + jt->prim[j].z);
	    dcomplex cumsum(0);
	    for(int io = 0; io < sto.size_basis(); io++ ) {
	      dcomplex sto_c(sto.basis(io).c);
	      int      sto_n(sto.basis(io).n);
	      dcomplex sto_z(sto.basis(io).z);
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
  void R1GTOs::CalcMatH(const R1STOs& sto, string label) {
    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg); 
      msg += ": coef is not set up.";
      throw runtime_error(msg);
    }

    static MultArray<dcomplex, 2> prim_(1000);

    int num = size_basis();    

    if(mat_.find(label) == mat_.end())
      mat_[label] = MatrixXcd::Zero(num, num);

    MatrixXcd& s = mat_[label];

    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {      
      for(ItCont jt = conts_.begin(); jt != end; ++jt) {

	prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	for(int i = 0; i < it->size_prim(); i++) {
	  for(int j = 0; j < jt->size_prim(); j++) {
	    int      n(it->prim[i].n + jt->prim[j].n);
	    dcomplex z(conj(it->prim[i].z) +jt->prim[j].z);
	    dcomplex cumsum(0);
	    for(int io = 0; io < sto.size_basis(); io++ ) {
	      dcomplex sto_c(sto.basis(io).c);
	      int      sto_n(sto.basis(io).n);
	      dcomplex sto_z(sto.basis(io).z);
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
		dcomplex cc(conj(it->coef(ib, i)) * jt->coef(jb, j));
		s0 += cc * prim_(i, j);
	      }
	    }
	    s(idx, jdx) = s0;
	  }
	}
      }
    }
  }
  dcomplex HmES(int ni, int nj, dcomplex zj, int L, double energy, dcomplex *gs) {
    return (TMat(ni,   nj, zj, L, gs)
	    -gs[ni+nj-1]
	    -energy * gs[ni+nj]);
  }
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

	/*
	     u  = (N  r^n                         ) Exp[-zr^2]
	     u' = (N' r^n - N r^(n+2)             ) Exp[-zr^2]
	     u''= (N''r^n -2N'r^(n+2) + N r^(n+4) ) Exp[-zr^2]
	 */
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
    /*
    if(!this->normalized_q())
      this->Normalize();

    int maxn = this->max_n() + o.max_n() + 1;

    static vector<dcomplex> buf(maxn);
    if((int)buf.capacity() < maxn)
      buf.reserve(maxn);
    dcomplex* gs = &buf[0];


    int num(this->size_basis());

    if(vec_.find("m") == vec_.end())
      vec_["m"] = VectorXcd::Zero(num);

    VectorXcd& m = vec_["m"];

    int i(0);
    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      if(it->basis.size() != 1) {
	throw runtime_error("contraction!");
      }
      for(vector<R1GTO>::iterator it_g = it->basis.begin(), end_g = it->basis.end();
	  it_g != end_g; ++it_g) {
	for(int io = 0; io < o.size_basis(); io++ ) {
	  dcomplex c(it_g->c * o.basis(io).c);
	  int      n(it_g->n + o.basis(io).n);
	  dcomplex z(it_g->z + o.basis(io).z);
	  CalcGTOInt(n, z, gs);
	  m(i) += c * gs[n];
	}
	++i;
      }
    }
    */
    /*
    for(int i = 0; i < num; i++) {
      m[i] = 0.0;
      for(int j = 0; j < o.size_basis(); j++) {
	int      n(this->basis(i).n + o.basis(j).n);
	dcomplex z(this->basis(i).z + o.basis(j).z);
	dcomplex c(this->basis(i).c * o.basis(j).c);

	CalcGTOInt(n, z, gs);
	m[i] += c * gs[n];
      }      
    }
    */
    
  }
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
  void R1GTOs::CalcDerivVec(const R1STOs& o) { 
    
    /*
    u  = (N  r^n                         ) Exp[-zr^2]
    u' = (N' r^n - N r^(n+2)             ) Exp[-zr^2]
    u''= (N''r^n -2N'r^(n+2) + N r^(n+4) ) Exp[-zr^2]
    */
    
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
  void R1GTOs::AtR(const VectorXcd& cs, const VectorXcd& rs,
		   VectorXcd* ys) {

    if(!this->coef_set_q_) {
      string msg; SUB_LOCATION(msg); 
      msg += ": coef is not set up.";
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
  VectorXcd* R1GTOs::AtR(const VectorXcd& cs, const VectorXcd& rs) {
    VectorXcd* ys = new VectorXcd();
    this->AtR(cs, rs, ys);
    return ys;
  }
  void R1GTOs::swap(R1GTOs& o) {
    ::swap(this->coef_set_q_, o.coef_set_q_);
    this->conts_.swap(o.conts_);
    ::swap(this->L_, o.L_);
  }

  void R1STOs::Add(dcomplex c, int n, dcomplex zeta) {
    stos_.push_back(R1STO(c, n, zeta));
  }
  void R1STOs::Add(int n, dcomplex zeta) {
    stos_.push_back(R1STO(1, n, zeta));
  }

  ostream& operator<<(ostream& out, const R1GTOs& us) {
    out << "R1GTOs" << endl;
    out << "coef_set_q : "
	<< (us.coef_set_q_ ? "Yes" : "No")
	<< endl;
    out << "coef type : " << us.coef_type_ << endl;
    
    for(R1GTOs::cItCont it = us.conts_.begin();
	it != us.conts_.end(); ++it) {
      if(it->size_prim() == 1 && it->size_basis() == 1) {
	cout << it->prim[0].n << ", " << it->prim[0].z << endl;
      } else {
	cout << "offset:" << it->offset << endl;
	cout << "prim:" << endl;
	for(R1GTOs::cItPrim iprim = it->prim.begin();
	    iprim != it->prim.end(); ++iprim) {

	  cout << iprim->n << ", " << iprim->z << endl;
	  
	}
	cout << "coef:" << endl;
	cout << it->coef << endl;
	
      }
    }


    for(int i = 0; i < us.size_prim(); ++i) {
      out << us.prim(i).n << us.prim(i).z << endl;
    }
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
