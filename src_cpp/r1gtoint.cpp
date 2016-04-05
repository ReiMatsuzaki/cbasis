#include <exception>
#include <algorithm>
#include "r1gtoint.hpp"
#include "macros.hpp"
#include "exp_int.hpp"
#include "mult_array.hpp"

using namespace Eigen;
using namespace std;

namespace l2func {

  R1GTO::R1GTO(dcomplex _c, int _n, dcomplex _z): c(_c), n(_n), z(_z) {
    if(n < 1) {
      string msg; SUB_LOCATION(msg);
      msg += "n must be positive.";
      throw runtime_error(msg);
    }
  }
  R1STO::R1STO(dcomplex _c, int _n, dcomplex _z): c(_c), n(_n), z(_z) {
    if(n < 1) {
      string msg; SUB_LOCATION(msg);
      msg += "n must be positive.";
      throw runtime_error(msg);
    }
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
    normalized_q_(false),
    calc_mat_q_(false),
    calc_vec_q_(false),
    L_(_L) { }
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
  const R1GTOs::Prim& R1GTOs::basis(int i) const {

    typedef vector<Contraction>::const_iterator ContIt;
    int cumsum(0);

    for(ContIt it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      if(cumsum <= i && i < cumsum + it->size_prim() )
	return it->basis[i-cumsum];
      cumsum += it->size_prim();
    }

    string msg; SUB_LOCATION(msg);
    ostringstream oss; 
    oss << msg << endl << ": Failed to find index i." << endl;
    oss << "i: " << i << endl;
    throw runtime_error(oss.str());
    
  }
  R1GTOs::Prim& R1GTOs::basis(int i) {
    typedef vector<Contraction>::iterator ContIt;
    int cumsum(0);

    for(ContIt it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      if(cumsum <= i && i < cumsum + it->size_prim() )
	return it->basis[i-cumsum];
      cumsum += it->size_prim();
    }

    string msg; SUB_LOCATION(msg);
    ostringstream oss; 
    oss << msg << endl << ": Failed to find index i." << endl;
    oss << "i: " << i << endl;
    throw runtime_error(oss.str());
  }
  void R1GTOs::Add(int _n, dcomplex _zeta) {
    this->normalized_q_ = false;
    this->calc_mat_q_ = false;
    this->calc_vec_q_ = false;
    Contraction cont;
    cont.basis.push_back(Prim(_n, _zeta));
    cont.coef = MatrixXcd::Ones(1, 1);
    cont.offset = this->size_basis();
    conts_.push_back(cont);    
  }
  void R1GTOs::Add(int n, const VectorXcd& zs) {

    for(int i = 0; i < zs.size(); i++) {
      this->Add(n, zs[i]);
    }
    
  }
  void R1GTOs::Add(int n, const VectorXcd& zs, const MatrixXcd& coef) {

    if(zs.size() != coef.cols()) {
      string msg; SUB_LOCATION(msg);
      throw runtime_error(msg);
    }

    this->normalized_q_ = false;
    this->calc_mat_q_ = false;
    this->calc_vec_q_ = false;

    Contraction cont;
    for(int i = 0; i < zs.size(); i++) 
      cont.basis.push_back(Prim(n, zs(i)));
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
	it->basis[ip].n = n;
	it->basis[ip].z = zs[i];
      }
    }

    this->normalized_q_ = false;
    this->calc_mat_q_ = false;
    this->calc_vec_q_ = false;

  }
  MatrixXcd& R1GTOs::mat(string label) {

    if(!this->calc_mat_q_) {
      string msg; SUB_LOCATION(msg);
      msg += ": matrix is calculated yet.";
      throw runtime_error(msg);
    }

    return mat_[label];
  }
  VectorXcd& R1GTOs::vec(string label) {

    if(!this->calc_vec_q_) {
      string msg; SUB_LOCATION(msg);
      msg += ": vector is calculated yet.";
      throw runtime_error(msg);
    }
    return vec_[label];
  }
  int  R1GTOs::max_n() const {

    int maxn(0);
    for(cItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      for(cItPrim it_p = it->basis.begin(), end_p = it->basis.end();
	  it_p != end_p; ++it_p) {
	if(maxn < it_p->n)
	  maxn = it_p->n;
      }
    }

    return maxn;    
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
	    int      n(it->basis[i].n + it->basis[j].n);
	    dcomplex z(it->basis[i].z + it->basis[j].z);
	    CalcGTOInt(n, z, gs);
	    norm2 += c*gs[n];
	  }
	}
	dcomplex scale(1.0/sqrt(norm2));
	for(int i = 0; i < it->size_prim(); i++) {
	  it->coef(ibasis, i) *= scale;
	}
      }
    }

    this->normalized_q_ = true;
  }
  void R1GTOs::CalcMat() {
    if(!this->normalized_q())
      this->Normalize();

    static MultArray<dcomplex, 2> s_prim_(1000);
    static MultArray<dcomplex, 2> t_prim_(1000);
    static MultArray<dcomplex, 2> v_prim_(1000);

    int maxn = 2*this->max_n() + 2 + 1;
    // this->Reserve(maxn);
    // dcomplex* gs = &buf_[0];
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

    typedef vector<Contraction>::const_iterator ContIt;
    typedef vector<R1GTO>::const_iterator BasisIt;
    for(ContIt it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      for(ContIt jt = conts_.begin(); jt != end; ++jt) {

	s_prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	t_prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	v_prim_.SetRange(0, it->size_prim(), 0, jt->size_prim());
	for(int i = 0; i < it->size_prim(); i++) {
	  for(int j = 0; j < jt->size_prim(); j++) {
	    
	    int      ni(it->basis[i].n);
	    int      nj(jt->basis[j].n);
	    dcomplex zi(it->basis[i].z);
	    dcomplex zj(jt->basis[j].z);
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

    this->calc_mat_q_ = true;
  }
  void R1GTOs::CalcVec(const R1GTOs& o) {
    string msg; SUB_LOCATION(msg);
    msg += ": Not implemented yet";
    throw runtime_error(msg);
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
    
    calc_vec_q_ = true;
    
  }
  void R1GTOs::CalcVec(const R1STOs& o) {

    static MultArray<dcomplex, 1> m_prim_(20);
    
    if(!this->normalized_q())
      this->Normalize();
    
    int num(this->size_basis());
    if(vec_.find("m") == vec_.end())
      vec_["m"] = VectorXcd::Zero(num);
    VectorXcd& m = vec_["m"];
    
    typedef vector<Contraction>::iterator ContIt;
    for(ContIt it = conts_.begin(), end = conts_.end(); it != end; ++it) {

      int num_prim(it->size_prim());
      m_prim_.SetRange(0, num_prim);
      for(int i = 0; i < num_prim; i++ ) {
	dcomplex cumsum(0);
	for(int io = 0; io < o.size_basis(); io++ ) {
	  dcomplex c(o.basis(io).c);
	  int      n(it->basis[i].n + o.basis(io).n);
	  dcomplex z(it->basis[i].z);
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
    calc_vec_q_ = true;
  }
  void R1GTOs::AtR(const VectorXcd& cs, const VectorXcd& rs,
		   VectorXcd* ys) {

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
	    int      n(it->basis[ip].n); 
	    dcomplex z(it->basis[ip].z); 
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
    ::swap(this->normalized_q_, o.normalized_q_);
    ::swap(this->calc_mat_q_, o.calc_mat_q_);
    ::swap(this->calc_vec_q_, o.calc_vec_q_);
    this->conts_.swap(o.conts_);
    ::swap(this->L_, o.L_);
  }

  void R1STOs::Add(dcomplex c, int n, dcomplex zeta) {
    this->normalized_q_ = false;
    stos_.push_back(R1STO(c, n, zeta));
  }
  void R1STOs::Add(int n, dcomplex zeta) {
    this->normalized_q_ = false;
    stos_.push_back(R1STO(1, n, zeta));
  }

  ostream& operator<<(ostream& out, const R1GTOs& us) {
    out << "R1GTOs" << endl;
    out << "normalized_q : "
	<< (us.normalized_q() ? "Yes" : "No")
	<< endl;
    for(int i = 0; i < us.size_basis(); ++i) {
      out << us.basis(i).n << us.basis(i).z << endl;
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
