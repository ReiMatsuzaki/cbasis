#include <exception>
#include <algorithm>
#include "r1gtoint.hpp"
#include "macros.hpp"
#include "cip_exp.hpp"
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
    L_(_L) {
    buf_.reserve(1);
  }
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
  const R1GTO& R1GTOs::basis(int i) const {
    typedef vector<Contraction>::const_iterator ContIt;
    for(ContIt it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      if(i < it->offset ) {
	ContIt it0 = it - 1;
	return it0->basis[i-it0->offset];
      }
    }

    if(i < this->size_basis()) {
      ContIt it = (conts_.end()-1);
      return it->basis[i-it->offset];
    } else {
      string msg; SUB_LOCATION(msg);
      ostringstream oss; 
      oss << msg << endl << ": Failed to find index i." << endl;
      oss << "i: " << i << endl;
      throw runtime_error(oss.str());
    }
    
  }
  R1GTO& R1GTOs::basis(int i) {
    typedef vector<Contraction>::iterator ContIt;
    for(ContIt it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      if(i < it->offset ) {
	ContIt it0 = it - 1;
	return it0->basis[i-it0->offset];
      }
    }

    if(i < this->size_basis()) {
      ContIt it = (conts_.end()-1);
      return it->basis[i-it->offset];
    } else {
      string msg; SUB_LOCATION(msg);
      ostringstream oss; 
      oss << msg << endl << ": Failed to find index i." << endl;
      oss << "i: " << i << endl;
      throw runtime_error(oss.str());
    }
  }
  void R1GTOs::Add(dcomplex c, int _n, dcomplex _zeta) {
    this->normalized_q_ = false;
    this->calc_mat_q_ = false;
    this->calc_vec_q_ = false;

    Contraction cont;
    cont.basis.push_back(R1GTO(c, _n, _zeta));
    cont.coef = MatrixXcd::Ones(1, 1);
    cont.offset = this->size_basis();
    conts_.push_back(cont);
  }
  void R1GTOs::Add(int _n, dcomplex _zeta) {
    this->Add(1.0, _n, _zeta);
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
      cont.basis.push_back(R1GTO(1, n, zs(i)));
    cont.coef = coef;
    cont.offset = this->size_basis();
    this->conts_.push_back(cont);
    
  }
  void R1GTOs::Set(int n, const Eigen::VectorXcd& zs) {
    
    if(zs.size() != this->size_basis()) {
      string msg; SUB_LOCATION(msg);
      msg += ": size mismatch.";
      throw runtime_error(msg);
    }

    for(vector<Contraction>::iterator it = conts_.begin(), end = conts_.end();
	it != end; ++it) {
      int i(it->offset);
      for(vector<R1GTO>::iterator it_g = it->basis.begin(), end_g = it->basis.end();
	  it_g != end_g; ++it_g, ++i) {
	it_g->n = n;
	it_g->z = zs[i];
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

    typedef vector<Contraction>::const_iterator ContIt;
    typedef vector<R1GTO>::const_iterator BasisIt;
    int maxn(0);
    for(ContIt it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      for(BasisIt it_g = it->basis.begin(), end_g = it->basis.end();
	  it_g != end_g; ++it_g) {
	if(maxn < it_g->n)
	  maxn = it_g->n;
      }
    }

    return maxn;    
  }
  void R1GTOs::Normalize() {

    int maxn = this->max_n() + this->max_n() + 2 + 1;
    this->Reserve(maxn);
    dcomplex* gs = &buf_[0];

    //    MultArray<dcomplex, 2> smat_prim(10*10);

    typedef vector<Contraction>::iterator ItCont;
    typedef vector<R1GTO>::iterator ItPrim;

    for(ItCont it = conts_.begin(), end = conts_.end(); it != end; ++it) {
      for(ItPrim it_g = it->basis.begin(), end_g = it->basis.end();
	  it_g != end_g; ++it_g) {
	  
	dcomplex z(it_g->z);
	int n(it_g->n);
	CalcGTOInt(2*n, 2.0*z, gs);
	dcomplex norm2 = gs[n*2];
	it_g->c = 1.0/sqrt(norm2);
      }
    }

    this->normalized_q_ = true;
  }
  void R1GTOs::Reserve(int n) {

    if((int)buf_.capacity() < n)
      buf_.reserve(n);

  }
  void R1GTOs::CalcMat() {
    if(!this->normalized_q())
      this->Normalize();

    int maxn = 2*this->max_n() + 2 + 1;
    this->Reserve(maxn);
    dcomplex* gs = &buf_[0];
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
      if(it->basis.size() != 1) {
	throw runtime_error("contraction!");
      }
      for(ContIt jt = conts_.begin(); jt != end; ++jt) {
	int i(it->offset);
	for(BasisIt ib = it->basis.begin(), end_ib = it->basis.end();
	    ib != end_ib; ++ib, ++i) {
	  int j(jt->offset);
	  for(BasisIt jb = jt->basis.begin(), end_jb = jt->basis.end();
	      jb != end_jb; ++jb, ++j) {
	    
	    dcomplex cc(ib->c * jb->c);
	    int      ni(ib->n);
	    int      nj(jb->n);
	    dcomplex zi(ib->z);
	    dcomplex zj(jb->z);

	    CalcGTOInt(ni+nj+2, zi+zj, gs);
	    
	    dcomplex sval = cc*gs[ni+nj];
	    dcomplex tval(0.0);
	    tval -= zj * dcomplex(4*nj+2) * gs[ni+nj];
	    if(ni+nj >= 2)
	      tval += dcomplex(nj*nj-nj-L_*L_-L_) * gs[ni+nj-2];
	    tval += 4.0*zj*zj*gs[ni+nj+2];
	    tval *= -0.5*cc;
	    dcomplex vval = -cc*gs[ni+nj-1];

	    
	    s(i, j) = sval;
	    t(i, j) = tval;
	    v(i, j) = vval;
	  
	    if(i != j) {
	      s(j, i) = sval;
	      t(j, i) = tval;
	      v(j, i) = vval;
	    }

	  }
	}
      }
    }


/*
    for(int i = 0; i < num; i++) {
      for(int j = i; j < num; j++) {
	int ni(this->gtos_[i].n);
	int nj(    this->gtos_[j].n);
	dcomplex zi(this->gtos_[i].z);
	dcomplex zj(this->gtos_[j].z);
	dcomplex cc = this->gtos_[i].c * this->gtos_[j].c;

	CalcGTOInt(ni+nj+2, zi+zj, gs);
	
	dcomplex sval = cc*gs[ni+nj];
	dcomplex tval(0.0);
	tval -= zj * dcomplex(4*nj+2) * gs[ni+nj];
	if(ni+nj >= 2)
	  tval += dcomplex(nj*nj-nj-L_*L_-L_) * gs[ni+nj-2];
	tval += 4.0*zj*zj*gs[ni+nj+2];
	tval *= -0.5*cc;
	dcomplex vval = -cc*gs[ni+nj-1];

	s(i, j) = sval;
	t(i, j) = tval;
	v(i, j) = vval;
	  
	if(i != j) {
	  s(j, i) = sval;
	  t(j, i) = tval;
	  v(j, i) = vval;
	}
      }
    }
*/    

    this->calc_mat_q_ = true;
  }
  void R1GTOs::CalcVec(const R1GTOs& o) {
    
    if(!this->normalized_q())
      this->Normalize();

    int maxn = this->max_n() + o.max_n() + 1;
    this->Reserve(maxn);
    dcomplex* gs = &buf_[0];

    int num(this->size_basis());

    if(vec_.find("m") == vec_.end())
      vec_["m"] = VectorXcd::Zero(num);

    VectorXcd& m = vec_["m"];

    int i(0);
    for(vector<Contraction>::iterator it = conts_.begin(), end = conts_.end();
	it != end; ++it) {
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

    if(!this->normalized_q())
      this->Normalize();

    int num(this->size_basis());
    if(vec_.find("m") == vec_.end())
      vec_["m"] = VectorXcd::Zero(num);
    VectorXcd& m = vec_["m"];

    for(vector<Contraction>::iterator it = conts_.begin(), end = conts_.end();
	it != end; ++it) {
      if(it->basis.size() != 1) {
	throw runtime_error("contraction!");
      }
      int i(it->offset);
      for(vector<R1GTO>::iterator it_g = it->basis.begin(), end_g = it->basis.end();
	  it_g != end_g; ++it_g, ++i) {
	m(i) = 0.0;
	for(int io = 0; io < o.size_basis(); io++ ) {
	  dcomplex c(it_g->c * o.basis(io).c);
	  int      n(it_g->n + o.basis(io).n);
	  m(i) += c * STO_GTO_Int(o.basis(io).z, it_g->z, n);
	}
      }
    }


    /*
    for(int i = 0; i < this->size_basis(); i++) {
      m(i) = 0.0;
      for(int j = 0; j < o.size_basis(); ++j) {
	int n(o.basis(j).n + this->basis(i).n);
	dcomplex c(o.basis(j).c * this->basis(i).c);
	m(i) += c * STO_GTO_Int(o.basis(j).z, this->basis(i).z, n);
      }
    }
    */

    calc_vec_q_ = true;

  }
  void R1GTOs::AtR(const VectorXcd& cs, const VectorXcd& rs,
		   VectorXcd* ys) {

    if(cs.size() != this->size_basis()) {
      string msg; SUB_LOCATION(msg);
      msg += ": \n size of cs must be equal to basis size";
      throw runtime_error(msg);
    }
    
    VectorXcd res = VectorXcd::Zero(rs.size());
    for(int i = 0; i < rs.size(); i++) {
      dcomplex r(rs[i]);
      dcomplex cumsum(0);

      for(vector<Contraction>::iterator it = conts_.begin(), end = conts_.end();
	it != end; ++it) {
	int ib(it->offset);
	for(vector<R1GTO>::iterator it_g = it->basis.begin(), end_g = it->basis.end();
	    it_g != end_g; ++it_g, ++ib) {
	  dcomplex c(it_g->c);
	  int      n(it_g->n);
	  dcomplex z(it_g->z);
	  cumsum += cs[ib] * c * pow(r, n) * exp(-z*r*r);
	}
      }

	  /*
      for(int ib = 0; ib < this->size_basis(); ib++) {
	R1GTO b(basis(ib));
	cumsum += cs[ib] * b.c * pow(r, b.n) * exp(-b.z*r*r);
      }
      */
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
      out << us.basis(i) << endl;
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
