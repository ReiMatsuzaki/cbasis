#include <exception>
#include "r1gtoint.hpp"
#include "macros.hpp"
#include "cip_exp.hpp"

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

  R1GTOs::R1GTOs(int _L): normalized_q_(false), L_(_L) {
    buf_.reserve(1);
  }
  void R1GTOs::Add(dcomplex c, int _n, dcomplex _zeta) {
    this->normalized_q_ = false;
    gtos_.push_back(R1GTO(c, _n, _zeta));
  }
  void R1GTOs::Add(int _n, dcomplex _zeta) {
    this->normalized_q_ = false;
    gtos_.push_back(R1GTO(1.0, _n, _zeta));
  }
  void R1GTOs::Add(int n, const VectorXcd& zs) {

    this->normalized_q_ = false;
    for(int i = 0; i < zs.size(); i++) {
      this->Add(n, zs[i]);
    }
    
  }
  void R1GTOs::Set(int n, const Eigen::VectorXcd& zs) {
    
    if(zs.size() != this->size_basis()) {
      string msg; SUB_LOCATION(msg);
      msg += ": size mismatch.";
      throw runtime_error(msg);
    }

    this->normalized_q_ = false;
    for(int i = 0; i < zs.size(); i++) {
      gtos_[i].n = n;
      gtos_[i].z = zs(i);
    }

    this->normalized_q_ = false;

  }
  int  R1GTOs::max_n() const {
    int maxn(0);
    for(CIt it = this->gtos_.begin(); it != this->gtos_.end(); ++it)
      if(maxn < it->n)
	maxn = it->n;
    return maxn;    
  }
  void R1GTOs::Normalize() {

    int maxn = this->max_n() + this->max_n() + 2 + 1;
    this->Reserve(maxn);
    dcomplex* gs = &buf_[0];

    for(int i = 0; i < size_basis(); i++) {
      dcomplex z(this->basis(i).z);
      int n(this->basis(i).n);
      CalcGTOInt(2*n, 2.0*z, gs);
      dcomplex norm2 = gs[n*2];
      this->gtos_[i].c = 1.0/sqrt(norm2);
    }    
    this->normalized_q_ = true;
  }
  void R1GTOs::Reserve(int n) {

    if((int)buf_.capacity() < n)
      buf_.reserve(n);

  }
  void R1GTOs::CalcMat(MatMap* res) {

    if(!this->normalized_q())
      this->Normalize();

    int maxn = 2*this->max_n() + 2 + 1;
    this->Reserve(maxn);
    dcomplex* gs = &buf_[0];
    int num = size_basis();    

    if(res->find("s") == res->end())
      (*res)["s"] = MatrixXcd::Zero(num, num);
    if(res->find("t") == res->end())
      (*res)["t"] = MatrixXcd::Zero(num, num);
    if(res->find("v") == res->end())
      (*res)["v"] = MatrixXcd::Zero(num, num);

    MatrixXcd& s = (*res)["s"];
    MatrixXcd& t = (*res)["t"];
    MatrixXcd& v = (*res)["v"];
    
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

    // delete[] gs;
  }
  MatMap* R1GTOs::CalcMatNew() {
    MatMap* res = new MatMap();
    this->CalcMat(res);
    return res;
  }
  void R1GTOs::CalcVec(R1GTOs& o, VecMap* res) {
    
    if(!this->normalized_q())
      this->Normalize();

    int maxn = this->max_n() + o.max_n() + 1;
    this->Reserve(maxn);
    dcomplex* gs = &buf_[0];

    int num(this->size_basis());

    if(res->find("m") == res->end())
      (*res)["m"] = VectorXcd::Zero(num);

    VectorXcd& m = (*res)["m"];

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

    (*res)["m"] = MatrixXcd::Zero(1, 1);
    (*res)["m"].swap(m);
  }
  void R1GTOs::CalcVecSTO(const R1STOs& o, VecMap* res) {

    if(!this->normalized_q())
      this->Normalize();

    int num(this->size_basis());
    if(res->find("m") == res->end())
      (*res)["m"] = VectorXcd::Zero(num);
    VectorXcd& m = (*res)["m"];

    for(int i = 0; i < this->size_basis(); i++) {
      m(i) = 0.0;
      for(int j = 0; j < o.size_basis(); ++j) {
	int n(o.basis(j).n + this->basis(i).n);
	dcomplex c(o.basis(j).c * this->basis(i).c);
	m(i) += c * STO_GTO_Int(o.basis(j).z, this->basis(i).z, n);
      }
    }
  }
  VecMap* R1GTOs::CalcVecSTONew(const R1STOs& vs) {
    VecMap* res = new VecMap();
    this->CalcVecSTO(vs, res);
    return res;
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
      for(int ib = 0; ib < this->size_basis(); ib++) {
	R1GTO b(basis(ib));
	cumsum += cs[ib] * b.c * pow(r, b.n) * exp(-b.z*r*r);
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
    this->gtos_.swap(o.gtos_);
    bool tmp(this->normalized_q_);
    this->normalized_q_ = o.normalized_q_;
    o.normalized_q_     = tmp;
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
