#include <stdexcept>
#include "../src_cpp/macros.hpp"
#include "r1basis.hpp"
#include "r1_lc.hpp"

using namespace std;
using namespace Eigen;

namespace cbasis {

  // ==== Calculation ====
  dcomplex GTOInt(int n, dcomplex a) {

    if(n < 0) {
      std::string msg; SUB_LOCATION(msg);
      msg += "maxn must be bigger than 0";
      throw std::runtime_error(msg);
    }

    if(n == 0)
      return sqrt(M_PI)/(2.0*sqrt(a));
    if(n == 1)
      return 0.5/a;

    return dcomplex(n-1)/(2.0*a) * GTOInt(n-2, a);

  }

  // ==== member field ====
  _GTOs::_GTOs() {
    this->setupq_ =false;
  }
  bool _GTOs::OnlyPrim() const {

    bool acc(true);

    for(vector<LC_GTOs>::const_iterator
	  it = this->basis_.begin();
	it != this->basis_.end();
	++it) {
      acc &= (*it)->size() == 1;
    }

    return acc;
    
  }
  VectorXcd _GTOs::AtR(const VectorXcd& rs, const VectorXcd& cs) const {

    if(cs.size() != this->size()) {
      string msg;
      SUB_LOCATION(msg);
      msg += " : size mismatch.";
      throw runtime_error(msg);
    }

    int num(this->size());
    VectorXcd ys = VectorXcd::Zero(rs.size());

    for(int i = 0; i < num; i++) {
      VectorXcd ys0 = this->basis(i)->AtR(rs);
      ys += ys0;
    }

    return ys;

  }
  VectorXcd _GTOs::DAtR(const VectorXcd& rs, const VectorXcd& cs) const {

    if(cs.size() != this->size()) {
      string msg;
      SUB_LOCATION(msg);
      msg += " : size mismatch.";
      throw runtime_error(msg);
    }

    int num(this->size());
    VectorXcd ys = VectorXcd::Zero(rs.size());

    for(int i = 0; i < num; i++) {
      VectorXcd ys0 = this->basis(i)->DAtR(rs);
      ys += ys0;
    }

    return ys;

  }
  string _GTOs::str() const {
    return "GTOs";
  }

  _GTOs* _GTOs::AddPrim(int n, dcomplex z) {
    this->setupq_ = false;
    LC_GTOs g = Create_LC_GTOs();
    g->Add(1.0, n, z);
    this->basis_.push_back(g);
    return this;
  }
  _GTOs* _GTOs::AddPrims(int n, Eigen::VectorXcd zs) {
    this->setupq_ = false;
    for(int i = 0; i < zs.size(); i++)
      this->AddPrim(n, zs[i]);
    return this;
  }
  _GTOs* _GTOs::AddLC(LC_GTOs lc) {
    this->setupq_ =false;
    this->basis_.push_back(lc);
    return this;
  }
  _GTOs* _GTOs::SetUp() {

    if(this->setupq_)
      return this;

    int num(this->size());
    
    for(int i = 0; i < num; i++) {

      LC_GTOs bi = this->basis(i);

      dcomplex acc(0);
      for(int ii = 0; ii < bi->size(); ii++) {
	for(int jj = 0; jj < bi->size(); jj++) {
	  dcomplex c(bi->c(ii) * bi->c(jj));
	  int      n(bi->n(ii) + bi->n(jj));
	  dcomplex z(bi->z(ii) * bi->z(jj));
	  acc +=  c * GTOInt(n, z);
	}
      }
      dcomplex nterm(1.0/sqrt(acc));
      for(int ii = 0; ii < bi->size(); ii++) {
	bi->c(ii) *= nterm;
      }
      
    }

    this->setupq_ = true;
    return this;
  }

  GTOs _GTOs::Clone() const {

        GTOs ptr = Create_GTOs();
    typedef vector<LC_GTOs>::const_iterator It;
    for(It it = this->basis_.begin();
	it != this->basis_.end();
	++it) {
      LC_GTOs lc = (*it)->Clone();
      ptr->AddLC(*it);
    }
    ptr->SetUp();
    return ptr;

  }
  GTOs _GTOs::Conj() const {

    GTOs ptr = Create_GTOs();
    typedef vector<LC_GTOs>::const_iterator It;
    for(It it = this->basis_.begin();
	it != this->basis_.end();
	++it) {
      LC_GTOs lc = (*it)->Conj();
      ptr->AddLC(lc);
    }
    ptr->SetUp();
    return ptr;

  }

  MatrixXcd _GTOs::CalcRmMat(int m) const {
    
    if(!this->setupq_) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "GTO is not setup.";
      throw runtime_error(msg);
    }
    
    int num(this->size());
    MatrixXcd mat(num, num);
    
    for(int i = 0; i < num; i++)
      for(int j = 0; j < num; j++) {

	LC_GTOs bi = this->basis(i);
	LC_GTOs bj = this->basis(j);
	dcomplex acc(0);
	for(int ii = 0; ii < bi->size(); ii++)
	  for(int jj = 0; jj < bj->size(); jj++) {
	    dcomplex c(bi->c(ii) * bj->c(jj));
	    int      n(bi->n(ii)+bj->n(jj));
	    dcomplex z(bi->z(ii) * bj->z(jj));
	    acc +=  c * GTOInt(n+m, z);
	  }
	mat(i, j) = acc;
      }

    return mat;

  }
  MatrixXcd _GTOs::CalcD2Mat() const {

    if(!this->setupq_) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "GTO is not setup.";
      throw runtime_error(msg);
    }


    int num(this->size());
    MatrixXcd mat(num, num);
    
    for(int i = 0; i < num; i++)
      for(int j = 0; j < num; j++) {

	LC_GTOs bi = this->basis(i);
	LC_GTOs bj = this->basis(j);
	dcomplex acc(0);
	for(int ii = 0; ii < bi->size(); ii++)
	  for(int jj = 0; jj < bj->size(); jj++) {
	    dcomplex c(bi->c(ii) * bj->c(jj));
	    int      ni(bi->n(ii));
	    int      nj(bj->n(jj));
	    dcomplex zi(bi->z(ii));
	    dcomplex zj(bj->z(jj));
	    acc += c * ( +4.0*zj*zj      * GTOInt(ni+nj+2, zi+zj)
			 -2.0*(2*nj+1)   * GTOInt(ni+nj, zi+zj));
	    if(nj > 1)
	      acc += 1.0*(nj*nj-nj) * c * GTOInt(ni+nj-2, zi+zj);
	  }
	mat(i, j) = acc;
      }

    return mat;
  }
  VectorXcd _GTOs::CalcVecSTO(LC_STOs stos) const {
    
    if(!this->setupq_) {
      string msg; SUB_LOCATION(msg);
      msg += " : GTO is not setup.";
      throw runtime_error(msg);
    }

    return VectorXcd::Zero(1);
  }
  VectorXcd _GTOs::CalcVecGTO(LC_GTOs gtos) const {

    if(!this->setupq_) {
      string msg; SUB_LOCATION(msg);
      msg += " : GTO is not setup.";
      throw runtime_error(msg);
    }
    
    int num(this->size());
    int numo(gtos->size());

    VectorXcd vec(num);
        
    for(int i = 0; i < num; i++) {
      LC_GTOs bi = this->basis(i);
      dcomplex acc;
      for(int oo = 0; oo < numo; oo++) {
	for(int ii = 0; ii < num; ii++) {
	  dcomplex c(bi->c(ii) * gtos->c(oo));
	  int      n(bi->n(ii) + gtos->n(oo));
	  dcomplex z(bi->z(ii) + gtos->z(oo));
	  acc += c * GTOInt(n, z);
	}
      }
      
      vec(i) = acc;
      
    }
    
    return vec;

  }

  // ==== create =====
  GTOs Create_GTOs() {

    GTOs ptr(new _GTOs());
    return ptr;
    
  }

}
