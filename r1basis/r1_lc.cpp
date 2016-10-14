#include <sstream>
#include "r1_lc.hpp"


using namespace std;
using namespace Eigen;

namespace cbasis {

  LC_STOs create_LC_STOs() {
    LC_STOs ptr(new _LC_STOs());
    return ptr;
  }
  int _LC_STOs::max_n() const {
    return *max_element(this->ns.begin(), this->ns.end());
  }
  int _LC_STOs::size() const {
    return this->ns.size();
  }
  VectorXcd _LC_STOs::AtR(const VectorXcd& rs) const {

    int num_rs(rs.size());
    VectorXcd ys = VectorXcd::Zero(num_rs);
    
    for(int i_r = 0; i_r < num_rs; i_r++) {
      
      dcomplex r(rs[i_r]);
      dcomplex y(0);
      for(int j_b = 0; j_b < this->size(); j_b++) {
	dcomplex c(this->cs[j_b]);
	int      n(this->ns[j_b]);
	dcomplex z(this->zs[j_b]);
	y += c * pow(r, n) * exp(-z * r);
      }
      ys[i_r] = y;

    }

    return ys;

  }
  VectorXcd _LC_STOs::DAtR(const VectorXcd& rs) const {
    
    int num_rs(rs.size());
    VectorXcd ys = VectorXcd::Zero(num_rs);
    
    for(int i_r = 0; i_r < num_rs; i_r++) {
      
      dcomplex r(rs[i_r]);
      dcomplex y(0);
      for(int j_b = 0; j_b < this->size(); j_b++) {
	dcomplex c(this->cs[j_b]);
	int      n(this->ns[j_b]);
	double  dn(n);
	dcomplex z(this->zs[j_b]);
	y += c * pow(r, n) * exp(-z * r) * (-z);
	if (n != 0) 
	  y += c * dn * pow(r, n-1) * exp(-z * r) * (-z);
	
      }
      ys[i_r] = y;

    }

    return ys;

  }
  //VectorXcd _LC_STOs::D2AtR(const VectorXcd& rs, VectorXcd& ys) const 

  //  void _LC_STOs::AddOne(dcomplex c, int n, dcomplex z) {

  
  _LC_STOs* _LC_STOs::Add(dcomplex c, int n, dcomplex z) {
    
    this->cs.push_back(c);
    this->ns.push_back(n);
    this->zs.push_back(z);
    
    return this;
  }
  LC_STOs _LC_STOs::Conj() const {
    LC_STOs stos = create_LC_STOs();

    for(int i = 0; i < this->size(); i++) {
      dcomplex c(this->cs[i]);
      int      n(this->ns[i]);
      dcomplex z(this->zs[i]);
      stos->Add(conj(c), n, conj(z));
    }

    return stos;
  }
  string _LC_STOs::str() const {

    ostringstream oss;
    oss << "==== LC_STOs ====" << endl;
    oss << "number : " << this->size() << endl;
    for(int i = 0; i < this->size(); i++) {
      oss << i << " : ";
      oss << this->cs[i] << ", ";
      oss << this->ns[i] << ", ";
      oss << this->zs[i] << ", ";
      oss << endl;
    }
    return oss.str();
  }

  LC_GTOs create_LC_GTOs() {
    LC_GTOs ptr(new _LC_GTOs());
    return ptr;
  }
  int _LC_GTOs::max_n() const {
    return *max_element(this->ns.begin(), this->ns.end());
  }
  int _LC_GTOs::size() const {
    return this->ns.size();
  }
  void _LC_GTOs::AtR(const VectorXcd& rs, VectorXcd& ys) const {

  }
  void _LC_GTOs::DAtR(const VectorXcd& rs, VectorXcd& ys) const {

  }
  void _LC_GTOs::D2AtR(const VectorXcd& rs, VectorXcd& ys) const {

  }
  _LC_GTOs* _LC_GTOs::Add(dcomplex c, int n, dcomplex z) {
    return this;
  }
  LC_GTOs _LC_GTOs::Conj() const {
    LC_GTOs gtos = create_LC_GTOs();
    return gtos;
  }

}
