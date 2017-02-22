#include <sstream>
#include <boost/shared_ptr.hpp>
#include "../utils/macros.hpp"
#include "r1_lc.hpp"


using namespace std;
using namespace Eigen;
using namespace boost;

namespace cbasis {
  
  // ==== Create ====
  template<int M>
  shared_ptr<_LC_EXPs<M> > Create_LC_EXPs() {
    shared_ptr<_LC_EXPs<M> > ptr(new _LC_EXPs<M>());
    return ptr;
  }
  LC_STOs Create_LC_STOs() {
    return Create_LC_EXPs<1>();
  }
  LC_GTOs Create_LC_GTOs() {
    return Create_LC_EXPs<2>();
  }
  template<int M>
  shared_ptr<_LC_EXPs<M> > Create_LC_EXPs(int n) {
    shared_ptr<_LC_EXPs<M> > ptr(new _LC_EXPs<M>());
    for(int i = 0; i < n; i++) {
      ptr->Add(1.0, 1, 1.0);
    }
    return ptr;
  }
  
  // ==== Member functions ====
  template<int M>
  int _LC_EXPs<M>::max_n() const {
    return *max_element(this->ns.begin(), this->ns.end());
  }
  template<int M>
  int _LC_EXPs<M>::size() const {
    return this->ns.size();
  }
  template<int M>
  VectorXcd _LC_EXPs<M>::AtR(const VectorXcd& rs) const {

    int num_rs(rs.size());
    VectorXcd ys = VectorXcd::Zero(num_rs);
    
    for(int i_r = 0; i_r < num_rs; i_r++) {
      
      dcomplex r(rs[i_r]);
      dcomplex rr(pow(r, M));
      dcomplex y(0);
      for(int j_b = 0; j_b < this->size(); j_b++) {
	dcomplex c(this->cs[j_b]);
	int      n(this->ns[j_b]);
	dcomplex z(this->zs[j_b]);
	y += c * pow(r, n) * exp(-z * rr);
      }
      ys[i_r] = y;

    }

    return ys;

  }
  template<int M> dcomplex Exp_DAtR(dcomplex r, dcomplex c, int n, dcomplex z);
  template<>
  dcomplex Exp_DAtR<1>(dcomplex r, dcomplex c, int n, dcomplex z) {
    dcomplex dn(n);
    dcomplex y = c * pow(r, n) * exp(-z * r) * (-z);
    if (n != 0) 
      y += c * dn * pow(r, n-1) * exp(-z * r) * (-z);
    return y;
  }
  template<>
  dcomplex Exp_DAtR<2>(dcomplex r, dcomplex c, int n, dcomplex z) {
    dcomplex dn(n);
    dcomplex y = c * pow(r, n+1) * exp(-z * r * r) * (-2.0*z);
    if (n != 0) 
      y += c * dn * pow(r, n-1) * exp(-z * r) * (-z);
    return y;
  }
  template<int M>
  VectorXcd _LC_EXPs<M>::DAtR(const VectorXcd& rs) const {
    
    int num_rs(rs.size());
    VectorXcd ys = VectorXcd::Zero(num_rs);
    
    for(int i_r = 0; i_r < num_rs; i_r++) {
      
      dcomplex r(rs[i_r]);
      dcomplex y(0);
      for(int j_b = 0; j_b < this->size(); j_b++) {
	dcomplex c(this->cs[j_b]);
	int      n(this->ns[j_b]);
	dcomplex z(this->zs[j_b]);
	y += Exp_DAtR<M>(r, c, n, z);
      }
      ys[i_r] = y;

    }

    return ys;

  }
  template<int M>
  _LC_EXPs<M>* _LC_EXPs<M>::Add(dcomplex c, int n, dcomplex z) {
    
    this->cs.push_back(c);
    this->ns.push_back(n);
    this->zs.push_back(z);
    
    return this;
  }

  template<int M>
  shared_ptr<_LC_EXPs<M> > _LC_EXPs<M>::Clone() const {
    LC_EXPs ptr;
    this->Clone(INITIAL, &ptr);
    return ptr;
  }
  template<int M>
  shared_ptr<_LC_EXPs<M> > _LC_EXPs<M>::Conj() const {
    LC_EXPs ptr;
    this->Conj(INITIAL, &ptr);
    return ptr;
  }
  template<int M>
  void _LC_EXPs<M>::Clone_Symbolic(shared_ptr<_LC_EXPs<M> > *o) const{
    *o = Create_LC_EXPs<M>(this->size());
  }
  template<int M>
  void _LC_EXPs<M>::Clone_Numeric(shared_ptr<_LC_EXPs<M> > o) const{
    for(int i = 0; i < this->size(); i++) {
      o->c(i) = this->cs[i];
      o->n(i) = this->ns[i];
      o->z(i) = this->zs[i];
    }
  }  
  template<int M>
  void _LC_EXPs<M>::Clone(int reuse, shared_ptr<_LC_EXPs<M> > *o) const{
    if(reuse == INITIAL) {
      this->Clone_Symbolic(o);
    } else if (reuse != REUSE) {
      THROW_ERROR("reuse <- {INITIAL, REUSE}");
    }
    this->Clone_Numeric(*o);
  }
  template<int M>
  void _LC_EXPs<M>::Conj_Numeric(shared_ptr<_LC_EXPs<M> > o) const{
    for(int i = 0; i < this->size(); i++) {
      o->c(i) = conj(this->cs[i]);
      o->n(i) = this->ns[i];
      o->z(i) = conj(this->zs[i]);
    }
  }  
  template<int M>
  void _LC_EXPs<M>::Conj(int reuse, shared_ptr<_LC_EXPs<M> > *o) const{
    if(reuse == INITIAL) {
      this->Clone_Symbolic(o);
    } else if (reuse != REUSE) {
      THROW_ERROR("reuse <- {INITIAL, REUSE}");
    }
    this->Conj_Numeric(*o);
  }  
  template<int M>  
  string _LC_EXPs<M>::str() const {
    ostringstream oss;
    string basis_name = (M==1 ? "STO" : "GTO");
    for(int i = 0; i < this->size(); i++) {
      oss << basis_name << "(";
      oss << this->cs[i] << ", ";
      oss << this->ns[i] << ", ";
      oss << this->zs[i] << ")";
      if (i != this->size()-1) {
	oss << " + ";
      }
    }
    return oss.str();
  }

  // ==== Realize ====
  template boost::shared_ptr<_LC_EXPs<1> > Create_LC_EXPs<1>();
  template boost::shared_ptr<_LC_EXPs<2> > Create_LC_EXPs<2>();
  template class _LC_EXPs<1>;
  template class _LC_EXPs<2>;


}
