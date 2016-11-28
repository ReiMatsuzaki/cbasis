#include <stdexcept>
#include "../utils/macros.hpp"
#include "../math/int_exp.hpp"
#include "opt_green.hpp"

using namespace std;
using namespace Eigen;

namespace cbasis {

  template<int m_s, int m, int m_r>
  OptGreen<m_s, m, m_r>::OptGreen(LC_EXPs_S _driv_S, EXPs _basis, LC_EXPs_R _driv_R,
				  int _L, double _Z, double _W): driv_S_(_driv_S), basis_(_basis), driv_R_(_driv_R), L_(_L), Z_(_Z), W_(_W) {
    int num = _basis->size();
    int num_opt = _basis->size();
    L00_ = Mat::Zero(num, num);
    L01_ = Mat::Zero(num, num_opt);
    L11_ = Mat::Zero(num_opt, num_opt);
    L02_ = Mat::Zero(num, num_opt);

    S0_ = Vec::Zero(num);
    S1_ = Vec::Zero(num_opt);
    S2_ = Vec::Zero(num_opt);
    R0_ = Vec::Zero(num);
    R1_ = Vec::Zero(num_opt);
    R2_ = Vec::Zero(num_opt);    
    
  }
  
  template<int m>
  void CalcL00(typename _EXPs<m>::EXPs us, vector<dcomplex>& buf, 
	       int L, double Z, double W, MatrixXcd *L00);
		  
  template<>
  void CalcL00<1>(STOs us, std::vector<dcomplex>& buf,
		  int L, double Z, double W, MatrixXcd *L00) {

    if(not us->HasCoefAll()) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "No coefficient.";
      throw runtime_error(msg);
    }

    int num(us->size());
    
    for(int i = 0; i < num; i++) {
      for(int j = 0; j < num; j++) {
	
	_EXPs<1>::LC_EXPs bi = us->basis(i);
	_EXPs<1>::LC_EXPs bj = us->basis(j);
	
	dcomplex acc(0);
	for(int ii = 0; ii < bi->size(); ii++) {
	  for(int jj = 0; jj < bj->size(); jj++) {
	    dcomplex c(bi->c(ii) * bj->c(jj));
	    int      ni(bi->n(ii));
	    int      nj(bj->n(jj));
	    dcomplex zi(bi->z(ii));
	    dcomplex zj(bj->z(jj));
	    STOInt_Rplus_array(ni+nj, zi+zj, &buf);
	    acc += c * (+W * buf[ni+nj]
			+0.5*(zj*zj*buf[ni+nj] 
			      -2.0*nj*zj*buf[ni+nj-1]
			      + (nj>1 ? dcomplex(nj*nj-nj) * buf[ni+nj-2] : 0.0))
			-0.5*L*(L+1) * buf[ni+nj-2]
			+Z * buf[ni+nj-1]);
	  }
	}
	(*L00).coeffRef(i, j) = acc;
      }
    }
  }

  
  template<int m_s, int m, int m_r>
  void OptGreen<m_s,m,m_r>::Calc_S0_L00_R0() {
    cbasis::CalcL00<m>(this->basis_, this->buf,
		       this->L_, this->Z_, this->W_, &this->L00_);
  }

  // ==== Explicit instance ====
  template class OptGreen<1,1,1>;
}
