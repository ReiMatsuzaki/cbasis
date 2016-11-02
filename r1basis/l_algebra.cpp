#include "l_algebra.hpp"
#include "typedef.hpp"
using namespace Eigen;

namespace cbasis {


  void Calc_a_A10_b(const VectorXcd& a,
		    const MatrixXcd& A10,
		    const VectorXcd& b,
		    VectorXcd* res) {

    VectorXcd A10_b = A10 * b;
    *res = a.array() * A10_b.array() ;

  }

  void Calc_a_b1(const VectorXcd& a, 
		 const VectorXcd& b1, 
		 VectorXcd* res) {
    *res = a.array() * b1.array();
  }

  VectorXcd Calc_a_Aj_b(const VectorXcd& a,
			const MatrixXcd& A10,
			const VectorXcd& b) {

    VectorXcd A10_b = A10 * b;
    VectorXcd res1  = a.array() * A10_b.array() ;

    VectorXcd A10_a = A10 * a;
    VectorXcd res2 = b.array() * A10_a.array();

    return res1 + res2;
  }

  MatrixXcd Calc_a_Ai_B_Aj_b(const VectorXcd& a,
			     const MatrixXcd& A10,
			     const MatrixXcd& B,
			     const VectorXcd& b) {

    typedef MatrixXcd M;

    M tmp1 = (a * (A10*b).transpose());
    M tmp2 = A10 * B;
    M t1 = tmp1.array() * tmp2.array();

    M t2 = ((A10*a) * (A10*b).transpose()).array() * B.array();
    M t3 = (A10*B*A10.transpose()).array() * (a*b.transpose()).array();
    M t4 = (A10*a * b.transpose() ).array()
	    * (B*A10.transpose()).array();
    return t1 + t2 + t3 + t4;

  }

  MatrixXcd Calc_a_Aij_a(const VectorXcd& a,
			 const MatrixXcd& A20,
			 const MatrixXcd& A11){
    MatrixXcd res;
    res = (a * a.transpose()).array() * A11.array();
    res +=  ((A20 * a).array() * a.array()).matrix().asDiagonal();
    res *= dcomplex(2);
    return res;
  }

  MatrixXcd Calc_a_Aij_b(const VectorXcd& a,
			 const MatrixXcd& A20,
			 const MatrixXcd& A11,
			 const VectorXcd& b) {
    MatrixXcd res;
    res = (a * b.transpose()).array() * A11.array();
    res+= ((b * a.transpose()).array() * A11.array()).matrix();
    res+= ((A20 * a).array() * b.array()).matrix().asDiagonal();
    res+= ((A20 * b).array() * a.array()).matrix().asDiagonal();
    return res;

  }

  MatrixXcd Calc_ai_A_Bj_b(const VectorXcd& a1,
			   const MatrixXcd& A,
			   const MatrixXcd& B10,
			   const VectorXcd& b) {

    return 
      (a1 * (B10 * b).transpose()).array() * A.array() +
      (a1 * b.transpose()).array() * (A * B10.transpose()).array();
  }

  VectorXcd Calc_ai_b(const VectorXcd& a1,
		      const VectorXcd& b) {

    return  a1.array() * b.array();
  }


}
