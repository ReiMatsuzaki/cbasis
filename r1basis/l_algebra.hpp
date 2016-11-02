#ifndef L_ALGEBRA_HPP
#define L_ALGEBRA_HPP

#include <Eigen/Core>

namespace cbasis {

  // compute aA^(1,0)b
  // aA^(10)b = { \sum(i)  a_kA^(10)_ki b_i   | k }
  void Calc_a_A10_b(const Eigen::VectorXcd& a,
		    const Eigen::MatrixXcd& A10,
		    const Eigen::VectorXcd& b,
		    Eigen::VectorXcd* res);
  
  // ab^(1) = { a_k b^(1)_k   | k }
  void Calc_a_b1(const Eigen::VectorXcd& a, 
		 const Eigen::VectorXcd& b1, 
		 Eigen::VectorXcd* res);

  void Calc_a_A10_B_A10_b(const Eigen::VectorXcd& a, 
			  const Eigen::MatrixXcd& A10,
			  const Eigen::MatrixXcd& B,
			  const Eigen::VectorXcd& b, 
			  Eigen::MatrixXcd* res);

  void Calc_a_A20_b(const Eigen::VectorXcd& a, 
		    const Eigen::MatrixXcd& A20,
		    const Eigen::VectorXcd& b, 
		    Eigen::MatrixXcd* res);

  void Calc_a_A20_b(const Eigen::VectorXcd& a, 
		    const Eigen::MatrixXcd& A20,
		    const Eigen::VectorXcd& b, 
		    Eigen::MatrixXcd* res);


  

  // compute aA^jb
  // where a and b are vector and A^j is derivative of j th basis:
  // (A^j)_rs = \delta(j,r) (r'|A|s) + \delta(j,s) (r|A|s')
  // input is A10_rs = (r'|A|s)
  Eigen::VectorXcd Calc_a_Aj_b(const Eigen::VectorXcd& a,
			       const Eigen::MatrixXcd& A10,
			       const Eigen::VectorXcd& b);

  
  // compute aA^iBA^jb
  // where a and b are vector and A^j is derivative of j th basis:
  // (A^j)_rs = \delta(j,r) (r'|A|s) + \delta(j,s) (r|A|s')
  // input is A10_rs = (r'|A|s)
  Eigen::MatrixXcd Calc_a_Ai_B_Aj_b(const Eigen::VectorXcd& a,
				    const Eigen::MatrixXcd& A10,
				    const Eigen::MatrixXcd& B,
				    const Eigen::VectorXcd& b);

  
  // compute (aA^{i,j}a)_ij
  // i == j
  // a_x A^{ii}_xy a_y = a_i A^{20}_iy a_y + a_x A^{02}_xi a_i
  //                   + 2 a_i A^{11}_ii a_i
  Eigen::MatrixXcd Calc_a_Aij_a(const Eigen::VectorXcd& a,
				const Eigen::MatrixXcd& A20,
				const Eigen::MatrixXcd& A11);

  // compute (aA^{i,j}b)_xy
  // i == j
  // a_x A^{ii}_xy b_y = a_i A^{20}_iy b_y + a_x A^{02}_xi b_i
  //                   + 2 a_i A^{11}_ii b_i
  Eigen::MatrixXcd Calc_a_Aij_b(const Eigen::VectorXcd& a,
				const Eigen::MatrixXcd& A20,
				const Eigen::MatrixXcd& A11,
				const Eigen::VectorXcd& b);


  // ai_A_Bj_b = d_i(a_x) A_xy d_j(B_yz) b_z
  //           = a'_i A_ij B^(1,0)(j, z) b_z +
  //             a'_i A_ix B^(0,1)(x, j) b_j
  Eigen::MatrixXcd Calc_ai_A_Bj_b(const Eigen::VectorXcd& a1,
				  const Eigen::MatrixXcd& A,
				  const Eigen::MatrixXcd& B10,
				  const Eigen::VectorXcd& b);


  // compute a^i b = a^i_i * b_i
  Eigen::VectorXcd Calc_ai_b(const Eigen::VectorXcd& a1,
			     const Eigen::VectorXcd& b);

}

#endif
