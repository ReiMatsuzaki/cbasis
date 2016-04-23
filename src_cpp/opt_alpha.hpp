#ifndef OPT_ALPHA_H
#define OPT_ALPHA_H

#include <Eigen/Core>
#include "r1gtoint.hpp"

namespace l2func {

  void AlphaGrad(const Eigen::VectorXcd& d, const Eigen::VectorXcd& c,
   		 const Eigen::MatrixXcd& D10,
		 const Eigen::VectorXcd& r1,
		 const Eigen::VectorXcd& s0, const Eigen::VectorXcd& s1,
		 dcomplex *a, Eigen::VectorXcd *g);
  
  void AlphaGradHess(const Eigen::VectorXcd& c,
		     const Eigen::VectorXcd& d,
		     const Eigen::MatrixXcd& U,
		     const Eigen::MatrixXcd& D00,
		     const Eigen::MatrixXcd& D10,
		     const Eigen::MatrixXcd& D20,
		     const Eigen::MatrixXcd& D11,
		     const Eigen::VectorXcd& r0,
		     const Eigen::VectorXcd& r1,
		     const Eigen::VectorXcd& r2,
		     const Eigen::VectorXcd& s0,
		     const Eigen::VectorXcd& s1,
		     const Eigen::VectorXcd& s2,
		     dcomplex* a, Eigen::VectorXcd* g, Eigen::MatrixXcd* h);


  /*
    Optimize orbital exponents of cGTO for alpha.    
   */
  void OptAlphaShift(const R1STOs& driv,
		     const Eigen::VectorXi& opt_idx,
		     int L,
		     dcomplex ene,
		     double h,
		     int max_iter,
		     double eps,
		     double eps_canonical,
		     R1GTOs& gtos,
		     dcomplex& z_shift,
		     bool* convq,
		     dcomplex* alpha);
  /*
    Compute coefficient for driven type equation.
    Notice that calculation speed is not fast.
   */
  Eigen::VectorXcd SolveAlpha(const R1STOs& driv, R1GTOs& gtos, dcomplex ene, double eps);
  Eigen::VectorXcd SolveAlpha(const R1STOs& driv, R1GTOs& gtos, dcomplex ene,
			      const R1STOs& sto_pot);

  dcomplex CalcAlphaFull(const R1STOs& driv,R1GTOs& gtos, int L, dcomplex ene);
  dcomplex CalcAlpha(const R1STOs& driv, R1GTOs& gtos, int L,  dcomplex ene,
		     double eps);
  dcomplex CalcAlpha(const R1STOs& driv, R1GTOs& gtos, int L, dcomplex ene, const R1STOs&);
		     
}
#endif
