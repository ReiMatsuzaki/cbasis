#ifndef OPT_ALPHA_H
#define OPT_ALPHA_H

#include <Eigen/Core>
#include "r1gtoint.hpp"

namespace l2func {

  /*
    Optimize orbital exponents of cGTO for alpha.    
   */
  void OptAlphaShift(const R1STOs& driv,
		     const Eigen::VectorXi& opt_idx,
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

  dcomplex CalcAlphaFull(const R1STOs& driv,R1GTOs& gtos, dcomplex ene);    
  dcomplex CalcAlpha(const R1STOs& driv, R1GTOs& gtos, dcomplex ene, double eps);
}
#endif
