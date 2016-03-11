#ifndef SPEC_FUNC_H
#define SPEC_FUNC_H

#include "math_utils.hpp"

namespace l2func {

  // ==== Incomplete Gamma ====
  // K.Ishida J.Comput.Chem. 25, (2004), 739
  // F1 and F2 algorithms 
  // Warning:
  // Eq. (24)(25)(26) are incorrect.
  void IncompleteGamma(int max_m, dcomplex z, dcomplex* res);

  dcomplex coef_d(dcomplex zetap,
		  dcomplex wPk, dcomplex wAk, dcomplex wBk,
		  int nAk, int nBk, int Nk);

  MultArray<dcomplex, 3> calc_d_coef(int max_ni, int max_nj, int max_n,
				   dcomplex zetaP, dcomplex wPx,
				   dcomplex xi, dcomplex xj, dcomplex* buf);


  dcomplex coef_d(dcomplex zetap,
		  dcomplex wPk, dcomplex wAk, dcomplex wBk,
		  int nAk, int nBk, int Nk);

  dcomplex coef_R(dcomplex zetaP,
		  dcomplex wPx, dcomplex wPy, dcomplex wPz,
		  dcomplex cx,  dcomplex cy,  dcomplex cz,
		  int mx, int my, int mz, int j, dcomplex* Fjs);
}

#endif
