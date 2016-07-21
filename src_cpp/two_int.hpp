#ifndef TWO_INT_H
#define TWO_INT_H

#include "mol_func.hpp"
#include "symmolint.hpp"

namespace l2func {

  // ==== Slow routines ====
  dcomplex ERIEle(CartGTO& a, CartGTO& b, CartGTO& c, CartGTO& d);

  // ==== aa ====
  void coef_R_eri_switch(dcomplex zarg,
			 dcomplex wPx, dcomplex wPy, dcomplex wPz,
			 dcomplex wPpx,dcomplex wPpy,dcomplex wPpz,
			 int max_n, dcomplex *Fjs, dcomplex mult_coef,
			 MultArray<dcomplex, 3>& res, ERIMethod method);
			 

  // ==== SymGTOs ====
  void CalcERI_Complex(SymGTOs& i, IB2EInt* eri, ERIMethod m);
  void CalcERI_Hermite(SymGTOs& i, IB2EInt* eri, ERIMethod m);
  void SymGTOs_CalcERI(SymGTOs& i, SymGTOs& j, SymGTOs& k, SymGTOs& l, IB2EInt* eri, ERIMethod method);  
	       
}

#endif
