#ifndef TWO_INT_H
#define TWO_INT_H

#include "mol_func.hpp"
#include "symmolint.hpp"

namespace l2func {

  // ==== Slow routines ====
  dcomplex ERIEle(CartGTO& a, CartGTO& b, CartGTO& c, CartGTO& d);

  // ==== SymGTOs ====
  void SymGTOs_CalcERI(SymGTOs& i, SymGTOs& j, SymGTOs& k, SymGTOs& l,
		       IB2EInt* eri);  
}

#endif
