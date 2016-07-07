#ifndef ONE_INT_H
#define ONE_INT_H

#include <Eigen/Core>
#include "typedef.hpp"
#include "mol_func.hpp"
#include "symmolint.hpp"

namespace l2func {

  // ==== Definition ====
  typedef vector<SubSymGTOs>::iterator SubIt;
  typedef vector<Reduction>::iterator RdsIt;
  typedef vector<SubSymGTOs>::const_iterator cSubIt;
  typedef vector<Reduction>::const_iterator cRdsIt;
  typedef MultArray<dcomplex, 1> A1dc;
  typedef MultArray<dcomplex, 2> A2dc;
  typedef MultArray<dcomplex, 3> A3dc;
  typedef MultArray<dcomplex, 4> A4dc;


  // ==== Slow routines ====
  dcomplex SMatEle(CartGTO& a, CartGTO& b);
  dcomplex TMatEle(CartGTO& a, CartGTO& b);
  dcomplex VMatEle(CartGTO& a, Eigen::Vector3cd at, CartGTO& b);

  // ==== SymGTOs ====
  void SymGTOs_CalcMatrix(SymGTOs& a, SymGTOs& b, bool calc_coulomb, BMatSet* res);
  //  void SymGTOs_CalcERI(SymGTOs& i, SymGTOs& j, SymGTOs& k, SymGTOs& l, IB2EInt* eri);

}

#endif
