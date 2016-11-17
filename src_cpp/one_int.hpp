#ifndef ONE_INT_H
#define ONE_INT_H

#include <Eigen/Core>
#include "../utils/typedef.hpp"
#include "mol_func.hpp"
#include "symmolint.hpp"

namespace cbasis {

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
  dcomplex DXMatEle(CartGTO& a, CartGTO& b);
  dcomplex DYMatEle(CartGTO& a, CartGTO& b);
  dcomplex DZMatEle(CartGTO& a, CartGTO& b);
  dcomplex PWVecEle(const Eigen::Vector3cd& k, CartGTO& a);

  // ==== SymGTOs ====
  BMatSet CalcMat(SymGTOs a, SymGTOs b, bool calc_coulomb);
  BMatSet CalcMat_Complex(SymGTOs g, bool calc_coulomb);
  BMatSet CalcMat_Hermite(SymGTOs g, bool calc_coulomb);
  //  BMatSet CalcMat_V(SymGTOs a, SymGTOs b, Eigen::Vector3cd xyz, dcomplex q);

  // ==== SymGTOs(new) ====
  void InitBVec(SymGTOs a, BVec *ptr_bvec);
  void InitBMat(SymGTOs a, Irrep krrep, SymGTOs b, BMat *ptr_mat);
  void CalcMat(SymGTOs a, SymGTOs b, 
	       BMat *S, BMat *T, BMat *V, BMat *X, BMat *Y, BMat *Z);
  void CalcPWVec(SymGTOs a, const Eigen::Vector3cd& k,
		 BVec *S, BVec *X, BVec *Y, BVec *Z);
  
  
}

#endif
