#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
#include <iomanip>
#include "macros.hpp"
#include "l_algebra.hpp"
#include "eigen_plus.hpp"
#include "opt_alpha.hpp"


namespace l2func {
  
  using namespace Eigen;
  using namespace std;

  void AlphaGrad(const VectorXcd& d, const VectorXcd& c,
   		 const MatrixXcd& D10,
		 const VectorXcd& r1,
		 const VectorXcd& s0, const VectorXcd& s1,
		 dcomplex *a, VectorXcd *g) {

    // alpha
    *a = (d.array() * s0.array()).sum();

    // alpha
    *a = (d.array() * s0.array()).sum();

    // grad
    VectorXcd g1 = Calc_a_Aj_b(d, D10, c);
    VectorXcd g2 = Calc_ai_b(r1, c);
    VectorXcd g3 = Calc_ai_b(s1, d);
      
    *g = g2 + g3 - g1;
  }
  
  void AlphaGradHess(const VectorXcd& c, const VectorXcd& d, const MatrixXcd& U,
		     const MatrixXcd& D00, const MatrixXcd& D10,
		     const MatrixXcd& D20, const MatrixXcd& D11,
		     const VectorXcd& r0,  const VectorXcd& r1,  const VectorXcd& r2,
		     const VectorXcd& s0,  const VectorXcd& s1,  const VectorXcd& s2,
		     dcomplex* a, VectorXcd* g, MatrixXcd* h) {

    // alpha
    *a = (d.array() * s0.array()).sum();

    // grad
    VectorXcd g1 = Calc_a_Aj_b(d, D10, c);
    VectorXcd g2 = Calc_ai_b(r1, c);
    VectorXcd g3 = Calc_ai_b(s1, d);
      
    *g = g2 + g3 - g1;

    // hess
    MatrixXcd rij_c       = (r2.array() * c.array()).matrix().asDiagonal();
    MatrixXcd d_sij       = (s2.array() * d.array()).matrix().asDiagonal();
    MatrixXcd d_Lij_c     = Calc_a_Aij_b(d, D20, D11, c);

    MatrixXcd tmp = Calc_ai_A_Bj_b(r1, U, D10, c);
    MatrixXcd ri_U_Lj_c   = tmp + tmp.transpose();    

    tmp = (r1 * s1.transpose()).array() * U.array();
    MatrixXcd ri_U_sj = tmp + tmp.transpose();

    tmp = Calc_a_Ai_B_Aj_b(d, D10, U, c);
    MatrixXcd d_Li_U_Lj_c = tmp + tmp.transpose();

    tmp = Calc_ai_A_Bj_b(s1, U, D10, d);
    MatrixXcd d_Li_U_sj   = tmp + tmp.transpose();

    *h = rij_c - d_Lij_c + d_sij
      - ri_U_Lj_c + d_Li_U_Lj_c - ri_U_Lj_c + ri_U_sj;

  }



  void VectorShift(VectorXcd& zs_in, VectorXcd& zs_out,
		   const VectorXi& opt_idx, dcomplex z_shift) {

    int num(zs_out.size());
    for(int i = 0; i < num; i++)
      zs_out[i] = zs_in[i];
    for(int i = 0; i < opt_idx.size(); i++) {
      zs_out[opt_idx[i]] += z_shift;
    }
  }

  VectorXcd SolveAlpha(const R1STOs& driv,  R1GTOs& gtos, int L, dcomplex ene,
		       MatVecMap& d, double eps) {

    gtos.CalcMatSTV(L, d, "s", "t", "v");
    gtos.CalcVec(driv, d, "m");

    VectorXcd& m = d.vec("m");
    MatrixXcd& t = d.mat("t");
    MatrixXcd& v = d.mat("v");
    MatrixXcd& s = d.mat("s");

    if(eps > pow(10.0, -14.0)) {
      MatrixXcd X;
      CanonicalMatrix(s, eps, &X);

      MatrixXcd Lp = (X.transpose() *
		      (t+v-ene*s) *
		      X);
      VectorXcd mp = X.transpose() * m;
      return X*Lp.fullPivLu().solve(mp);
    } else {
      return (t+v-ene*s).fullPivLu().solve(m);
    }
  }

  void OptAlphaShiftFull(const R1STOs& driv,
			 const VectorXi& opt_idx,
			 int L,
			 dcomplex ene,
			 double h,
			 int max_iter,
			 double eps,
			 R1GTOs& gtos,
			 dcomplex& z_shift,
			 bool* convq,
			 dcomplex* alpha) {

    if(opt_idx.array().minCoeff() < 0 ||
       opt_idx.array().maxCoeff() >= gtos.size_prim()) {
      string msg; SUB_LOCATION(msg);
      ostringstream oss;
      oss << msg << ": invalid opt_idx" << endl;
      oss << "min(opt_idx): " << opt_idx.array().minCoeff() << endl;
      oss << "max(opt_idx): " << opt_idx.array().maxCoeff() << endl;
      oss << "size_prim   : " << gtos.size_prim() << endl;
      throw runtime_error(oss.str());
    }

    VectorXcd zs0(gtos.size_prim());
    VectorXcd zs1(gtos.size_prim());
    dcomplex ih(0.0, h);
    dcomplex ii(0.0, 1.0);
    dcomplex a00;

    for(int i = 0; i < gtos.size_prim(); i++) {
      zs0[i] = gtos.prim(i).z;
    }

    int i(0);
    for(i = 0; i < max_iter; i++) {
      dcomplex ap0, am0, a0p, a0m;
      VectorShift(zs0, zs1, opt_idx, z_shift + h);
      gtos.Set(2, zs1); gtos.Normalize();
      ap0 = CalcAlphaFull(driv, gtos, L, ene);

      VectorShift(zs0, zs1, opt_idx, z_shift - h);
      gtos.Set(2, zs1); gtos.Normalize();
      am0 = CalcAlphaFull(driv, gtos, L, ene);

      VectorShift(zs0, zs1, opt_idx, z_shift + ih);
      gtos.Set(2, zs1); gtos.Normalize();
      a0p = CalcAlphaFull(driv, gtos, L, ene);

      VectorShift(zs0, zs1, opt_idx, z_shift - ih);
      gtos.Set(2, zs1);       gtos.Normalize();
      a0m = CalcAlphaFull(driv, gtos, L, ene);
      
      dcomplex grad = (ap0 - am0 + ii*a0m - ii*a0p) / (4.0 * h);
      dcomplex hess = (ap0 + am0 - a0p - a0m) / (2.0*h*h);

      dcomplex dz(-grad/hess);
      z_shift += dz;
      if(abs(grad) < eps && abs(dz) < eps)
	break;
    }

    *convq = (i < max_iter-2);
    cout << i << endl;

    VectorShift(zs0, zs1, opt_idx, z_shift);
    gtos.Set(2, zs1); gtos.Normalize();
    a00 = CalcAlphaFull(driv, gtos, L, ene);

    *alpha = a00;

  }
/*
  void OptAlphaShiftCanonical(const R1STOs& driv,
			      const VectorXi& opt_idx,
			      dcomplex ene,
			      double h,
			      int max_iter,
			      double eps,
			      double eps_canonical,
			      R1GTOs& gtos,
			      dcomplex& z_shift,
			      bool* convq,
			      dcomplex* alpha) {

    VectorXcd zs0(gtos.size_basis());
    VectorXcd zs1(gtos.size_basis());
    dcomplex ih(0.0, h);
    dcomplex ii(0.0, 1.0);
    
    MatrixXcd L;

    for(int i = 0; i < gtos.size_basis(); i++) {
      zs0[i] = gtos.prim(i).z;
    }

    int i(0);
    for(i = 0; i < max_iter; i++) {
      dcomplex ap0, am0, a0p, a0m;
      VectorShift(zs0, zs1, opt_idx, z_shift + h);
      gtos.Set(2, zs1);
      ap0 = CalcAlpha(driv, gtos, ene, eps_canonical);

      VectorShift(zs0, zs1, opt_idx, z_shift - h);
      gtos.Set(2, zs1);
      am0 = CalcAlpha(driv, gtos, ene, eps_canonical);

      VectorShift(zs0, zs1, opt_idx, z_shift + ih);
      gtos.Set(2, zs1);
      a0p = CalcAlpha(driv, gtos, ene, eps_canonical);	

      VectorShift(zs0, zs1, opt_idx, z_shift - ih);
      gtos.Set(2, zs1);      
      a0m = CalcAlpha(driv, gtos, ene, eps_canonical);	
      
      dcomplex grad = (ap0 - am0 + ii*a0m - ii*a0p) / (4.0 * h);
      dcomplex hess = (ap0 + am0 - a0p - a0m) / (2.0*h*h);

      dcomplex dz(-grad/hess);
      z_shift += dz;
      if(abs(grad) < eps && abs(dz) < eps)
	break;
    }

    cout << i << endl;
    *convq = (i < max_iter-2);
    VectorShift(zs0, zs1, opt_idx, z_shift);
    gtos.Set(2, zs1);
    *alpha = CalcAlpha(driv, gtos, ene, 0.0);;

  }
*/      
  void OptAlphaShift(const R1STOs& driv,
		     const VectorXi& opt_idx,
		     int L,
		     dcomplex ene,
		     double h,
		     int max_iter,
		     double eps,
		     double eps_canonical,
		     R1GTOs& gtos,
		     dcomplex& z_shift,
		     bool* convq,
		     dcomplex* alpha) {

    double eps_canonical_th(pow(10.0, -14.0));
    if(eps_canonical < eps_canonical_th)
      OptAlphaShiftFull(driv, opt_idx, L, ene, h, max_iter, eps,
			gtos, z_shift, convq, alpha);
    else {
      string msg; SUB_LOCATION(msg);
      msg += ": canonical calculation is not supported";
      throw runtime_error(msg);
      //OptAlphaShiftCanonical(driv, opt_idx, ene, h, max_iter, eps, eps_canonical,
      //gtos, z_shift, convq, alpha);
		    }
  }
  VectorXcd SolveAlpha(const R1STOs& driv, R1GTOs& gtos,
		       int L, dcomplex ene, const R1STOs& sto_pot) {

    static MatVecMap d;
    
    gtos.CalcMatSTV(L, d, "s", "t", "v");
    gtos.CalcMatSTO(sto_pot, d, "v2");
    gtos.CalcVec(driv, d, "m");

    VectorXcd& m = d.vec("m");
    MatrixXcd& t = d.mat("t");
    MatrixXcd& v = d.mat("v");
    MatrixXcd& v2= d.mat("v2");
    MatrixXcd& s = d.mat("s");
    return (t+v+v2-ene*s).fullPivLu().solve(m);

  }
  dcomplex CalcAlpha(const R1STOs& driv, R1GTOs& gtos,
		     int L, dcomplex ene, const R1STOs& sto_pot) {

    static MatVecMap d;

    gtos.CalcMatSTV(L, d, "s", "t", "v");
    gtos.CalcMatSTO(sto_pot, d, "v2");
    gtos.CalcVec(driv, d, "m");

    VectorXcd& m = d.vec("m");
    MatrixXcd& t = d.mat("t");
    MatrixXcd& v = d.mat("v");
    MatrixXcd& v2= d.mat("v2");
    MatrixXcd& s = d.mat("s");
    return (m.array() * 
	    ((t+v+v2-ene*s).fullPivLu().solve(m)).array()).sum();
  }

  dcomplex CalcAlpha(const R1STOs& driv, R1GTOs& gtos,
		     int L, dcomplex ene, double eps) {

    
    static MatVecMap mat_vec;
    
    gtos.CalcMatSTV(L, mat_vec, "s", "t", "v");
    gtos.CalcVec(driv, mat_vec, "m");
    MatrixXcd X;
    CanonicalMatrix(mat_vec.mat("s"), eps, &X);
    MatrixXcd Lp = (X.transpose() *
		    (mat_vec.mat("t") + mat_vec.mat("v") - ene * mat_vec.mat("s")) *
		    X);
    VectorXcd m = X.transpose()*mat_vec.vec("m");
    dcomplex a = (m.array() *
		  Lp.fullPivLu().solve(m).array()).sum();
  
    return a;

  }  
  dcomplex CalcAlphaFull(const R1STOs& driv,R1GTOs& gtos, int L, dcomplex ene) {

    static MatVecMap mat_vec;
    
    gtos.CalcMatSTV(L, mat_vec, "s", "t", "v");
    gtos.CalcVec(driv, mat_vec, "m");

    VectorXcd& m = mat_vec.vec("m");
    MatrixXcd& t = mat_vec.mat("t");
    MatrixXcd& v = mat_vec.mat("v");
    MatrixXcd& s = mat_vec.mat("s");
    dcomplex a = (m.array() * 
		  ((t+v-ene*s).fullPivLu().solve(m)).array()).sum();
    return a;
  }

}
