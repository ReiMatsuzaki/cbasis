#include <iostream>
#include <Eigen/Dense>
#include <iomanip>
#include "eigen_plus.hpp"
#include "opt_alpha.hpp"


namespace l2func {
  
  using namespace Eigen;
  using namespace std;
  void VectorShift(VectorXcd& zs_in, VectorXcd& zs_out,
		   const VectorXi& opt_idx, dcomplex z_shift) {

    int num(zs_out.size());
    for(int i = 0; i < num; i++)
      zs_out[i] = zs_in[i];
    for(int i = 0; i < opt_idx.size(); i++) {
      zs_out[opt_idx[i]] += z_shift;
    }
  }

  dcomplex CalcAlpha(const R1STOs& driv,
		     R1GTOs& gtos, dcomplex ene,
		     MatMap& mat, VecMap& vec, double eps) {
    
    gtos.CalcMat(&mat);
    gtos.CalcVecSTO(driv, &vec);
    
    MatrixXcd X;
    CanonicalMatrix(mat["s"], eps, &X);
    MatrixXcd Lp = X.transpose()*(mat["t"] + mat["v"] - ene * mat["s"])*X;
    VectorXcd m = X.transpose()*vec["m"];
    dcomplex a = (m.array() *
		  Lp.fullPivLu().solve(m).array()).sum();
  
    return a;
  }
  dcomplex CalcAlpha(const R1STOs& driv, R1GTOs& gtos,
		     dcomplex ene, double eps) {
    MatMap mat;
    VecMap vec;
    return CalcAlpha(driv, gtos, ene, mat, vec, eps);
  }  
  dcomplex CalcAlphaFull(const R1STOs& driv,R1GTOs& gtos,
			 dcomplex ene,
			 MatMap& mat, VecMap& vec) {
    
    gtos.CalcMat(&mat);
    gtos.CalcVecSTO(driv, &vec);

    dcomplex a = (vec["m"].array() *
		  ((mat["t"] + mat["v"] - ene * mat["s"]).fullPivLu().solve(vec["m"])).array()).sum();
    return a;
  }

  VectorXcd SolveAlpha(const R1STOs& driv,  R1GTOs& gtos, dcomplex ene, double eps) {

    MatMap mat;
    VecMap vec;
    gtos.CalcMat(&mat);
    gtos.CalcVecSTO(driv, &vec);

    if(eps > pow(10.0, -14.0)) {
      MatrixXcd X;
      CanonicalMatrix(mat["s"], eps, &X);

      MatrixXcd Lp = (X.transpose() *
		      (mat["t"] + mat["v"] - ene * mat["s"]) *
		      X);
      VectorXcd m = X.transpose() * vec["m"];
      return X*Lp.fullPivLu().solve(m);
    } else {
      return (mat["t"] + mat["v"] - ene * mat["s"]).fullPivLu().solve(vec["m"]);
    }
  }

  void OptAlphaShiftFull(const R1STOs& driv,
			 const VectorXi& opt_idx,
			 dcomplex ene,
			 double h,
			 int max_iter,
			 double eps,
			 R1GTOs& gtos,
			 dcomplex& z_shift,
			 bool* convq,
			 dcomplex* alpha) {

    MatMap mat;
    VecMap vec;
    //    MatrixXcd L;
    VectorXcd zs0(gtos.size_basis());
    VectorXcd zs1(gtos.size_basis());
    dcomplex ih(0.0, h);
    dcomplex ii(0.0, 1.0);
    dcomplex a00;

    for(int i = 0; i < gtos.size_basis(); i++) {
      zs0[i] = gtos.basis(i).z;
    }

    int i(0);
    for(i = 0; i < max_iter; i++) {
      dcomplex ap0, am0, a0p, a0m;
      VectorShift(zs0, zs1, opt_idx, z_shift + h);
      gtos.Set(2, zs1); gtos.Normalize();
      ap0 = CalcAlphaFull(driv, gtos, ene, mat, vec);      

      VectorShift(zs0, zs1, opt_idx, z_shift - h);
      gtos.Set(2, zs1); gtos.Normalize();
      am0 = CalcAlphaFull(driv, gtos, ene, mat, vec);	

      VectorShift(zs0, zs1, opt_idx, z_shift + ih);
      gtos.Set(2, zs1); gtos.Normalize();
      a0p = CalcAlphaFull(driv, gtos, ene, mat, vec);      

      VectorShift(zs0, zs1, opt_idx, z_shift - ih);
      gtos.Set(2, zs1);       gtos.Normalize();
      a0m = CalcAlphaFull(driv, gtos, ene, mat, vec);      

      /*
      VectorShift(zs0, zs1, opt_idx, z_shift);
      gtos.Set(2, zs1); gtos.Normalize();
      a00 = CalcAlpha(driv, gtos, mat, vec, 0.0);
      */
      
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
    a00 = CalcAlphaFull(driv, gtos, ene, mat, vec);

    *alpha = a00;

  }
      
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

    MatMap mat;
    VecMap vec;

    VectorXcd zs0(gtos.size_basis());
    VectorXcd zs1(gtos.size_basis());
    dcomplex ih(0.0, h);
    dcomplex ii(0.0, 1.0);
    
    MatrixXcd L;

    for(int i = 0; i < gtos.size_basis(); i++) {
      zs0[i] = gtos.basis(i).z;
    }

    int i(0);
    for(i = 0; i < max_iter; i++) {
      dcomplex ap0, am0, a0p, a0m;
      VectorShift(zs0, zs1, opt_idx, z_shift + h);
      gtos.Set(2, zs1);
      ap0 = CalcAlpha(driv, gtos, ene, mat, vec, eps_canonical);

      VectorShift(zs0, zs1, opt_idx, z_shift - h);
      gtos.Set(2, zs1);
      am0 = CalcAlpha(driv, gtos, ene, mat, vec, eps_canonical);

      VectorShift(zs0, zs1, opt_idx, z_shift + ih);
      gtos.Set(2, zs1);
      a0p = CalcAlpha(driv, gtos, ene, mat, vec, eps_canonical);	

      VectorShift(zs0, zs1, opt_idx, z_shift - ih);
      gtos.Set(2, zs1);      
      a0m = CalcAlpha(driv, gtos, ene, mat, vec, eps_canonical);	
      
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
    *alpha = CalcAlpha(driv, gtos, ene, mat, vec, 0.0);;


  }

  void OptAlphaShift(const R1STOs& driv,
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

    double eps_canonical_th(pow(10.0, -14.0));
    if(eps_canonical < eps_canonical_th)
      OptAlphaShiftFull(driv, opt_idx, ene, h, max_iter, eps,
			gtos, z_shift, convq, alpha);
    else
      OptAlphaShiftCanonical(driv, opt_idx, ene, h, max_iter, eps, eps_canonical,
			     gtos, z_shift, convq, alpha);
  }

}
