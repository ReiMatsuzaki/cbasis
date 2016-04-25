#include <Eigen/Dense>
#include <gtest/gtest.h>

#include "typedef.hpp"
#include "eigen_plus.hpp"
#include "gtest_plus.hpp"

#include "r1gtoint.hpp"
#include "op_driv.hpp"
#include "opt_alpha.hpp"
#include "l_algebra.hpp"

using namespace std;
using namespace l2func;
using namespace Eigen;

class OptTargetCos :public IOptTarget {
public:
  OptTargetCos() {}
  ~OptTargetCos() {}
  int dim() { return 1; }
  void ValGradHess(const VectorXcd& zs, dcomplex* a, VectorXcd& g, MatrixXcd& h) {
    g = VectorXcd::Zero(1);
    h = MatrixXcd::Zero(1, 1);
    *a      = +cos(zs[0]);
    g(0)    = -sin(zs[0]);
    h(0, 0) = -cos(zs[0]);
  }
};

dcomplex NormalTerm(int n, dcomplex z0) {
  R1GTOs gtos;
  gtos.Add(n, z0);
  gtos.Normalize();
  return gtos.conts_[0].coef(0, 0);
}
MatrixXcd CalcL(dcomplex z0, dcomplex z1) {
  R1GTOs gtos;
  gtos.Add(2, z0);
  gtos.Add(2, z1);
  gtos.Normalize();

  MatVecMap res;
  gtos.CalcMatSTV(1, res, "s", "t", "v");
  dcomplex energy(0.5);

  MatrixXcd hmes00 = res.mat("t") + res.mat("v") -energy*res.mat("s");
  return hmes00;
}
VectorXcd CalcM(dcomplex z0, dcomplex z1, const R1STOs& driv) {
  R1GTOs gtos;
  gtos.Add(2, z0);
  gtos.Add(2, z1);
  gtos.Normalize();

  MatVecMap res;
  gtos.CalcVec(driv, res, "m");

  return res.vec("m");

}

TEST(DerivBasis, NormalTerm) {

  dcomplex z(1.1, 0.2);
  dcomplex h(0.0001, 0.0);
  dcomplex hi(0.0, 0.0001);
  dcomplex ii(0.0, 1.0);
  int pn(2);

  R1STOs driv; driv.Add(1.0, 2, 1.0);

  dcomplex vp0, vm0, v0p, v0m, grad, hess;  
  vp0 = NormalTerm(pn, z+h);
  vm0 = NormalTerm(pn, z-h);
  v0p = NormalTerm(pn, z+hi);
  v0m = NormalTerm(pn, z-hi);

  grad = (vp0 - vm0 - ii*v0p +ii*v0m)/(4.0*h);
  hess = (vp0 + vm0 - v0p - v0m)/(2.0*h*h);

  dcomplex grad_calc = NPrimeGTO(      NormalTerm(pn, z), pn, z);
  dcomplex hess_calc = NDoublePrimeGTO(NormalTerm(pn, z), pn, z);
  EXPECT_C_NEAR(grad, grad_calc, pow(10.0, -8.0));
  EXPECT_C_NEAR(hess, hess_calc, pow(10.0, -8.0));

}
TEST(DerivBasis, DrivVec) {
  dcomplex z0(1.1, 0.2);
  dcomplex z1(1.3, 1.2);
  dcomplex h(0.0001, 0.0);
  dcomplex hi(0.0, 0.0001);
  dcomplex ii(0.0, 1.0);

  R1STOs driv; driv.Add(1.0, 2, 1.0);

  dcomplex vp0, vm0, v0p, v0m, grad, hess;  
  vp0 = CalcM(z0+h, z1, driv)[0];
  vm0 = CalcM(z0-h, z1, driv)[0];
  v0p = CalcM(z0+hi,z1, driv)[0];
  v0m = CalcM(z0-hi,z1, driv)[0];
  grad = (vp0 - vm0 - ii*v0p +ii*v0m)/(4.0*h);
  hess = (vp0 + vm0 - v0p - v0m)/(2.0*h*h);

  R1GTOs gtos;
  gtos.Add(2, z0); gtos.Add(2, z1);
  gtos.Normalize();

  R1GTOs d1_gtos;
  d1_gtos.SetOneDeriv(gtos);

  R1GTOs d2_gtos;
  d2_gtos.SetTwoDeriv(gtos);
  
  MatVecMap res;
  d1_gtos.CalcVec(driv, res, "d1_m");
  d2_gtos.CalcVec(driv, res, "d2_m");
  
  double eps(pow(10.0, -7.0));
  EXPECT_C_NEAR(grad, res.vec("d1_m")(0), eps);
  EXPECT_C_NEAR(hess, res.vec("d2_m")(0), eps);
  
}
TEST(DerivBasis, HMES) {
  dcomplex z0(1.1, 0.2);
  dcomplex z1(1.3, 1.2);
  dcomplex h(0.0001, 0.0);
  dcomplex hi(0.0, 0.0001);
  dcomplex ii(0.0, 1.0);
  
  dcomplex vp0, vm0, v0p, v0m, grad, hess;

  // -- for 01 element --
  vp0 = CalcL(z0+h,  z1)(0, 1);
  vm0 = CalcL(z0-h,  z1)(0, 1);
  v0p = CalcL(z0+hi, z1)(0, 1);
  v0m = CalcL(z0-hi, z1)(0, 1);

  grad = (vp0 - vm0 - ii*v0p +ii*v0m)/(4.0*h);
  hess = (vp0 + vm0 -    v0p -   v0m)/(2.0*h*h);

  R1GTOs gtos, d1_gs, d2_gs;
  gtos.Add(2, z0); gtos.Add(2, z1);
  gtos.Normalize();
  d1_gs.SetOneDeriv(gtos);
  d2_gs.SetTwoDeriv(gtos);
  
  MatVecMap res;
  d1_gs.CalcMatSTV(gtos,  1, res, "d10_s", "d10_t", "d10_v");
  d1_gs.CalcMatSTV(d1_gs, 1, res, "d11_s", "d11_t", "d11_v");
  d2_gs.CalcMatSTV(gtos, 1, res, "d20_s", "d20_t", "d20_v");
  //  gtos.CalcDerivMat(0.5);

  double eps(pow(10.0, -7.0));
  double energy(0.5);
  MatrixXcd d10_hmes = res.mat("d10_t") + res.mat("d10_v") - energy*res.mat("d10_s");
  MatrixXcd d11_hmes = res.mat("d11_t") + res.mat("d11_v") - energy*res.mat("d11_s");
  MatrixXcd d20_hmes = res.mat("d20_t") + res.mat("d20_v") - energy*res.mat("d20_s");
  EXPECT_C_NEAR(grad, d10_hmes(0, 1), eps);
  EXPECT_C_NEAR(hess, d20_hmes(0, 1), eps);
  EXPECT_C_NEAR(dcomplex(0.8990314915292628, -1.6007590026235783),
    d11_hmes(0, 1), eps);

}

TEST(LinearAlgebra, a_Aj_b) {

  // A^0 = (2.0, 2.0, 
  //        2.0, 0.0 )
  // A^1 = (0.0,   1.5)
  //       (1.5,   0.2)
  // aA^0b = (0.1,0.2) (1.6, 0.6) = 0.16 + 0.12 = 0.28
  // aA^1b = (0.1,0.2) (0.75, 0.45+0.1) = 0.075+0.11=0.185

  MatrixXcd A10(2,2); A10 << 1.0, 2.0, 1.5, 0.1;
  VectorXcd a(2), b(2);  a << 0.1, 0.2;  b << 0.3, 0.5;
  VectorXcd res(2);
  res = Calc_a_Aj_b(a, A10, b);

  EXPECT_C_EQ(0.28, res(0));
  EXPECT_C_EQ(0.185, res(1));
}
TEST(LinearAlgebra, a_A01_b) {

  MatrixXcd A10(2,2); A10 << 1.0, 2.0, 1.5, 0.1;
  VectorXcd a(2), b(2);  a << 0.1, 0.2;  b << 0.3, 0.5;

  VectorXcd res(2);
  VectorXcd res0(2), res1(2);

  res = Calc_a_Aj_b(a, A10, b);
  Calc_a_A10_b(a, A10, b, &res0);
  Calc_a_A10_b(b, A10, a, &res1);

  EXPECT_C_EQ(res(0), res1(0) + res0(0));
  EXPECT_C_EQ(res(1), res1(1) + res0(1));
}
TEST(LinearAlgebra, a_b1) {

  VectorXcd a(2); a << 1.0, 2.0;
  VectorXcd b1(2); b1 << 3.0, -2.0;
  VectorXcd res(2);
  Calc_a_b1(a, b1, &res);
  EXPECT_C_EQ(3.0, res(0));
  EXPECT_C_EQ(-4.0, res(1));

}
TEST(LinearAlgebra, sym) {

  MatrixXcd A10(2,2); A10 << 1.0, 2.0, 1.5, 0.1;
  MatrixXcd A20(2,2); A20 << 1.0, 2.0, 1.5, 0.1;
  MatrixXcd A11(2,2); A11 << 1.0, 2.0, 2.0, 0.1;
  VectorXcd a(2), b(2);  a << 0.1, 0.2;  b << 0.3, 0.5;

  MatrixXcd res = Calc_a_Aij_b(a, A20, A11, b);

  EXPECT_C_EQ(0.0, (res - res.transpose())(0, 1));

  res = Calc_a_Aij_b(a, A20, A11, b);
  EXPECT_C_EQ(0.0, (res - res.transpose())(0, 1));

  res = Calc_a_Ai_B_Aj_b(a, A10, A11, b) + Calc_a_Ai_B_Aj_b(b, A10, A11, a);
  EXPECT_C_EQ(0.0, (res - res.transpose())(0, 1));
  

}

dcomplex CalcAlpha(dcomplex z_shift) {

  R1GTOs gtos;
  VectorXcd zs(15);  
  zs << 0.463925,
    1.202518,
    3.379649,
    10.6072,
    38.65163,
    173.5822,
    1170.498,
    0.16934112166516593 + z_shift,
    0.08989389391311804 + z_shift,
    0.055610873913491725 + z_shift,
    0.03776599632952126 + z_shift,
    0.02731159914174668 + z_shift,
    0.020665855224060142 + z_shift,
    0.016180602421004654 + z_shift,
    0.013011569667967734 + z_shift;
  gtos.Add(2, zs);
  gtos.Normalize();
  

  R1STOs driv;
  driv.Add(2.0, 2, 1.0);

  MatVecMap res;
  gtos.CalcMatSTV(1, res, "s", "t", "v");
  gtos.CalcVec(driv, res, "m");

  MatrixXcd L = res.mat("t") + res.mat("v") -0.5* res.mat("s");
  VectorXcd& m = res.vec("m");
  dcomplex a = (m.array() *
		(L.fullPivLu().solve(m)).array()).sum();
  //dcomplex a = vec["m"].dot(L.fullPivLu().solve(vec["m"]));
  return a;

}
TEST(Optimization, model) {
  
  
  IOptTarget* target = new OptTargetCos();
  IOptimizer* optimizer = new OptNewton(100, pow(10.0, -5.0), target, 0);

  VectorXcd z(1); z << 0.1;
  optimizer->Optimize(z);
  EXPECT_TRUE(optimizer->conv_q);
  EXPECT_C_EQ(0.0, optimizer->zs[0]);
  EXPECT_C_EQ(1.0, optimizer->val);

  delete target;
  delete optimizer;

}
TEST(OptAlphaInd, grad_hess) {
  
  R1STOs driv; driv.Add(2.0, 2, 1.0);
  IDriv* p_driv = new DrivSTO(driv);

  IOp*   p_op   = new OpCoulomb(1, 0.5);

  R1GTOs gs;
  VectorXcd zs(1); zs << dcomplex(1.1, -0.4);
  dcomplex ii(0, 1);
  VectorXcd h(1);  h << 0.0001;
  VectorXcd ih(1); ih << 0.0001*ii;
  
  gs.Add(2, zs);

  OptAlpha opt(p_driv, p_op, gs);

  dcomplex vp0, vm0, v0p, v0m;
  try {
    vp0 = opt.Val(zs+h);
  } catch(const exception& e){
    cout << "+h" << endl;
    cout << e.what();
  }
  try{
    vm0 = opt.Val(zs-h);
  } catch(const exception& e){
    cout << "-h" << endl;
    cout << e.what();
  }
  try{
    v0p = opt.Val(zs+ih);
  } catch(const exception& e){
    cout << "+ih" << endl;
    cout << e.what();
  }
  try{
    v0m = opt.Val(zs-ih);
    cout << v0m << endl;
  } catch(const exception& e){
    cout << "-ih" << endl;
    cout << e.what();
  }
  dcomplex v00;
  VectorXcd grad(1);
  MatrixXcd hess(1, 1);
  try {
    opt.ValGradHess(zs, &v00, grad, hess);
  } catch(const exception& e) {
    cout << "failed to call ValGradHess:" << e.what() << endl;
  } catch(...) {
    throw runtime_error("failed to call ValGradHess:");
  }
  
  dcomplex g0 = (vp0-vm0-ii*v0p+ii*v0m)/(4.0*h[0]);
  dcomplex h0 = (vp0+vm0-v0p-v0m)/(2.0*h[0]*h[0]);
  EXPECT_C_EQ(g0, grad[0]);
  EXPECT_C_NEAR(h0, hess(0, 0), pow(10.0, -6.0));

  delete p_driv;
  delete p_op;
}
TEST(OptAlphaInd, two) {

  R1STOs driv; driv.Add(2.0, 2, 1.0);
  IDriv* p_driv = new DrivSTO(driv);

  IOp*   p_op   = new OpCoulomb(1, 0.5);

  R1GTOs gs;
  VectorXcd zs(2); zs << dcomplex(1.1, -0.4), dcomplex(0.3, 0.1);
  gs.Add(2, zs);

  IOptTarget* opt = new OptAlpha(p_driv, p_op, gs);

  dcomplex a;
  VectorXcd g(2);
  MatrixXcd h(2, 2);
  try {
    opt->ValGradHess(zs, &a, g, h);
  } catch(const exception& e) {
    cout << "failed to calc\n";
    cout << e.what() << endl;
  }

  delete p_driv;
  delete p_op;
  delete opt;
}
TEST(OptAlphaInd, val) {

  R1GTOs gtos;
  VectorXcd zs(15);  
  dcomplex zshift(0, -0.001);
  zs << 0.463925,
    1.202518,
    3.379649,
    10.6072,
    38.65163,
    173.5822,
    1170.498,
    0.16934112166516593+zshift,
    0.08989389391311804+zshift,
    0.05561087391349172+zshift,
    0.03776599632952126+zshift,
    0.02731159914174668+zshift,
    0.020665855224060142+zshift,
    0.016180602421004654+zshift,
    0.013011569667967734+zshift;
  gtos.Add(2, zs);
  gtos.Normalize();

  R1STOs driv; driv.Add(2.0, 2, 1.0);
  IDriv* p_driv = new DrivSTO(driv);
  IOp*   p_op   = new OpCoulomb(1, 0.5);
  OptAlpha opt(p_driv, p_op, gtos);
  
  dcomplex v0 = CalcAlpha(zshift);
  dcomplex v1 = opt.Val(zs);

  EXPECT_C_NEAR(v0, v1, pow(10.0, -8.0));
  delete p_driv;
  delete p_op;

}
TEST(OptAlphaInd, opt_one) {
  
  R1STOs driv; driv.Add(2.0, 2, 1.0);
  IDriv* p_driv = new DrivSTO(driv);

  IOp*   p_op   = new OpCoulomb(1, 0.5);

  R1GTOs gs;
  VectorXcd zs(1);
  zs << dcomplex(0.2, -0.1);
  gs.Add(2, zs);
  gs.Normalize();
  
  IOptTarget* target = new OptAlpha(p_driv, p_op, gs);
  IOptimizer* optimizer = new OptNewton(100, pow(10.0, -5.0), target);
  try {
    optimizer->Optimize(zs);
  } catch(const exception& e) {
  cout << "exception in optimization:" << endl << e.what() << endl;
  }

  EXPECT_TRUE(optimizer->conv_q);
  cout << optimizer->zs << endl;
  cout << optimizer->val << endl;

  delete p_driv;
  delete p_op;
  delete target;
  delete optimizer;
}
TEST(OptAlphaInd, opt_two) {
  
  R1STOs driv; driv.Add(2.0, 2, 1.0);
  IDriv* p_driv = new DrivSTO(driv);
  IOp*   p_op   = new OpCoulomb(1, 0.5);

  R1GTOs gs;
  VectorXcd zs(2);
  zs << dcomplex(0.2, -0.1), dcomplex(0.8, -0.5);
  gs.Add(2, zs);
  gs.Normalize();
  
  IOptTarget* target = new OptAlpha(p_driv, p_op, gs);
  IOptimizer* optimizer = new OptNewton(100, pow(10.0, -5.0), target);
  try {
    optimizer->Optimize(zs);
  } catch(const exception& e) {
  cout << "exception in optimization:" << endl << e.what() << endl;
  }

  EXPECT_TRUE(optimizer->conv_q);
  cout << optimizer->zs << endl;
  //  cout << optimizer->val << endl;

  //  EXPECT_C_NEAR(optimizer->zs[0], dcomplex(0.964095, -0.0600633), eps);
  //  EXPECT_C_NEAR(optimizer->zs[1], dcomplex(0.664185, -1.11116), eps);

  CheckOptTarget(target, zs, 0.0001, 2.0*pow(10.0, -5.0));

  delete p_driv;
  delete p_op;
  delete target;
  delete optimizer;
}
TEST(OptAlphaInd, shift) {

  R1GTOs gtos;
  VectorXcd zs(15);  
  zs << 0.463925,
    1.202518,
    3.379649,
    10.6072,
    38.65163,
    173.5822,
    1170.498,
    0.16934112166516593,
    0.08989389391311804,
    0.055610873913491725,
    0.03776599632952126,
    0.02731159914174668,
    0.020665855224060142,
    0.016180602421004654,
    0.013011569667967734;
  gtos.Add(2, zs);
  gtos.Normalize();

  R1STOs driv; driv.Add(2.0, 2, 1.0);
  IDriv* p_driv = new DrivSTO(driv);
  IOp*   p_op   = new OpCoulomb(1, 0.5);

  VectorXi index(8); index << 7, 8, 9, 10, 11, 12, 13, 14;
  IOptTarget *target = new OptAlphaShift(p_driv, p_op, gtos, index);
  IOptimizer *optimizer = new OptNewton(10, pow(10.0, -5.0), target);

  double shift0(0.02);
  CheckOptTarget(target, dcomplex(0.0, -shift0), shift0/300.0, 1.0);
  
  optimizer->Optimize(dcomplex(0.0, -shift0));
  EXPECT_TRUE(optimizer->conv_q);

  double eps(0.000003);
  dcomplex ref(-0.00293368, -0.0204361);
  EXPECT_C_NEAR(ref, optimizer->zs[0], eps*10);
  dcomplex ref_alpha(dcomplex(-5.6568937518988989, 1.0882823480377297));
  EXPECT_C_NEAR(ref_alpha, optimizer->val, pow(10.0, -9.0));  

  delete p_driv;
  delete p_op;
  delete target;
  delete optimizer;
}
TEST(OptAlphaInd, part) {

  R1GTOs gtos;
  VectorXcd zs(3);
  zs << dcomplex(0.18,-0.1), dcomplex(0.8, -0.6), dcomplex(2.9, -3.2);
  gtos.Add(2, zs);
  gtos.Normalize();

  R1STOs driv_sto; driv_sto.Add(2.0, 2, 1.0);
  IDriv* driv = new DrivSTO(driv_sto);
  IOp*   op   = new OpCoulomb(1, 0.5);

  VectorXi index(2); index << 1, 2;
  IOptTarget* target = new OptAlphaPartial(driv, op, gtos, index);
  IOptimizer* optimizer = new OptNewton(10, pow(10.0, -5.0), target, 0);

  VectorXcd z0(2); z0 << dcomplex(0.8, -0.6), dcomplex(2.9, -3.2);
  CheckOptTarget(target, z0, 0.001, pow(10.0, -5.0));

  optimizer->Optimize(z0);
  EXPECT_TRUE(optimizer->conv_q);
  EXPECT_C_EQ(optimizer->zs[0], gtos.prim(1).z);
  EXPECT_C_EQ(optimizer->zs[1], gtos.prim(2).z);
  cout << optimizer->zs << endl;

  delete driv;
  delete op;
  delete target;
  delete optimizer;
}
 
TEST(OptAlpha, WithContraction) {
  cout << "not impl" << endl;
  /*
  VectorXcd zs_h(7);  
  zs_h <<
    0.463925,
    1.202518,
    3.379649,
    10.6072,
    38.65163,
    173.5822,
    1170.498;

  VectorXcd zs_k(10);
  zs_k << 0.16934112166516593, 
    0.08989389391311804, 
    0.055610873913491725,
    0.03776599632952126, 
    0.02731159914174668, 
    0.020665855224060142,
    0.016180602421004654,
    0.013011569667967734,
    0.01068994588552496, 
    0.008938476331399546;
  
  R1GTOs gtos;
  gtos.Add(2, zs_k);
  gtos.Normalize();
  
  MatVecMap res;
  gtos.CalcMatSTV(1, res, "s", "t", "v");
  MatrixXcd xmat;
  CanonicalMatrix(res.mat("s"), pow(10.0, -5.0), &xmat);

  R1GTOs gtos2;
  gtos2.Add(2, zs_h);
  gtos2.Add(2, zs_k, xmat.transpose());
  gtos2.Normalize();

  EXPECT_EQ(10+7, gtos2.size_prim());

  R1STOs driv; driv.Add(2.0, 2, 1.0);
  VectorXi opt_idx(10);
  opt_idx << 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

  double h(0.0001);
  double eps(0.000003);
  dcomplex z_shift(0.0, -0.02);
  bool convq;
  dcomplex alpha;
  OptAlphaShift(driv, opt_idx, 1, 0.5,
		h, 20, eps, 0.0,
		gtos2, z_shift, &convq, &alpha);

  EXPECT_TRUE(convq);
  dcomplex ref_alpha(dcomplex(-5.6568937518988989, 1.0882823480377297));
  EXPECT_C_NEAR(ref_alpha, alpha, pow(10.0, -9.0));  
  */

}
TEST(OptAlpha, ShiftKaufmann) {

  // 2016/3/26
  // ARCH=fast
  // 258ms

  R1GTOs gtos;
  
  int max_iter(100);
  dcomplex h(0.0001);
  dcomplex ih(0.0, 0.0001);
  dcomplex ii(0.0, 1.0);
  double eps(0.000003);

  dcomplex z_shift(0.0, -0.02);
  
  int i(0);
  dcomplex alpha;
  for(i = 0; i < max_iter; i++) {

    alpha = CalcAlpha(z_shift);
    dcomplex ap0 = CalcAlpha(z_shift + h);
    dcomplex am0 = CalcAlpha(z_shift - h);
    dcomplex a0p = CalcAlpha(z_shift + ih);
    dcomplex a0m = CalcAlpha(z_shift - ih);

    dcomplex grad = (ap0 - am0 + ii*a0m - ii*a0p) / (4.0 * h);
    dcomplex hess = (ap0 + am0 - a0p - a0m) / (2.0*h*h);

    dcomplex dz(grad/hess);
    z_shift -= dz;

    if(abs(grad) < eps && abs(dz) < eps) {
      break;
    }
  }
  cout << i << endl;

  // -- from 2016/3/25
  // -- -0.00293368-0.0204361j 
  
  EXPECT_TRUE(i < 90) << i;
  dcomplex ref(-0.00293368, -0.0204361);
  EXPECT_C_NEAR(ref, z_shift, eps*10);
  dcomplex ref_alpha(dcomplex(-5.6568937518988989, 1.0882823480377297));
  EXPECT_C_NEAR(ref_alpha, alpha, pow(10.0, -9.0));
}
TEST(OptAlpha, CalcAlpha) {

  R1GTOs gtos;
  VectorXcd zs(15);  
  zs << 0.463925,
    1.202518,
    3.379649,
    10.6072,
    38.65163,
    173.5822,
    1170.498,
    0.16934112166516593,
    0.08989389391311804,
    0.055610873913491725,
    0.03776599632952126,
    0.02731159914174668,
    0.020665855224060142,
    0.016180602421004654,
    0.013011569667967734;
  gtos.Add(2, zs);
  gtos.Normalize();

  R1STOs driv; driv.Add(2.0, 2, 1.0);
  IDriv* p_driv = new DrivSTO(driv);
  IOp*   p_op   = new OpCoulomb(1, 0.5);

  dcomplex v0 = CalcAlpha(0.0);
  dcomplex v1 = CalcAlpha(p_driv, p_op, gtos);

  EXPECT_C_NEAR(v0, v1, pow(10.0, -8.0));

  delete p_driv;
  delete p_op;

}
TEST(OptAlpha, CalcAlphaSTO) {
  R1GTOs gtos;
  VectorXcd zs(15);  
  zs << 0.463925,
    1.202518,
    3.379649,
    10.6072,
    38.65163,
    173.5822,
    1170.498,
    0.16934112166516593 -dcomplex(0, 0.02),
    0.08989389391311804 -dcomplex(0, 0.02),
    0.055610873913491725 -dcomplex(0, 0.02),
    0.03776599632952126 -dcomplex(0, 0.02),
    0.02731159914174668 -dcomplex(0, 0.02),
    0.020665855224060142 -dcomplex(0, 0.02),
    0.016180602421004654 -dcomplex(0, 0.02),
    0.013011569667967734 -dcomplex(0, 0.02);
  gtos.Add(2, zs);
  gtos.Normalize();

  R1STOs sto; sto.Add(2.0, 2, 1.0);
  IDriv *driv = new DrivSTO(sto);

  R1STOs pot;  pot.Add(-1.0, 0, 1.0);
  IOp *op = new OpCoulombShort(1, 0.5, pot);

  // 2016/4/12
  // alpha = -4.770452, -0.235177

  EXPECT_C_NEAR(
    CalcAlpha(driv, op, gtos),
    dcomplex(-4.770452, +0.235177),
    pow(10.0, -3.0));

  delete driv;
  delete op;

}

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
