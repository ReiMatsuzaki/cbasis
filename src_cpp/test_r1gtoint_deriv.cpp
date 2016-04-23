#include <Eigen/Dense>
#include <gtest/gtest.h>

#include "typedef.hpp"
#include "gtest_plus.hpp"
#include "r1gtoint.hpp"
#include "l_algebra.hpp"

using namespace l2func;
using namespace Eigen;

dcomplex NormalTerm(int n, dcomplex z0) {
  R1GTOs gtos(0);
  gtos.Add(n, z0);
  gtos.Normalize();
  return gtos.conts_[0].coef(0, 0);
}
MatrixXcd CalcL(dcomplex z0, dcomplex z1) {
  R1GTOs gtos(1);
  gtos.Add(2, z0);
  gtos.Add(2, z1);
  gtos.Normalize();
  gtos.CalcMat();
  dcomplex energy(0.5);

  MatrixXcd hmes00 = gtos.mat("t") + gtos.mat("v") -energy*gtos.mat("s");
  return hmes00;
}
VectorXcd CalcM(dcomplex z0, dcomplex z1, const R1STOs& driv) {
  R1GTOs gtos(1);
  gtos.Add(2, z0);
  gtos.Add(2, z1);
  gtos.Normalize();
  gtos.CalcVec(driv);

  return gtos.vec("m");

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

  R1GTOs gtos(1);
  gtos.Add(2, z0); gtos.Add(2, z1);
  gtos.Normalize();
  gtos.CalcDerivVec(driv);
  
  double eps(pow(10.0, -7.0));
  EXPECT_C_NEAR(grad, gtos.vec("d1_m")(0), eps);
  EXPECT_C_NEAR(hess, gtos.vec("d2_m")(0), eps);
  
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

  R1GTOs gtos(1);
  gtos.Add(2, z0); gtos.Add(2, z1);
  gtos.Normalize();
  gtos.CalcDerivMat(0.5);

  double eps(pow(10.0, -7.0));
  EXPECT_C_NEAR(grad, gtos.mat("d10_hmes")(0, 1), eps);
  EXPECT_C_NEAR(hess, gtos.mat("d20_hmes")(0, 1), eps);
  EXPECT_C_NEAR(dcomplex(0.8990314915292628, -1.6007590026235783),
    gtos.mat("d11_hmes")(0, 1), eps);

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

TEST(OptAlpha, Grad) {

 }

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
