#include <iostream>
#include <gtest/gtest.h>
#include <Eigen/Core>

#include "../src_cpp/gtest_plus.hpp"
#include "../src_cpp/eigen_plus.hpp"
#include "erfc.hpp"
#include "r1_lc.hpp"
#include "r1basis.hpp"

using namespace std;
using namespace Eigen;
using namespace cbasis;

TEST(Erfc, real_erfc) {

  using namespace erfc_mori;

  // this function is forbidden
  //  erfc_add_Eh_q<int>(1, 2);
  
  double y, x, expect;
  ErfcCalcData calc_data;

  x = 1.0;
  expect =0.157299207050285130658779364917390740703933002034;
  //erfc_d(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);
  
  x = 1.0;
  x = 1.0 / 100.0;
  expect = 0.98871658444415038308409047645193078905089904517;
  //  erfc_d(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);

  x = 3.0;
  expect = 0.0000220904969985854413727761295823203798477070873992;
  //  erfc_d(x, y, calc_data);
  // erfc(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);

}
TEST(Erfc, complex_erfc) {

  using namespace erfc_mori;

  dcomplex y;
  ErfcCalcData calc_data;
  double eps = 10.0 * machine_eps();

  dcomplex x(1, -1);
  dcomplex y_expect( -0.31615128169794764488027108024367,
	       +0.190453469237834686284108861969162);

  Erfc(x, y, calc_data);

  EXPECT_DOUBLE_EQ(y.real(), y_expect.real());
  EXPECT_NEAR(y.real(), y_expect.real(), eps);
  EXPECT_DOUBLE_EQ(y.imag(), y_expect.imag());
  EXPECT_NEAR(y.imag(), y_expect.imag(), eps);

  x = dcomplex(0.0157073173118206757532953533099,
	 0.9998766324816605986389071277312);
  y_expect = dcomplex(0.95184545524179913420177473658805,
		-1.64929108965086517748934245403332);

  Erfc(x, y, calc_data);
  EXPECT_TRUE(calc_data.convergence);
  EXPECT_NEAR(     y.real(), y_expect.real(), eps);
  EXPECT_NEAR(     y.imag(), y_expect.imag(), eps);

  x = dcomplex(0.001564344650402308690101053194671668923139,
	 0.009876883405951377261900402476934372607584);
  y_expect = dcomplex(0.9982346553205423153337357292658472915601,
		-0.0111452046101524188315708507537751407281);

  Erfc(x, y, calc_data);
  EXPECT_TRUE(calc_data.convergence);
  EXPECT_NEAR(y.real(), y_expect.real(), eps);
  EXPECT_NEAR(y.imag(), y_expect.imag(), eps);

}

TEST(TestR1LC, TestSTO) {

  LC_STOs s = Create_LC_STOs();
  s->Add(1.1, 2, 1.1);
  s->Add(1.1, 3, 1.3);
  
  cout << s->str() << endl;

}
TEST(TestR1Basis, TestSTO) {

  STOs s = Create_STOs();
  s->AddPrim(2, 1.1);
  s->AddPrim(3, 1.3);
  s->SetUp();
  
  cout << s->str() << endl;

}
TEST(TestGTO, TestSTO) {

  GTOs g = Create_GTOs();
  for(int n = -5; n < 5; n++) 
    g->AddPrim(1, pow(2.5, n));
  g->SetUp();

  LC_STOs s = Create_LC_STOs();
  s->Add(1.1, 2, 0.5);
  s->Add(0.1, 2, 1.5);
  
  VectorXcd x;
  g->InitVec(x);
  g->CalcVec(s, x);
  EXPECT_EQ(10, x.size());
  cout << x << endl;
}
TEST(TestGTO, H_atom) {

  GTOs g = Create_GTOs();
  for(int n = -5; n < 5; n++) 
    g->AddPrim(1, pow(2.0, n));
  g->SetUp();
  
  MatrixXcd s, d2, v, h;
  g->InitMat(s);
  g->InitMat(h);
  g->InitMat(d2);
  g->InitMat(v);

  g->CalcRmMat(0, s);
  g->CalcD2Mat(d2);
  g->CalcRmMat(-1, v);

  h = -0.5 * d2 - v;

  double eps = pow(10.0, -2);
  SymGenComplexEigenSolver solver(h, s);
  EXPECT_C_NEAR(-0.5, solver.eigenvalues()[0], eps);

}


int main(int argc, char **args) {
  cout << "wa:" << endl;
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}

