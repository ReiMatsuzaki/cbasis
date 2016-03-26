#include <iostream>
#include <gtest/gtest.h>

#include "eigen_plus.hpp"
#include "gtest_plus.hpp"
#include "r1gtoint.hpp"
#include "cip_exp.hpp"
#include "linspace.hpp"

using namespace l2func;
using namespace Eigen;
using namespace std;

//TEST(calc_gto_int, exception)
TEST(calc_gto_int, value) {

  int maxn(10);
  dcomplex a(1.1, -0.2);
  dcomplex* calc = new dcomplex[maxn+1];
  CalcGTOInt(maxn, a, calc);

  EXPECT_C_EQ(calc[0], sqrt(M_PI)/(2.0*sqrt(a)));
  for(int n = 1; n <= maxn; n++) 
    EXPECT_C_EQ(calc[n], GTO_Int(a, n)) << n;
  delete[] calc;

}

class TestR1GTOs : public ::testing::Test {
public:
  R1GTOs gtos;
  vector<CGTO> gs;
  TestR1GTOs() : gtos(0) {
    gtos.Add(2, dcomplex(1.1, 0.2));
    gtos.Add(3, dcomplex(1.3, 1.2));
    VectorXcd zs(2); zs << dcomplex(0.1, 0.2), dcomplex(0.2, 0.3);
    gtos.Add(4, zs);
    gtos.Normalize();

    gs.push_back(CGTO(1.0, 2, dcomplex(1.1, 0.2)));
    gs.push_back(CGTO(1.0, 3, dcomplex(1.3, 1.2)));
    gs.push_back(CGTO(1.0, 4, dcomplex(0.1, 0.2)));
    gs.push_back(CGTO(1.0, 4, dcomplex(0.2, 0.3)));

  }
};
TEST_F(TestR1GTOs, add) {

  EXPECT_EQ(4, gtos.size_basis());
  EXPECT_C_EQ(dcomplex(1.1, 0.2), gtos.basis(0).z);
  EXPECT_EQ(  3                 , gtos.basis(1).n);
  EXPECT_C_EQ(dcomplex(0.1, 0.2), gtos.basis(2).z);
  EXPECT_C_EQ(dcomplex(0.2, 0.3), gtos.basis(3).z);

}
TEST_F(TestR1GTOs, max_n) {
  EXPECT_EQ(4, gtos.max_n());
}
TEST_F(TestR1GTOs, normalized) {

  for(int i = 0; i < gtos.size_basis(); i++)
    EXPECT_C_EQ(gtos.basis(i).c,
		1.0/sqrt(GTO_Int(2.0*gtos.basis(i).z, 2*gtos.basis(i).n)))
      << i;

  MatMap res;
  gtos.CalcMat(&res);
  for(int i = 0; i < 4; i++)
    EXPECT_C_EQ(1.0, res["s"](i, i)) << i;
}
TEST_F(TestR1GTOs, matrix) {

  MatMap res;
  gtos.CalcMat(&res);

  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++) {
      EXPECT_C_EQ(CIP(CNormalize(gs[i]), CNormalize(gs[j])),
		  res["s"](i, j)) << i << j;
      dcomplex tref = -0.5*CIP(CNormalize(gs[i]),
			       OpD<dcomplex, double, 2>(),
			       CNormalize(gs[j]));
      EXPECT_C_EQ(tref, res["t"](i, j)) << i << j;
      dcomplex vref = -CIP(CNormalize(gs[i]),
			   OpRm<dcomplex, double>(-1),
			   CNormalize(gs[j]));
      EXPECT_C_EQ(vref, res["v"](i, j)) << i << j;
    }
}
TEST_F(TestR1GTOs, vector_sto) {

  vector<CSTO> ss;
  R1STOs stos;

  ss.push_back(  CSTO( 1.1, 2, dcomplex(0.1, 0.01)));
  stos.Add(1.1, 2, dcomplex(0.1, 0.01));

  ss.push_back(  CSTO( 0.11, 3, 0.3));
  stos.Add(0.11, 3, 0.3);

  VecMap res;
  gtos.CalcVecSTO(stos, &res);
  for(int i = 0; i < 4; i++) {
    dcomplex ref = CIP(CNormalize(gs[i]),
		       func_add_func(ss[0], ss[1]));
    EXPECT_C_EQ(ref, res["m"](i));
  }
  
}
TEST(TestHAtom, s_state) {
  
  R1GTOs gtos(0);
  for(int n = -5; n <= 5; n++)
    gtos.Add(1, pow(2.0, n));
  gtos.Normalize();
  MatMap mat;
  gtos.CalcMat(&mat);
  
  MatrixXcd H = mat["t"] + mat["v"];
  MatrixXcd S = mat["s"];
  
  MatrixXcd c;
  VectorXcd eig;
  generalizedComplexEigenSolve(H, S, &c, &eig);
  EXPECT_C_NEAR(-0.5, eig(0), 0.0001);
  
}
TEST(TestHAtom, p_state) {

  R1GTOs gtos(1);
  for(int n = -10; n <= 10; n++)
    gtos.Add(2, pow(2.0, n));
  gtos.Normalize();
  MatMap mat;
  gtos.CalcMat(&mat);
  
  MatrixXcd H = mat["t"] + mat["v"];
  MatrixXcd S = mat["s"];
  
  MatrixXcd c;
  VectorXcd eig;
  generalizedComplexEigenSolve(H, S, &c, &eig);
  EXPECT_C_NEAR(-0.125, eig(0), 0.0001);

}
TEST(TestHAtom, dense) {

  R1GTOs gtos(0);
  for(int n = -20; n <= 20; n++)
    gtos.Add(1, pow(1.2, n));
  gtos.Normalize();
  MatMap mat;
  gtos.CalcMat(&mat);
  
  MatrixXcd H = mat["t"] + mat["v"];
  MatrixXcd S = mat["s"];
  
  MatrixXcd c;
  VectorXcd eig;
  
  CEigenSolveCanonical(H, S, pow(10.0, -7.0), &c, &eig);
  EXPECT_C_NEAR(-0.5, eig(0), 0.0001);
  
  // -- below test will fail.
  // CEigenSolveCanonical(H, S, 0.0, &c, &eig);
  // EXPECT_C_NEAR(-0.5, eig(0), 0.0001);

}

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
