#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <gtest/gtest.h>

#include "eigen_plus.hpp"
#include "gtest_plus.hpp"
#include "r1gtoint.hpp"
#include "cip_exp.hpp"
#include "linspace.hpp"
#include "opt_alpha.hpp"

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

  gtos.CalcMat();
  for(int i = 0; i < 4; i++)
    EXPECT_C_EQ(1.0, gtos.mat("s")(i, i)) << i;
}
TEST_F(TestR1GTOs, matrix) {

  gtos.CalcMat();

  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++) {
      EXPECT_C_EQ(CIP(CNormalize(gs[i]), CNormalize(gs[j])),
		  gtos.mat("s")(i, j)) << i << j;
      dcomplex tref = -0.5*CIP(CNormalize(gs[i]),
			       OpD<dcomplex, double, 2>(),
			       CNormalize(gs[j]));
      EXPECT_C_EQ(tref, gtos.mat("t")(i, j)) << i << j;
      dcomplex vref = -CIP(CNormalize(gs[i]),
			   OpRm<dcomplex, double>(-1),
			   CNormalize(gs[j]));
      EXPECT_C_EQ(vref, gtos.mat("v")(i, j)) << i << j;
    }
}
TEST_F(TestR1GTOs, vector_sto) {

  vector<CSTO> ss;
  R1STOs stos;

  ss.push_back(  CSTO( 1.1, 2, dcomplex(0.1, 0.01)));
  stos.Add(1.1, 2, dcomplex(0.1, 0.01));

  ss.push_back(  CSTO( 0.11, 3, 0.3));
  stos.Add(0.11, 3, 0.3);

  gtos.CalcVecSTO(stos);
  for(int i = 0; i < 4; i++) {
    dcomplex ref = CIP(CNormalize(gs[i]),
		       func_add_func(ss[0], ss[1]));
    EXPECT_C_EQ(ref, gtos.vec("m")(i));
  }
  
}
dcomplex CalcAlpha(dcomplex z_shift) {

  R1GTOs gtos(1);
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

  gtos.CalcMat();
  gtos.CalcVecSTO(driv);

  MatrixXcd L = gtos.mat("t") + gtos.mat("v") -0.5* gtos.mat("s");
  VectorXcd& m = gtos.vec("m");
  dcomplex a = (m.array() *
		(L.fullPivLu().solve(m)).array()).sum();
  //dcomplex a = vec["m"].dot(L.fullPivLu().solve(vec["m"]));
  return a;

}
TEST(OptAlpha, ShiftKaufmann) {

  // 2016/3/26
  // ARCH=fast
  // 258ms

  R1GTOs gtos(1);
  
  int max_iter(100);
  dcomplex h(0.0001);
  dcomplex ih(0.0, 0.0001);
  dcomplex ii(0.0, 1.0);
  double eps(0.000001);

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

  // -- from 2016/3/25
  // -- -0.00293368-0.0204361j 
  
  EXPECT_TRUE(i < 90) << i;
  dcomplex ref(-0.00293368, -0.0204361);
  EXPECT_C_NEAR(ref, z_shift, eps*10);
  dcomplex ref_alpha(dcomplex(-5.6568937518988989, 1.0882823480377297));
  EXPECT_C_NEAR(ref_alpha, alpha, pow(10.0, -9.0));
}
TEST(OptAlpha, WithoutCanonical) {

  R1GTOs gtos(1);
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

  R1STOs driv; driv.Add(2.0, 2, 1.0);
  VectorXi opt_idx(8);
  opt_idx << 7, 8, 9, 10, 11, 12, 13, 14;

  double h(0.0001);
  double eps(0.000001);
  dcomplex z_shift(0.0, -0.02);
  bool convq;
  dcomplex alpha;
  OptAlphaShift(driv, opt_idx, 0.5,
		h, 100, eps, /*pow(10.0, -8.0)*/ 0.0,
		gtos, z_shift, &convq, &alpha);

  // -- from 2016/3/25
  // -- -0.00293368-0.0204361j 
  EXPECT_TRUE(convq);
  dcomplex ref(-0.00293368, -0.0204361);
  EXPECT_C_NEAR(ref, z_shift, eps*10);
  dcomplex ref_alpha(dcomplex(-5.6568937518988989, 1.0882823480377297));
  EXPECT_C_NEAR(ref_alpha, alpha, pow(10.0, -9.0));  

}
TEST(OptAlpha, WithCanonical) {

  R1GTOs gtos(1);
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

  R1STOs driv; driv.Add(2.0, 2, 1.0);
  VectorXi opt_idx(8);
  opt_idx << 7, 8, 9, 10, 11, 12, 13, 14;

  double h(0.0001);
  double eps(0.000001);
  dcomplex z_shift(0.0, -0.02);
  bool convq;
  dcomplex alpha;
  OptAlphaShift(driv, opt_idx, 0.5,
		h, 100, eps, pow(10.0, -8.0),
		gtos, z_shift, &convq, &alpha);

  // -- from 2016/3/25
  // -- -0.00293368-0.0204361j 
  EXPECT_TRUE(convq);
  dcomplex ref(-0.00293368, -0.0204361);
  EXPECT_C_NEAR(ref, z_shift, eps*10);
  dcomplex ref_alpha(dcomplex(-5.6568937518988989, 1.0882823480377297));
  EXPECT_C_NEAR(ref_alpha, alpha, pow(10.0, -9.0));  

}
TEST(TestHAtom, s_state) {
  
  R1GTOs gtos(0);
  for(int n = -5; n <= 5; n++)
    gtos.Add(1, pow(2.0, n));
  gtos.Normalize();
  gtos.CalcMat();
  
  MatrixXcd H = gtos.mat("t") + gtos.mat("v");
  MatrixXcd S = gtos.mat("s");
  
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
  gtos.CalcMat();
  
  MatrixXcd H = gtos.mat("t") + gtos.mat("v");
  MatrixXcd S = gtos.mat("s");
  
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
  gtos.CalcMat();
  
  MatrixXcd H = gtos.mat("t") + gtos.mat("v");
  MatrixXcd S = gtos.mat("s");
  
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
