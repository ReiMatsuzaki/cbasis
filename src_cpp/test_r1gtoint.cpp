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
#include "exp_int.hpp"

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
class TestR1GTOsCont : public ::testing::Test {
public:
  R1GTOs gs0;
  R1GTOs gs1;
  TestR1GTOsCont(): gs0(1), gs1(1) {
    gs0.Add(2, dcomplex(1.1, 0.2)); gs1.Add(2, dcomplex(1.1, 0.2));
    VectorXcd zs(3); zs << 2.1, dcomplex(0.1, 0.2), dcomplex(0.2, 0.3);
    MatrixXcd cs(2, 3); 
    cs <<
      1, 1, 0,
      1, 0, 2;
    gs0.Add(2, zs, cs); gs1.Add(2, zs);
    
    gs0.Add(3, 0.3); gs1.Add(3, 0.3);
    gs0.Add(3, 0.4); gs1.Add(3, 0.4);
  }
};
TEST_F(TestR1GTOsCont, offset) {

  EXPECT_EQ(4, gs0.conts_.size());

  EXPECT_EQ(0, gs0.conts_[0].offset);
  EXPECT_EQ(1, gs0.conts_[1].offset);
  EXPECT_EQ(3, gs0.conts_[2].offset);
  EXPECT_EQ(4, gs0.conts_[3].offset);

}
TEST_F(TestR1GTOsCont, size) {
  EXPECT_EQ(1+3+2, gs0.size_prim());
  EXPECT_EQ(1+2+2, gs0.size_basis());
}
TEST_F(TestR1GTOsCont, prim) {

  for(int i = 0; i < 5; i++) {
    EXPECT_C_EQ(gs0.prim(i).z, gs1.prim(i).z)
      << "i = " << i;
  }

}
TEST_F(TestR1GTOsCont, set) {
  
  VectorXcd zs(6); zs << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
  gs0.Set(2, zs);
  EXPECT_C_EQ(0.1, gs0.prim(0).z);
  EXPECT_C_EQ(0.2, gs0.prim(1).z);
  EXPECT_C_EQ(0.3, gs0.prim(2).z);
  EXPECT_C_EQ(0.4, gs0.prim(3).z);
  EXPECT_C_EQ(0.5, gs0.prim(4).z);
  EXPECT_C_EQ(0.6, gs0.prim(5).z);

}
TEST_F(TestR1GTOsCont, normalized) {

  gs0.Normalize();
  gs0.CalcMat();
  for(int i = 0; i < gs0.size_basis(); i++) {
    EXPECT_C_EQ(1.0, gs0.mat("s")(i, i))
      << "i = " << i;
  }
}
TEST_F(TestR1GTOsCont, vec) {

  R1STOs driv; driv.Add(2.0, 2, 1.0);
  gs0.CalcVec(driv); gs1.CalcVec(driv);
  EXPECT_EQ(5, gs0.vec("m").size());

  EXPECT_C_EQ(gs0.vec("m")(0), gs1.vec("m")(0));

  /*
  dcomplex mcalc = gs0.vec("m")(1);
  dcomplex mref(0);
  for(int i = 0; i < 3; i++)
    mref += gs1.vec("m")(1+i);
    */  

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
TEST_F(TestR1GTOs, offset) {

  EXPECT_EQ(4, gtos.conts_.size());
  EXPECT_EQ(0, gtos.conts_[0].offset);
  EXPECT_EQ(1, gtos.conts_[1].offset);
  EXPECT_EQ(2, gtos.conts_[2].offset);
  EXPECT_EQ(3, gtos.conts_[3].offset);

}
TEST_F(TestR1GTOs, add) {

  EXPECT_EQ(4, gtos.size_basis());
  EXPECT_C_EQ(dcomplex(1.1, 0.2), gtos.prim(0).z);
  EXPECT_EQ(  3                 , gtos.prim(1).n);
  EXPECT_C_EQ(dcomplex(0.1, 0.2), gtos.prim(2).z);
  EXPECT_C_EQ(dcomplex(0.2, 0.3), gtos.prim(3).z);

}
TEST_F(TestR1GTOs, set) {

  EXPECT_EQ(4, gtos.size_basis());
  VectorXcd zs(4); zs << 1.1, 1.2, 1.3, 1.4;
  gtos.Set(5, zs);

  EXPECT_EQ(5, gtos.prim(1).n);
  EXPECT_C_EQ(1.1, gtos.prim(0).z);
  EXPECT_C_EQ(1.2, gtos.prim(1).z);
  EXPECT_C_EQ(1.3, gtos.prim(2).z);
  EXPECT_C_EQ(1.4, gtos.prim(3).z);

  VectorXcd zs1(5); zs1 << 1.1, 1.2, 1.3, 1.4, 1.2;
  EXPECT_ANY_THROW(gtos.Set(4, zs1));


}
TEST_F(TestR1GTOs, Exception) {

  EXPECT_ANY_THROW(gtos.mat("s"));
  EXPECT_ANY_THROW(gtos.vec("m"));

  gtos.CalcMat();
  R1STOs stos;
  stos.Add(1.1, 2, dcomplex(0.1, 0.01));
  gtos.CalcVec(stos);

  EXPECT_NO_THROW(gtos.mat("s"));
  EXPECT_NO_THROW(gtos.vec("m"));

}
TEST_F(TestR1GTOs, max_n) {
  EXPECT_EQ(4, gtos.max_n());
}
TEST_F(TestR1GTOs, normalized) {

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
TEST_F(TestR1GTOs, matrix_sto) {

  R1STOs sto; sto.Add(1.1, 1, 0.3); sto.Add(1.2, 1, 0.35);
  cout << "calc" << endl;
  gtos.CalcMat(sto, "vexp");
  cout << "end" << endl;

  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++) {
      dcomplex nconst = 1.0/(CNorm(gs[i]) * CNorm(gs[j]));
      dcomplex val1 = STO_GTO_Int(dcomplex(0.3, 0.0),
				 gs[i].z() + gs[j].z(),
				  gs[i].n()+gs[j].n()+1);
      dcomplex val2 = STO_GTO_Int(dcomplex(0.35, 0.0),
				 gs[i].z() + gs[j].z(),
				 gs[i].n()+gs[j].n()+1);
      EXPECT_C_EQ(nconst*(1.1*val1+1.2*val2), gtos.mat("vexp")(i, j))
	<< i << j;

    }
  

}
TEST_F(TestR1GTOs, vector_sto) {

  vector<CSTO> ss;
  R1STOs stos;

  ss.push_back(  CSTO( 1.1, 2, dcomplex(0.1, 0.01)));
  stos.Add(1.1, 2, dcomplex(0.1, 0.01));

  ss.push_back(  CSTO( 0.11, 3, 0.3));
  stos.Add(0.11, 3, 0.3);

  gtos.CalcVec(stos);
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
  gtos.CalcVec(driv);

  MatrixXcd L = gtos.mat("t") + gtos.mat("v") -0.5* gtos.mat("s");
  VectorXcd& m = gtos.vec("m");
  dcomplex a = (m.array() *
		(L.fullPivLu().solve(m)).array()).sum();
  //dcomplex a = vec["m"].dot(L.fullPivLu().solve(vec["m"]));
  return a;

}
TEST(OptAlpha, WithContraction) {
  
  
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
  
  R1GTOs gtos(1);
  gtos.Add(2, zs_k);
  gtos.CalcMat();
  MatrixXcd xmat;
  CanonicalMatrix(gtos.mat("s"), pow(10.0, -5.0), &xmat);

  R1GTOs gtos2(1);
  gtos2.Add(2, zs_h);
  gtos2.Add(2, zs_k, xmat.transpose());

  EXPECT_EQ(10+7, gtos2.size_prim());

  R1STOs driv; driv.Add(2.0, 2, 1.0);
  VectorXi opt_idx(10);
  opt_idx << 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

  double h(0.0001);
  double eps(0.000003);
  dcomplex z_shift(0.0, -0.02);
  bool convq;
  dcomplex alpha;
  OptAlphaShift(driv, opt_idx, 0.5,
		h, 20, eps, /*pow(10.0, -8.0)*/ 0.0,
		gtos2, z_shift, &convq, &alpha);

  EXPECT_TRUE(convq);
  dcomplex ref_alpha(dcomplex(-5.6568937518988989, 1.0882823480377297));
  EXPECT_C_NEAR(ref_alpha, alpha, pow(10.0, -9.0));  

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
  
  EXPECT_C_EQ(CalcAlpha(0.0),
	      CalcAlphaFull(driv, gtos, 0.5));
  

}
TEST(OptAlpha, CalcAlphaSTO) {
  R1GTOs gtos(1);
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

  R1STOs driv; driv.Add(2.0, 2, 1.0);
  R1STOs pot;  pot.Add(-1.0, 0, 1.0);

  // 2016/4/12
  // alpha = -4.770452, -0.235177

  EXPECT_C_NEAR(CalcAlpha(driv, gtos, 0.5, pot),
		dcomplex(-4.770452, -0.235177),
		pow(10.0, -4.0));

}
TEST(OptAlpha, WithoutCanonical) {

  R1GTOs gtos(1);
  VectorXcd zs(15);  
  zs <<
    0.463925,
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
  double eps(0.000003);
  dcomplex z_shift(0.0, -0.02);
  bool convq;
  dcomplex alpha;
  OptAlphaShift(driv, opt_idx, 0.5,
		h, 20, eps, /*pow(10.0, -8.0)*/ 0.0,
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
  double eps(0.000003);
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
  for(int n = -10; n <= 10; n++)
    gtos.Add(1, pow(2.0, n));
  gtos.Normalize();
  gtos.CalcMat();
  
  MatrixXcd H = gtos.mat("t") + gtos.mat("v");
  MatrixXcd S = gtos.mat("s");
  
  MatrixXcd c;
  VectorXcd eig;
  generalizedComplexEigenSolve(H, S, &c, &eig);
  EXPECT_C_NEAR(-0.5, eig(0), 0.0001);

  VectorXcd rs(1); rs << 1.1;
  VectorXcd ys;
  gtos.AtR(c.col(0), rs, &ys);
  dcomplex r(rs[0]);
  EXPECT_C_NEAR(2.0*r*exp(-r), ys[0], pow(10.0, -6.0));
  
}
TEST(TestHAtom, contraction) {

  R1GTOs gtos(0);
  VectorXcd zs(11);
  for(int n = -5; n <= 5; n++)
    zs(n + 5) = pow(2.0, n);

  MatrixXcd coef = MatrixXcd::Zero(10, 11);
  for(int i = 0; i < 10; i++) {
    coef(i, i)   = 1;
    coef(i, i+1) = 2;
  }
  
  gtos.Add(1, zs, coef);
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
