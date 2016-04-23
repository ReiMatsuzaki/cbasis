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
    gs0.Normalize(); gs1.Normalize();
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
  MatVecMap res;
  gs0.CalcMatSTV(res);
  for(int i = 0; i < gs0.size_basis(); i++) {
    EXPECT_C_EQ(1.0, res.mat["s"](i, i))
      << "i = " << i;
  }
}
TEST_F(TestR1GTOsCont, vec) {

  MatVecMap mat_vec;

  R1STOs driv; driv.Add(2.0, 2, 1.0);
  gs0.CalcVec(driv, mat_vec, "m"); gs1.CalcVec(driv);

  EXPECT_EQ(5, mat_vec.vec["m"].size());
  //EXPECT_C_EQ(mat_vec.vec["m"](0), vec["m"](0));

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

  MatVecMap res;
  EXPECT_ANY_THROW(res.mat["s"]);
  EXPECT_ANY_THROW(res.vec["m"]);

  gtos.CalcMat();
  R1STOs stos;
  stos.Add(1.1, 2, dcomplex(0.1, 0.01));
  gtos.CalcVec(stos);

  EXPECT_NO_THROW(res.mat["s"]);
  EXPECT_NO_THROW(res.vec["m"]);

}
TEST_F(TestR1GTOs, max_n) {
  EXPECT_EQ(4, gtos.max_n());
}
TEST_F(TestR1GTOs, normalized) {

  MatVecMap res;
  gtos.CalcMatSTV(res);
  for(int i = 0; i < 4; i++)
    EXPECT_C_EQ(1.0, res.mat["s"](i, i)) << i;

}
TEST_F(TestR1GTOs, matrix) {

  MatVecMap res;
  gtos.CalcMatSTV(res);

  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++) {
      EXPECT_C_EQ(CIP(CNormalize(gs[i]), CNormalize(gs[j])),
		  res.mat["s"](i, j)) << i << j;
      dcomplex tref = -0.5*CIP(CNormalize(gs[i]),
			       OpD<dcomplex, double, 2>(),
			       CNormalize(gs[j]));
      EXPECT_C_EQ(tref, res.mat["t"](i, j)) << i << j;
      dcomplex vref = -CIP(CNormalize(gs[i]),
			   OpRm<dcomplex, double>(-1),
			   CNormalize(gs[j]));
      EXPECT_C_EQ(vref, res.mat["v"](i, j)) << i << j;
    }

}
TEST_F(TestR1GTOs, matrix_sto) {

  R1STOs sto; sto.Add(1.1, 1, 0.3); sto.Add(1.2, 1, 0.35);
  cout << "calc" << endl;
  
  MatVecMap res;

  gtos.CalcMatSTO(sto, res, "vexp");
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
      EXPECT_C_EQ(nconst*(1.1*val1+1.2*val2), res.mat["vexp"](i, j))
	<< i << j;

    }
  

}
TEST_F(TestR1GTOs, matrix_h) {

  gtos.CalcMatH();
  R1STOs sto; sto.Add(1.1, 1, 0.3); sto.Add(1.2, 1, 0.35);
  gtos.CalcMatH(sto, "hv2");

  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++) {
      CGTO c_gs_i(gs[i]); c_gs_i.SetComplexConjugate();
      EXPECT_C_EQ(CIP(CNormalize(c_gs_i), CNormalize(gs[j])),
		  gtos.mat("hs")(i, j)) << i << j;
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
