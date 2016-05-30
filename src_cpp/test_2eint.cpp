#include <iostream>
#include <gtest/gtest.h>
#include "gtest_plus.hpp"
#include "b2eint.hpp"
#include "symmolint.hpp"

using namespace std;
using namespace Eigen;
using namespace l2func;

TEST(first, first) {
  EXPECT_EQ(2, 1+1);
}
TEST(B2EInt, Simple) {
  IB2EInt *eri = new B2EIntMem(100);
  eri->Set(1, 2, 3, 4,
	   5, 6, 7, 8, 1.1);
  eri->Set(0, 0, 0, 1,
	   0, 0, 3, 0, 1.2);
  eri->Set(0, 2, 0, 0,
	   0, 0, 1, 0, 1.3);
  eri->Set(1, 0, 0, 0,
	   0, 1, 0, 0, 1.4);
  int ib,jb,kb,lb,i,j,k,l,t;
  dcomplex v;
  EXPECT_TRUE(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v));
  EXPECT_EQ(1, ib);
  EXPECT_EQ(2, jb);
  EXPECT_EQ(3, kb);
  EXPECT_EQ(4, lb);
  EXPECT_EQ(5, i);
  EXPECT_EQ(6, j);
  EXPECT_EQ(7, k);
  EXPECT_EQ(8, l);
  EXPECT_EQ(0, t);
  EXPECT_C_EQ(1.1, v);
  
  EXPECT_TRUE(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v));
  EXPECT_EQ(k, 3);
  EXPECT_TRUE(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v));
  EXPECT_EQ(jb, 2);
  EXPECT_TRUE(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v));
  EXPECT_EQ(v, 1.4);
  EXPECT_FALSE(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v));
  
  eri->Reset();
  while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
    cout << ib << jb << kb << lb << endl;
  }
  delete eri;
}
TEST(SymGTOs, CalcERI) {

  Irrep Ap = Cs_Ap();
  SymGTOs gtos(SymmetryGroup_Cs());

  // -- A' --
  VectorXcd zeta0(3); zeta0 << 0.1, 0.4, 1.1;
  gtos.AddSub(Sub_s(Ap, Vector3cd(0, 0, 0), zeta0));

  // -- potential --
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  // -- Calculation --
  IB2EInt *eri = new B2EIntMem(100);
  gtos.CalcERI(eri);

  // -- print --
  eri->Reset();

  int ib,jb,kb,lb,i,j,k,l,t;
  dcomplex v;
  while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
    cout << ib << jb << kb << lb << v << endl;
  }
  delete eri;
}

int main (int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
