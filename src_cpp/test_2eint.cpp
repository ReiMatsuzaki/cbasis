#include <iostream>
#include <gtest/gtest.h>
#include "gtest_plus.hpp"
#include "b2eint.hpp"
#include "symmolint.hpp"
#include "timer.hpp"

using namespace std;
using namespace Eigen;
using namespace l2func;

TEST(first, first) {
  EXPECT_EQ(2, 1+1);
}
class TestB2EInt :public ::testing::Test {
public:
  IB2EInt *eri;
public:
  TestB2EInt() {
    eri = new B2EIntMem(100);
    eri->Set(1, 2, 3, 4,
	     5, 6, 7, 8, 1.1);
    eri->Set(0, 0, 0, 1,
	     0, 0, 3, 0, 1.2);
    eri->Set(0, 2, 0, 0,
	     0, 0, 1, 0, 1.3);
    eri->Set(1, 0, 0, 0,
	     0, 1, 0, 0, 1.4);
  }
  ~TestB2EInt() {
    delete eri;
  }
};
TEST_F(TestB2EInt, Get) {
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
}
TEST_F(TestB2EInt, At) {

  dcomplex v;
  v = eri->At(1, 2, 3, 4, 5, 6, 7, 8);
  EXPECT_C_EQ(1.1, v);

  EXPECT_ANY_THROW(eri->At(0, 0, 0, 0, 1, 1, 1, 1));
}
TEST_F(TestB2EInt, size) {
  EXPECT_EQ(4, eri->size());
  EXPECT_EQ(100, eri->capacity());
}
TEST_F(TestB2EInt, IO) {

  string fn("eri.bin");
  eri->Write(fn);
  
  IB2EInt *eri2 = new B2EIntMem(fn);

  EXPECT_C_EQ(eri->At( 1, 2, 3, 4, 5, 6, 7, 8),
	      eri2->At(1, 2, 3, 4, 5, 6, 7, 8));
  EXPECT_C_EQ(eri->At( 0, 0, 0, 1, 0, 0, 3, 0),
	      eri2->At(0, 0, 0, 1, 0, 0, 3, 0));
  EXPECT_C_EQ(eri->At( 0, 2, 0, 0, 0, 0, 1, 0),
	      eri2->At(0, 2, 0, 0, 0, 0, 1, 0));
  EXPECT_C_EQ(eri->At( 1, 0, 0, 0, 0, 1, 0, 0),
	      eri2->At(1, 0, 0, 0, 0, 1, 0, 0));

  delete eri2;
  
}
VectorXcd OneVec(dcomplex z) {
  VectorXcd zs(1); zs << z;
  return zs;
}
TEST(SymGTOs, CalcERI) {

  pSymmetryGroup Cs = SymmetryGroup::Cs();
  SymGTOs gtos(Cs);

  // -- A --
  gtos.AddSub(Sub_s(Cs, 0, Vector3cd(0.0, 0.0,  0.4), OneVec(1.2)));
  gtos.AddSub(Sub_s(Cs, 0, Vector3cd(0.0, 0.0,  0.0), OneVec(1.4)));
  gtos.AddSub(Sub_s(Cs, 0, Vector3cd(0.0, -0.2, 0.0), OneVec(1.1)));
  gtos.AddSub(Sub_s(Cs, 0, Vector3cd(0.2, 0.0,  0.1), OneVec(1.0)));
  

  // -- potential --
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos.SetAtoms(xyzq);
  try {
    gtos.SetUp();
  } catch(exception& e) {
    cout << "in setup" << endl;
    cout << e.what() << endl;
    throw e;
  }

  // -- Calculation --
  IB2EInt *eri = new B2EIntMem(pow(4, 4));
  gtos.CalcERI(eri, 1);
  // -- size check --
  EXPECT_EQ(pow(4, 4), eri->size());

  // -- Symmetry check --
  // eri is chemist's notation.
  EXPECT_C_EQ(eri->At(0, 0, 0, 0,  0, 1, 2, 3),
	      eri->At(0, 0, 0, 0,  1, 0, 2, 3));
  EXPECT_C_EQ(eri->At(0, 0, 0, 0,  0, 1, 2, 3),
	      eri->At(0, 0, 0, 0,  0, 1, 3, 2));
  EXPECT_C_EQ(eri->At(0, 0, 0, 0,  0, 1, 2, 3),
	      eri->At(0, 0, 0, 0,  2, 3, 0, 1));

  // -- copied from neo_ccol/aoint/gto_utest.cpp --
  // (00|00): (1.23607744235395,0)
  // (01|23): (1.01774464000324,0)
  // (00|11): (1.03695966484363,0)
  // (00|01): (1.11023051353687,0)
  // (00|12): (0.992414842562312,0)

  double eps(pow(10.0, -8.0));
  dcomplex ref0000(1.23607744235395, 0);
  dcomplex ref0123(1.01774464000324, 0);
  dcomplex ref0011(1.03695966484363,0);
  dcomplex ref0001(1.11023051353687,0);
  dcomplex ref0012(0.992414842562312,0);
  
  EXPECT_C_NEAR(ref0000, eri->At(0, 0, 0, 0, 0, 0, 0, 0), eps);
  EXPECT_C_NEAR(ref0123, eri->At(0, 0, 0, 0, 0, 2, 1, 3), eps);
  EXPECT_C_NEAR(ref0011, eri->At(0, 0, 0, 0, 0, 1, 0, 1), eps);
  EXPECT_C_NEAR(ref0001, eri->At(0, 0, 0, 0, 0, 0, 0, 1), eps);
  EXPECT_C_NEAR(ref0012, eri->At(0, 0, 0, 0, 0, 1, 0, 2), eps);
  delete eri;

}
SubSymGTOs SubC1(Vector3i ns, Vector3cd xyz, dcomplex z) {
  pSymmetryGroup C1 = SymmetryGroup::C1();
  SubSymGTOs sub(C1);
  sub.AddNs(ns);
  sub.AddXyz(xyz);
  sub.AddRds(Reduction(0, MatrixXcd::Ones(1, 1)));
  VectorXcd zeta(1); zeta << z;
  sub.AddZeta(zeta);
  sub.SetUp();
  return sub;
}
TEST(SymGTOs, CalcERI2) {

  pSymmetryGroup Cs = SymmetryGroup::Cs();
  SymGTOs gtos(Cs);
  SubSymGTOs s1(Cs);
  s1.AddXyz(Vector3cd(0.0, 0.0, 0.4));
  s1.AddNs( Vector3i( 0,   1,   0));
  VectorXcd z1(1); z1 << 1.2; s1.AddZeta(z1);
  s1.AddRds(Reduction(0, MatrixXcd::Ones(1, 1)));

  gtos.AddSub(Sub_mono(Cs, 0, Vector3cd(0.0, 0.0,  0.4),
		       Vector3i(0, 1, 0), OneVec(1.2)));
  gtos.AddSub(Sub_mono(Cs, 0, Vector3cd(0.0, 0.0,  0.0),
		       Vector3i(1, 1, 0), OneVec(1.4)));
  gtos.AddSub(Sub_mono(Cs, 0, Vector3cd(0.0, -0.2, 0.0),
		       Vector3i(1, 1, 1), OneVec(1.1)));
  gtos.AddSub(Sub_mono(Cs, 0, Vector3cd(0.2, 0.0,  0.1),
		       Vector3i(0, 3, 0), OneVec(1.0)));

  // -- potential --
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  // -- Calculation --
  IB2EInt *eri = new B2EIntMem(pow(4, 4));
  gtos.CalcERI(eri, 1);

  // -- size check --
  //  cout << 1 << endl;
  //  EXPECT_EQ(pow(4, 4), eri->size());

  // -- reference --
  // (00|00): (1.00946324498146,0)
  // (00|11): (0.133478467219949,0)
  // (00|01): (0,0)
  // (00|12): (0.0481552620386253,3.14588102745493e-19)
  // (01|23): (0.0240232215271162,1.00002719131888e-18)

  double eps(pow(10.0, -8.0));
  EXPECT_C_NEAR(1.00946324498146,  eri->At(0, 0, 0, 0, 0, 0, 0, 0), eps);
  EXPECT_C_NEAR(0.133478467219949, eri->At(0, 0, 0, 0, 0, 1, 0, 1), eps);
  EXPECT_C_NEAR(0.0, eri->At(0, 0, 0, 0, 0, 0, 0, 1), eps);
  EXPECT_C_NEAR(0.0481552620386253, eri->At(0, 0, 0, 0, 0, 1, 0, 2), eps);
  EXPECT_C_NEAR(0.0240232215271162, eri->At(0, 0, 0, 0, 0, 2, 1, 3), eps);

  delete eri;

}
TEST(SymGTOs, CalcERI_time) {
/*
  SymGTOs gtos(SymmetryGroup_C1());

  Vector3cd xyz1(0.0, 0.0, 0.0);
  Vector3cd xyz2(0.0, 0.4, 0.0);

  SubSymGTOs sub1;
  
  int num_at(2);
  sub1.AddXyz(Vector3cd(0.0, 0.4, 0.0));
  sub1.AddXyz(Vector3cd(0.0, 0.0, 0.5));

  int  num_ns(3);
  sub1.AddNs(Vector3i(1, 0, 0));
  sub1.AddNs(Vector3i(0, 1, 0));
  sub1.AddNs(Vector3i(0, 0, 1));


  int num_z(5);
  VectorXcd zs(num_z); zs << 1.1, 1.2, 1.3, 1.4, 1.5;
  sub1.AddZeta(zs);
  MatrixXcd cs = MatrixXcd::Ones(num_at, num_ns);
  sub1.AddRds(Reduction(0, cs));
  sub1.SetSym(MatrixXi::Ones(1, 2*num_ns),
	      MatrixXi::Ones(1, 2*num_ns));
  sub1.SetUp();
  
  gtos.AddSub(sub1);
  gtos.SetUp();

  // -- potential --
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  IB2EInt *eri = new B2EIntMem(pow(num_z, 4));
  gtos.CalcERI(eri);
  delete eri;
*/
}
TEST(SymGTOs, CalcERI_sym_p) {

  pSymmetryGroup Cs = SymmetryGroup::Cs();
  SymGTOs gtos(Cs);

  Vector3cd xyz0(0.0, 0.0, 0.0);

  // ---- p orbital ----
  SubSymGTOs sub_p(Cs);
  sub_p.AddXyz(Vector3cd(0, 0, 0));
  sub_p.AddNs( Vector3i( 1, 0, 0));
  sub_p.AddNs( Vector3i( 0, 1, 0));

  int num_z_p(1);
  VectorXcd zs_p(num_z_p); zs_p << 1.1;
  sub_p.AddZeta(zs_p);
  MatrixXcd c(1, 2);
  c(0, 0) = 1.0; c(0, 1) = 0.0;
  sub_p.AddRds(Reduction(0, c));
  c(0, 0) = 0.0; c(0, 1) = 1.0;
  sub_p.AddRds(Reduction(1, c));

  sub_p.SetUp();
  gtos.AddSub(sub_p);

  // -- potential --
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();
  //  cout << gtos.str() << endl;

  IB2EInt *eri1 = new B2EIntMem(pow(3*num_z_p, 4));
  IB2EInt *eri2 = new B2EIntMem(pow(3*num_z_p, 4));
  gtos.CalcERI(eri1, 1);
  gtos.CalcERI(eri2, 2);

  int ib,jb,kb,lb,i,j,k,l,t;
  dcomplex v;
  eri1->Reset();
  while(eri2->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
    if(abs(v) > 0.000001)
      EXPECT_C_EQ(v, eri1->At(ib, jb, kb, lb, i, j, k, l)) <<
	ib << jb << kb << lb << i << j << k << l;
  }

  delete eri1;
  delete eri2;

}
TEST(SymGTOs, CalcERI_sym) {

  pSymmetryGroup sym = SymmetryGroup::D2h();
  Irrep irrep_s  = sym->GetIrrep("Ag");
  Irrep irrep_x = sym->GetIrrep("B1u");
  Irrep irrep_y = sym->GetIrrep("B2u");
  Irrep irrep_z = sym->GetIrrep("B3u");

  /*
  pSymmetryGroup sym = SymmetryGroup::Cs();
  Irrep irrep_s  = sym->GetIrrep("A'");
  Irrep irrep_x = sym->GetIrrep("A'");
  Irrep irrep_y = sym->GetIrrep("A'");
  Irrep irrep_z = sym->GetIrrep("A''");
  */

  SymGTOs gtos(sym);

  Vector3cd xyz0(0.0, 0.0, 0.0);

  // ---- s orbital ----
  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs( Vector3i( 0, 0, 0));
  int num_z(1);
  VectorXcd zs(num_z); zs << 1.1;
  sub_s.AddZeta(zs);
  sub_s.AddRds(Reduction(irrep_s, MatrixXcd::Ones(1, 1)));
  sub_s.SetUp();
  gtos.AddSub(sub_s);

  // ---- p orbital ----
  SubSymGTOs sub_p(sym);
  sub_p.AddXyz(Vector3cd(0, 0, 0));
  sub_p.AddNs( Vector3i( 1, 0, 0));
  sub_p.AddNs( Vector3i( 0, 1, 0));
  sub_p.AddNs( Vector3i( 0, 0, 1));
  VectorXcd zs_p(1); zs_p << 1.1; sub_p.AddZeta(zs_p);
  MatrixXcd c1(1, 3); c1 << 1.0, 0.0, 0.0; sub_p.AddRds(Reduction(irrep_x, c1));
  MatrixXcd c2(1, 3); c2 << 0.0, 1.0, 0.0; sub_p.AddRds(Reduction(irrep_y, c2));
  MatrixXcd c3(1, 3); c3 << 0.0, 0.0, 1.0; sub_p.AddRds(Reduction(irrep_z, c3));
  sub_p.SetUp();
  gtos.AddSub(sub_p);

  // -- potential --
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();  

  BMatSet mat;
  gtos.CalcMat(&mat);

  int n(gtos.size_basis());
  IB2EInt *eri1 = new B2EIntMem(pow(n, 4));
  IB2EInt *eri2 = new B2EIntMem(pow(n, 4));
  gtos.CalcERI(eri1, 1);
  gtos.CalcERI(eri2, 2);

  int ib,jb,kb,lb,i,j,k,l,t;
  dcomplex v;
  eri1->Reset();
  while(eri2->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
    if(eri1->Exist(ib, jb, kb, lb, i, j, k, l)) {
      EXPECT_C_EQ(v, eri1->At(ib, jb, kb, lb, i, j, k, l)) <<
	ib << jb << kb << lb << " : " << i << j << k << l;
    } else {
      EXPECT_TRUE(abs(v) < 0.000001) <<
	v << " : " <<
	ib << jb << kb << lb << " : " <<
	i << j << k << l;
    }

  }

  delete eri1;
  delete eri2;
}
TEST(SymGTOs, CalcERI_time2) {

  Timer timer;
  pSymmetryGroup sym = SymmetryGroup::D2h();
  Irrep irrep_s  = sym->GetIrrep("Ag");
  Irrep irrep_x = sym->GetIrrep("B1u");
  Irrep irrep_y = sym->GetIrrep("B2u");
  Irrep irrep_z = sym->GetIrrep("B3u");
  /*
  pSymmetryGroup sym = SymmetryGroup::Cs();
  Irrep irrep_s  = sym->GetIrrep("A'");
  Irrep irrep_x = sym->GetIrrep("A'");
  Irrep irrep_y = sym->GetIrrep("A'");
  Irrep irrep_z = sym->GetIrrep("A''");
  */
  /*
  pSymmetryGroup sym = SymmetryGroup::C4();
  Irrep irrep_s  = sym->GetIrrep("A");
  Irrep irrep_x = sym->GetIrrep("E");
  Irrep irrep_y = sym->GetIrrep("E");
  Irrep irrep_z = sym->GetIrrep("A");
  */

  SymGTOs gtos(sym);

  Vector3cd xyz0(0.0, 0.0, 0.0);

  // ---- s orbital ----
  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs( Vector3i( 0, 0, 0));

  int num_z(10);
  VectorXcd zs(num_z); zs << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0;

  sub_s.AddZeta(zs);
  sub_s.AddRds(Reduction(irrep_s, MatrixXcd::Ones(1, 1)));
  sub_s.SetUp();
  gtos.AddSub(sub_s);

  // ---- p orbital ----
  SubSymGTOs sub_p(sym);
  sub_p.AddXyz(Vector3cd(0, 0, 0));
  sub_p.AddNs( Vector3i( 1, 0, 0));
  sub_p.AddNs( Vector3i( 0, 1, 0));
  sub_p.AddNs( Vector3i( 0, 0, 1));
  int num_z_p(8);
  VectorXcd zs_p(num_z_p); zs_p << 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8;
  sub_p.AddZeta(zs_p);
  MatrixXcd c1(1, 3); c1 << 1.0, 0.0, 0.0; sub_p.AddRds(Reduction(irrep_x, c1));
  MatrixXcd c2(1, 3); c2 << 0.0, 1.0, 0.0; sub_p.AddRds(Reduction(irrep_y, c2));
  MatrixXcd c3(1, 3); c3 << 0.0, 0.0, 1.0; sub_p.AddRds(Reduction(irrep_z, c3));
  sub_p.SetUp();
  gtos.AddSub(sub_p);

  // -- potential --
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  BMatSet mat;
  gtos.CalcMat(&mat);

  IB2EInt *eri1 = new B2EIntMem(pow(num_z+3*num_z_p, 4)); 
  IB2EInt *eri2 = new B2EIntMem(pow(num_z+3*num_z_p, 4)); 

  timer.Start("method1"); gtos.CalcERI(eri1, 1); timer.End("method1");
  timer.Start("method2"); gtos.CalcERI(eri2, 2); timer.End("method2");
    
  EXPECT_C_EQ(eri1->At(0, 0, 0, 0, 2, 1, 5, 0),
	      eri2->At(0, 0, 0, 0, 2, 1, 5, 0));

  timer.Display();
  cout << "size(eri1): " << eri1->size() << endl;
  cout << "size(eri2): " << eri2->size() << endl;
  delete eri1;
  delete eri2;
}
TEST(Time, MatrixAccess) {
  int n(1000);
  MatrixXcd a(n, n);
  dcomplex *b = new dcomplex[n*n];
  
  
  Timer timer;

  timer.Start("Eigen");
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      a(i, j) = 1.2;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      dcomplex x = a(i, j);
  timer.End("Eigen");

  timer.Start("Array");
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      b[i+n*j] = 1.2;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      dcomplex x = b[i+n*j];
  timer.End("Array");

  timer.Start("Array2");
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      b[j+n*i] = 1.2;
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      dcomplex x = b[j+n*i];
  timer.End("Array2");
  timer.Display();

  Timer timer1;
  MatrixXcd xyz(3, 2);
  timer1.Start("MatrixXyz");
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      xyz(0, 0) = 1.0;
      xyz(0, 1) = 1.2;
    }
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      dcomplex x= xyz(0, 0);
      x =  xyz(0, 1);
    }
  timer1.End("MatrixXyz");

  dcomplex xs[2];
  timer1.Start("ArrayXyz");
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      xs[0]= 1.0*i*j;
      xs[1] = 1.2*i*j;
    }
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      dcomplex x= xs[0];
      x =  xs[1];
    }
  timer1.End("ArrayXyz");

  vector<dcomplex> cs(2);
  timer1.Start("VectorXyz");
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      cs[0]= 1.0*i*j;
      cs[1] = 1.2*i*j;
    }
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++) {
      dcomplex x= cs[0];
      x =  cs[1];
    }
  timer1.End("VectorXyz");  
  
  timer1.Display();

}
  /*
TEST(H2mole, element) {

  // ==== Symmetry ====
  pSymmetryGroup D2h = SymmetryGroup::D2h();

  // ==== Sub ====
  SubSymGTOs sub1(D2h);
  sub1.AddXyz(Vector3cd(0, 0, +0.7));
  sub1.AddXyz(Vector3cd(0, 0, -0.7));
  sub1.AddNs( Vector3i( 0, 0, 0));
  VectorXcd z1(4); z1 << 2.013, 0.1233, 0.0411, 0.0137; sub1.AddZeta(z1);
  // VectorXcd z1(2); z1 << 0.1233, 0.0411; sub1.AddZeta(z1);
  MatrixXcd c1_1(2, 1); c1_1 <<+1.0,+1.0; sub1.AddRds(Reduction(0, c1_1));
  MatrixXcd c1_2(2, 1); c1_2 <<+1.0,-1.0; sub1.AddRds(Reduction(1, c1_2));
  sub1.SetUp();

  SubSymGTOs sub2(D2h);
  sub2.AddXyz(Vector3cd(0, 0, +0.7));
  sub2.AddXyz(Vector3cd(0, 0, -0.7));
  sub2.AddNs( Vector3i( 0, 0, 1));
  VectorXcd z2(1); z2 << 1.0; sub2.AddZeta(z2);
  MatrixXcd C2_1(2, 1); C2_1 << +1,-1; sub2.AddRds(Reduction(0, C2_1));
  MatrixXcd C2_2(2, 1); C2_2 << +1,+1; sub2.AddRds(Reduction(1, C2_2));
  sub2.SetUp();

  SubSymGTOs sub3(D2h);
  sub3.AddXyz(Vector3cd(0, 0, 0));
  sub3.AddNs( Vector3i( 0, 0, 0));
  VectorXcd z3(1); z3 << dcomplex(0.011389, -0.002197); sub3.AddZeta(z3); 
  MatrixXcd C3_1(1, 1); C3_1 << 1; sub3.AddRds(Reduction(0, C3_1));
  sub3.SetUp();
  
  SubSymGTOs sub4(D2h);
  sub4.AddXyz(Vector3cd(0, 0, 0));
  sub4.AddNs( Vector3i( 2, 0, 0));
  sub4.AddNs( Vector3i( 0, 2, 0));
  sub4.AddNs( Vector3i( 0, 0, 2));
  VectorXcd z4(1); z4 << dcomplex(5.063464, -0.024632); sub4.AddZeta(z4);
  MatrixXcd C4_1(1, 3); C4_1 << -1,-1,+2; sub4.AddRds(Reduction(0, C4_1 ));
  sub4.SetUp();

  // ==== GTOs ====
  SymGTOs gtos(D2h);
  gtos.AddSub(sub1); gtos.AddSub(sub2); gtos.AddSub(sub3); gtos.AddSub(sub4);
  MatrixXcd xyzq(4, 2); xyzq <<
			  0,    0,
			  0,    0,
			  +0.7, -0.7,
			  1.0,  1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  // ==== matrix evaluation ====
  //  BMatSet mat;
  //  gtos.CalcMat(&mat);
  IB2EInt *eri = new B2EIntMem(pow(gtos.size_basis(), 4));
  gtos.CalcERI(eri, 1);

  // copied from ~/calc/cCOLUMBUS
  //1  1  1  1  6  2  3  1        0.18344       -0.02309
  //  EXPECT_C_EQ(dcomplex(0.18344, -0.02309),
  //	      eri->At(0, 0, 0, 0, 5, 1, 2, 0));
}
*/

int main (int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
