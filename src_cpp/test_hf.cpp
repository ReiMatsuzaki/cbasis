#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include <gtest/gtest.h>
#include <Eigen/Core>
#include "gtest_plus.hpp"

#include "one_int.hpp"
#include "two_int.hpp"
#include "symmolint.hpp"
#include "b2eint.hpp"
#include "symgroup.hpp"
#include "eigen_plus.hpp"
#include "mo.hpp"
#include "timer.hpp"
#include "trans_eri.hpp"

using namespace std;
using namespace l2func;
using namespace Eigen;
TEST(Matrix, J) {

  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::Cs();
  
  // GTO
  SymGTOs g_i(new _SymGTOs); g_i->SetSym(sym);
  SymGTOs g_1(new _SymGTOs); g_1->SetSym(sym);
  SymGTOs g_0(new _SymGTOs); g_0->SetSym(sym);
  SymGTOs g_full(new _SymGTOs); g_full->SetSym(sym);

  VectorXcd zeta_i(2); zeta_i << 0.4, 1.0;
  SubSymGTOs sub_i(Sub_mono(0,
			    Vector3cd(0, 0, 0), Vector3i(0, 0, 0), zeta_i));
  g_i->AddSub(     sub_i);
  g_full->AddSub(  sub_i);

  VectorXcd zeta0(2); zeta0 << dcomplex(0.5, 0.0), dcomplex(0.4, 0.1);
  SubSymGTOs sub_0(Sub_mono(1, Vector3cd(0, 0, 0), Vector3i(0, 0, 1), zeta0));
  g_0->AddSub(sub_0);
  g_full->AddSub(sub_0);

  VectorXcd zeta1(2); zeta1 << dcomplex(1.0, 0.4), dcomplex(0.4, 0.1);
  SubSymGTOs sub_1(Sub_mono(1, Vector3cd(0, 0, 0), Vector3i(0, 0, 1), zeta1));
  g_1->AddSub(sub_1);
  g_full->AddSub(sub_1);

  g_i->SetUp(); g_1->SetUp(); g_0->SetUp(); g_full->SetUp();
  
  // Set coefficient
  BMat C;
  C[make_pair(0, 0)] = MatrixXcd::Zero(2, 2);
  C[make_pair(0, 0)] <<
    1.0, 1.5,
    0.3, 0.2;
  C[make_pair(1, 1)] = MatrixXcd::Zero(2, 2);
  C[make_pair(1, 1)] <<
    0.2, 0.5,
    1.1, 1.4;

  
  pair<int, int> i00(0, 0), i11(1, 1);

  // compute J/K
  BMat J_full;
  
  ERIMethod method_sym; method_sym.symmetry = 1;
  ERIMethod method; method.coef_R_memo = 1;

  B2EInt eri_full = CalcERI_Complex(g_full, method_sym);
  J_full[i00] = MatrixXcd::Zero(2, 2); J_full[i11] = MatrixXcd::Zero(4, 4);
  AddJK(eri_full, C, 0, 0, 1.0, 0.0, J_full);
  
  BMat J; 
  J[i00] = MatrixXcd::Zero(2, 2); J[i11] = MatrixXcd::Zero(2, 2);
  B2EInt eri_J = CalcERI(g_0, g_1, g_i, g_i, method);
  cout << eri_J->size() << endl;
  VectorXcd C0 = C[make_pair(0, 0)].col(0);
  AddJ(eri_J, C0, 0, J);
  cout << J[i11](0, 0) << J_full[i11](0, 0+2) << endl;

  // Check values
  EXPECT_C_EQ(J[i11](0, 0), J_full[i11](0, 0+2));
  EXPECT_C_EQ(J[i11](0, 1), J_full[i11](0, 1+2));
  EXPECT_C_EQ(J[i11](1, 0), J_full[i11](1, 0+2));
  EXPECT_C_EQ(J[i11](1, 1), J_full[i11](1, 1+2));

}
TEST(Matrix, JK) {

  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::Cs();
  
  // GTO
  SymGTOs gtos(new _SymGTOs); gtos->SetSym(sym);
  SymGTOs gtos_cc = gtos->Clone();
  SymGTOs gtos_full = gtos->Clone();

  VectorXcd zeta1(2); zeta1 << 0.4, 1.0;
  SubSymGTOs sub_s(Sub_mono(0, Vector3cd(0, 0, 0), Vector3i(0, 0, 0), zeta1));
  gtos->AddSub(     sub_s);
  gtos_cc->AddSub(  sub_s);
  gtos_full->AddSub(sub_s);

  VectorXcd zeta2(2); zeta2 << dcomplex(1.0, 0.4), dcomplex(0.4, 0.1);
  SubSymGTOs sub_z(Sub_mono(1, Vector3cd(0, 0, 0), Vector3i(0, 0, 1), zeta2));
  SubSymGTOs sub_zc(Sub_mono(1, Vector3cd(0, 0, 0), Vector3i(0, 0, 1),
			     zeta2.conjugate()));
  gtos->AddSub(   sub_z);
  gtos_cc->AddSub(sub_zc);
  gtos_full->AddSub(sub_z);
  gtos_full->AddSub(sub_zc);

  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 2.0;
  gtos->SetAtoms(xyzq); gtos->SetUp();  
  gtos_cc->SetAtoms(xyzq); gtos_cc->SetUp();  
  gtos_full->SetAtoms(xyzq); gtos_full->SetUp();  
  
  // Set coefficient
  BMat C;
  C[make_pair(0, 0)] = MatrixXcd::Zero(2, 2);
  C[make_pair(0, 0)] <<
    1.0, 1.5,
    0.3, 0.2;
  C[make_pair(1, 1)] = MatrixXcd::Zero(2, 2);
  C[make_pair(1, 1)] <<
    0.2, 0.5,
    1.1, 1.4;

  // compute J/K
  BMat JK, JK_cc, JK_h, JK_full;
  pair<int, int> i00(0, 0), i11(1, 1);
  ERIMethod method; method.symmetry = 1;
  B2EInt eri = CalcERI_Complex(gtos, method);
  B2EInt eri_cc = CalcERI_Complex(gtos_cc, method);
  B2EInt eri_h = CalcERI_Hermite(gtos, method);
  B2EInt eri_full = CalcERI_Complex(gtos_full, method);
  
  JK[i00] = MatrixXcd::Zero(2, 2); JK[i11] = MatrixXcd::Zero(2, 2);
  AddJK(eri, C, 0, 0, 1.1, 1.2, JK);
  
  JK_cc[i00] = MatrixXcd::Zero(2, 2); JK_cc[i11] = MatrixXcd::Zero(2, 2);
  AddJK(eri_cc, C, 0, 0, 1.1, 1.2, JK_cc);

  JK_h[i00] = MatrixXcd::Zero(2, 2); JK_h[i11] = MatrixXcd::Zero(2, 2);
  AddJK(eri_h,    C, 0, 0, 1.1, 1.2, JK_h);
  
  JK_full[i00] = MatrixXcd::Zero(4, 4); JK_full[i11] = MatrixXcd::Zero(4, 4);
  AddJK(eri_full, C, 0, 0, 1.1, 1.2, JK_full);

  // Check values
  EXPECT_C_EQ(JK[i11](0, 0), JK_full[i11](0, 0));
  EXPECT_C_EQ(JK[i11](0, 1), JK_full[i11](0, 1));
  EXPECT_C_EQ(JK[i11](1, 0), JK_full[i11](1, 0));
  EXPECT_C_EQ(JK[i11](1, 1), JK_full[i11](1, 1));

  EXPECT_C_EQ(JK_h[i11](0, 0), JK_full[i11](2, 0));
  EXPECT_C_EQ(JK_h[i11](0, 1), JK_full[i11](2, 1));
  EXPECT_C_EQ(JK_h[i11](1, 0), JK_full[i11](3, 0));
  EXPECT_C_EQ(JK_h[i11](1, 1), JK_full[i11](3, 1));

}
TEST(Matrix, JK2) {

  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::Cs();
  
  // GTO
  SymGTOs gtos(new _SymGTOs);
  gtos->SetSym(sym);

  VectorXcd zeta1(2); zeta1 << 0.4, 1.0;
  SubSymGTOs sub_s(Sub_mono(0, Vector3cd(0, 0, 0), Vector3i(0, 0, 0), zeta1));
  gtos->AddSub(     sub_s);

  VectorXcd zeta2(2); zeta2 << dcomplex(1.0, 0.4), dcomplex(0.4, 0.1);
  SubSymGTOs sub_z(Sub_mono(1, Vector3cd(0, 0, 0), Vector3i(0, 0, 1), zeta2));
  gtos->AddSub(   sub_z);

  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 2.0;
  gtos->SetAtoms(xyzq);
  gtos->SetUp();  

  BMatSet mat = CalcMat_Complex(gtos, false);

  ERIMethod m; m.set_symmetry(1);
  B2EInt eri = CalcERI_Complex(gtos, m);
  
  BMat H_slow = mat->GetBlockMatrix("t");
  BMat H_fast = mat->GetBlockMatrix("t");

  BMat C;
  MatrixXcd c00(2, 2); c00 << 1.1, 1.2, 1.3, 1.4;
  MatrixXcd c11(2, 2); c11 << 2.1, 2.2, 2.3, 2.4;
  C[make_pair(0, 0)] = c00;
  C[make_pair(1, 1)] = c11;

  AddJK(     eri, C, 0, 0, 1.1, 1.2, H_fast);
  AddJK_Slow(eri, C, 0, 0, 1.1, 1.2, H_slow);

  pair<Irrep, Irrep> ii(0, 0);
  pair<Irrep, Irrep> jj(1, 1);
  EXPECT_MATXCD_EQ(H_slow[ii], H_fast[ii]);
  EXPECT_MATXCD_EQ(H_slow[jj], H_fast[jj]);
  
}
TEST(Trans, Slow) {

  cout << "Trans test is comment out now" << endl;

  /*
  IB2EInt *eri_ao = new B2EIntMem();
  IB2EInt *eri_mo = new B2EIntMem();
  
  cout << 0 << endl;

  MO mo(new _MO());
  MatrixXcd C1(1, 1); C1 << 1.2;
  MatrixXcd C2(2, 2); C1 << 2.0, 2.1, 2.2, 2.3;
  mo->C[make_pair(0, 0)] = C1;
  mo->C[make_pair(1, 1)] = C2;

  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::D2h();

  // sub set (S orbital)
  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  VectorXcd zeta_s(1);
  zeta_s << 0.107951;
  sub_s.AddZeta(zeta_s);
  sub_s.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));
  sub_s.SetUp();

  // GTO set
  SymGTOs gtos(sym);
  gtos.AddSub(sub_s);
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 2.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  cout << 1 << endl;
  CalcERI_Complex(gtos, eri_ao, ERIMethod());
  cout << 2 << endl;
  TransformERI_Slow(eri_ao, mo, eri_mo);
  cout << 3 << endl;
  BMat jk_mo;
  CalcJK_MO(eri_ao, mo->C, 0, 0, 0, 0, 1.0, 0.0, &jk_mo);
  
  pair<Irrep, Irrep> II(0, 0);
  EXPECT_C_EQ(jk_mo[II](0, 0),
	      eri_mo->At(0, 0, 0, 0, 0, 0, 0, 0));
  */

}
TEST(HF, first) {

  /*
  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::D2h();
  Irrep irrep_s  = sym->GetIrrep("Ag");
  Irrep irrep_x = sym->GetIrrep("B1u");
  Irrep irrep_y = sym->GetIrrep("B2u");
  Irrep irrep_z = sym->GetIrrep("B3u");

  // sub set (S orbital)
  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  int num_zeta(5);
  VectorXcd zetas(num_zeta);
  for(int i = 0; i < num_zeta; i++) {
    zetas[i] = pow(2.5, num_zeta/2-i);
  }
  sub_s.AddZeta(zetas);
  sub_s.AddRds(Reduction(irrep_s, MatrixXcd::Ones(1, 1)));
  sub_s.SetUp();

  // sub set (P orbital)
  SubSymGTOs sub_p(sym);
  sub_p.AddXyz(Vector3cd(0, 0, 0));
  sub_p.AddNs(Vector3i(1, 0, 0));
  sub_p.AddNs(Vector3i(0, 1, 0));
  sub_p.AddNs(Vector3i(0, 0, 1));
  sub_p.AddZeta(zetas);
  MatrixXcd cx(1, 3); cx << 1, 0, 0; sub_p.AddRds(Reduction(irrep_x, cx));
  MatrixXcd cy(1, 3); cy << 0, 1, 0; sub_p.AddRds(Reduction(irrep_y, cy));
  MatrixXcd cz(1, 3); cz << 0, 0, 1; sub_p.AddRds(Reduction(irrep_z, cz));
  sub_p.SetUp();

  // GTO set
  SymGTOs gtos(sym);
  gtos.AddSub(sub_s);
  gtos.AddSub(sub_p);
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 2.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();
  
  // Compute basic matrix
  int n(gtos.size_basis());
  cout << "size_basis: " << n << endl;
  IB2EInt *eri = new B2EIntMem(pow(n, 4)); 
  BMatSet mat;
  gtos.CalcMat(&mat);
  gtos.CalcERI(eri, 2);

  // -- set non0 symmetry --
  vector<Irrep> irrep_list;
  for(Irrep irrep = 0; irrep < sym->num_class(); irrep++) 
    if(mat.Exist("t", irrep, irrep))
      irrep_list.push_back(irrep);

  // -- initialize --
  BMat H, S, F;   // -- Core hamiltonian, Overlap, Fock --
  BMat coef, rho; // -- MO coefficient, density matrix --
  BVec eig;       // -- orbital energies --
  typedef vector<Irrep>::iterator It;
  for(It it = irrep_list.begin(); it != irrep_list.end(); ++it) {
    Irrep irrep = *it;
    pair<Irrep, Irrep> ii(make_pair(irrep, irrep));
    int n_irrep(gtos.size_basis_isym(irrep));
    coef[ii] = MatrixXcd::Zero(n_irrep, n_irrep);
    rho[ii]  = MatrixXcd::Zero(n_irrep, n_irrep);
    H[ii]  = mat.GetMatrix("t", irrep, irrep);
    H[ii] += mat.GetMatrix("v", irrep, irrep);
    S[ii]  = mat.GetMatrix("s", irrep, irrep);	      
    F[ii] = H[ii];
    eig[irrep] = VectorXcd::Zero(n_irrep);
  }

  int iter_max(10);
  for(int iter = 0; iter < iter_max; iter++) {
    // -- solve eigen value problem --
    for(It it = irrep_list.begin(); it != irrep_list.end(); ++it) {
      pair<Irrep, Irrep> ii(make_pair(*it, *it));
      generalizedComplexEigenSolve(F[ii], S[ii], &coef[ii], &eig[*it]);
    }

    // -- number of occupied orbitals --
    vector<int> num_occ_orb(sym->num_class(), 0);
    num_occ_orb[irrep_s] = 1;

    // -- density matrix --
    for(It it = irrep_list.begin(); it != irrep_list.end(); ++it) {
      pair<Irrep, Irrep> ii(make_pair(*it, *it));
      int n_irrep(gtos.size_basis_isym(*it));
      MatrixXcd& rho_ii = rho[ii];
      for(int i = 0; i < num_occ_orb[*it]; i++) {
	for(int k = 0; k < n_irrep; k++)
	  for(int l = 0; l < n_irrep; l++)
	    rho_ii(k, l) = 2.0 * coef[ii](k, i) * coef[ii](l, i);
      }
    }

    // -- update Fock matrix --
    for(It it = irrep_list.begin(); it != irrep_list.end(); ++it) {
      pair<Irrep, Irrep> ii(make_pair(*it, *it));
      F[ii] = H[ii];
    }
    int ib,jb,kb,lb,i,j,k,l,t;
    dcomplex v;
    eri->Reset();
    while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
      if(ib == jb && kb == lb) {
	pair<Irrep, Irrep> ii(make_pair(ib, jb));
	pair<Irrep, Irrep> kk(make_pair(kb, lb));
	F[ii](i, j) += rho[kk](k, l) * v;
      }
      if(ib == lb && kb == jb) {
	pair<Irrep, Irrep> il(make_pair(ib, lb));
	pair<Irrep, Irrep> jk(make_pair(jb, kb));
	F[il](i, l) -= 0.5 * rho[jk](k, j) * v;
      }
    }

    // -- convergence check --
    
    // -- calculate total energy --
    dcomplex ene(0);
    for(It it = irrep_list.begin(); it != irrep_list.end(); ++it) {
      pair<Irrep, Irrep> ii(make_pair(*it, *it));
      MatrixXcd& rho_ii(rho[ii]);
      MatrixXcd& H_ii(H[ii]);
      MatrixXcd& F_ii(F[ii]);
      int n(gtos.size_basis_isym(*it));
      for(int i = 0; i < n; i++)
	for(int j = 0; j < n; j++) 
	  ene += 0.5 * rho_ii(i, j) * (H_ii(i, j) + F_ii(i, j));
    }

    // -- finalize --
    cout << ene;
    for(It it = irrep_list.begin(); it != irrep_list.end(); ++it) {
      for(int i = 0; i < num_occ_orb[*it]; i++) {
	cout << eig[*it][i];
      }
    }
    cout << endl;
  }

  */

}
TEST(HF, OccupiedNumber) {

  BVec eigs;
  
  VectorXcd e1(3); e1 << dcomplex(0.01, 0.3), dcomplex(1.1, 0.1), 2.1;
  VectorXcd e2(4); e2 << dcomplex(0.12, 0.3), 10.0, 11.0, 12.0;
  
  eigs[0] = e1;
  eigs[2] = e2;
  
  vector<int> num_occ = CalcOccNum(eigs, 3, 3);
  
  EXPECT_EQ(2, num_occ[0]);
  EXPECT_EQ(0, num_occ[1]);
  EXPECT_EQ(1, num_occ[2]);

}
TEST(HF, He) {

  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::D2h();

  // sub set (S orbital)
  SubSymGTOs sub_s;
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
    VectorXcd zeta_s(10);
  zeta_s << 0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053, 30.17990, 108.7723, 488.8941, 3293.694;
  sub_s.AddZeta(zeta_s);
  sub_s.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));

  // GTO set
  SymGTOs gtos(new _SymGTOs);
  gtos->SetSym(sym);
  gtos->AddSub(sub_s);
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 2.0;
  gtos->SetAtoms(xyzq);
  gtos->SetUp();

  bool conv;
  MO mo = CalcRHF(gtos, 2, 10, 0.0000001, &conv);
  EXPECT_TRUE(conv);
  EXPECT_C_NEAR(dcomplex(-2.8617,0.0), mo->energy,
		0.0001);

  // Coef check
  VectorXcd c = mo->C[make_pair(0, 0)].col(0);
  double eps(0.00001);
  EXPECT_C_NEAR(-0.91795 , mo->eigs[0][0], eps);

  EXPECT_C_NEAR(-0.05243, c(0), eps);
  EXPECT_C_NEAR(-0.24887, c(1), eps);
  EXPECT_C_NEAR(-0.36002, c(2), eps);
  EXPECT_C_NEAR(-0.28403, c(3), eps);
  EXPECT_C_NEAR(-0.14909, c(4), eps);
  EXPECT_C_NEAR(-0.05709, c(5), eps);
  EXPECT_C_NEAR(-0.01721, c(6), eps);
  EXPECT_C_NEAR(-0.00412, c(7), eps);
  EXPECT_C_NEAR(-0.00076, c(8), eps);
  EXPECT_C_NEAR(-0.00010, c(9), eps);

}
TEST(HF, H2) {

  // structure
  dcomplex R0(1.4);

  pSymmetryGroup sym = SymmetryGroup::C1();
  Irrep irrep_s  = 0;

  // see calculation result in ylcls
  SubSymGTOs sub_s;
  sub_s.AddXyz(Vector3cd(0, 0,+R0/2.0));
  sub_s.AddXyz(Vector3cd(0, 0,-R0/2.0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  int num_zeta(6);
  VectorXcd zetas(num_zeta);
  zetas << 1.336, 2.013, 0.4538, 0.1233, 0.0411, 0.0137;
  sub_s.AddZeta(zetas); 
  MatrixXcd c(2, 1); c << 1, 1;
  sub_s.AddRds(Reduction(irrep_s, c));
  
  SubSymGTOs sub_p;
  sub_p.AddXyz(Vector3cd(0, 0, +R0/2.0));
  sub_p.AddXyz(Vector3cd(0, 0, -R0/2.0));
  sub_p.AddNs(Vector3i(  0, 0, 1));
  VectorXcd zeta_p(1); zeta_p << 1.1;
  sub_p.AddZeta(zeta_p);
  MatrixXcd cp(2, 1); cp << 1, -1;
  sub_p.AddRds(Reduction(irrep_s, cp));
  
  SymGTOs gtos(new _SymGTOs);
  gtos->SetSym(sym);
  gtos->AddSub(sub_s); 
  gtos->AddSub(sub_p);
  MatrixXcd xyzq(4, 2); xyzq <<
			  0.0,      0.0,
			  0.0,      0.0,
			  +R0/2.0, -R0/2.0,
			  1.0,      1.0;
  gtos->SetAtoms(xyzq);
  gtos->SetUp();

  bool conv;
  double eps(pow(10.0, -5.0));
  BMatSet mat_set = CalcMat_Complex(gtos, true);
  ERIMethod method; method.symmetry = 1;
  B2EInt eri = CalcERI_Complex(gtos, method);

  MO mo = CalcRHF(sym, mat_set, eri, 2, 50, eps, &conv, 0);
  EXPECT_TRUE(conv);
  EXPECT_C_NEAR(mo->energy + 1.0/R0, -1.11881323240, 0.00001);
  EXPECT_C_NEAR(mo->eigs[0](0), -0.59012, 0.00002);
  EXPECT_C_NEAR(mo->eigs[0](1), +0.02339, 0.00002);
}
TEST(HF, H2_eig) {

/*
  // ==== Symmetry ====
  pSymmetryGroup D2h = SymmetryGroup::D2h();

  // ==== Sub ====
  SubSymGTOs sub1(D2h);
  sub1.AddXyz(Vector3cd(0, 0, +0.7));
  sub1.AddXyz(Vector3cd(0, 0, -0.7));
  sub1.AddNs( Vector3i( 0, 0, 0));
  //  VectorXcd z1(4); z1 << 2.013, 0.1233, 0.0411, 0.0137; sub1.AddZeta(z1);
  //  VectorXcd z1(2); z1 << 2.013, 0.0411; sub1.AddZeta(z1);
  VectorXcd z1(2); z1 << 0.1233, 0.0411; sub1.AddZeta(z1); // raise Error
  // VectorXcd z1(1); z1 << 0.1233; sub1.AddZeta(z1);
  //  VectorXcd z1(1); z1 << 0.0411; sub1.AddZeta(z1);
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

  // ==== RHF =====
  BMatSet mat_set;  
  IB2EInt *eri = new B2EIntMem(pow(gtos.size_basis(), 4));
  gtos.CalcMat(&mat_set);
  gtos.CalcERI(eri, 1);
  
  bool conv;
  MO mo = CalcRHF(gtos, 2, 10, 0.0000001, &conv);
  EXPECT_TRUE(conv);
  
  // copied from ~/calc/ccolumbus/
  // EIGENVALUES   -1.162672  0.001574   -0.351197  0.000572   -0.170858  0.000287   -0.100557 -0.001012    0.420850 -0.009795
  EXPECT_C_EQ(dcomplex(-1.162672, 0.001574), mo->eigs[0](0));
  EXPECT_C_EQ(dcomplex(-0.351197, 0.000572), mo->eigs[0](1));  
*/


/*
  // ==== Symmetry ====
  pSymmetryGroup D2h = SymmetryGroup::D2h();

  // ==== Sub ====
  SubSymGTOs sub1(D2h);
  sub1.AddXyz(Vector3cd(0, 0, +0.7));
  sub1.AddXyz(Vector3cd(0, 0, -0.7));
  sub1.AddNs( Vector3i( 0, 0, 0));
  VectorXcd z1(4); z1 << 2.013, 0.1233, 0.0411, 0.0137; sub1.AddZeta(z1);
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

  // ==== RHF =====
  bool conv;
  MO mo = CalcRHF(gtos, 2, 10, 0.0000001, &conv);
  EXPECT_TRUE(conv);
  
  // copied from ~/calc/ccolumbus/
  // EIGENVALUES   -1.162672  0.001574   -0.351197  0.000572   -0.170858  0.000287   -0.100557 -0.001012    0.420850 -0.009795
  EXPECT_C_EQ(dcomplex(-1.162672, 0.001574), mo->eigs[0](0));
  EXPECT_C_EQ(dcomplex(-0.351197, 0.000572), mo->eigs[0](1));
  */  

}
TEST(HF, H2O) {

  /*

  // Symmetry
  pSymmetryGroup C2v = SymmetryGroup::C2v();

  // Molecule
  dcomplex H_z(0.2);
  dcomplex H_x(1.2);

  // sub set (on H)
  SubSymGTOs sub_h(C2v);
  sub_h.AddXyz(Vector3cd(+H_x, 0, H_z));
  sub_h.AddXyz(Vector3cd(-H_x, 0, H_z));
  sub_h.AddNs(Vector3i(0, 0, 0));
  // -- uncontracted 3-21
  //VectorXcd zs_h(3); zs_h << 5.4471780, 0.8245470, 0.1831920;
  VectorXcd zs_h(1); zs_h << 0.8245470;
  sub_h.AddZeta(zs_h);
  MatrixXcd ch1(2,1); ch1<<1,1;  sub_h.AddRds(Reduction(C2v->irrep_s, ch1));
  MatrixXcd ch2(2,1); ch2<<1,-1; sub_h.AddRds(Reduction(C2v->irrep_x, ch2));
  sub_h.SetUp();

  // sub set (on O, s only)
  SubSymGTOs sub_o_1(C2v);
  sub_o_1.AddXyz(Vector3cd::Zero());
  sub_o_1.AddNs(Vector3i::Zero());  
  //VectorXcd zs_o_1(3); zs_o_1 << 322.0370000, 48.4308000, 10.4206000;
  VectorXcd zs_o_1(1); zs_o_1 << 48.4308000;
  sub_o_1.AddZeta(zs_o_1);
  sub_o_1.AddRds(Reduction(C2v->irrep_s, MatrixXcd::Ones(1, 1)));
  sub_o_1.SetUp();

  // sub set (on O, sp)
  SubSymGTOs sub_o_2(C2v);
  sub_o_2.AddXyz(Vector3cd::Zero());
  sub_o_2.AddNs(Vector3i(0, 0, 0));
  sub_o_2.AddNs(Vector3i(1, 0, 0));
  sub_o_2.AddNs(Vector3i(0, 1, 0));
  sub_o_2.AddNs(Vector3i(0, 0, 1));
  //VectorXcd zs_o_2(3); zs_o_2 << 7.4029400, 1.5762000, 0.3736840;
  VectorXcd zs_o_2(1); zs_o_2 << 1.5762000;
  sub_o_2.AddZeta(zs_o_2);
  MatrixXcd co1(1,4); co1 << 1,0,0,0; sub_o_2.AddRds(Reduction(C2v->irrep_s, co1));
  MatrixXcd co2(1,4); co2 << 0,1,0,0; sub_o_2.AddRds(Reduction(C2v->irrep_x, co2));
  MatrixXcd co3(1,4); co3 << 0,0,1,0; sub_o_2.AddRds(Reduction(C2v->irrep_y, co3));
  MatrixXcd co4(1,4); co4 << 0,0,0,1; sub_o_2.AddRds(Reduction(C2v->irrep_z, co4));
  sub_o_2.SetUp();

  // GTO set
  SymGTOs gtos(C2v);
  gtos.AddSub(sub_h); gtos.AddSub(sub_o_1); gtos.AddSub(sub_o_2);
  MatrixXcd xyzq(4, 3); xyzq <<
			  0, +H_x, -H_x,
			  0, 0,     0,
			  0, +H_z, -H_z,
			  8, 1,     1;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  //cout << gtos.str() << endl;

  // Calculation
  bool conv;
  MO mo = CalcRHF(gtos, 10, 20, 0.000001, &conv);
  EXPECT_TRUE(conv);
  cout << mo->energy << endl;
  */

}
TEST(HF, H) {

  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::D2h();

  // sub set (S orbital)
  SubSymGTOs sub_s;
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  int num_zeta(10);
  VectorXcd zetas(num_zeta);
  for(int i = 0; i < num_zeta; i++) {
    zetas[i] = pow(2.5, 5-i);
  }
  sub_s.AddZeta(zetas);
  sub_s.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));

  // sub set (Pz orbital)
  SubSymGTOs sub_pz;
  sub_pz.AddXyz(Vector3cd(0, 0, 0));
  sub_pz.AddNs(Vector3i(0, 0, 1));
  int num_zeta_pz(12);
  VectorXcd zetas_pz(num_zeta_pz);
  double theta(30.0*M_PI/180.0);
  for(int i = 0; i < num_zeta_pz; i++)
    zetas_pz[i] = pow(2.0, num_zeta_pz/2-i) * dcomplex(cos(theta), -sin(theta));
  sub_pz.AddZeta(zetas_pz);
  sub_pz.AddRds(Reduction(sym->irrep_z, MatrixXcd::Ones(1, 1)));

  // GTO set
  SymGTOs gtos(new _SymGTOs);
  gtos->SetSym(sym);
  gtos->AddSub(sub_s);
  gtos->AddSub(sub_pz);
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos->SetAtoms(xyzq);
  gtos->SetUp();

  // Calculate matrix
  BMatSet bmat_set = CalcMat_Complex(gtos, true);
  
  // Solve eigen value problem
  MO mo = CalcOneEle(sym, bmat_set);;
  EXPECT_C_NEAR(-0.5, mo->energy, pow(10.0, -4));

  //  cout << "z >>> " << endl;
  //  cout << bmat_set.GetMatrix("z", sym->irrep_z, sym->irrep_s) << endl;
  //  cout << "z <<< " << endl;

  // Solve driven type equation
  double w10(1.0), w14(1.4);
  dcomplex alpha10 = CalcAlpha(mo, bmat_set, sym->irrep_s, 0, mo->H, w10, CoordZ);
  dcomplex alpha14 = CalcAlpha(mo, bmat_set, sym->irrep_s, 0, mo->H, w14, CoordZ);
  double cs10 = PITotalCrossSection(alpha10, w10, 1);
  double cs14 = PITotalCrossSection(alpha14, w14, 1);
  cout << "w : calc : exact" << endl;
  cout << w10 << ": " << cs10 << " : " << 0.954 << endl ;
  cout << w14 << ": " << cs14 << " : " << 0.353 << endl;
}
TEST(STEX, small) {
  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::D2h();

  // sub set (S orbital)
  SubSymGTOs sub_s;
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  VectorXcd zeta_s(2);
  zeta_s << 0.1, 0.5;
  sub_s.AddZeta(zeta_s);
  sub_s.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));

  // sub set (P orbital)
  SubSymGTOs sub_z;
  sub_z.AddXyz(Vector3cd(0, 0, 0));
  sub_z.AddNs(Vector3i(0, 0, 1));
  VectorXcd zeta_z(2);
  zeta_z << dcomplex(1.0, 0.1), dcomplex(3.0, 0.2);
  sub_z.AddZeta(zeta_z);
  sub_z.AddRds(Reduction(sym->irrep_z, MatrixXcd::Ones(1, 1)));

  // GTO set
  SymGTOs gtos(new _SymGTOs);
  gtos->SetSym(sym);
  gtos->AddSub(sub_s);
  gtos->AddSub(sub_z);
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 2.0;
  gtos->SetAtoms(xyzq);
  gtos->SetUp();

  // compute basic matrix
  BMatSet mat_set = CalcMat_Complex(gtos, true);
  ERIMethod method; method.symmetry = 1;
  B2EInt eri = CalcERI_Complex(gtos, method);
  
  // RHF
  MO mo = CalcOneEle(sym, mat_set, 0);

  // build Static Exchange Hamiltonian
  BMat hmat0, hmat1;
  CalcSEHamiltonian(mo, eri, 0, 0, &hmat0, 0);
  CalcSEHamiltonian(mo, eri, 0, 0, &hmat1, 1);
  int i(sym->irrep_z);
  cout << hmat0[make_pair(i, i)] << endl;
  cout << "-------" << endl;
  cout << hmat1[make_pair(i, i)] << endl;

  MatrixXcd& S00 = mo->S[make_pair(0, 0)];
  VectorXcd C0 = mo->C[make_pair(0, 0)].col(0);
  dcomplex cumsum(0);
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < 2; j++)
      cumsum += S00(i, j) * C0(j) * C0(i);
  cout << "1? " << cumsum << endl;

}
TEST(STEX, He) {

  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::D2h();

  // sub set (S orbital)
  SubSymGTOs sub_s;
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  VectorXcd zeta_s(10);
  zeta_s << 0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053, 30.17990, 108.7723, 488.8941, 3293.694;
  sub_s.AddZeta(zeta_s);
  sub_s.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));

  // sub set (P orbital)
  SubSymGTOs sub_z;
  sub_z.AddXyz(Vector3cd(0, 0, 0));
  sub_z.AddNs(Vector3i(0, 0, 1));
  int num_zeta(19);
  VectorXcd zetas(num_zeta);
  zetas << 
    dcomplex(0.00256226, -0.01559939),
    dcomplex(0.00389597, -0.02240632),
    dcomplex(0.01229986, -0.03080238),
    dcomplex(0.03010506, -0.04378147),
    0.01,
    0.0177160054,
    0.0313856847,
    0.0556028960,
    0.0985061205,
    0.174513496,
    0.309168204,
    0.547722558,
    0.970345578,
    1.71906475,
    3.04549604,
    5.39540243,
    9.55849785,
    16.9338400,
    30.0000000;

  sub_z.AddZeta(zetas);
  sub_z.AddRds(Reduction(sym->irrep_z, MatrixXcd::Ones(1, 1)));

  // GTO set
  SymGTOs gtos(new _SymGTOs);
  gtos->SetSym(sym);
  gtos->AddSub(sub_s);
  gtos->AddSub(sub_z);
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 2.0;
  gtos->SetAtoms(xyzq);
  gtos->SetUp();

  // compute basic matrix
  BMatSet mat_set = CalcMat_Complex(gtos, true);
  ERIMethod method; method.symmetry = 1;
  B2EInt eri;
  
  //  string fn("eri_hf.bin");
  //  eri= CalcERI_Complex(gtos, method);
  //  eri->Write(fn);
  string fn("eri_hf.bin");
  FILE *fp;  
  fp = fopen(fn.c_str(), "r");
  if(fp == NULL) {
    // If ERI file is not exist:
    cout << "calculate ERI" << endl;
    eri= CalcERI_Complex(gtos, method);
    eri->Write(fn);
  } else {
    // If ERI file is exist:
    cout << "find eri.bin" << endl;
    eri = ERIRead(fn);
    fclose(fp);
  }

  cout << 1 << endl;
  // RHF
  bool conv;
  MO mo = CalcRHF(sym, mat_set, eri, 2, 10, 0.0000001, &conv);
  EXPECT_TRUE(conv);
  cout << 2 << endl;

  // build Static Exchange Hamiltonian
  BMat hmat;
  CalcSEHamiltonian(mo, eri, 0, 0, &hmat, 0);
  cout << 3 << endl;
  
  // Compute alpha
  double au2ev(27.2114);
  double w25(25/au2ev);
  double w35(35/au2ev);
  double w45(45/au2ev);
  
  dcomplex alpha25 = CalcAlpha(mo, mat_set, 0, 0, hmat, w25, CoordZ);
  dcomplex alpha35 = CalcAlpha(mo, mat_set, 0, 0, hmat, w35, CoordZ);
  dcomplex alpha45 = CalcAlpha(mo, mat_set, 0, 0, hmat, w45, CoordZ);

  cout << PITotalCrossSection(alpha25, w25, 2) << " " << 7.21 << endl;
  cout << PITotalCrossSection(alpha35, w35, 2) << " " << 4.09 << endl;
  cout << PITotalCrossSection(alpha45, w45, 2) << " " << 2.48 << endl;
  cout << 4 << endl;
  
  /*
  EXPECT_C_NEAR(7.21,
  PITotalCrossSection(alpha25, w25, 2),
		0.01);
  EXPECT_C_NEAR(4.09,
  PITotalCrossSection(alpha35, w35, 2),
		0.01);
  EXPECT_C_NEAR(2.48,
		PITotalCrossSection(alpha45, w45, 2),
		0.01);
		*/
  
  //  [[ 25.     7.21]]
  //  [[ 35.     4.09]]
  //  [[ 45.     2.48]]
}
TEST(STEX, He_one_gto) {

  dcomplex z1(0.09154356, -0.24865707);
  double w_ev = 24.98;
  double sig(5.53533505);
  //  dcomplex z1(0.42258809, -0.93849382);
  //  double w_ev = 99.92;
  //  double sig(3.66579556);
  
  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::Cs();

  // sub set (S orbital)
  SubSymGTOs sub_s;
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  VectorXcd zeta_s(10);
  zeta_s << 0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053, 30.17990, 108.7723, 488.8941, 3293.694;
  sub_s.AddZeta(zeta_s);
  sub_s.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));

  // sub set (P orbital)
  SubSymGTOs sub_z;
  sub_z.AddXyz(Vector3cd(0, 0, 0));
  sub_z.AddNs(Vector3i(0, 0, 1));
  VectorXcd zeta_z(1); zeta_z << z1;
  sub_z.AddZeta(zeta_z);
  sub_z.AddRds(Reduction(sym->irrep_z, MatrixXcd::Ones(1, 1)));

  // GTO set
  SymGTOs gtos(new _SymGTOs);
  gtos->SetSym(sym);
  gtos->AddSub(sub_s);
  gtos->AddSub(sub_z);
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 2.0;
  gtos->SetAtoms(xyzq);
  gtos->SetUp();

  // compute basic matrix
  BMatSet mat_set = CalcMat_Complex(gtos, true);
  ERIMethod method; method.symmetry = 1;
  B2EInt eri = CalcERI_Complex(gtos, method);

  // RHF
  bool conv;
  MO mo = CalcRHF(sym, mat_set, eri, 2, 20, 0.0000001, &conv);
  EXPECT_TRUE(conv);
  
  // coef check
  VectorXcd c = mo->C[make_pair(0, 0)].col(0);
  double eps(0.00001);
  EXPECT_C_NEAR(-0.91795 , mo->eigs[0][0], eps);

  EXPECT_C_NEAR(-0.05243, c(0), eps);
  EXPECT_C_NEAR(-0.24887, c(1), eps);
  EXPECT_C_NEAR(-0.36002, c(2), eps);
  EXPECT_C_NEAR(-0.28403, c(3), eps);
  EXPECT_C_NEAR(-0.14909, c(4), eps);
  EXPECT_C_NEAR(-0.05709, c(5), eps);
  EXPECT_C_NEAR(-0.01721, c(6), eps);
  EXPECT_C_NEAR(-0.00412, c(7), eps);
  EXPECT_C_NEAR(-0.00076, c(8), eps);
  EXPECT_C_NEAR(-0.00010, c(9), eps);
  
  // build Static Exchange Hamiltonian
  BMat hmat0, hmat2, hmat3;
  CalcSEHamiltonian(mo, eri, 0, 0, &hmat0, 0);
  //  CalcSEHamiltonian(mo, eri, 0, 0, &hmat1, 1);
  CalcSEHamiltonian(mo, eri, 0, 0, &hmat2, 2);
  CalcSEHamiltonian(mo, eri, 0, 0, &hmat3, 3);
  
  // Compute alpha
  double au2ev(27.2114);
  double w = w_ev/au2ev;

  cout << PITotalCrossSection(CalcAlpha(mo, mat_set, 0, 0, hmat0, w, CoordZ, 0), w, 2)
       << endl;
  cout << PITotalCrossSection(CalcAlpha(mo, mat_set, 0, 0, hmat0, w, CoordZ, 1), w, 2)
       << endl;
  cout << PITotalCrossSection(CalcAlpha(mo, mat_set, 0, 0, hmat0, w, CoordZ, 2), w, 2)
       << endl;

  cout << sig << endl;
	   
}
SymGTOs CreateGTOs_H2() {

  double R0 = 1.4;
  pSymmetryGroup sym = SymmetryGroup::D2h();
  SymGTOs gtos(new _SymGTOs);
  gtos->SetSym(sym);

  SubSymGTOs sub_s;
  sub_s.AddXyz(Vector3cd(0, 0,+R0/2.0));
  sub_s.AddXyz(Vector3cd(0, 0,-R0/2.0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  int num_zeta(6);
  VectorXcd zetas(num_zeta);
  zetas << 1.336, 2.013, 0.4538, 0.1233, 0.0411, 0.0137;
  sub_s.AddZeta(zetas); 
  MatrixXcd c(2, 1); c << 1, 1;
  sub_s.AddRds(Reduction(sym->irrep_s, c));
  gtos->AddSub(sub_s); 

  SubSymGTOs sub_p;
  sub_p.AddXyz(Vector3cd(0, 0, +R0/2.0));
  sub_p.AddXyz(Vector3cd(0, 0, -R0/2.0));
  sub_p.AddNs(Vector3i(  0, 0, 1));
  VectorXcd zeta_p(1); zeta_p << 1.1;
  sub_p.AddZeta(zeta_p);
  MatrixXcd cp(2, 1); cp << 1, -1;
  sub_p.AddRds(Reduction(sym->irrep_s, cp));
  gtos->AddSub(sub_p);

  SubSymGTOs sub_p_cen;
  sub_p_cen.AddXyz(Vector3cd(0, 0, 0));
  sub_p_cen.AddNs( Vector3i( 1, 0, 0));
  sub_p_cen.AddNs( Vector3i( 0, 1, 0));
  sub_p_cen.AddNs( Vector3i( 0, 0, 1));
  VectorXcd zeta_cen(19);
  zeta_cen << 
    dcomplex(0.00256226, -0.01559939),
    dcomplex(0.00389597, -0.02240632),
    dcomplex(0.01229986, -0.03080238),
    dcomplex(0.03010506, -0.04378147),
    0.01,
    0.0177160054,
    0.0313856847,
    0.0556028960,
    0.0985061205,
    0.174513496,
    0.309168204,
    0.547722558,
    0.970345578,
    1.71906475,
    3.04549604,
    5.39540243,
    9.55849785,
    16.9338400,
    30.0000000;
  sub_p_cen.AddZeta(zeta_cen);
  MatrixXcd cx(1, 3); cx << 1, 0, 0;
  MatrixXcd cy(1, 3); cy << 0, 1, 0;
  MatrixXcd cz(1, 3); cz << 0, 0, 1;
  sub_p_cen.AddRds(Reduction(sym->irrep_x, cx));
  sub_p_cen.AddRds(Reduction(sym->irrep_y, cy));
  sub_p_cen.AddRds(Reduction(sym->irrep_z, cz));
  gtos->AddSub(sub_p_cen);

  MatrixXcd xyzq(4, 2); xyzq <<
			  0.0,      0.0,
			  0.0,      0.0,
			  +R0/2.0, -R0/2.0,
			  1.0,      1.0;
  gtos->SetAtoms(xyzq);
  gtos->SetUp();

  return gtos;

}
TEST(STEX, H2) {

  Timer timer;

  SymGTOs gtos = CreateGTOs_H2();
  double R0 = 1.4;

  bool conv;
  double eps(pow(10.0, -5.0));
  BMatSet mat_set = CalcMat_Complex(gtos, true);
  ERIMethod method; method.symmetry = 1;
  B2EInt eri;

  timer.Start("ERI");
  eri = CalcERI_Complex(gtos, method);
  timer.End("ERI");

  pSymmetryGroup sym = gtos->sym_group;
  MO mo;
  timer.Start("HF");
  mo = CalcRHF(sym, mat_set, eri, 2, 50, eps, &conv, 0);
  timer.End("HF");

  EXPECT_TRUE(conv);
  EXPECT_C_NEAR(mo->energy + 1.0/R0, -1.11881323240, 0.00001);
  EXPECT_C_NEAR(mo->eigs[0](0), -0.59012, 0.00004);
  EXPECT_C_NEAR(mo->eigs[0](1), +0.02339, 0.00004);

  // Build Static exchange Hamiltonian
  timer.Start("STEX");
  BMat hmat;
  CalcSEHamiltonian(mo, eri, 0, 0, &hmat);
  
  // Compute alpha
  double au2ev(27.2114);
  double w_ev = 20.0;
  double w = w_ev/au2ev;

  dcomplex alpha_x = CalcAlpha(mo, mat_set, 0, 0, hmat, w, CoordX, 0);
  dcomplex alpha_y = CalcAlpha(mo, mat_set, 0, 0, hmat, w, CoordY, 0);
  dcomplex alpha_z = CalcAlpha(mo, mat_set, 0, 0, hmat, w, CoordZ, 0);
  cout << PITotalCrossSection(alpha_x + alpha_y + alpha_z, w, 2)/3.0 << endl;
  timer.End("STEX");
  
  timer.Display();
  
}

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
