#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

#include <gtest/gtest.h>
#include <Eigen/Core>
#include "gtest_plus.hpp"

#include "symmolint.hpp"
#include "b2eint.hpp"
#include "symgroup.hpp"
#include "eigen_plus.hpp"
#include "mo.hpp"

using namespace std;
using namespace l2func;
using namespace Eigen;

TEST(Matrix, HeSmall) {
  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::D2h();

  // sub set (S orbital)
  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  VectorXcd zeta_s(2);
  zeta_s << 0.5, 1.5;
  sub_s.AddZeta(zeta_s);
  sub_s.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));
  sub_s.SetUp();

  // sub set (P orbital)
  SubSymGTOs sub_z(sym);
  sub_z.AddXyz(Vector3cd(0, 0, 0));
  sub_z.AddNs(Vector3i(0, 0, 1));
  VectorXcd zetas(2); zetas <<
			dcomplex(0.0025, -0.016),
			dcomplex(0.004,  -0.022),


  sub_z.AddZeta(zetas);
  sub_z.AddRds(Reduction(sym->irrep_z, MatrixXcd::Ones(1, 1)));
  sub_z.SetUp();

  // GTO set
  SymGTOs gtos(sym);
  gtos.AddSub(sub_z);
  gtos.AddSub(sub_s);
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 2.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  // compute basic matrix
  BMatSet mat_set;  
  gtos.CalcMat(&mat_set);
  IB2EInt *eri;
  eri = new B2EIntMem(pow(gtos.size_basis(), 4));
  gtos.CalcERI(eri, 1);

  cout << "S(00)" << endl;
  cout << mat_set.GetMatrix("s", 0, 0)<< endl;
  cout << "S(11)" << endl;
  cout << mat_set.GetMatrix("s", sym->irrep_z, sym->irrep_z)<< endl;

  cout << "T(00)" << endl;
  cout << mat_set.GetMatrix("t", 0, 0)<< endl;
  cout << "T(11)" << endl;
  cout << mat_set.GetMatrix("t", sym->irrep_z, sym->irrep_z)<< endl;

  cout << "V(00)" << endl;
  cout << mat_set.GetMatrix("v", 0, 0)<< endl;
  cout << "V(11)" << endl;
  cout << mat_set.GetMatrix("v", sym->irrep_z, sym->irrep_z)<< endl;

  int ib,jb,kb,lb,i,j,k,l,t;
  dcomplex v;
  eri->Reset();
    while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
      cout <<  ib << jb << kb << lb << i << j << k << l << v << endl;
    }
  
  delete eri;
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
  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
    VectorXcd zeta_s(10);
  zeta_s << 0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053, 30.17990, 108.7723, 488.8941, 3293.694;
  sub_s.AddZeta(zeta_s);
  sub_s.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));
  sub_s.SetUp();

  // GTO set
  SymGTOs gtos(sym);
  gtos.AddSub(sub_s);
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 2.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  bool conv;
  MO mo = CalcRHF(gtos, 2, 10, 0.0000001, &conv);
  EXPECT_TRUE(conv);
  EXPECT_C_NEAR(dcomplex(-2.8617,0.0), mo->energy,
		0.0001);
}
TEST(HF, H2) {

  // structure
  dcomplex R0(1.4);
  //dcomplex R0(0.0);

  // set symmetry
  /*
  pSymmetryGroup sym = SymmetryGroup::D2h();
  Irrep irrep_s  = sym->GetIrrep("Ag");

  // sub set (S orbital)
  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0,+R0/2.0));
  sub_s.AddXyz(Vector3cd(0, 0,-R0/2.0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  int num_zeta(6);
  VectorXcd zetas(num_zeta);
  zetas << 1.336, 2.013, 0.4538, 0.1233, 0.0411, 0.0137;
  sub_s.AddZeta(zetas);
  sub_s.AddRds(Reduction(irrep_s, MatrixXcd::Ones(2, 1)));
  sub_s.SetUp();

  // GTO set
  SymGTOs gtos(sym);
  gtos.AddSub(sub_s);
  //  gtos.AddSub(sub_p);
  */

  pSymmetryGroup sym = SymmetryGroup::C1();
  Irrep irrep_s  = 0;

  // sub set (S orbital)
  SubSymGTOs sub_s(sym), sub_s2(sym);
  sub_s.AddXyz(Vector3cd(0, 0,+R0/2.0)); sub_s2.AddXyz(Vector3cd(0, 0,-R0/2.0));
  sub_s.AddNs(Vector3i(0, 0, 0));        sub_s2.AddNs(Vector3i(0, 0, 0));
  int num_zeta(6);
  VectorXcd zetas(num_zeta);
  zetas << 1.336, 2.013, 0.4538, 0.1233, 0.0411, 0.0137;
  sub_s.AddZeta(zetas); sub_s2.AddZeta(zetas);
  sub_s.AddRds(Reduction(irrep_s, MatrixXcd::Ones(1, 1)));
  sub_s2.AddRds(Reduction(irrep_s, MatrixXcd::Ones(1, 1)));
  sub_s.SetUp(); sub_s2.SetUp();

  // GTO set
  SymGTOs gtos(sym);
  gtos.AddSub(sub_s); gtos.AddSub(sub_s2);
  MatrixXcd xyzq(4, 2); xyzq <<
			  0.0,      0.0,
			  0.0,      0.0,
			  +R0/2.0, -R0/2.0,
			  1.0,      1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  bool conv;
  double eps(pow(10.0, -7.0));
  int num(gtos.size_basis());
  BMatSet mat_set; gtos.CalcMat(&mat_set);
  IB2EInt *eri = new B2EIntMem(pow(num, 4)); gtos.CalcERI(eri, 1);

  //  MatrixXcd H = mat_set.GetMatrix("v", 0, 0) + mat_set.GetMatrix("t", 0, 0);
  //  MatrixXcd S = mat_set.GetMatrix("s", 0, 0);
  //  MatrixXcd C; VectorXcd eig;
  //  cout << H << endl;
  //  generalizedComplexEigenSolve(H, S, &C, &eig);
  //  cout << eig << endl;

  MO mo = CalcRHF(sym, mat_set, eri, 2, 30, eps, &conv, 0);
  EXPECT_TRUE(conv);
  cout << mo->energy << endl;
  cout << "Eigs: " << endl;
  cout << mo->eigs[0] << endl;

}
TEST(HF, H2_eig) {

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
  //  BMatSet mat_set;  
  IB2EInt *eri = new B2EIntMem(pow(gtos.size_basis(), 4));
  //  gtos.CalcMat(&mat_set);
  gtos.CalcERI(eri, 1);
  
  //  bool conv;
  //  MO mo = CalcRHF(gtos, 2, 10, 0.0000001, &conv);
  //  EXPECT_TRUE(conv);
  
  // copied from ~/calc/ccolumbus/
  // EIGENVALUES   -1.162672  0.001574   -0.351197  0.000572   -0.170858  0.000287   -0.100557 -0.001012    0.420850 -0.009795
  //  EXPECT_C_EQ(dcomplex(-1.162672, 0.001574), mo->eigs[0](0));
  //  EXPECT_C_EQ(dcomplex(-0.351197, 0.000572), mo->eigs[0](1));  

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
  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  int num_zeta(10);
  VectorXcd zetas(num_zeta);
  for(int i = 0; i < num_zeta; i++) {
    zetas[i] = pow(2.5, 5-i);
  }
  sub_s.AddZeta(zetas);
  sub_s.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));
  sub_s.SetUp();

  // sub set (Pz orbital)
  SubSymGTOs sub_pz(sym);
  sub_pz.AddXyz(Vector3cd(0, 0, 0));
  sub_pz.AddNs(Vector3i(0, 0, 1));
  int num_zeta_pz(12);
  VectorXcd zetas_pz(num_zeta_pz);
  double theta(30.0*M_PI/180.0);
  for(int i = 0; i < num_zeta_pz; i++)
    zetas_pz[i] = pow(2.0, num_zeta_pz/2-i) * dcomplex(cos(theta), -sin(theta));
  sub_pz.AddZeta(zetas_pz);
  sub_pz.AddRds(Reduction(sym->irrep_z, MatrixXcd::Ones(1, 1)));
  sub_pz.SetUp();

  // GTO set
  SymGTOs gtos(sym);
  gtos.AddSub(sub_s);
  gtos.AddSub(sub_pz);
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  // Calculate matrix
  BMatSet bmat_set;
  gtos.CalcMat(&bmat_set);
  
  // Solve eigen value problem
  MO mo = CalcOneEle(sym, bmat_set);;
  EXPECT_C_NEAR(-0.5, mo->energy, pow(10.0, -4));

  //  cout << "z >>> " << endl;
  //  cout << bmat_set.GetMatrix("z", sym->irrep_z, sym->irrep_s) << endl;
  //  cout << "z <<< " << endl;

  // Solve driven type equation
  double w10(1.0), w14(1.4);
  dcomplex alpha10 = CalcAlpha(mo, bmat_set, sym->irrep_s, 0, mo->H, w10);
  dcomplex alpha14 = CalcAlpha(mo, bmat_set, sym->irrep_s, 0, mo->H, w14);
  double cs10 = PITotalCrossSection(alpha10, w10, 1);
  double cs14 = PITotalCrossSection(alpha14, w14, 1);
  cout << "w : calc : exact" << endl;
  cout << w10 << ": " << cs10 << " : " << 0.954 << endl ;
  cout << w14 << ": " << cs14 << " : " << 0.353 << endl;
}
TEST(STEX, He) {

  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::D2h();

  // sub set (S orbital)
  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  VectorXcd zeta_s(10);
  zeta_s << 0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053, 30.17990, 108.7723, 488.8941, 3293.694;
  sub_s.AddZeta(zeta_s);
  sub_s.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));
  sub_s.SetUp();

  // sub set (P orbital)
  SubSymGTOs sub_z(sym);
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
  sub_z.SetUp();

  // GTO set
  SymGTOs gtos(sym);
  gtos.AddSub(sub_s);
  gtos.AddSub(sub_z);
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 2.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  // compute basic matrix
  BMatSet mat_set;  
  IB2EInt *eri;
  gtos.CalcMat(&mat_set);

  string fn("eri_hf.bin");
  FILE *fp;  
  fp = fopen(fn.c_str(), "r");
  if(fp == NULL) {
    // If ERI file is not exist:
    cout << "calculate ERI" << endl;
    eri = new B2EIntMem(pow(gtos.size_basis(), 4));
    gtos.CalcERI(eri, 1);
    eri->Write(fn);
  } else {
    // If ERI file is exist:
    cout << "find eri.bin" << endl;
    fclose(fp);
    eri = new B2EIntMem(fn);
  }

  // RHF
  bool conv;
  MO mo = CalcRHF(sym, mat_set, eri, 2, 10, 0.0000001, &conv);
  EXPECT_TRUE(conv);

  // build Static Exchange Hamiltonian
  BMat hmat;
  CalcSEHamiltonian(mo, eri, 0, 0, &hmat);
  
  // Compute alpha
  double au2ev(27.2114);
  double w30(30/au2ev);
  double w40(40/au2ev);
  double w50(50/au2ev);
  
  dcomplex alpha30 = CalcAlpha(mo, mat_set, 0, 0, hmat, w30);
  dcomplex alpha40 = CalcAlpha(mo, mat_set, 0, 0, hmat, w40);
  dcomplex alpha50 = CalcAlpha(mo, mat_set, 0, 0, hmat, w50);
  cout << 30 << " : " << PITotalCrossSection(alpha30, w30, 2) << " : " << 5.0 << endl;
  cout << 40 << " : " << PITotalCrossSection(alpha40, w40, 2) << " : " << 3.0 << endl;
  cout << 50 << " : " << PITotalCrossSection(alpha50, w50, 2) << " : " << 2.0 << endl;
  
  delete eri;

}

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
