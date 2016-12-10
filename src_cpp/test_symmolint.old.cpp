
TEST(SymGTOsMatrix, OneIntNewOld) {

  // ==== Symmetry ====
  SymmetryGroup D2h = SymmetryGroup_D2h();

  // ==== Molecule ====
  Molecule mole = NewMolecule(D2h);
  mole
    ->Add(NewAtom("H", 1.0)->Add(0,0,0.7))
    ->Add(NewAtom("Cen", 0.0)->Add(0,0,0))
    ->SetSymPos();
  EXPECT_EQ(3, mole->size());

  // ==== Sub ====
// ==== GTOs ====
  SymGTOs gtos = NewSymGTOs(mole);
  VectorXcd z1(4); z1 << 2.013, 0.1233, 0.0411, 0.0137;
  MatrixXcd c1_1(2, 1); c1_1 <<+1.0,+1.0;
  MatrixXcd c1_2(2, 1); c1_2 <<+1.0,-1.0;  
  gtos->NewSub("H")
    .AddNs( 0, 0, 0)
    .AddRds(Reduction(D2h->irrep_s(), c1_1))
    .AddRds(Reduction(D2h->irrep_z(), c1_2))
    .AddConts_Mono(z1);
  VectorXcd z2(1); z2 << 1.0;
  MatrixXcd C2_1(2, 1); C2_1 << +1,-1;
  MatrixXcd C2_2(2, 1); C2_2 << +1,+1;
  gtos->NewSub("H")
    .AddNs( 0, 0, 1)
    .AddRds(Reduction(D2h->irrep_s(), C2_1))
    .AddRds(Reduction(D2h->irrep_z(), C2_2))
    .AddConts_Mono(z2);
  VectorXcd z3(1); z3 << dcomplex(0.011389, -0.002197);
  gtos->NewSub("Cen")
    .SolidSH_M(0, 0, z3);
  VectorXcd z4(1); z4 << dcomplex(5.063464, -0.024632);
  MatrixXcd C4_1(1, 3); C4_1 << -1,-1,+2; 
  gtos->NewSub("Cen")
    .AddNs(2, 0, 0)
    .AddNs(0, 2, 0)
    .AddNs(0, 0, 2)
    .AddRds(Reduction(D2h->irrep_s(), C4_1))
    .AddConts_Mono(z4);
  gtos->SetUp();

  BMatSet mat = CalcMat(gtos, gtos, true);
  BMat S,T,V,X,Y,Z,DX,DY,DZ;
  CalcSTVMat(gtos, gtos, &S, &T, &V);
  CalcDipMat(gtos, gtos, &X, &Y, &Z, &DX, &DY, &DZ);
  EXPECT_MATXCD_EQ(S(0,0), mat->GetMatrix("s", 0, 0));
  EXPECT_MATXCD_EQ(Z(0,D2h->irrep_z()),
		   mat->GetMatrix("z", 0, D2h->irrep_z()));

}
void test_SymGTOsOneInt(CartGTO a, Vector3cd at, CartGTO b) {
  /*  
  SymmetryGroup sym = SymmetryGroup_C1();

  Atom h = NewAtom("H", 1.0); h->Add(at);
  Molecule mole(new _Molecule(sym));
  mole->Add(h);

  SymGTOs gtos(new _SymGTOs(mole));

  SubSymGTOs sub_a(sym, h);
  sub_a.AddNs( Vector3i( a.nx,a.ny,a.nz));
  sub_a.AddRds(Reduction(sym->irrep_s(), MatrixXcd::Ones(1, 1)));
  VectorXcd zeta_a(1); zeta_a << a.zeta;
  sub_a.AddConts_Mono(zeta_a);
  gtos->AddSub(sub_a);
  
  SubSymGTOs sub_b(sym, h);
  sub_b.AddNs( b.nx,b.ny,b.nz);
  sub_b.AddRds(Reduction(sym->irrep_s(), MatrixXcd::Ones(1, 1)));
  VectorXcd zeta_b(1); zeta_b << b.zeta;
  sub_b.AddConts_Mono(zeta_b);
  gtos->AddSub(sub_b);

  gtos->SetUp();
  BMatSet mat = CalcMat_Complex(gtos, true);

  const MatrixXcd& S_sym  = mat->GetMatrix("s", 0, 0);
  MatrixXcd S_cart(2, 2);
  S_cart(0, 0) = SMatEle(a, a); S_cart(0, 1) = SMatEle(a, b);
  S_cart(1, 0) = SMatEle(b, a); S_cart(1, 1) = SMatEle(b, b);
  dcomplex c_sym = 1.0/sqrt(S_sym(0, 0) *S_sym( 1, 1));
  dcomplex c_cart= 1.0/sqrt(S_cart(0, 0)*S_cart(1, 1));
  dcomplex s_sym  = S_sym(0, 1)  * c_sym;
  dcomplex s_cart = S_cart(0, 1) * c_cart;
  EXPECT_C_EQ(s_cart, s_sym) << endl
			     << "S matrix" << endl
			     << "a: " << a.str() << endl
			     << "b: " << b.str() << endl;
  dcomplex T_sym = mat->GetMatrix("t", 0, 0)(0, 1)*c_sym;
  dcomplex T_cart= TMatEle(a, b)*c_cart;
  EXPECT_C_EQ(T_cart, T_sym) << endl
			     << "T matrix" << endl
			     << "a: " << a.str() << endl
			     << "b: " << b.str() << endl;
  dcomplex V_sym = mat->GetMatrix("v", 0, 0)(0, 1)*c_sym;
  dcomplex V_cart= VMatEle(a, at, b)*c_cart;
  EXPECT_C_EQ(V_cart, V_sym) << endl
			     << "V matrix" << endl
			     << "a: " << a.str() << endl
			     << "b: " << b.str() << endl;

  dcomplex DX_sym = mat->GetMatrix("dx", 0, 0)(0, 1)*c_sym;
  dcomplex DX_cart= DXMatEle(a, b)*c_cart;
  EXPECT_C_EQ(DX_cart, DX_sym) << endl
			     << "Dx matrix" << endl
			     << "a: " << a.str() << endl
			     << "b: " << b.str() << endl;
  dcomplex DY_sym = mat->GetMatrix("dy", 0, 0)(0, 1)*c_sym;
  dcomplex DY_cart= DYMatEle(a, b)*c_cart;
  EXPECT_C_EQ(DY_cart, DY_sym) << endl
			     << "Dy matrix" << endl
			     << "a: " << a.str() << endl
			     << "b: " << b.str() << endl;
  dcomplex DZ_sym = mat->GetMatrix("dz", 0, 0)(0, 1)*c_sym;
  dcomplex DZ_cart= DZMatEle(a, b)*c_cart;
  EXPECT_C_EQ(DZ_cart, DZ_sym) << endl
			     << "Dz matrix" << endl
			     << "a: " << a.str() << endl
			     << "b: " << b.str() << endl;
  */
  
}
TEST(SymGTOsMatrix, OneInt) {
  
  /*
  CartGTO s0(0, 0, 0, 0.0, 0.0, +0.0, 1.336);
  CartGTO s1(0, 0, 0, 0.0, 0.0, -0.7, 1.336);
  CartGTO p0(0, 0, 1, 0.0, 0.0, +0.7, 1.0);
  CartGTO p1(0, 0, 1, 0.0, 0.0, -0.7, 1.0);
  dcomplex zeta_d(0.00256226, -0.01559939);
  CartGTO dx(2, 0, 0, 0.0, 0.0, 0.0, zeta_d);
  CartGTO dy(0, 2, 0, 0.0, 0.0, 0.0, zeta_d);
  CartGTO dz(0, 0, 2, 0.0, 0.0, 0.0, zeta_d);
  */
  
  /*
  try {
    test_SymGTOsOneInt(s0, Vector3cd(0, 0, 0), s0);
  } catch(runtime_error& e) {
    cout << "s0,s0" << endl;
    cout << e.what() << endl;
    throw runtime_error("exception");
  }    
  
  try {
    test_SymGTOsOneInt(s0, Vector3cd(0, 0, 0.35), s1);
  } catch(runtime_error& e) {
    cout << "s,s" << endl;
    cout << e.what() << endl;
    throw runtime_error("exception");
  }  
  try {
    test_SymGTOsOneInt(s0, Vector3cd(0, 0, 0.35), dz);
  } catch(runtime_error& e) {
    cout << "s,zz" << endl;
    cout << e.what() << endl;
    throw runtime_error("exception");
  }
  try {
    test_SymGTOsOneInt(p0, Vector3cd(0, 0, 0.35), dz);
  } catch(runtime_error& e) {
    cout << "z,zz" << endl;
    cout << e.what() << endl;
  }
  test_SymGTOsOneInt(CartGTO(2, 1, 3, 0.1, 0.2, 0.3, dcomplex(1.0, -0.4)),
		     Vector3cd(-0.1, 0, 0.35),
		     CartGTO(0, 2, 2, 0.4, 0.3, 0.0, dcomplex(0.1, -0.1)));
  test_SymGTOsOneInt(p0, Vector3cd(0, 0, 0.7), dz);
 
 */
  
}

TEST(CompareCColumbus, small_he) {

  dcomplex z1(0.09154356, -0.24865707);
  
  SymmetryGroup sym = SymmetryGroup_Cs();
  Molecule mole = NewMolecule(sym);
  Atom he = NewAtom("He", 2.0); he->Add(0,0,0);
  mole->Add(he);

  // sub set (S orbital)
  SubSymGTOs sub_s(sym, he);
  sub_s.AddNs(Vector3i(0, 0, 0));
  VectorXcd zeta_s(2); zeta_s << 0.107951, 3293.694;
  sub_s.AddConts_Mono(zeta_s);
  sub_s.AddRds(Reduction(sym->irrep_s(), MatrixXcd::Ones(1, 1)));

  // sub set (P orbital)
  SubSymGTOs sub_z(sym, he);
  sub_z.AddNs(Vector3i(0, 0, 1));
  VectorXcd zeta_z(1); zeta_z << z1;
  sub_z.AddConts_Mono(zeta_z);
  sub_z.AddRds(Reduction(sym->irrep_z(), MatrixXcd::Ones(1, 1)));

  // GTO set
  SymGTOs gtos(new _SymGTOs(mole));
  gtos->AddSub(sub_s);
  gtos->AddSub(sub_z);
  gtos->SetUp();

  // compute basic matrix
  BMatSet mat_set = CalcMat_Complex(gtos, true);
  ERIMethod m; m.symmetry = 1;
  B2EInt eri = CalcERI_Complex(gtos, m);
  
  const MatrixXcd& S00 = mat_set->GetMatrix("s", 0, 0);
  EXPECT_C_EQ(dcomplex(1.00000000000000,0.000000000000000), S00(0, 0));
  EXPECT_C_EQ(dcomplex(0.001225127274895041, 0.000000000000000), S00(0, 1));
  EXPECT_C_EQ(dcomplex(0.001225127274895041, 0.000000000000000), S00(1, 0));
  EXPECT_C_EQ(1.0, S00(1, 1));
  const MatrixXcd& S11 = mat_set->GetMatrix("s", 1, 1);
  EXPECT_C_EQ(1.0, S11(0, 0));
  
  const MatrixXcd& T00 = mat_set->GetMatrix("t", 0, 0);
  EXPECT_C_EQ(0.1619265, T00(0, 0));
  EXPECT_C_EQ(0.0003967481399147181, T00(0, 1));
  EXPECT_C_EQ(0.0003967481399147181, T00(1, 0));
  EXPECT_C_EQ(4940.54100000000, T00(1, 1));
  const MatrixXcd& T11 = mat_set->GetMatrix("t", 1, 1);  
  EXPECT_C_EQ(dcomplex(0.2288589, -0.621642675), T11(0, 0));

  const MatrixXcd& V00 = mat_set->GetMatrix("v", 0, 0);
  EXPECT_C_EQ(-1.04860853360520,  V00(0, 0));
  EXPECT_C_EQ(-0.158677374090625, V00(0, 1));
  EXPECT_C_EQ(0.0003967481399147181, T00(1, 0));
  EXPECT_C_EQ(-183.164657050577, V00(1, 1));
  const MatrixXcd& V11 = mat_set->GetMatrix("v", 1, 1);  
  EXPECT_C_EQ(dcomplex(-0.898325031872102,0.626548799637203), V11(0, 0));

  EXPECT_C_EQ(eri->At(1,  1,  1,  1,  0,  0,  0,  0),
	      dcomplex(0.389067179569661,  -0.271360104292752));
  EXPECT_C_EQ(eri->At(1, 1, 0, 0,  0, 0, 0, 0),
	      dcomplex(0.431669280913818,  -0.113052122000587));
  EXPECT_C_EQ(eri->At(1, 0, 1, 0,  0, 0, 0, 0),
	      dcomplex(0.148545492311628,  0.012837791008275));
  EXPECT_C_EQ(eri->At(0, 0, 0, 0,  0, 0, 0, 0),
	      dcomplex(0.370739102461159,   0.000000000000000));
  EXPECT_C_EQ(eri->At(1, 1, 0, 0, 0, 0, 1, 0),
	      dcomplex(0.000550281255615,  -0.000383801011115));
  EXPECT_C_EQ(eri->At(1, 0, 1, 0,  0, 1, 0, 0),
	      dcomplex(-0.000000049115720,  -0.000000075073399));
  EXPECT_C_EQ(eri->At(0, 0, 0, 0, 1, 0, 0, 0),
	      dcomplex(0.000642318406618,   0.000000000000000));
  EXPECT_C_EQ(eri->At(1, 1, 0, 0, 0, 0, 1, 1),
	      dcomplex( 0.449162517258858,  -0.313274399690404));
  EXPECT_C_EQ(eri->At(1, 0, 1, 0, 0, 1, 0, 1),
	      dcomplex(-0.000000007054518,  -0.000000000684658));
  EXPECT_C_EQ(eri->At(0, 0, 0, 0, 1, 1, 0, 0),
	      dcomplex( 0.524295674963366,   0.000000000000000));
  EXPECT_C_EQ(eri->At(0, 0, 0, 0, 1, 0, 1, 0),
	      dcomplex(0.000068730771674,   0.000000000000000));
  EXPECT_C_EQ(eri->At(0, 0, 0, 0, 1, 1, 1, 0),
	      dcomplex( 0.064779412850629,   0.000000000000000));
  EXPECT_C_EQ(eri->At(0, 0, 0, 0, 1, 1, 1, 1),
	      dcomplex(64.758485537085718,   0.000000000000000));
	   
}
