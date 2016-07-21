#include <gtest/gtest.h>
#include <Eigen/Core>

#include "gtest_plus.hpp"
#include "eigen_plus.hpp"

#include "mol_func.hpp"
#include "one_int.hpp"
#include "two_int.hpp"
#include "symmolint.hpp"


using namespace std;
using namespace l2func;
using namespace Eigen;

// other calculation:
// 2016/3/22

class TestValue : public ::testing::Test {
public:
  MatrixXcd S_ref, T_ref, V_ref;
  CartGTO *cart_gtos;
  double R0;
  IB2EInt *eri;
  VectorXcd N;
  TestValue() {

    // -- copied from 
    // ~/calc/ccolumbus/h2/look_matrix/prim/
    // Warnining basis number 2 is not normalized.

    S_ref = MatrixXcd::Zero(3, 3);
    S_ref(1, 0) = S_ref(0, 1) = dcomplex(0.513873194611748, 0.0);
    S_ref(2, 0) = S_ref(0, 2) = dcomplex(-0.004134825663905873,-0.003633459594216835);
    S_ref(2, 1) = S_ref(1, 2) = dcomplex(0.008021777897592569, 0.007549817884085504);
    S_ref(0, 0) = 1.0;
    S_ref(1, 1) = 1.0;
    S_ref(2, 2) = 3.0;

    T_ref = MatrixXcd::Zero(3, 3);
    T_ref(0, 0) = dcomplex(2.00400000000000,0.0);
    T_ref(1, 1) = 2.5;
    T_ref(2, 2) = dcomplex(0.01665468999999993,-0.101396035000000);
    T_ref(1, 0) = T_ref(0, 1) = dcomplex(0.810581687160939,0.0);
    T_ref(2, 0) = T_ref(0, 2) = dcomplex(0.004400790809944475, 0.004508684316644799);
    T_ref(2, 1) = T_ref(1, 2) = dcomplex(0.001116507965258936,-0.0008268252225118336);

    V_ref = MatrixXcd::Zero(3, 3);
    V_ref(0, 0) = -2.55789824622471;
    V_ref(1, 0) = V_ref(0, 1) = -1.28506460356472;
    V_ref(2, 0) = V_ref(0, 2) = dcomplex(0.006725376894693752, 0.005806630112951551);
    V_ref(1, 1) = -1.91602738492149;
    V_ref(2, 1) = V_ref(1, 2) = dcomplex(-0.008460539669617330, -0.007855139694167059);
    V_ref(2, 2) = dcomplex(-0.487930940839014, 0.418008125164049);
    
    eri = new B2EIntMem(21);
    eri->Set(0, 0, 0, 0, 0,0, 0, 0, dcomplex( 1.304242320953500, 0.000000000000000));
    eri->Set(0, 0, 0, 0, 1,0, 0, 0, dcomplex( 0.628781702685092, 0.000000000000000));
    eri->Set(0, 0, 0, 0, 1,1, 0, 0, dcomplex( 0.759778545089558, 0.000000000000000));
    eri->Set(0, 0, 0, 0, 1,0, 1, 0, dcomplex( 0.339943011162754, 0.000000000000000));
    eri->Set(0, 0, 0, 0, 1,1, 1, 0, dcomplex( 0.446428664582801, 0.000000000000000));
    eri->Set(0, 0, 0, 0, 1,1, 1, 1, dcomplex(0.921509653127999,	 0.000000000000000)); 
    eri->Set(0, 0, 0, 0, 2,0, 0, 0, dcomplex(-0.004008078843243,-0.003473830310264)); 
    eri->Set(0, 0, 0, 0, 2,1, 0, 0, dcomplex( 0.002518384782399, 0.002384106992764)); 
    eri->Set(0, 0, 0, 0, 2,0, 1, 0, dcomplex(-0.001736294955885,-0.001503863734779)); 
    eri->Set(0, 0, 0, 0, 2,1, 1, 0, dcomplex( 0.001535095096629, 0.001443367598624)); 
    eri->Set(0, 0, 0, 0, 2,0, 1, 1, dcomplex(-0.002169086387108,-0.001883459993939)); 
    eri->Set(0, 0, 0, 0, 2,1, 1, 1, dcomplex( 0.005972593235128, 0.005519274508113)); 
    eri->Set(0, 0, 0, 0, 2,2, 0, 0, dcomplex( 0.243965040060797,-0.209006616246956)); 
    eri->Set(0, 0, 0, 0, 2,0, 2, 0, dcomplex( 0.000003926040252, 0.000029032254199)); 
    eri->Set(0, 0, 0, 0, 2,2, 1, 0, dcomplex( 0.125696676700414,-0.106858517235463)); 
    eri->Set(0, 0, 0, 0, 2,1, 2, 0, dcomplex(-0.000001456206062,-0.000015419411288)); 
    eri->Set(0, 0, 0, 0, 2,2, 1, 1, dcomplex( 0.243188051473405,-0.210254097028053)); 
    eri->Set(0, 0, 0, 0, 2,1, 2, 1, dcomplex( 0.000009097897649, 0.000121343795672)); 
    eri->Set(0, 0, 0, 0, 2,2, 2, 0, dcomplex(-0.001771869914966,-0.000003726441273)); 
    eri->Set(0, 0, 0, 0, 2,2, 2, 1, dcomplex( 0.003556365185446, 0.000070319128288)); 
    eri->Set(0, 0, 0, 0, 2,2, 2, 2, dcomplex( 0.736988908548374 ,-0.625811204311571));

    cart_gtos = new CartGTO[3];
    R0 = 1.4;
    dcomplex zeta_s(1.336);
    CartGTO s0(0, 0, 0, 0.0, 0.0, +R0/2.0, +zeta_s);
    
    dcomplex zeta_p(1.0);
    CartGTO p1(0, 0, 1, 0.0, 0.0, -R0/2.0, +zeta_p);
    
    dcomplex zeta_d(0.00256226, -0.01559939);
    CartGTO dz(0, 0, 2, 0.0, 0.0, 0.0, zeta_d);

    cart_gtos[0] = s0; cart_gtos[1] = p1; cart_gtos[2] = dz; 

    // Adjust cCoulombus normalization.
    // I dont knoe the origion of this difference of sign for basis 2.
    N = VectorXcd::Zero(3);
    N[0] = +1.0/sqrt(SMatEle(cart_gtos[0], cart_gtos[0]));
    N[1] = +1.0/sqrt(SMatEle(cart_gtos[1], cart_gtos[1]));
    N[2] = -1.0/sqrt(SMatEle(cart_gtos[2], cart_gtos[2]));
    
  }
};
TEST_F(TestValue, OneInt) {

  MatrixXcd S(3, 3), T(3, 3), V(3, 3);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      S(i, j) = SMatEle(cart_gtos[i], cart_gtos[j]);
      T(i, j) = TMatEle(cart_gtos[i], cart_gtos[j]);
      V(i, j) = (VMatEle(cart_gtos[i], Vector3cd(0,0,+R0/2.0), cart_gtos[j]) +
		 VMatEle(cart_gtos[i], Vector3cd(0,0,-R0/2.0), cart_gtos[j]));
    }
  }

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j <= i; j++) {
      // dcomplex c = sqrt(S(i, i)) * sqrt(S(j, j)) / (coef(i) * coef(j));
      dcomplex c = N[i] * N[j];
      dcomplex cr = 1.0/sqrt(S_ref(i, i) * S_ref(j, j));
      EXPECT_C_EQ(S_ref(i, j)*cr, S(i, j)*c) << "S" << i << j;
      EXPECT_C_EQ(T_ref(i, j)*cr, T(i, j)*c) << "T" << i << j;
      EXPECT_C_EQ(V_ref(i, j)*cr, V(i, j)*c) << "V" << i << j;
    }
  }

}
TEST_F(TestValue, TwoInt) {
  
  int ib,jb,kb,lb,i,j,k,l,t;
  dcomplex v;
  eri->Reset();
  while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
    dcomplex cc = N[i]*N[j]*N[k]*N[l];
    dcomplex cr = 1.0/sqrt(S_ref(i, i) * S_ref(j, j) * S_ref(k, k) * S_ref(l, l));
    dcomplex eri_c = ERIEle(cart_gtos[i], cart_gtos[j], cart_gtos[k], cart_gtos[l]);
    EXPECT_C_EQ(cr*v, cc*eri_c) << i << j << k << l;
  }
    
}

void SymGTOs_AtR_Ylm_NDeriv(SymGTOs* gtos, int L, int M, int irrep,
				const VectorXcd& cs, dcomplex r,
				dcomplex* res_dv, dcomplex* res_dv_nd) {
  
  VectorXcd v, dv;
  double h(0.001);
  dcomplex hr(h, 0), hi(0, h);
  VectorXcd r_in(1);
  dcomplex ii(0, 1);

  r_in[0] = r + hr;
  gtos->AtR_Ylm(L, M, irrep, cs, r_in, &v, &dv);
  dcomplex vp0 = v[0];

  r_in[0] = r - hr;
  gtos->AtR_Ylm(L, M, irrep, cs, r_in, &v, &dv);
  dcomplex vm0 = v[0];

  r_in[0] = r + hi;
  gtos->AtR_Ylm(L, M, irrep, cs, r_in, &v, &dv);
  dcomplex iv0p = ii*v[0];

  r_in[0] = r - hi;
  gtos->AtR_Ylm(L, M, irrep, cs, r_in, &v, &dv);
  dcomplex iv0m = ii*v[0];

  r_in[0] = r;
  gtos->AtR_Ylm(L, M, irrep, cs, r_in, &v, &dv);

  *res_dv = dv[0];
  *res_dv_nd = (vp0 + iv0m - vm0 - iv0p) / (4.0 * h);
  
}

TEST(SubSymGTOs, AddZeta) {

  SubSymGTOs sub(SymmetryGroup::C1());
  VectorXcd z1(2); z1 << 1.0, 1.1;
  VectorXcd z2(3); z2 << 2.0, 2.1, 2.2;
  sub.AddZeta(z1);
  sub.AddZeta(z2);

  EXPECT_C_EQ(1.0, sub.get_zeta_iz()(0));
  EXPECT_C_EQ(1.1, sub.get_zeta_iz()(1));
  EXPECT_C_EQ(2.0, sub.get_zeta_iz()(2));
  EXPECT_C_EQ(2.1, sub.get_zeta_iz()(3));
  EXPECT_C_EQ(2.2, sub.get_zeta_iz()(4));

}

void test_SymGTOsOneInt(CartGTO a, Vector3cd at, CartGTO b) {
  
  pSymmetryGroup sym = SymmetryGroup::C1();


  // Build SymGTOs
  SubSymGTOs sub_a(sym);
  sub_a.AddXyz(Vector3cd(a.x, a.y, a.z));
  sub_a.AddNs( Vector3i( a.nx,a.ny,a.nz));
  sub_a.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));
  VectorXcd zeta_a(1); zeta_a << a.zeta;
  sub_a.AddZeta(zeta_a);
  sub_a.SetUp();

  SubSymGTOs sub_b(sym);
  sub_b.AddXyz(Vector3cd(b.x, b.y, b.z));
  sub_b.AddNs( Vector3i( b.nx,b.ny,b.nz));
  sub_b.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));
  VectorXcd zeta_b(1); zeta_b << b.zeta;
  sub_b.AddZeta(zeta_b);
  sub_b.SetUp();

  SymGTOs gtos(sym);
  gtos.AddSub(sub_a);
  gtos.AddSub(sub_b);
  MatrixXcd xyzq(4, 1); 
  xyzq << at[0], at[1], at[2], 1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();
  BMatSet mat; gtos.CalcMat(&mat);

  // Check Matrix
  const MatrixXcd& S_sym  = mat.GetMatrix("s", 0, 0);
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
  dcomplex T_sym = mat.GetMatrix("t", 0, 0)(0, 1)*c_sym;
  dcomplex T_cart= TMatEle(a, b)*c_cart;
  EXPECT_C_EQ(T_cart, T_sym) << endl
			     << "T matrix" << endl
			     << "a: " << a.str() << endl
			     << "b: " << b.str() << endl;
  dcomplex V_sym = mat.GetMatrix("v", 0, 0)(0, 1)*c_sym;
  dcomplex V_cart= VMatEle(a, at, b)*c_cart;
  EXPECT_C_EQ(V_cart, V_sym) << endl
			     << "V matrix" << endl
			     << "a: " << a.str() << endl
			     << "b: " << b.str() << endl;
  
}
TEST(SymGTOsMatrix, OneInt) {

  CartGTO s0(0, 0, 0, 0.0, 0.0, +0.7, 1.336);
  CartGTO s1(0, 0, 0, 0.0, 0.0, -0.7, 1.336);
  CartGTO p0(0, 0, 1, 0.0, 0.0, +0.7, 1.0);
  CartGTO p1(0, 0, 1, 0.0, 0.0, -0.7, 1.0);
  dcomplex zeta_d(0.00256226, -0.01559939);
  CartGTO dx(2, 0, 0, 0.0, 0.0, 0.0, zeta_d);
  CartGTO dy(0, 2, 0, 0.0, 0.0, 0.0, zeta_d);
  CartGTO dz(0, 0, 2, 0.0, 0.0, 0.0, zeta_d);

  test_SymGTOsOneInt(s0, Vector3cd(0, 0, 0.35), s1);
  test_SymGTOsOneInt(s0, Vector3cd(0, 0, 0.35), dz);
  test_SymGTOsOneInt(p0, Vector3cd(0, 0, 0.35), dz);
  test_SymGTOsOneInt(CartGTO(2, 1, 3, 0.1, 0.2, 0.3, dcomplex(1.0, -0.4)),
		     Vector3cd(-0.1, 0, 0.35),
		     CartGTO(0, 2, 2, 0.4, 0.3, 0.0, dcomplex(0.1, -0.1)));
  test_SymGTOsOneInt(p0, Vector3cd(0, 0, 0.7), dz);

}
void test_SymGTOsTwoInt(CartGTO a, CartGTO b, CartGTO c, CartGTO d) {
  
  // == Cart GTO ==
  CartGTO *cart_gtos[4];
  cart_gtos[0] = &a;
  cart_gtos[1] = &b;
  cart_gtos[2] = &c;
  cart_gtos[3] = &d;

  dcomplex c_cart(1);
  for(int i = 0; i < 4; i++) {
    CartGTO *o = cart_gtos[i];
    c_cart *= 1.0/sqrt(SMatEle(*o, *o));
  }

  dcomplex eri_cart = ERIEle(a, b, c, d) * c_cart;
  
  
  // == SymGTOs ==
  pSymmetryGroup sym = SymmetryGroup::C1();
  SymGTOs gtos(sym);
  for(int i = 0; i < 4; i++) {
    SubSymGTOs sub(sym);
    CartGTO *o = cart_gtos[i];
    sub.AddXyz(Vector3cd(o->x, o->y, o->z));
    sub.AddNs( Vector3i( o->nx,o->ny,o->nz));
    sub.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));
    VectorXcd zeta(1); zeta << o->zeta;
    sub.AddZeta(zeta);
    sub.SetUp();
    gtos.AddSub(sub);
  }
  gtos.SetUp();

  BMatSet mat; gtos.CalcMat(&mat);
  const MatrixXcd& S_sym  = mat.GetMatrix("s", 0, 0);
  dcomplex c_sym(1);
  for(int i = 0; i < 4; i++) {
    c_sym *= 1.0/sqrt(S_sym(i, i));
  }

  IB2EInt *eri;
  eri = new B2EIntMem(100);
  ERIMethod m; m.symmetry = 1;
  gtos.CalcERI(eri, m);

  dcomplex eri_sym  = eri->At(0, 0, 0, 0, 0, 1, 2, 3) * c_sym;  

  // == Compare ==
  EXPECT_C_EQ(eri_cart, eri_sym)
    << "a: " << a.str() << endl
    << "b: " << b.str() << endl
    << "c: " << c.str() << endl
    << "d: " << d.str() << endl;
  
}
TEST(SymGTOsMatrix, TwoInt) {

  CartGTO s0(0, 0, 0, 0.0, 0.0, +0.7, 1.336);
  CartGTO s1(0, 0, 0, 0.0, 0.0, -0.7, 1.336);
  CartGTO p0(0, 0, 1, 0.0, 0.0, +0.7, 1.0);
  CartGTO p1(0, 0, 1, 0.0, 0.0, -0.7, 1.0);
  dcomplex zeta_d(0.00256226, -0.01559939);
  CartGTO dx(2, 0, 0, 0.0, 0.0, 0.0, zeta_d);
  CartGTO dy(0, 2, 0, 0.0, 0.0, 0.0, zeta_d);
  CartGTO dz(0, 0, 2, 0.0, 0.0, 0.0, zeta_d);
  CartGTO fff(1, 1, 2, 0.1, 0.2, 0.3, dcomplex(1.0, -0.4));

  test_SymGTOsTwoInt(s0, s0, s0, s0);
  test_SymGTOsTwoInt(s0, s1, p0, p1);
  test_SymGTOsTwoInt(s0, dz, dx, p1);
  test_SymGTOsTwoInt(s0, dz, fff, p1);
  
}

// -- to be removed
SymGTOs* GTOs_h2_small() {

  double R0 = 1.4;
  pSymmetryGroup sym = SymmetryGroup::D2h();
  SymGTOs* gtos = new SymGTOs(sym);

  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0,+R0/2.0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  VectorXcd zetas(1);
  zetas << 1.336;
  sub_s.AddZeta(zetas); 
  MatrixXcd c(1, 1); c << 1;
  sub_s.AddRds(Reduction(sym->irrep_s, c));
  sub_s.SetUp(); 
  gtos->AddSub(sub_s); 

  SubSymGTOs sub_p(sym);
  sub_p.AddXyz(Vector3cd(0, 0, -R0/2.0));
  sub_p.AddNs(Vector3i(  0, 0, 1));
  VectorXcd zeta_p(1); zeta_p << 1.0;
  sub_p.AddZeta(zeta_p);
  MatrixXcd cp(1, 1); cp << 1;
  sub_p.AddRds(Reduction(sym->irrep_s, cp));
  sub_p.SetUp();
  gtos->AddSub(sub_p);

  SubSymGTOs sub_p_cen(sym);
  sub_p_cen.AddXyz(Vector3cd(0, 0, 0));
  sub_p_cen.AddNs( Vector3i( 0, 0, 2));
  VectorXcd zeta_cen(1);
  zeta_cen << dcomplex(0.00256226, -0.01559939);
  sub_p_cen.AddZeta(zeta_cen);
  MatrixXcd cd(1, 1); cd << 1;
  sub_p_cen.AddRds(Reduction(sym->irrep_s, cd));
  sub_p_cen.SetUp();
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
// -- to be removed
void CalcMatByCart(MatrixXcd& S, MatrixXcd& T, MatrixXcd& V) {
  CartGTO s0(0, 0, 0, 0.0, 0.0, +0.7, 1.336);
  CartGTO s1(0, 0, 0, 0.0, 0.0, -0.7, 1.336);
  CartGTO p0(0, 0, 1, 0.0, 0.0, +0.7, 1.0);
  CartGTO p1(0, 0, 1, 0.0, 0.0, -0.7, 1.0);
  CartGTO dx(2, 0, 0, 0.0, 0.0, 0.0, 1.336);
  CartGTO dy(0, 2, 0, 0.0, 0.0, 0.0, 1.336);
  CartGTO dz(0, 0, 2, 0.0, 0.0, 0.0, 1.336);

  // 0 vs 2
  S = MatrixXcd::Zero(3, 3);
  S(0, 0) = SMatEle(s0, s0) + 2.0*SMatEle(s0, s1) + SMatEle(s1, s1);
  S(1, 1) = SMatEle(p0, p0) - 2.0*SMatEle(p0, p1) + SMatEle(p1, p1);
  S(2, 2) = (SMatEle(dx, dx)        + SMatEle(dx, dy)        -2.0 * SMatEle(dx, dz) +
	     SMatEle(dy, dx)        + SMatEle(dy, dy)        -2.0 * SMatEle(dy, dz) +
	     (-2.0)*SMatEle(dz, dx) + (-2.0)*SMatEle(dz, dy) + 4.0*SMatEle(dz, dz));
  S(0, 1) = S(1, 0) = (+SMatEle(s0, p0) - SMatEle(s0, p1) 
		       +SMatEle(s1, p0) - SMatEle(s1, p1));
  S(0, 2) = S(2, 0) = (SMatEle(s0, dx) + SMatEle(s0, dy) -2.0*SMatEle(s0, dz) + 
		       SMatEle(s1, dx) + SMatEle(s1, dy) -2.0*SMatEle(s1, dz));
  S(1, 2) = S(2, 1) = (+SMatEle(p0, dx) + SMatEle(p0, dy) -2.0*SMatEle(p0, dz)
		       -SMatEle(p1, dx) - SMatEle(p1, dy) +2.0*SMatEle(p1, dz));
  T = MatrixXcd::Zero(3, 3);
  T(0, 0) = TMatEle(s0, s0) + 2.0*TMatEle(s0, s1) + TMatEle(s1, s1);
  T(1, 1) = TMatEle(p0, p0) - 2.0*TMatEle(p0, p1) + TMatEle(p1, p1);
  T(2, 2) = (TMatEle(dx, dx)        + TMatEle(dx, dy)        -2.0 * TMatEle(dx, dz) +
	     TMatEle(dy, dx)        + TMatEle(dy, dy)        -2.0 * TMatEle(dy, dz) +
	     (-2.0)*TMatEle(dz, dx) + (-2.0)*TMatEle(dz, dy) + 4.0*TMatEle(dz, dz));
  T(0, 1) = T(1, 0) = (+TMatEle(s0, p0) - TMatEle(s0, p1) 
		       +TMatEle(s1, p0) - TMatEle(s1, p1));
  T(0, 2) = T(2, 0) = (TMatEle(s0, dx) + TMatEle(s0, dy) -2.0*TMatEle(s0, dz) + 
		       TMatEle(s1, dx) + TMatEle(s1, dy) -2.0*TMatEle(s1, dz));
  T(1, 2) = T(2, 1) = (+TMatEle(p0, dx) + TMatEle(p0, dy) -2.0*TMatEle(p0, dz)
		       -TMatEle(p1, dx) - TMatEle(p1, dy) +2.0*TMatEle(p1, dz));
  
  V = MatrixXcd::Zero(3, 3);
  Vector3cd kat0(0, 0, +0.7);
  Vector3cd kat1(0, 0, -0.7);
  V(0, 0) = (VMatEle(s0, kat0, s0) + 2.0*VMatEle(s0, kat0, s1) + VMatEle(s1, kat0, s1)+
	     VMatEle(s0, kat1, s0) + 2.0*VMatEle(s0, kat1, s1) + VMatEle(s1, kat1, s1));
  V(1, 1) = (VMatEle(p0, kat0, p0) - 2.0*VMatEle(p0, kat0, p1) + VMatEle(p1, kat0, p1) +
	     VMatEle(p0, kat1, p0) - 2.0*VMatEle(p0, kat1, p1) + VMatEle(p1, kat1, p1));
  V(2, 2) = (VMatEle(dx, kat0, dx) + VMatEle(dx, kat0, dy) -2.0 * VMatEle(dx, kat0, dz) +
	     VMatEle(dy, kat0, dx) + VMatEle(dy, kat0, dy) -2.0 * VMatEle(dy, kat0, dz) +
	     (-2.0)*VMatEle(dz, kat0, dx) + (-2.0)*VMatEle(dz, kat0, dy) + 4.0*VMatEle(dz, kat0, dz)+
	     VMatEle(dx, kat1, dx) + VMatEle(dx, kat1, dy) -2.0 * VMatEle(dx, kat1, dz) +
	     VMatEle(dy, kat1, dx) + VMatEle(dy, kat1, dy) -2.0 * VMatEle(dy, kat1, dz) +
	     (-2.0)*VMatEle(dz, kat1, dx) + (-2.0)*VMatEle(dz, kat1, dy) + 4.0*VMatEle(dz, kat1, dz));

  V(0, 1) = V(1, 0) = (+VMatEle(s0, kat0, p0) - VMatEle(s0, kat0, p1) 
		       +VMatEle(s1, kat0, p0) - VMatEle(s1, kat0, p1)
		       +VMatEle(s0, kat1, p0) - VMatEle(s0, kat1, p1) 
		       +VMatEle(s1, kat1, p0) - VMatEle(s1, kat1, p1));
  V(0, 2) = V(2, 0) = (VMatEle(s0, kat0, dx) + VMatEle(s0, kat0, dy)
		       -2.0*VMatEle(s0, kat0, dz) + 
		       VMatEle(s1, kat0, dx) + VMatEle(s1, kat0, dy)
		       -2.0*VMatEle(s1, kat0, dz) +
		       VMatEle(s0, kat1, dx) + VMatEle(s0, kat1, dy)
		       -2.0*VMatEle(s0, kat1, dz) + 
		       VMatEle(s1, kat1, dx) + VMatEle(s1, kat1, dy)
		       -2.0*VMatEle(s1, kat1, dz));
  V(1, 2) = V(2, 1) = (+VMatEle(p0, kat0, dx) + VMatEle(p0, kat0, dy)
		       -2.0*VMatEle(p0, kat0, dz)
		       -VMatEle(p1, kat0, dx) - VMatEle(p1, kat0, dy)
		       +2.0*VMatEle(p1, kat0, dz)
		       +VMatEle(p0, kat1, dx) + VMatEle(p0, kat1, dy)
		       -2.0*VMatEle(p0, kat1, dz)
		       -VMatEle(p1, kat1, dx) - VMatEle(p1, kat1, dy)
		       +2.0*VMatEle(p1, kat1, dz));
  
}

TEST(SymGTOs, at_r_ylm_lin) {

  pSymmetryGroup C1 = SymmetryGroup::C1();
  SymGTOs gtos(C1);
  VectorXcd zeta(1); zeta << 1.1;
  gtos.AddSub(Sub_s(C1, 0, Vector3cd(0, 0, 2), zeta));
  gtos.SetUp();
  VectorXcd cs1(1); cs1 << 1.1;  
  VectorXcd cs2(1); cs2 << 2.2;  
  VectorXcd rs(10); 
  for(int i = 0; i < 10; i++)
    rs(i) = i * 0.1 + 0.1;
  VectorXcd vs1, dvs1;
  VectorXcd vs2, dvs2;

  gtos.AtR_Ylm(0, 0, 0, cs1, rs, &vs1, &dvs1);
  gtos.AtR_Ylm(0, 0, 0, cs2, rs, &vs2, &dvs2);
  EXPECT_C_EQ(2.0 * vs1(2), vs2(2));

  gtos.AtR_Ylm(1, 0, 0, cs1, rs, &vs1, &dvs1);
  gtos.AtR_Ylm(1, 0, 0, cs2, rs, &vs2, &dvs2);
  EXPECT_TRUE(abs(vs1(2)) > 0.001);
  EXPECT_TRUE(abs(dvs1(2)) > 0.001);
  EXPECT_C_EQ(2.0 * vs1(2), vs2(2));

}
void test_SymGTOs_at_r_ylm_nd(SymGTOs& gtos, dcomplex r, VectorXcd& cs, string label) {

  dcomplex nd, dv;
  SymGTOs_AtR_Ylm_NDeriv(&gtos, 0, 0, 0,
			 cs, r, &dv, &nd);
  EXPECT_C_EQ(nd, dv) << label << endl
		      << "r = " << r << endl;
		      

}
TEST(SymGTOs, at_r_ylm_nd_s) {

  pSymmetryGroup C1 = SymmetryGroup::C1();
  SymGTOs gtos(C1);

  /*
  VectorXcd zeta(1); zeta << 1.2;
  gtos.AddSub(Sub_s(0, Vector3cd(0.0, 0.0, 0.0), zeta));
  VectorXcd cs(1); cs << 0.2;
  */

  VectorXcd zeta(2); zeta << 1.2, 1.4;
  gtos.AddSub(Sub_s(C1, 0, Vector3cd(0.0, 0.0, 0.0), zeta));
  VectorXcd cs(2); cs << 0.2, 1.1;

  gtos.SetUp();
  test_SymGTOs_at_r_ylm_nd(gtos, 1.1, cs, "s-GTO on center");

}
TEST(SymGTOs, at_r_ylm_nd) {

  pSymmetryGroup C1 = SymmetryGroup::C1();
  SymGTOs gtos(C1);
  // VectorXcd zeta(2); zeta << 1.2, 1.1;
  // gtos.AddSub(Sub_s(0, Vector3cd(0.1, 0.2, 0.3), zeta));
  // VectorXcd cs(2); cs << 0.2, 0.4;

  VectorXcd zeta(1); zeta << 1.2;
  gtos.AddSub(Sub_s(C1, 0, Vector3cd(0.0, 0.0, 2.1), zeta));
  VectorXcd cs(1); cs << 1.0;

  gtos.SetUp();
  test_SymGTOs_at_r_ylm_nd(gtos, 1.1, cs, "s-GTO on (0.1, 0.2, 0.3)");

}
TEST(SymGTOs, Create) {

  pSymmetryGroup sym_group = SymmetryGroup::Cs();
  cout << sym_group->num_class() << endl;
  SymGTOs gtos(sym_group);

  // -- Symmetries --
  Irrep Ap = sym_group->GetIrrep("A'");
  //  Irrep App= sym_group->GetIrrep("A''");

  // -- s-GTO at Center, Sym = (0, 0)
  SubSymGTOs sgto_cen(sym_group);
  VectorXcd zetas(10);
  for(int n = 0; n < 10; n++)
    zetas(n) = pow(2.0, n-5);
  sgto_cen.AddXyz( Vector3cd(0, 0, 0));
  sgto_cen.AddNs(  Vector3i(0, 0, 0));
  sgto_cen.AddZeta(zetas);
  sgto_cen.AddRds( Reduction(Ap, MatrixXcd::Ones(1, 1)));
  sgto_cen.SetUp();
   gtos.AddSub(sgto_cen);

   // -- s-GTO at z-axis, sym = (1, 0)
   SubSymGTOs pgto_z(sym_group);
   VectorXcd zetas2(6);
   for(int n = 0; n < 6; n++)
     zetas2[n] = pow(1.6, n-3);
   pgto_z.AddXyz(Vector3cd(0, 0, +1));
   pgto_z.AddXyz(Vector3cd(0, 0, -1));
   pgto_z.AddNs( Vector3i( 0, 0, +1));	
   pgto_z.AddZeta(zetas2);
   MatrixXcd c1(2, 1); c1 << 1.0, -1.0;
   pgto_z.AddRds(Reduction(Ap, c1));
   MatrixXcd c2(2, 1); c2 << 1.0, +1.0;
   pgto_z.AddRds(Reduction(Ap, c2));
   pgto_z.SetUp();
   gtos.AddSub(pgto_z);

   EXPECT_EQ(1,	 gtos.subs[0].size_at());
   EXPECT_EQ(1,	 gtos.subs[0].size_pn());
   EXPECT_EQ(1,	 gtos.subs[0].size_cont());
   EXPECT_EQ(10, gtos.subs[0].size_zeta());

   EXPECT_EQ(2,	 gtos.subs[1].size_at());
   EXPECT_EQ(1,	 gtos.subs[1].size_pn());
   EXPECT_EQ(2,	 gtos.subs[1].size_cont());
   EXPECT_EQ(6,	 gtos.subs[1].size_zeta());

   BMatSet mat;
   gtos.CalcMat(&mat);
   EXPECT_C_EQ(1.0, mat.GetValue("s", Ap, Ap, 0, 0));
 }
TEST(SymGTOs, add_atom) {

   SymGTOs gtos(SymmetryGroup::C1()); 
   gtos.AddAtom(Vector3cd(1.2, 1.3, 1.4), 1.5);
   gtos.AddAtom(Vector3cd(2.2, 2.3, 2.4), 2.5);

   EXPECT_C_EQ(1.2, gtos.xyzq_iat(0, 0));
   EXPECT_C_EQ(1.3, gtos.xyzq_iat(1, 0));
   EXPECT_C_EQ(1.4, gtos.xyzq_iat(2, 0));
   EXPECT_C_EQ(1.5, gtos.xyzq_iat(3, 0));

   EXPECT_C_EQ(2.2, gtos.xyzq_iat(0, 1));
   EXPECT_C_EQ(2.3, gtos.xyzq_iat(1, 1));
   EXPECT_C_EQ(2.4, gtos.xyzq_iat(2, 1));
   EXPECT_C_EQ(2.5, gtos.xyzq_iat(3, 1));

 }
TEST(SymGTOs, CalcMatOther) {

  pSymmetryGroup Cs = SymmetryGroup::Cs();
   SymGTOs gtos_full(Cs);
   SymGTOs gtos_1(   Cs);
   SymGTOs gtos_2(   Cs);
   Irrep Ap = Cs->GetIrrep("A'");
   Irrep App= Cs->GetIrrep("A''");
   dcomplex z_gh(1.1);

   // -- A' symmetry --
   VectorXcd zeta_h = VectorXcd::Zero(3); zeta_h << 0.4, 1.0, 2.0;
   gtos_full.AddSub(Sub_s(Cs, Ap, Vector3cd(0, 0, 0), zeta_h));
   gtos_1.AddSub(Sub_s(Cs, Ap, Vector3cd(0, 0, 0), zeta_h));


   // -- A'' symmetry, Center --
   VectorXcd zeta_gh = VectorXcd::Zero(3); zeta_gh << 0.6, 1.2, 2.2;
   gtos_full.AddSub(Sub_TwoSGTO(Cs, App,
			       Vector3cd(0.0, 0.0, z_gh), zeta_gh));
   VectorXcd zeta_cen(2); zeta_cen << 0.2, 1.2;
   gtos_full.AddSub(Sub_pz(Cs, App, Vector3cd(0, 0, 0), zeta_cen));

   gtos_2.AddSub(Sub_TwoSGTO(Cs, App,
			    Vector3cd(0.0, 0.0, z_gh), zeta_gh));
   gtos_2.AddSub(Sub_pz(Cs, App, Vector3cd(0, 0, 0), zeta_cen));


   // -- potential --
   MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
   gtos_full.SetAtoms(xyzq);
   gtos_1.SetAtoms(xyzq);
   gtos_2.SetAtoms(xyzq);


   // -- calculate matrix --
   BMatSet mat_full(2), mat_12(2);
   gtos_full.SetUp();
   gtos_1.SetUp();
   gtos_2.SetUp();
   gtos_full.CalcMat(&mat_full);
   gtos_1.CalcMatOther(gtos_2, false, &mat_12);


   // -- compare results --
   EXPECT_MATXCD_EQ(mat_full.GetMatrix("z", Ap, App),
		   mat_12.GetMatrix(  "z", Ap, App));

 }
TEST(SymGTOs, conjugate) {

  pSymmetryGroup Cs = SymmetryGroup::Cs();
  SymGTOs gtos_1(Cs);
  SymGTOs gtos_3(Cs);
  SymGTOs gtos_2(Cs);
  Irrep Ap = Cs->GetIrrep("A'");
  //  Irrep App= Cs->GetIrrep("A''");
  dcomplex z_gh(1.1);
  
  // -- A' symmetry --
  VectorXcd zeta_h(3);
  zeta_h << dcomplex(0.4, 0.3), dcomplex(1.0, 0.4), dcomplex(2.0, 0.1);
  gtos_1.AddSub(Sub_s(Cs, Ap, Vector3cd(0, 0, 0), zeta_h));
  gtos_2.AddSub(Sub_s(Cs, Ap, Vector3cd(0, 0, 0), zeta_h.conjugate()));

  // -- potential --
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos_1.SetAtoms(xyzq);
  gtos_2.SetAtoms(xyzq);
  
  gtos_1.SetUp();
  gtos_2.SetUp();
  gtos_3.SetComplexConj(gtos_1);
  
  BMatSet mat2, mat3;
  gtos_2.CalcMat(&mat2);
  gtos_3.CalcMat(&mat3);
  
  EXPECT_C_EQ(mat2.GetMatrix("t", 0, 0)(0, 1),
	      mat3.GetMatrix("t", 0, 0)(0, 1));

}

TEST(CompareCColumbus, small_h2) {

  double R0 = 1.4;
  pSymmetryGroup sym = SymmetryGroup::D2h();
  SymGTOs gtos(sym);

  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0,+R0/2.0));
  sub_s.AddXyz(Vector3cd(0, 0,-R0/2.0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  VectorXcd zetas(1);
  zetas << 1.336;
  sub_s.AddZeta(zetas); 
  MatrixXcd c(2, 1); c << 1, 1;
  sub_s.AddRds(Reduction(sym->irrep_s, c));
  sub_s.SetUp(); 
  gtos.AddSub(sub_s); 

  SubSymGTOs sub_p(sym);
  sub_p.AddXyz(Vector3cd(0, 0, +R0/2.0));
  sub_p.AddXyz(Vector3cd(0, 0, -R0/2.0));
  sub_p.AddNs(Vector3i(  0, 0, 1));
  VectorXcd zeta_p(1); zeta_p << 1.0;
  sub_p.AddZeta(zeta_p);
  MatrixXcd cp(2, 1); cp << 1, -1;
  sub_p.AddRds(Reduction(sym->irrep_s, cp));
  sub_p.SetUp();
  gtos.AddSub(sub_p);

  SubSymGTOs sub_p_cen(sym);
  sub_p_cen.AddXyz(Vector3cd(0, 0, 0));
  sub_p_cen.AddNs( Vector3i( 2, 0, 0));
  sub_p_cen.AddNs( Vector3i( 0, 2, 0));
  sub_p_cen.AddNs( Vector3i( 0, 0, 2));
  VectorXcd zeta_cen(1);
  zeta_cen << dcomplex(0.00256226, -0.01559939);
  //zeta_cen << dcomplex(0.00256226, -0.01559939);
  sub_p_cen.AddZeta(zeta_cen);
  MatrixXcd cd(1, 3); cd << 1, 1, -2;
  sub_p_cen.AddRds(Reduction(sym->irrep_s, cd));
  sub_p_cen.SetUp();
  gtos.AddSub(sub_p_cen);
  
  MatrixXcd xyzq(4, 2); xyzq <<
			  0.0,      0.0,
			  0.0,      0.0,
			  +R0/2.0, -R0/2.0,
			  1.0,      1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  BMatSet mat; gtos.CalcMat(&mat);

  // -- look ylcls:~/calc/ccolumbus/h2/look_matrix
  MatrixXcd S_ref00(3, 3);
  S_ref00(0, 0) = dcomplex(2.54002879355886,0.00000000000);
  S_ref00(0, 1) = dcomplex(-1.02774638922350,0.00000000000);
  S_ref00(0, 2) = dcomplex(-0.009327584786309644,-0.00828045703077020);
  S_ref00(1, 1) = dcomplex(2.72059730979469,0.0000000000000);
  S_ref00(1, 2) = dcomplex(-0.03236115842797820,-0.02998289626140810);
  S_ref00(2, 2) = dcomplex(12.0);
  S_ref00(2, 1) = S_ref00(1, 2);
  S_ref00(1, 0) = S_ref00(0, 1); 
  S_ref00(2, 0) = S_ref00(0, 2); 

  MatrixXcd T_ref00(3, 3);
  T_ref00(0, 0) = dcomplex(4.14560037345408,0.0);
  T_ref00(0, 1) = dcomplex(-1.62116337432188,0.0);
  T_ref00(0, 2) = dcomplex(-0.001080858125006526,0.000853094902854401);
  T_ref00(1, 1) = dcomplex(7.56652741838541,0.0);
  T_ref00(1, 2) = dcomplex(-0.003899884013199560,0.002908577808236390);
  T_ref00(2, 2) = dcomplex(0.107614920000000,-0.655174379999998);
  T_ref00(2, 1) = T_ref00(1, 2);
  T_ref00(1, 0) = T_ref00(0, 1); 
  T_ref00(2, 0) = T_ref00(0, 2); 

  MatrixXcd V_ref00(3, 3);
  V_ref00(0, 0) = dcomplex(-6.49577025469466,0.0);
  V_ref00(0, 1) = dcomplex(2.98842408078067,0.0);
  V_ref00(0, 2) = dcomplex(0.01636122032740119,0.01419384999265892);
  V_ref00(1, 1) = dcomplex(-5.53999525981359,0.0);
  V_ref00(1, 2) = dcomplex(0.03653639231763511, 0.03327467147633245);
  V_ref00(2, 2) = dcomplex(-1.95460635098110,1.66714310112737);
  V_ref00(2, 1) = V_ref00(1, 2);
  V_ref00(1, 0) = V_ref00(0, 1); 
  V_ref00(2, 0) = V_ref00(0, 2);   

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++) {
      dcomplex c = 1.0/sqrt(S_ref00(i, i) * S_ref00(j, j));
      EXPECT_C_EQ(S_ref00(i, j)*c, mat.GetMatrix("s", 0, 0)(i, j)) << i << j;
      EXPECT_C_EQ(T_ref00(i, j)*c, mat.GetMatrix("t", 0, 0)(i, j)) << i << j;
      EXPECT_C_EQ(V_ref00(i, j)*c, mat.GetMatrix("v", 0, 0)(i, j)) << i << j;
    }
  

}
TEST(CompareCColumbus, small_he) {

  dcomplex z1(0.09154356, -0.24865707);
  
  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::Cs();

  // sub set (S orbital)
  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0, 0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  VectorXcd zeta_s(2); zeta_s << 0.107951, 3293.694;
  sub_s.AddZeta(zeta_s);
  sub_s.AddRds(Reduction(sym->irrep_s, MatrixXcd::Ones(1, 1)));
  sub_s.SetUp();

  // sub set (P orbital)
  SubSymGTOs sub_z(sym);
  sub_z.AddXyz(Vector3cd(0, 0, 0));
  sub_z.AddNs(Vector3i(0, 0, 1));
  VectorXcd zeta_z(1); zeta_z << z1;
  sub_z.AddZeta(zeta_z);
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
  BMatSet mat_set; gtos.CalcMat(&mat_set);
  IB2EInt *eri = new B2EIntMem();
  ERIMethod m; m.symmetry = 1;
  gtos.CalcERI(eri, m);
  
  const MatrixXcd& S00 = mat_set.GetMatrix("s", 0, 0);
  EXPECT_C_EQ(dcomplex(1.00000000000000,0.000000000000000), S00(0, 0));
  EXPECT_C_EQ(dcomplex(0.001225127274895041, 0.000000000000000), S00(0, 1));
  EXPECT_C_EQ(dcomplex(0.001225127274895041, 0.000000000000000), S00(1, 0));
  EXPECT_C_EQ(1.0, S00(1, 1));
  const MatrixXcd& S11 = mat_set.GetMatrix("s", 1, 1);
  EXPECT_C_EQ(1.0, S11(0, 0));
  
  const MatrixXcd& T00 = mat_set.GetMatrix("t", 0, 0);
  EXPECT_C_EQ(0.1619265, T00(0, 0));
  EXPECT_C_EQ(0.0003967481399147181, T00(0, 1));
  EXPECT_C_EQ(0.0003967481399147181, T00(1, 0));
  EXPECT_C_EQ(4940.54100000000, T00(1, 1));
  const MatrixXcd& T11 = mat_set.GetMatrix("t", 1, 1);  
  EXPECT_C_EQ(dcomplex(0.2288589, -0.621642675), T11(0, 0));

  const MatrixXcd& V00 = mat_set.GetMatrix("v", 0, 0);
  EXPECT_C_EQ(-1.04860853360520,  V00(0, 0));
  EXPECT_C_EQ(-0.158677374090625, V00(0, 1));
  EXPECT_C_EQ(0.0003967481399147181, T00(1, 0));
  EXPECT_C_EQ(-183.164657050577, V00(1, 1));
  const MatrixXcd& V11 = mat_set.GetMatrix("v", 1, 1);  
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
TEST(CompareCColumbus, small_h2_2) {

  // set symmetry
  pSymmetryGroup sym = SymmetryGroup::D2h();
  double R0 = 1.4;

  // sub set (S orbital)
  SubSymGTOs sub_s(sym);
  sub_s.AddXyz(Vector3cd(0, 0, +R0/2.0));
  sub_s.AddXyz(Vector3cd(0, 0, -R0/2.0));
  sub_s.AddNs(Vector3i(0, 0, 0));
  VectorXcd zeta_s(1); zeta_s <<2.013;
  sub_s.AddZeta(zeta_s);
  MatrixXcd cs(2, 1); cs << 1.0, 1.0;
  sub_s.AddRds(Reduction(sym->irrep_s, cs));
  sub_s.SetUp();

  // sub set (P orbital)
  SubSymGTOs sub_z(sym);
  sub_z.AddXyz(Vector3cd(0, 0, +R0/2.0));
  sub_z.AddXyz(Vector3cd(0, 0, -R0/2.0));
  sub_z.AddNs(Vector3i(0, 0, 1));
  VectorXcd zeta_z(1); zeta_z << 1.0;
  MatrixXcd cz(2, 1); cz << 1.0, -1.0;
  sub_z.AddZeta(zeta_z);
  sub_z.AddRds(Reduction(sym->irrep_s, cz));
  sub_z.SetUp();

  // sub set (D orbital)
  SubSymGTOs sub_d(sym);
  sub_d.AddXyz(Vector3cd(0, 0, 0));
  sub_d.AddNs(Vector3i(2, 0, 0));
  sub_d.AddNs(Vector3i(0, 2, 0));
  sub_d.AddNs(Vector3i(0, 0, 2));
  VectorXcd zeta_d(1); zeta_d << dcomplex(5.063464, -0.024632);
  MatrixXcd cd(1, 3); cd << -1.0, -1.0, 2.0;
  sub_d.AddZeta(zeta_d);
  sub_d.AddRds(Reduction(sym->irrep_s, cd));
  sub_d.SetUp();

  // GTO set
  SymGTOs gtos(sym);
  gtos.AddSub(sub_s);
  gtos.AddSub(sub_z);
  gtos.AddSub(sub_d);
  MatrixXcd xyzq(4, 2);
  xyzq <<
    0.0, 0.0,
    0.0, 0.0,
    +R0/2.0, -R0/2.0,
    1.0, 1.0;
  gtos.SetAtoms(xyzq);
  gtos.SetUp();

  //  BMatSet mat_set; gtos.CalcMat(&mat_set);
  BMatSet mat_set; CalcMatrix_Complex(gtos, true, &mat_set);
  IB2EInt *eri = new B2EIntMem();
  ERIMethod m; m.symmetry = 1;
  gtos.CalcERI(eri, m);
  
  const MatrixXcd& S00 = mat_set.GetMatrix("s", 0, 0);
  //  const MatrixXcd& S11 = mat_set.GetMatrix("s", 1, 1);
  MatrixXcd S00_ref(3, 3);
  S00_ref(1-1, 1-1) = dcomplex (2.27815053488874,0.000000000000000E+000);
  S00_ref(1-1, 2-1) = dcomplex (-0.923122727265142,0.000000000000000E+000);
  S00_ref(1-1, 3-1) = dcomplex (1.35935784944244,6.316245846956289E-003);
  S00_ref(2-1, 1-1) = dcomplex (-0.923122727265142,0.000000000000000E+000);
  S00_ref(2-1, 2-1) = dcomplex (2.72059730979469,0.000000000000000E+000);
  S00_ref(2-1, 3-1) = dcomplex (0.774070190590069,5.100459746874277E-003);
  S00_ref(3-1, 1-1) = dcomplex (1.35935784944244,6.316245846956289E-003);
  S00_ref(3-1, 2-1) = dcomplex (0.774070190590069,5.100459746874277E-003);
  S00_ref(3-1, 3-1) = dcomplex (12.0000000000000,-1.214306433183765E-017);
  MatrixXcd S11_ref(2, 2);
  S11_ref(0, 0) = dcomplex(1.72184946511126,0.0);
  S11_ref(0, 1) = dcomplex(0.923122727265142,0.0);
  S11_ref(1, 0) = dcomplex(0.923122727265142,0.0);
  S11_ref(1, 1) = dcomplex(1.27940269020531,0.0);
  const MatrixXcd& T00 = mat_set.GetMatrix("t", 0, 0);
  MatrixXcd T00_ref(3, 3);
  T00_ref(1-1, 1-1) = dcomplex(5.77430482478317,0.0);
  T00_ref(1-1, 2-1) =  dcomplex(-1.46848241010191,0.0);
  T00_ref(1-1, 3-1) =  dcomplex(10.9421575645574,0.03952537783003958);
  T00_ref(2-1, 1-1) =  dcomplex(-1.46848241010191,0.0);
  T00_ref(2-1, 2-1) =  dcomplex(7.56652741838541,0.0);
  T00_ref(2-1, 3-1) =  dcomplex(3.10047956467259,0.01958232632537312);
  T00_ref(3-1, 1-1) =  dcomplex(10.9421575645574,0.03952537783003958);
  T00_ref(3-1, 2-1) =  dcomplex(3.10047956467259,0.01958232632537312);
  T00_ref(3-1, 3-1) =  dcomplex(212.665488000000,-1.03454400000000);
  const MatrixXcd& V00 = mat_set.GetMatrix("v", 0, 0);
  MatrixXcd V00_ref(3, 3);
  V00_ref(1-1, 1-1) = dcomplex(-6.71399792745968,0.0);
  V00_ref(1-1, 2-1) = dcomplex(2.80168590393072,0.0);
  V00_ref(1-1, 3-1) = dcomplex(-7.76431478876783,-0.02566801940514894);
  V00_ref(2-1, 1-1) = dcomplex(2.80168590393072,0.0);
  V00_ref(2-1, 2-1) = dcomplex(-5.53999525981359,0.0);
  V00_ref(2-1, 3-1) = dcomplex(0.338975375012722,-0.009191430579939701);
  V00_ref(3-1, 1-1) = dcomplex(-7.76431478876783,-0.02566801940514894);
  V00_ref(3-1, 2-1) = dcomplex(0.338975375012722,-0.009191430579939701);
  V00_ref(3-1, 3-1) = dcomplex(-43.3853460578071,-0.01186566895962939);

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      dcomplex c(1.0/sqrt(S00_ref(i, i)*S00_ref(j, j)));
      EXPECT_C_EQ(S00_ref(i, j)*c, S00(i, j)) << i << j;
      EXPECT_C_EQ(T00_ref(i, j)*c, T00(i, j)) << i << j;
      EXPECT_C_EQ(V00_ref(i, j)*c, V00(i, j)) << i << j;
    }
  }
  dcomplex ref =  dcomplex(249.194759129673912, -0.606119545098925)/(S00_ref(2, 2) * S00_ref(2, 2));
  EXPECT_C_NEAR(ref, eri->At(0, 0, 0, 0, 2, 2, 2, 2),
		250 * pow(10.0, -10.0));
  EXPECT_C_EQ(eri->At(0, 0, 0, 0, 0, 0, 0, 0)
	      , 6.082101996924533/(S00_ref(0, 0) * S00_ref(0, 0)));
  ref = dcomplex(4.713103999979870, 0.017955195142080)/
    (pow(S00_ref(0,0), 1.5) * sqrt(S00_ref(2, 2)));
  EXPECT_C_EQ(ref, eri->At(0, 0, 0, 0, 2, 0, 0, 0));
	      
//  1  1  1  1  3  1  1  1   4.713103999979870   0.017955195142080


  //  EXPECT_C_EQ(eri->At(1, 1, 0, 0, 0, 0, 0, 0)/
  //	      (S11_ref(0, 0) * S00_ref(0, 0)), 4.499506255091513);

  //    1  1  1  1  3  3  3  3 249.194759129673912  -0.606119545098925
  //1  1  1  1  1  1  1  1   6.082101996924533   0.000000000000000
  //2  2  1  1  1  1  1  1   4.499506255091513   0.000000000000000
  //2  2  2  2  1  1  1  1   3.412356981545473   0.000000000000000
  //2  1  2  1  1  1  1  1   1.780420011334297   0.000000000000000

}

TEST(H2Plus, energy) {
  pSymmetryGroup Cs = SymmetryGroup::Cs();
  SymGTOs gtos(Cs);
  int Ap = Cs->GetIrrep("A'");
  //  int App= 1;

  int num_zeta(10);
  Vector3cd pos1(0, 0, +0.7);
  Vector3cd pos2(0, 0, -0.7);
  VectorXcd zeta(num_zeta);
  for(int n = 0; n < num_zeta; n++)
    zeta(n) = pow(2.0, n-5);

  SubSymGTOs sub(Cs);
  sub.AddXyz( pos1);
  sub.AddXyz( pos2);
  sub.AddNs(  Vector3i(0, 0, 0));
  sub.AddZeta(zeta);  
  sub.AddRds( Reduction(Ap, MatrixXcd::Ones(2, 1)));
  sub.SetUp();

  gtos.AddSub(sub);
  gtos.AddAtom(pos1, 1.0);
  gtos.AddAtom(pos2, 1.0);
  
  gtos.SetUp();

  BMatSet mat;
  gtos.CalcMat(&mat);
  const MatrixXcd& t = mat.GetMatrix("t", Ap, Ap);
  const MatrixXcd& v = mat.GetMatrix("v", Ap, Ap);
  const MatrixXcd& s = mat.GetMatrix("s", Ap, Ap);
  MatrixXcd h(t+v);

  MatrixXcd c;
  VectorXcd eig;
  generalizedComplexEigenSolve(h, s, &c, &eig);
  cout << eig << endl;

}
TEST(H2Plus, matrix) {
  
  // ==== Symmetry ====
  pSymmetryGroup D2h = SymmetryGroup::D2h();

  // ==== Sub ====
  cout << 1 << endl;
  SubSymGTOs sub1(D2h);
  sub1.AddXyz(Vector3cd(0, 0, +0.7));
  sub1.AddXyz(Vector3cd(0, 0, -0.7));
  sub1.AddNs( Vector3i( 0, 0, 0));
  VectorXcd z1(4); z1 << 2.013, 0.1233, 0.0411, 0.0137; sub1.AddZeta(z1);
  MatrixXcd c1_1(2, 1); c1_1 <<+1.0,+1.0; sub1.AddRds(Reduction(0, c1_1));
  MatrixXcd c1_2(2, 1); c1_2 <<+1.0,-1.0; sub1.AddRds(Reduction(1, c1_2));
  sub1.SetUp();

  cout << 2 << endl;
  SubSymGTOs sub2(D2h);
  sub2.AddXyz(Vector3cd(0, 0, +0.7));
  sub2.AddXyz(Vector3cd(0, 0, -0.7));
  sub2.AddNs( Vector3i( 0, 0, 1));
  VectorXcd z2(1); z2 << 1.0; sub2.AddZeta(z2);
  MatrixXcd C2_1(2, 1); C2_1 << +1,-1; sub2.AddRds(Reduction(0, C2_1));
  MatrixXcd C2_2(2, 1); C2_2 << +1,+1; sub2.AddRds(Reduction(1, C2_2));
  sub2.SetUp();

  cout << 3 << endl;
  SubSymGTOs sub3(D2h);
  sub3.AddXyz(Vector3cd(0, 0, 0));
  sub3.AddNs( Vector3i( 0, 0, 0));
  VectorXcd z3(1); z3 << dcomplex(0.011389, -0.002197); sub3.AddZeta(z3); 
  MatrixXcd C3_1(1, 1); C3_1 << 1; sub3.AddRds(Reduction(0, C3_1));
  sub3.SetUp();
  
  cout << 4 << endl;
  SubSymGTOs sub4(D2h);
  sub4.AddXyz(Vector3cd(0, 0, 0));
  sub4.AddNs( Vector3i( 2, 0, 0));
  sub4.AddNs( Vector3i( 0, 2, 0));
  sub4.AddNs( Vector3i( 0, 0, 2));
  VectorXcd z4(1); z4 << dcomplex(5.063464, -0.024632); sub4.AddZeta(z4);
  MatrixXcd C4_1(1, 3); C4_1 << -1,-1,+2; sub4.AddRds(Reduction(0, C4_1 ));
  sub4.SetUp();

  // ==== GTOs ====
  cout << "GTOs" << endl;
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
  BMatSet mat;
  gtos.CalcMat(&mat);
  //  IB2EInt *eri = new B2EIntMem(pow(gtos.size_basis(), 4));
  //  gtos.CalcERI(eri, 1);

  // copied from ~/calc/ccolumbus
  dcomplex s00(2.2781505348887450);
  dcomplex s11(3.7723621068772371);
  dcomplex s01(1.1443980248362513);
  dcomplex s44(2.7205973097946861);
  dcomplex s55(1);
  dcomplex s54(-2.9723198989659500*0.001,  1.0081573799329419*0.001);
  dcomplex s16(3.5564593409758150*0.001,  2.8861860781924722 * 0.00001);
  dcomplex s66(11.999999999999995);
  EXPECT_C_EQ(s01/(sqrt(s00)*sqrt(s11)), mat.GetMatrix("s", 0, 0)(0, 1));
  EXPECT_C_EQ(s54/(sqrt(s55)*sqrt(s44)), mat.GetMatrix("s", 0, 0)(5, 4));
  EXPECT_C_EQ(s16/(sqrt(s11)*sqrt(s66)), mat.GetMatrix("s", 0, 0)(1, 6));

  s00 = 1.7218494651112561;
  s11 = 0.22763789312276095;
  s01 = 0.12974083877243192;
  EXPECT_C_EQ(s01/sqrt(s11*s00), mat.GetMatrix("s", 1, 1)(0, 1));

  dcomplex v01 = -0.28905219384317249;
  EXPECT_C_EQ(v01/sqrt(s11*s00), mat.GetMatrix("v", 1, 1)(0, 1));
  
}


int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}

