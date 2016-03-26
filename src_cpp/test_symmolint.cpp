#include <gtest/gtest.h>
#include <Eigen/Core>

#include "gtest_plus.hpp"
#include "eigen_plus.hpp"

#include "mol_func.hpp"

#include "symmolint.hpp"

using namespace std;
using namespace l2func;
using namespace Eigen;

// other calculation:
// 2016/3/22


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

TEST(SymmetryGroup, C1) {

  /*
    C1 Symmetry
    +--------+
    | C1 | E |
    +--------+
    | A  |+1 |
    +--------+
   */

  SymmetryGroup C1(SymmetryGroup_C1());
  Irrep A = 0; 
  
  EXPECT_EQ(1, C1.order());
  EXPECT_TRUE( C1.IncludeScalar_2(A, A));
  
  EXPECT_TRUE( C1.IncludeZ_2(A, A));

  EXPECT_ANY_THROW(C1.IncludeScalar_2(0, 1));
  EXPECT_ANY_THROW(C1.IncludeScalar_2(1, 0));
  EXPECT_ANY_THROW(C1.IncludeScalar_2(-1, 0));

}
TEST(SymmetryGroup, Cs) {

  /*
    Cs Symmetry
    +---------------+
    | Cs |  E | C_2 |
    +----+----+-----+----------+
    | A' | +1 | +1  | x,y,Rz   |
    +----+----+-----+----------+
    | A''| +1 | -1  | z,Rx,Rz  |
    +----+----+-----+----------+
   */

  SymmetryGroup Cs(SymmetryGroup_Cs());
  Irrep Ap = 0;   //Cs.GetIrrep(0);
  Irrep App = 1; //Cs.GetIrrep(1);
  
  EXPECT_EQ(2, Cs.order());
  EXPECT_TRUE( Cs.IncludeScalar_2(Ap, Ap));
  EXPECT_FALSE(Cs.IncludeScalar_2(Ap, App));
  EXPECT_FALSE(Cs.IncludeScalar_2(App, Ap));
  EXPECT_TRUE( Cs.IncludeScalar_2(App, App));

  EXPECT_FALSE( Cs.IncludeZ_2(Ap, Ap));
  EXPECT_TRUE(Cs.IncludeZ_2(Ap, App));
  EXPECT_TRUE(Cs.IncludeZ_2(App, Ap));
  EXPECT_FALSE( Cs.IncludeZ_2(App, App));

  EXPECT_ANY_THROW(Cs.IncludeScalar_2(0, 3));
  EXPECT_ANY_THROW(Cs.IncludeScalar_2(2, 0));
  EXPECT_ANY_THROW(Cs.IncludeScalar_2(-1, 0));
}
TEST(SymmetryGroup, C2h) {

  /*
    C2h Symmetry
    +---------------+----+----+
    |    |  E | C_2 | i  | sh |      
    +----+----+-----+----+----+----------+
    | Ag | +1 | +1  | +1 | +1 |   Rz     |
    +----+----+-----+----+----+----------+
    | Bg | +1 | -1  | +1 | -1 |  Rx, Rz  |
    +----+----+-----+---------+----------+
    | Au | +1 | +1  | -1 | -1 |    z     |
    +----+----+-----+---------+----------+
    | Bu | +1 | -1  | -1 | +1 |   x, y   |
    +----+----+-----+---------+----------+

    
    |    | Ag | Bg | Au | Bu |
    --------------------------
    | Ag | Ag | Bg | Au | Bu |
    | Bg | Bg | Ag | Bu | Au |
    | Au | Au | Bu | Ag | Bg |
    | Bu | Bu | Au | Bg | Ag | 
   */
/*
  SymmetryGroup Cs(SymmetryGroup_C2h());
  Irrep Ap = 0;   //Cs.GetIrrep(0);
  Irrep App = 1; //Cs.GetIrrep(1);
  
  EXPECT_EQ(2, Cs.order());
  EXPECT_TRUE( Cs.IncludeScalar_2(Ap, Ap));
  EXPECT_FALSE(Cs.IncludeScalar_2(Ap, App));
  EXPECT_FALSE(Cs.IncludeScalar_2(App, Ap));
  EXPECT_TRUE( Cs.IncludeScalar_2(App, App));

  EXPECT_FALSE( Cs.IncludeZ_2(Ap, Ap));
  EXPECT_TRUE(Cs.IncludeZ_2(Ap, App));
  EXPECT_TRUE(Cs.IncludeZ_2(App, Ap));
  EXPECT_FALSE( Cs.IncludeZ_2(App, App));

  EXPECT_ANY_THROW(Cs.IncludeScalar_2(0, 3));
  EXPECT_ANY_THROW(Cs.IncludeScalar_2(2, 0));
  EXPECT_ANY_THROW(Cs.IncludeScalar_2(-1, 0));
*/
}

TEST(SubSymGTOs, AddZeta) {

  SubSymGTOs sub;
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
TEST(SymGTOs, at_r_ylm_lin) {

  SymGTOs gtos(SymmetryGroup_C1());
  VectorXcd zeta(1); zeta << 1.1;
  gtos.AddSub(Sub_s(0, Vector3cd(0, 0, 2), zeta));
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

  SymGTOs gtos(SymmetryGroup_C1());

  /*
  VectorXcd zeta(1); zeta << 1.2;
  gtos.AddSub(Sub_s(0, Vector3cd(0.0, 0.0, 0.0), zeta));
  VectorXcd cs(1); cs << 0.2;
  */

  VectorXcd zeta(2); zeta << 1.2, 1.4;
  gtos.AddSub(Sub_s(0, Vector3cd(0.0, 0.0, 0.0), zeta));
  VectorXcd cs(2); cs << 0.2, 1.1;

  gtos.SetUp();
  test_SymGTOs_at_r_ylm_nd(gtos, 1.1, cs, "s-GTO on center");

}
TEST(SymGTOs, at_r_ylm_nd) {

  SymGTOs gtos(SymmetryGroup_C1());
  // VectorXcd zeta(2); zeta << 1.2, 1.1;
  // gtos.AddSub(Sub_s(0, Vector3cd(0.1, 0.2, 0.3), zeta));
  // VectorXcd cs(2); cs << 0.2, 0.4;

  VectorXcd zeta(1); zeta << 1.2;
  gtos.AddSub(Sub_s(0, Vector3cd(0.0, 0.0, 2.1), zeta));
  VectorXcd cs(1); cs << 1.0;

  gtos.SetUp();
  test_SymGTOs_at_r_ylm_nd(gtos, 1.1, cs, "s-GTO on (0.1, 0.2, 0.3)");

}
TEST(SymGTOs, Create) {

  SymmetryGroup sym_group(SymmetryGroup_Cs());
  SymGTOs gtos(sym_group);

  // -- Symmetries
  Irrep lm00(0);
  Irrep lm10(1);

  // -- s-GTO at Center, Sym = (0, 0)
  SubSymGTOs sgto_cen;
  VectorXcd zetas(10);
  for(int n = 0; n < 10; n++)
    zetas(n) = pow(2.0, n-5);
  sgto_cen.AddXyz( Vector3cd(0, 0, 0));
  sgto_cen.AddNs(  Vector3i(0, 0, 0));
  sgto_cen.AddZeta(zetas);
  sgto_cen.AddRds( Reduction(lm00, MatrixXcd::Ones(1, 1)));
  gtos.AddSub(sgto_cen);

  // -- s-GTO at z-axis, sym = (1, 0)
  SubSymGTOs pgto_z;
  VectorXcd zetas2(6);
  for(int n = 0; n < 6; n++)
    zetas2[n] = pow(1.6, n-3);
  pgto_z.AddXyz(Vector3cd(0, 0, +1));
  pgto_z.AddXyz(Vector3cd(0, 0, -1));
  pgto_z.AddNs( Vector3i( 0, 0, +1));  
  pgto_z.AddZeta(zetas2);
  MatrixXcd c1(2, 1); c1 << 1.0, -1.0;
  pgto_z.AddRds(Reduction(lm10, c1));
  MatrixXcd c2(2, 1); c2 << 1.0, +1.0;
  pgto_z.AddRds(Reduction(lm00, c2));
  gtos.AddSub(pgto_z);

  EXPECT_EQ(1,  gtos.subs[0].size_at());
  EXPECT_EQ(1,  gtos.subs[0].size_pn());
  EXPECT_EQ(1,  gtos.subs[0].size_cont());
  EXPECT_EQ(10, gtos.subs[0].size_zeta());

  EXPECT_EQ(2,  gtos.subs[1].size_at());
  EXPECT_EQ(1,  gtos.subs[1].size_pn());
  EXPECT_EQ(2,  gtos.subs[1].size_cont());
  EXPECT_EQ(6,  gtos.subs[1].size_zeta());

  BMatSet mat;
  gtos.CalcMat(&mat);
  EXPECT_C_EQ(1.0, mat.GetValue("s", lm00, lm00, 0, 0));
}
TEST(SymGTOs, one_center_S_gto) {

  dcomplex zeta1(1.1, -0.2);
  dcomplex zeta2(0.3, -0.01);

  SymGTOs gtos(SymmetryGroup_C1()); 
  VectorXcd zetas(2); zetas << zeta1, zeta2;
  gtos.AddSub(Sub_s(0, Vector3cd(0, 0, 0), zetas));
  gtos.AddAtom(Vector3cd(0, 0, 0), 1.0);
  gtos.SetUp();

  BMatSet res;
  gtos.CalcMat(&res);

  dcomplex c0 = GTONormConst(0, 0, 0, zeta1);
  dcomplex c1 = GTONormConst(0, 0, 0, zeta2);
  dcomplex s01 = GTOOverlap(0, 0, 0, 0.0, 0.0, 0.0, zeta1,
			    0, 0, 0, 0.0, 0.0, 0.0, zeta2);
  dcomplex t00 = GTOKinetic(0, 0, 0, 0.0, 0.0, 0.0, zeta1,
			    0, 0, 0, 0.0, 0.0, 0.0, zeta1);
  dcomplex t01 = GTOKinetic(0, 0, 0, 0.0, 0.0, 0.0, zeta1,
			    0, 0, 0, 0.0, 0.0, 0.0, zeta2);
  dcomplex t11 = GTOKinetic(0, 0, 0, 0.0, 0.0, 0.0, zeta2,
			    0, 0, 0, 0.0, 0.0, 0.0, zeta2);
  dcomplex v00 = GTONuclearAttraction(0, 0, 0, 0.0, 0.0, 0.0, zeta1,
				      0, 0, 0, 0.0, 0.0, 0.0, zeta1,
				      0.0, 0.0, 0.0);
  dcomplex v11 = GTONuclearAttraction(0, 0, 0, 0.0, 0.0, 0.0, zeta2,
				      0, 0, 0, 0.0, 0.0, 0.0, zeta2,
				      0.0, 0.0, 0.0);
  dcomplex v01 = GTONuclearAttraction(0, 0, 0, 0.0, 0.0, 0.0, zeta1,
				      0, 0, 0, 0.0, 0.0, 0.0, zeta2,
				      0.0, 0.0, 0.0);

  const MatrixXcd& S = res.GetMatrix("s", 0, 0);
  EXPECT_C_EQ(1.0, S(0, 0));
  EXPECT_C_EQ(1.0, S(1, 1));
  EXPECT_C_EQ(S(1, 0), S(0, 1));
  EXPECT_C_EQ(s01*c0*c1, S(0, 1));

  const MatrixXcd& T = res.GetMatrix("t", 0, 0);
  EXPECT_C_EQ(t00*c0*c0, T(0, 0));
  EXPECT_C_EQ(t11*c1*c1, T(1, 1));
  EXPECT_C_EQ(T(1, 0), T(0, 1));
  EXPECT_C_EQ(t01*c0*c1, T(0, 1));  

  const MatrixXcd& V = res.GetMatrix("v", 0, 0);
  EXPECT_C_EQ(v00*c0*c0, V(0, 0));
  EXPECT_C_EQ(v11*c1*c1, V(1, 1));
  EXPECT_C_EQ(V(1, 0), V(0, 1));
  EXPECT_C_EQ(v01*c0*c1, V(0, 1));    

}
TEST(SymGTOs, two_center_SP) {

  dcomplex zeta1(1.1, -0.2);
  dcomplex x1(0), y1(0), z1(0);
  dcomplex zeta2(0.3, -0.01);
  dcomplex x2(1.1), y2(1.2), z2(1.3);

  SymGTOs gtos(SymmetryGroup_C1()); 

  VectorXcd zeta_cen(1); zeta_cen << zeta1;
  gtos.AddSub(Sub_pz(0, Vector3cd(x1, y1, z1), zeta_cen));

  VectorXcd zeta_gh(1); zeta_gh << zeta2;
  gtos.AddSub(Sub_s(0, Vector3cd(x2, y2, z2), zeta_gh));

  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos.SetAtoms(xyzq);

  BMatSet res;
  gtos.CalcMat(&res);

  dcomplex s_c =    GTOOverlap(0, 0, 1, x1, y1, z1,  zeta1,
			       0, 0, 1, x1, y1, z1,  zeta1);
  dcomplex s_g00  = GTOOverlap(0, 0, 0, x2, y2, z2, zeta2,
			       0, 0, 0, x2, y2, z2, zeta2);
  dcomplex s_cg0  = GTOOverlap(0, 0, 1, x1, y1, z1, zeta1,
			       0, 0, 0, x2, y2, z2, zeta2);

  const MatrixXcd& S = res.GetMatrix("s", 0, 0);
  EXPECT_EQ(2, S.rows());
  EXPECT_EQ(2, S.cols());

  dcomplex c_g = 1.0 / sqrt(s_g00);
  dcomplex c_c = 1.0 / sqrt(s_c);
  EXPECT_C_EQ(1.0, c_c * c_c * s_c);
  EXPECT_C_EQ(1.0, c_g * c_g * s_g00);
  EXPECT_C_EQ(1.0, S(0, 0));
  EXPECT_C_EQ(1.0, S(1, 1));
  EXPECT_C_EQ(S(1, 0), S(0, 1));
  EXPECT_C_EQ(c_c, gtos.subs[0].rds[0].coef_iz(0));
  EXPECT_C_EQ(c_g, gtos.subs[1].rds[0].coef_iz(0));
  EXPECT_C_EQ(c_g * c_c * s_cg0, S(0, 1));

}
TEST(SymGTOs, add_atom) {

  SymGTOs gtos(SymmetryGroup_C1()); 
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

  SymGTOs gtos_full(SymmetryGroup_Cs());
  SymGTOs gtos_1(   SymmetryGroup_Cs());
  SymGTOs gtos_2(   SymmetryGroup_Cs());
  Irrep Ap = Cs_Ap();
  Irrep App= Cs_App();
  dcomplex z_gh(1.1);

  // -- A' symmetry --
  VectorXcd zeta_h = VectorXcd::Zero(3); zeta_h << 0.4, 1.0, 2.0;
  gtos_full.AddSub(Sub_s(Ap, Vector3cd(0, 0, 0), zeta_h));
  gtos_1.AddSub(Sub_s(Ap, Vector3cd(0, 0, 0), zeta_h));


  // -- A'' symmetry, Center --
  VectorXcd zeta_gh = VectorXcd::Zero(3); zeta_gh << 0.6, 1.2, 2.2;
  gtos_full.AddSub(Sub_TwoSGTO(SymmetryGroup_Cs(), App,
			   Vector3cd(0.0, 0.0, z_gh), zeta_gh));
  VectorXcd zeta_cen(2); zeta_cen << 0.2, 1.2;
  gtos_full.AddSub(Sub_pz(App, Vector3cd(0, 0, 0), zeta_cen));
  
  gtos_2.AddSub(Sub_TwoSGTO(SymmetryGroup_Cs(), App,
			    Vector3cd(0.0, 0.0, z_gh), zeta_gh));
  gtos_2.AddSub(Sub_pz(App, Vector3cd(0, 0, 0), zeta_cen));


  // -- potential --
  MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
  gtos_full.SetAtoms(xyzq);
  gtos_1.SetAtoms(xyzq);
  gtos_2.SetAtoms(xyzq);


  // -- calculate matrix --
  BMatSet mat_full, mat_12;
  gtos_full.SetUp();
  gtos_1.SetUp();
  gtos_2.SetUp();
  gtos_full.CalcMat(&mat_full);
  gtos_1.CalcMatOther(gtos_2, false, &mat_12);

  
  // -- compare results --
  EXPECT_MATXCD_EQ(mat_full.GetMatrix("z", Ap, App),
		   mat_12.GetMatrix(  "z", Ap, App));
  
}

class SP_GTO : public ::testing::Test {
protected:
  SymGTOs* gtos;
  dcomplex zeta1, zeta2;
  dcomplex z_gh;
  VectorXcd zeta_h;
  BMatSet res;

  dcomplex c_g, c_c, c_h0;

public:
  SP_GTO() {
    gtos = new SymGTOs(SymmetryGroup_Cs());
    zeta1 = dcomplex(1.1, -0.2);
    zeta2 = dcomplex(0.3, -0.01);

    // -- A' symmetry --
    zeta_h = VectorXcd::Zero(3); zeta_h << 0.4, 1.0, 2.0;
    gtos->AddSub(Sub_s(Cs_Ap(), Vector3cd(0, 0, 0), zeta_h));

    // -- A'' symmetry, Center --
    VectorXcd zeta_cen(1); zeta_cen << zeta1;
    gtos->AddSub(Sub_pz(Cs_App(), Vector3cd(0, 0, 0), zeta_cen));

    // -- A'' symmetry, Ghost --
    z_gh = dcomplex(0.8);
    VectorXcd zeta_gh(1); zeta_gh << zeta2;
    gtos->AddSub(Sub_TwoSGTO(SymmetryGroup_Cs(), Cs_App(),
			     Vector3cd(0.0, 0.0, z_gh), zeta_gh));

    // -- potential --
    MatrixXcd xyzq(4, 1); xyzq << 0.0, 0.0, 0.0, 1.0;
    gtos->SetAtoms(xyzq);

    // -- calc --
    gtos->CalcMat(&res);

    // -- reference --
    dcomplex s_g00  = GTOOverlap(0, 0, 0, 0.0, 0.0, z_gh, zeta2,
				 0, 0, 0, 0.0, 0.0, z_gh, zeta2);
    dcomplex s_g01  = GTOOverlap(0, 0, 0, 0.0, 0.0, z_gh, zeta2,
				 0, 0, 0, 0.0, 0.0, -z_gh, zeta2);
    c_h0 = GTONormConst(0, 0, 0, zeta_h(0));
    c_c = GTONormConst(0, 0, 1, zeta1);
    c_g  = 1.0 / sqrt(s_g00*2.0 - s_g01*2.0);
    
  }
  ~SP_GTO() {
    delete gtos;
  }
};
TEST_F(SP_GTO, matrix) {

  dcomplex s_c = GTOOverlap(0, 0, 1, 0.0, 0.0, 0.0,  zeta1,
			    0, 0, 1, 0.0, 0.0, 0.0,  zeta1);
  dcomplex s_cg0  = GTOOverlap(0, 0, 1, 0.0, 0.0, 0.0,  zeta1,
			       0, 0, 0, 0.0, 0.0, z_gh, zeta2);

  dcomplex t_c = GTOKinetic(0, 0, 1, 0.0, 0.0, 0.0, zeta1,
			    0, 0, 1, 0.0, 0.0, 0.0, zeta1);
  dcomplex t_cg = GTOKinetic(0, 0, 1, 0.0, 0.0, 0.0,  zeta1,
			    0, 0, 0, 0.0, 0.0, z_gh, zeta2);
  dcomplex t_g00 = GTOKinetic(0, 0, 0, 0.0, 0.0, z_gh, zeta2,
			      0, 0, 0, 0.0, 0.0, z_gh, zeta2);
  dcomplex t_g01 = GTOKinetic(0, 0, 0, 0.0, 0.0, z_gh, zeta2,
			      0, 0, 0, 0.0, 0.0,-z_gh, zeta2);

  dcomplex z_h0g = GTODipZ(0, 0, 0, 0.0, 0.0, 0.0, zeta_h(0),
			   0, 0, 0, 0.0, 0.0, z_gh, zeta2);
  dcomplex z_h0c = GTODipZ(0, 0, 0, 0.0, 0.0, 0.0, zeta_h(0),
			   0, 0, 1, 0.0, 0.0, 0.0, zeta1);  

  dcomplex v_c = GTONuclearAttraction(0, 0, 1, 0.0, 0.0, 0.0, zeta1,
					0, 0, 1, 0.0, 0.0, 0.0, zeta1,
					0.0, 0.0, 0.0);
  dcomplex v_cg = GTONuclearAttraction(0, 0, 1, 0.0, 0.0, 0.0, zeta1,
				       0, 0, 0, 0.0, 0.0, z_gh, zeta2,
				      0.0, 0.0, 0.0);
  dcomplex v_g00 = GTONuclearAttraction(0, 0, 0, 0.0, 0.0, z_gh, zeta2,
					0, 0, 0, 0.0, 0.0, z_gh, zeta2,
					0.0, 0.0, 0.0);
  dcomplex v_g01 = GTONuclearAttraction(0, 0, 0, 0.0, 0.0, z_gh, zeta2,
					0, 0, 0, 0.0, 0.0,-z_gh, zeta2,
					0.0, 0.0, 0.0);

  const MatrixXcd& S = res.GetMatrix("s", 1, 1);

  EXPECT_EQ(2, S.rows());
  EXPECT_EQ(2, S.cols());
  
  EXPECT_C_EQ(1.0, c_c * c_c * s_c);
  EXPECT_C_EQ(1.0, S(0, 0));
  EXPECT_C_EQ(1.0, S(1, 1));
  EXPECT_C_EQ(S(1, 0), S(0, 1));
  EXPECT_C_EQ(c_c, gtos->subs[1].rds[0].coef_iz(0));
  EXPECT_C_EQ(c_g, gtos->subs[2].rds[0].coef_iz(0));
  EXPECT_C_EQ(c_g * c_c * s_cg0 * 2.0, S(0, 1));

  const MatrixXcd& T = res.GetMatrix("t", 1, 1);
  EXPECT_C_EQ(t_c*c_c*c_c, T(0, 0));
  EXPECT_C_EQ(2.0*(t_g00-t_g01)*c_g*c_g, T(1, 1));
  EXPECT_C_EQ(c_g * c_c * t_cg * 2.0, T(0, 1));
  EXPECT_C_EQ(T(1, 0), T(0, 1));  

  const MatrixXcd& V = res.GetMatrix("v", 1, 1);
  EXPECT_C_EQ(v_c*c_c*c_c, V(0, 0));
  EXPECT_C_EQ(2.0*(v_g00-v_g01)*c_g*c_g, V(1, 1));
  EXPECT_C_EQ(c_g * c_c * v_cg * 2.0, V(0, 1));
  EXPECT_C_EQ(V(1, 0), V(0, 1));  

  const MatrixXcd& Z = res.GetMatrix("z", 0, 0);
  EXPECT_MATXCD_EQ(MatrixXcd::Zero(3, 3), Z);
  const MatrixXcd& Z11 = res.GetMatrix("z", 1, 1);
  EXPECT_MATXCD_EQ(MatrixXcd::Zero(2, 2), Z11);
  const MatrixXcd& Z01 = res.GetMatrix("z", 0, 1);
  EXPECT_C_EQ(z_h0c * c_h0 * c_c, Z01(0, 0));
  EXPECT_C_EQ(z_h0g * c_h0 * c_g * 2.0, Z01(0, 1));

}
TEST_F(SP_GTO, at_r_0) {

  VectorXcd cs0_ibasis(3);
  cs0_ibasis << 1.1, 0.0, 0.0;

  VectorXcd rs(2); rs << 1.1, 1.4;
  VectorXcd vs, dvs; gtos->AtR_Ylm(0, 0, 0, cs0_ibasis, rs, &vs, &dvs);
  
  for(int i = 0; i < 2; i++) {
    dcomplex r(rs[i]);
    dcomplex c_expz(cs0_ibasis(0) * c_h0 * exp(-zeta_h[0]*r*r)*sqrt(4.0*M_PI));
    EXPECT_C_EQ(c_expz * r, vs[i]) << i;
    EXPECT_C_EQ(c_expz * (1.0 -2.0*zeta_h[0]*r*r), dvs[i]) << i;
  }  

}

TEST(H2Plus, energy) {
  SymGTOs gtos(SymmetryGroup_Cs());
  int num_zeta(10);
  Vector3cd pos1(0, 0, +0.7);
  Vector3cd pos2(0, 0, -0.7);
  VectorXcd zeta(num_zeta);
  for(int n = 0; n < num_zeta; n++)
    zeta(n) = pow(2.0, n-5);
  SubSymGTOs sub;
  sub.AddXyz( pos1);
  sub.AddXyz( pos2);
  sub.AddNs(  Vector3i(0, 0, 0));
  sub.AddZeta(zeta);  
  sub.AddRds( Reduction(Cs_Ap(), MatrixXcd::Ones(2, 1)));
  gtos.AddSub(sub);
  gtos.AddAtom(pos1, 1.0);
  gtos.AddAtom(pos2, 1.0);
  
  gtos.SetUp();

  BMatSet mat;
  gtos.CalcMat(&mat);
  Irrep Ap = Cs_Ap();
  const MatrixXcd& t = mat.GetMatrix("t", Ap, Ap);
  const MatrixXcd& v = mat.GetMatrix("v", Ap, Ap);
  const MatrixXcd& s = mat.GetMatrix("s", Ap, Ap);
  MatrixXcd h(t+v);

  MatrixXcd c;
  VectorXcd eig;
  generalizedComplexEigenSolve(h, s, &c, &eig);
  cout << eig << endl;

}

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
