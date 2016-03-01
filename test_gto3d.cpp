/*
  Unit test for cints.cpp, 
*/

#include "utils.hpp"
#include "cints.hpp"
#include "angmoment.hpp"
#include "gto3d.hpp"
#include "gto3dset.hpp"
#include <gtest/gtest.h>


using namespace std;
using namespace l2func;

int lm(int l, int m) {
  return lm_index(l, m);
}

TEST(CINTS, First) {
  EXPECT_TRUE(true);
}

TEST(CINTS, lgamma) {
  C val(-0.008197780565406, -0.05732294041672);
  C calc(lgamma(C(1.0, 0.1)));
  double eps(0.00000001);
  EXPECT_C_NEAR(val, calc, eps);
}
TEST(CINTS,  overlap) {
  double eps(0.000000001);
  C s = overlap(C(1.1,0.2), 0, 0, 0,
		0.1,  0.2,   0.3,
		C(1.3,0.5), 0, 0, 0,
		0.1, -0.22, -0.3);
  EXPECT_C_NEAR(s, C(0.88910293518825, -0.5004029971544), eps);
  //  EXPECT_C_EQ(s, C(0.88910293518825, -0.5004029971544));

  C s1 = overlap(C(1.1,0.2), 2, 1, 0,
		 0.1,  C(0.2, 0.03),   0.3,
		 C(1.3,0.02), 0, 1, 3,
		 0.1, -0.22, C(-0.3,-0.7));
  EXPECT_C_NEAR(s1, C(0.012403718856620667, 0.003134556944411788), eps);
  
  C s2 = overlap(C(0.2,0.3), 2, 1, 1,
		 0.0, 0.0, 0.0,
		 C(0.2,-0.7), 2,1,2,
		 0.0, C(0.2,-0.1), 0.3);
  EXPECT_C_NEAR(s2, C(-14.3092, -3.67256), 0.00005);
  
}
TEST(CINTS,  kinetic) {
  double eps(0.000001);
  C s = kinetic(C(1.1,0.2), 0, 0, 0,
		0.1,  0.2,   0.3,
		C(1.3,0.5), 0, 0, 0,
		0.1, -0.22, -0.3);
  C se(1.422908003563741695, -0.476445156336799588);
  EXPECT_C_NEAR(s, se, eps);


/*
  C s1 = kinetic(C(1.1,0.2), 2, 1, 0,
		 0.1,  C(0.2, 0.03),   0.3,
		 C(1.3,0.02), 0, 1, 3,
		 0.1, -0.22, C(-0.3,-0.7));
  //  EXPECT_NEAR_COMPLEX(s1, C(0.0143520521644942813,0.0108902064709103336), eps);
  
  C s2 = kinetic( C(0.1+0.3, 0.3), 1, 0, 1,
		  C(0.3, 0.1), 0.7, 0.0,
		  C(0.3, -0.1), 1, 1, 0,
		  C(0.3, 0.1), C(0.0,-0.0), 0.3);
  //  EXPECT_NEAR_COMPLEX(s2, C(0.404777, -0.50371), 0.000005);
  
  C s3 = kinetic( C(0.1+0.3, 0.3), 2, 1, 1,
		  C(0.3, 0.1), 0.7, 0.0,
		  C(0.3, -0.1), 2, 2, 1,
		  C(0.3, 0.1), C(0.0,-0.0), 0.3);
  //  EXPECT_NEAR_COMPLEX(s3, C(0.404777, -0.50371), 0.000005);
  
  C s3inv = kinetic( C(0.3, -0.1), 2, 2, 1,
		     C(0.3, 0.1), C(0.0,-0.0), 0.3,
		     C(0.1+0.3, 0.3), 2, 1, 1,
		     C(0.3, 0.1), 0.7, 0.0);

  //  EXPECT_NEAR_COMPLEX(s3inv, C(0.404777, -0.50371), eps);
*/
}
TEST(CINTS, product_center_1D) {


  /*
  C a1(1.1,0.2), a2(1.3,0.5), a3(1.1,0.2), a4(1.3,0.02);
  
  C rys   = product_center_1D_rys(1.1, 1.2,1.3,1.4);
  C cints = product_center_1D(1.1, 1.2,1.3,1.4);
  EXPECT_NEAR_COMPLEX(rys, cints, eps);

  
  rys   = product_center_1D_Rys(a1, 0.1, a2, 0.3);
  cints = product_center_1D(a1, 0.1, a2, 0.3);
  
  EXPECT_NEAR_COMPLEX(rys, cints, eps);
  */
}
/*
TEST(CINTS, ee_time) {

  C one(1.0, 0.0);
  C a1(1.1,0.2), a2(1.3,0.5), a3(1.1,0.2), a4(1.3,0.02);

  UniqueTimer timer;
  C ti;
  int stepNum = 60 * 60 * 60 * 60 / 8;

  timer.Start("Rys");
  for(int i = 0; i < stepNum; i++) {
    ti = coulomb_repulsion_Rys
      (0.1,  0.2,   0.3, one,
       0,1,0, a1,
       0.1, -0.22, -0.3, one,
       0, 0, 0, a2, 
       0.1,  C(0.2, 0.03), 0.3, one,
       1, 0, 0, a3,
       0.1, -0.22, C(-0.3,-0.7), one, 
       0, 1, 0, a4);
  }  
  timer.End("Rys");
  timer.Start("Huzinaga");
  for(int i = 0; i < stepNum; i++) {
    ti = coulomb_repulsion
      (0.1,  0.2,   0.3, one,
       0,1,0, a1,
       0.1, -0.22, -0.3, one,
       0, 0, 0, a2, 
       0.1,  C(0.2, 0.03), 0.3, one,
       1, 0, 0, a3,
       0.1, -0.22, C(-0.3,-0.7), one, 
       0, 1, 0, a4);
  }
  timer.End("Huzinaga");

  timer.Display();

}
TEST(CINTS,  ee ) {

  C one(1.0, 0.0);
  C a1(1.1,0.2), a2(1.3,0.5), a3(1.1,0.2), a4(1.3,0.02);
  
  C sc = coulomb_repulsion(
		    0.1,  0.2,   0.3, one,
		    0,0,0, a1,
		    0.1, -0.22, -0.3, one,
		    0, 0, 0, a2, 
		    0.1,  C(0.2, 0.03), 0.3, one,
		    2, 1, 0, a3,
		    0.1, -0.22, C(-0.3,-0.7), one, 
		    0, 1, 3, a4);

  C sc_Rys = coulomb_repulsion_Rys
    (0.1,  0.2,   0.3, one,
     0,0,0, a1,
     0.1, -0.22, -0.3, one,
     0, 0, 0, a2, 
     0.1,  C(0.2, 0.03), 0.3, one,
     2, 1, 0, a3,
     0.1, -0.22, C(-0.3,-0.7), one, 
     0, 1, 3, a4);

  
  C se(0.010125449275747, -0.000826749400705);
  EXPECT_NEAR_COMPLEX(sc, se, eps);
  EXPECT_NEAR_COMPLEX(sc_Rys, se, eps);

  sc = coulomb_repulsion(
			 0.1,  0.2,   0.3, one,
			 0,0,0, a1,
			 0.1,  0.2,   0.3, one,
			 0,0,0, a1,
			 0.1,  0.2,   0.3, one,
			 0,0,0, a1,
			 0.1,  0.2,   0.3, one,
			 0,0,0, a1);
  sc_Rys = coulomb_repulsion_Rys(
			 0.1,  0.2,   0.3, one,
			 0,0,0, a1,
			 0.1,  0.2,   0.3, one,
			 0,0,0, a1,
			 0.1,  0.2,   0.3, one,
			 0,0,0, a1,
			 0.1,  0.2,   0.3, one,
			 0,0,0, a1);			 

  se = C(2.979971963342084066,-1.438144700585544071);
  
  EXPECT_NEAR_COMPLEX(sc, se, 0.000005);
  EXPECT_NEAR_COMPLEX(sc_Rys, se, 0.000005);

}
TEST(CINTS,  ee_property) {

  C x=0.0; C y = 0.0; C z=0.0;
  C x1=0.0; C y1=0.0;
  C oe = 1.0;

  cout << endl;
  cout << "coulomb repulsion test for property " << endl;
  cout << "r is distance between u and v;" << endl;
  cout << "r, (uu|vv), (uv|uv)" << endl;

  for(C z1 = 0.0; real(z1) < 15.0; z1 += 1.0) {
    C uuvv = coulomb_repulsion(
			       x, y, z, 1.0, 0, 0, 0, oe,
			       x, y, z, 1.0, 0, 0, 0, oe,
			       x1,y1,z1,1.0, 0, 0, 0, oe,
			       x1,y1,z1,1.0, 0, 0, 0, oe);
    C uvuv = coulomb_repulsion(
			       x, y, z, 1.0, 0, 0, 0, oe,
			       x1,y1,z1,1.0, 0, 0, 0, oe,
			       x, y, z, 1.0, 0, 0, 0, oe,			   
			       x1,y1,z1,1.0, 0, 0, 0, oe);

    cout << z1 << "  " << uuvv << "  " << uvuv << endl;


    //      Seeing the result, (uu|vv) decay slowly while (uv|uv) decay rapidly.
    // So, (pq|rs) = Int p(1)q(1)r(2)s(2) .

  }

}
*/
TEST(CINTS,  nuclear_attraction) {
  double eps(0.000000001);
  C sc = nuclear_attraction(0.1,  C(0.2, 0.03), 0.3, 1.0, 
			    2, 1, 0, C(1.1,0.2),
			    0.1, -0.22, C(-0.3,-0.7), 1.0,
			    0, 1, 3, C(1.3,0.02),
			    -0.1, C(-0.1,0.1), 0.3);
  C se(-0.0111712963403, -0.0039848461450);
  EXPECT_C_NEAR(sc, se, eps);
}
TEST(CartGTO, SetScalarProd) {

  CartGTO<C, d3> g(1.0, i3(1, 1, 2), d3(0.0, 0.0, 0.0), C(1.0, 0.2));
  g.SetScalarProd(1.2);
  EXPECT_C_EQ(1.2, g.c());

}
TEST(CartGTO, SetNormalize) {
  
  double eps(0.000000001);
  CartGTO<C, d3> g(1.0, i3(1, 1, 2), d3(0.3, 0.1, 0.2), C(1.0, 0.2));
  g.SetNormalize();
  EXPECT_C_NEAR(1.0, CIP(g, g), eps);

}
TEST(CartGTO, Overlap) {
  double eps(0.00005);
  double c1(1.1);
  double c2(1.2);
  CartGTO<C, c3> g1(c1, i3(2, 1, 1), c3(0.0, 0.0, 0.0),        C(0.2, 0.3));
  CartGTO<C, c3> g2(c2, i3(2, 1, 2), c3(0.0, C(0.2,-0.1), 0.3), C(0.2, -0.7));
  C s = CIP(g1, g2);
  
  EXPECT_C_NEAR(s, c1*c2*C(-14.3092, -3.67256), eps);  
}
TEST(CargGTO, Kinetic) {

  double eps(0.000000001);
  double c1(1.1);
  double c2(1.2);
  CartGTO<C, c3> g1(c1,
		   i3(0, 0, 0),
		   c3(0.1, 0.2, 0.3),
		   C(1.1, 0.2));
  CartGTO<C, c3> g2(c2, 
		   i3(0, 0, 0),
		   c3(0.1, -0.22, -0.3),
		   C(1.3, 0.5));

  C s = CIP(g1, OpKE<C, c3>(), g2);
  C se(1.422908003563741695, -0.476445156336799588);
  EXPECT_C_NEAR(s, c1*c2*se, eps);    
  

}
TEST(CartGTO, NuclearAttraction) {
  double eps(0.000000001);
  double c1(0.33);
  double c2(0.11);
  double q1(1.1);
  CartGTO<C, c3> g1(c1,
		   i3(2, 1, 0),
		   c3(0.1,  C(0.2, 0.03), 0.3), 
		   C(1.1,0.2));
  CartGTO<C, c3> g2(c2, 
		   i3(0, 1, 3),
		   c3(0.1, -0.22, C(-0.3,-0.7)),
		   C(1.3,0.02));
  OpNA<C, c3> op_na(q1, c3(-0.1, C(-0.1,0.1), 0.3));
    
  C s = CIP(g1, op_na, g2);
  C se(-0.0111712963403, -0.0039848461450);
  
  EXPECT_C_NEAR(s, q1*c1*c2*se, eps);    
}
TEST(CartGTO, Dip) {

  double c1(1.1);
  double c2(1.2);
  CartGTO<C, c3> g1(c1, i3(2, 1, 1), c3(0.0, 0.0, 0.0),        C(0.2, 0.3));
  CartGTO<C, c3> g2(c2, i3(2, 1, 2), c3(0.0, C(0.2,-0.1), 0.3), C(0.2, -0.7));

  OpXyz<C, c3> op(1, 2, 3);

  CartGTO<C, c3> g3(c1, i3(3, 2, 2), c3(0.0, 0.0, 0.0),        C(0.2, 0.3));
  CartGTO<C, c3> g4(c2, i3(2, 2, 4), c3(0.0, C(0.2,-0.1), 0.3), C(0.2, -0.7));
  
  EXPECT_C_EQ(CIP(g1, op, g2), CIP(g3, g4));

}

TEST(Angmoment, cg_coef) {
  int j1 = 2;
  int m1 = 1;
  EXPECT_DOUBLE_EQ(-1.0/sqrt(2*j1+1),
		   cg_coef(j1,j1,m1,-m1, 0,0));
  EXPECT_NEAR(pow(fact(2*j1),2)/(6 * sqrt(fact(8)*1.0)),
	      cg_coef(j1,j1,m1,-m1,2*j1,0),
	      0.0000000001);
  std::cout << cg_coef(1,1,1,0,0,0) << std::endl;
}
TEST(Angmoment, lm_pair) {

  EXPECT_EQ(0, lm_index(0, 0));

  EXPECT_EQ(1, lm_index(1, -1));
  EXPECT_EQ(2, lm_index(1,  0));
  EXPECT_EQ(3, lm_index(1, +1));

  EXPECT_EQ(4, lm_index(2, -2));
  EXPECT_EQ(5, lm_index(2, -1));
  EXPECT_EQ(6, lm_index(2,  0));
  EXPECT_EQ(7, lm_index(2, +1));
  EXPECT_EQ(8, lm_index(2, +2));

  EXPECT_EQ(9,  lm_index(3, -3));
  EXPECT_EQ(10, lm_index(3, -2));
  EXPECT_EQ(11, lm_index(3, -1));
  EXPECT_EQ(12, lm_index(3,  0));
  EXPECT_EQ(13, lm_index(3, +1));
  EXPECT_EQ(14, lm_index(3, +2));
  EXPECT_EQ(15, lm_index(3, +3));
}
TEST(Angmoment, num_lm_pair) {

  EXPECT_EQ(1, num_lm_pair(0));
  EXPECT_EQ(1+3, num_lm_pair(1));
  EXPECT_EQ(1+3+5, num_lm_pair(2));
  EXPECT_EQ(1+3+5+7, num_lm_pair(3));

}
TEST(Angmoment, mod_spherical_bessel) {

  C xs[3];
  xs[0] = C(1.1, 0.3);
  xs[1] = C(0.1, 0.002);
  xs[2] = C(100.0, 20.0);
  
  for(int i = 0; i < 3; i++) {
    C x(xs[i]);
    C* vs = ModSphericalBessel(x, 10);
  
    EXPECT_C_EQ(sinh(x)/x,                 vs[0]);
    EXPECT_C_EQ((x*cosh(x)-sinh(x))/(x*x), vs[1]);
    EXPECT_C_EQ(((x*x+3.0)*sinh(x)-3.0*x*cosh(x))/(x*x*x), vs[2]);
    EXPECT_C_EQ(((x*x*x+15.0*x)*cosh(x)-(6.0*x*x+15.0)*sinh(x))/pow(x, 4.0), vs[3]);
    
    delete vs;
  }

  xs[0] = C(0.000001, 0.000002);
  xs[1] = C(0.00005, 0.002);
  xs[2] = C(0.0001, 0.002);
  for(int i = 0; i < 3; i++) {
    C* vs = ModSphericalBessel(xs[i], 10);
    for(int n = 2; n <= 8; n++) {
      EXPECT_C_EQ(vs[n-1]-vs[n+1], (2*n+1.0)*vs[n]/xs[i]);    
    }
    delete vs;
  }
    
  C x(0.0, 0.0);
  C* vs = ModSphericalBessel(x, 10);
  EXPECT_C_EQ(1.0, vs[0]);
  for(int i = 1; i < 10; i++)
    EXPECT_C_EQ(0.0, vs[i]);
  delete vs;

}
TEST(Angmoment, associated_legendre) {

  C x( 0.6, 0.0);
  C* vs = AssociatedLegendre(x, 10);  
  EXPECT_C_EQ(1.0, vs[lm(0, +0)]);

  EXPECT_C_EQ(-sqrt(1.0-x*x), vs[lm(1, +1)]);
  EXPECT_C_EQ(x,              vs[lm(1, +0)]);  

  EXPECT_C_EQ(0.5*(3.0*x*x-1.0),  vs[lm(2, 0)]);
  EXPECT_C_EQ(-3.0*x*sqrt(1.0-x*x), vs[lm(2, 1)]);
  EXPECT_C_EQ(3.0*(1.0-x*x),      vs[lm(2,+2)]);

  for(int L = 0; L < 4; L++) 
    for(int M = 1; M <= L; M++) 
      EXPECT_C_EQ(vs[lm(L, -M)], pow(-1.0, M)*fact(L-M)*1.0/fact(L+M)*vs[lm(L, M)]);

}
TEST(Angmoment, real_spherical_harm) {

  C t = C(1.1, 0.2);
  C p = C(1.1, 0.1);

  C* vs = RealSphericalHarmonics(t, p, 10);
  
  C a = 1.0/sqrt(4.0*M_PI);
  EXPECT_C_EQ(a,                         vs[lm(0, 0)]);
  EXPECT_C_EQ(a*sqrt(3.0)*cos(t),        vs[lm(1, 0)]);
  EXPECT_C_EQ(a*sqrt(3.0)*sin(t)*cos(p), vs[lm(1,+1)]);
  EXPECT_C_EQ(a*sqrt(3.0)*sin(t)*sin(p), vs[lm(1,-1)]);
  
  EXPECT_C_EQ(a*sqrt(5.0/4.0)*(3.0*cos(t)*cos(t)-1.0),   vs[lm(2, 0)]);
  EXPECT_C_EQ(a*sqrt(15.0)*sin(t)*cos(t)*cos(p),         vs[lm(2,+1)]);
  EXPECT_C_EQ(a*sqrt(15.0)*sin(t)*cos(t)*sin(p),         vs[lm(2,-1)]);
  EXPECT_C_EQ(a*sqrt(15.0/4.0)*sin(t)*sin(t)*cos(2.0*p), vs[lm(2,+2)]);
  EXPECT_C_EQ(a*sqrt(15.0/4.0)*sin(t)*sin(t)*sin(2.0*p), vs[lm(2,-2)]);
}

TEST(Angmoment, gto_00_r) {

  double pi(M_PI);  
  int num_r(5);
  C* rs = new C[num_r];
  for(int i = 0; i < num_r; i++) {
    rs[i] = 0.1 + 0.2*i;
  }
  C zeta(1.0, 0.0);
  int L(0);
  int M(0);
  C* vs = gto_00_r(0.0, 0.0, 0.0, L, M, rs, num_r, zeta);
  for(int i = 0; i < num_r; i++) {
    C r(rs[i]);
    C calc(sqrt(4.0*pi)*exp(-zeta*r*r));    
    EXPECT_C_EQ(vs[i], calc) << i;    
  }
}

TEST(SphericalGTO, OrthNormality) {
  double eps(0.000000001);
  c3 xyz(0.1, 0.2, 0.3);
  C zeta(1.1, -0.3);
  
  int num = 9;
  SphericalGTO<C, C>** gs = new SphericalGTO<C, C>*[num];
  
  // LinFunc<CartGTO<C, c3> >* gs;
  //gs = new LinFunc<CartGTO<C, c3> >[num];
  gs[0] = new SphericalGTO<C, C>(0, 0, xyz, zeta);
  gs[1] = new SphericalGTO<C, C>(1, 0, xyz, zeta);
  gs[2] = new SphericalGTO<C, C>(1, 1, xyz, zeta);
  gs[3] = new SphericalGTO<C, C>(1,-1, xyz, zeta);
  gs[4] = new SphericalGTO<C, C>(2, 0, xyz, zeta);
  gs[5] = new SphericalGTO<C, C>(2, 1, xyz, zeta);
  gs[6] = new SphericalGTO<C, C>(2,-1, xyz, zeta);
  gs[7] = new SphericalGTO<C, C>(2, 2, xyz, zeta);
  gs[8] = new SphericalGTO<C, C>(2,-2, xyz, zeta);

  
  for(int i = 0; i < num; i++) {
    EXPECT_C_NEAR(C(1), CIP(*gs[i], *gs[i]), eps) << i;
    for(int j = 0; j < i; j++) {
      EXPECT_C_NEAR(C(0), CIP(*gs[i], *gs[j]), eps) << i << j;
      EXPECT_C_NEAR(C(0), CIP(*gs[j], *gs[i]), eps) << i << j;
    }
  }
}
TEST(SphericalGTO, at) {

  c3 xyz(0.1, 0.2, 0.3);
  C zeta(1.1, -0.3);

  SphericalGTO<C, C> g(1, 0, xyz, zeta);
  EXPECT_C_EQ(0.0, g.at(c3(0.33, 0.44, 0.3)));

  SphericalGTO<C, C> g11(1, 1, xyz, zeta);
  EXPECT_C_EQ(0.0, g11.at(c3(0.1, 0.44, 0.44)));

  SphericalGTO<C, C> g1m1(1, -1, xyz, zeta);
  EXPECT_C_EQ(0.0, g1m1.at(c3(0.11, 0.2, 0.44)));
  
}
TEST(SphericalGTO, Exception) {
  
  LinFunc<CartGTO<C, c3> > func;
  try {
    SetSphericalGTO<C, C>(1, 2, c3(1.1, 1.2, 1.3), 1.2, &func);
  } catch(const std::exception& e) {
    std::cerr << e.what() << std::endl;    
  }

}

TEST(GTOSet, Create) {

  SphericalGTOSet gtos;
  for(int L = 0; L <= 2; L++)
    gtos.AddBasis(L,
		  0.1, 0.2, 0.3,
		  1.1);
  
  std::complex<double>* vs;
  vs = gtos.SMat(gtos);
  EXPECT_C_EQ(1, vs[0]);
  EXPECT_C_EQ(0, vs[1]);
  EXPECT_C_EQ(0, vs[2]);
}


int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
