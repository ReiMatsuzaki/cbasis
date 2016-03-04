/*
  Unit test for cints.cpp, 
*/

#include <gtest/gtest.h>
#include "utils.hpp"
#include "cints.hpp"
#include "angmoment.hpp"
#include "gto3d.hpp"
#include "gto3dset.hpp"
#include "molint.hpp"


using namespace std;
using namespace l2func;

int lm(int l, int m) {
  return lm_index(l, m);
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
TEST(GTOSet, Time) {

  #define USE_MOLINT
  //#undef USE_MOLINT
  // molint => 3318 ms (in macbook air)
  // cints  => 3260 ms (in macbook air)     
  std::cout << 0 << std::endl;
  SphericalGTOSet gtos;
  for(int L = 0; L <= 2; L++)
    for(int n = -5; n < 5; n++)
      gtos.AddBasis(L, 0.1, 0.2, 0.3, pow(2.0, n));
  std::cout << 0.5 << std::endl;
  dcomplex* smat = gtos.SMat(gtos);
  std::cout << 0.6 << std::endl;
  dcomplex* vmat = gtos.VMat(1.0, 0.0, 0.0, 0.0, gtos);
  std::cout << 0.7 << std::endl;
  dcomplex* tmat = gtos.TMat(gtos);
  std::cout << 1 << std::endl;
  delete smat;
  std::cout << 2 << std::endl;
  delete vmat;
  std::cout << 3 << std::endl;
  delete tmat;

}

TEST(GTOs, offset) {

  GTOs gtos;
  gtos.AddSphericalGTO(1, 0.0, 0.0, 0.0, 1.2);
  gtos.AddSphericalGTO(2, 0.0, 0.0, 0.0, 1.2);
  gtos.AddSphericalGTO(0, 0.0, 0.0, 0.0, 1.2);
  EXPECT_EQ(0, gtos.offset_ish[0]);
  EXPECT_EQ(3, gtos.offset_ish[1]);
  EXPECT_EQ(8, gtos.offset_ish[2]);

}
TEST(GTOs, d_coef) {

/*
  dcomplex ref[4*5*6];
  
  for(int ni = 0; ni <= 2; ni++)
    for(int nj = 0; nj <= 3; nj++)
      for(int n = 0; n <= 4; n++) {
	ref[ni + nj*3 + n*4*3] = coef_d(zetaP, 0.5, 0.1, 0.3, ni, nj, n);
      }
*/
  dcomplex zetaP(1.1, 0.3);
  dcomplex wPk(0.666);
  dcomplex xi(0.111);
  dcomplex xj(0.125);
  dcomplex* buff = new dcomplex[4*5*6];
  for(int i = 0; i < 4*5*6; i++)
    buff[i] = 1234.00;
  
  MultArray3<dcomplex> ary = calc_d_coef(2, 3, 4, zetaP, wPk, xi, xj, buff);

  int idx(0);
  for(int ni = 0; ni <= 2; ni++)
    for(int nj = 0; nj <= 3; nj++)
      for(int n = 0; n <=4; n++) {
	dcomplex calc = ary.get(ni, nj, n);
	dcomplex ref0 = coef_d(zetaP, wPk, xi, xj, ni, nj, n);
	EXPECT_C_EQ(ref0, calc)<< ni << nj << n;	  
	idx++;
      }
  

}
TEST(GTOs, size) {

  GTOs gtos;
  
  for(int L = 2; L < 3; L++)
    for(int n = 0; n < 1; n++ ) {

      dcomplex zeta(0.3);
      dcomplex x(0.0);
      dcomplex y(0.0);
      dcomplex z(0.0);

      gtos.AddSphericalGTO(L, x, y, z, zeta);
  }

  EXPECT_EQ(6, gtos.size_prim());
  EXPECT_EQ(5, gtos.size_basis());

}
TEST(GTOs, Create) {

  //  SphericalGTOSet gto_ref;
  typedef SphericalGTO<dcomplex, dcomplex> GY;
  GTOs gtos;
  vector<GY*> gs;
  for(int L = 0; L < 3; L++)
    for(int n = 0; n < 1; n++ ) {
      
      dcomplex zeta(0.1+n*0.1, 0.0);
      dcomplex x(0.4*(n+1));
      dcomplex y(-0.1+0.1*(n+1));
      dcomplex z(0.3+0.2*(n+1));

      gtos.AddSphericalGTO(L, x, y, z, zeta);
      //      for(int M = -L; M <= L; M++)
      //	gs.push_back(new GY(L, M, c3(x, y, z), zeta));
      //      gto_ref.AddBasis(    L, x, y, z, zeta);
  }
  
  //  dcomplex* smat_ref  = gto_ref.SMat(gto_ref);
  //  dcomplex* tmat_ref  = gto_ref.TMat(gto_ref);
  //  dcomplex* zmat_ref  = gto_ref.XyzMat(0, 0, 1, gto_ref);
  //  dcomplex* vmat_ref = gto_ref.VMat(1.0, 0.0, 0.0, 0.0, gto_ref);
  //  dcomplex* vmat_ref  = gto_ref.VMat(1.1, -0.1, 0.22, 0.31, gto_ref);
  dcomplex* smat_calc;
  dcomplex* tmat_calc;
  dcomplex* zmat_calc;
  dcomplex* vmat_calc;
  //gtos.Show();
  int num = gtos.size_basis();
  gtos.Normalize();
  gtos.CalcMat(&smat_calc, &tmat_calc, &zmat_calc, &vmat_calc);    

  double eps(pow(10.0, -10.0));

  //  EXPECT_EQ(gtos.size_basis(), gs.size());
  for(int i = 0; i < num; i++)
    for(int j = 0; j < num; j++) {
      int idx = i * num + j;
      std::cout << smat_calc[idx];
      //      std::cout << i << j << std::endl;
      //      EXPECT_C_NEAR(CIP(*gs[i], *gs[j]), smat_calc[idx], eps) 
      //	<< "(" << i << ", " << j << ")" << std::endl;
      /*
      EXPECT_C_NEAR(tmat_ref[idx], tmat_calc[idx], eps) 
	<< "(" << i << ", " << j << ")" << std::endl;
      EXPECT_C_NEAR(zmat_ref[idx], zmat_calc[idx], eps) 
	<< "(" << i << ", " << j << ")" << std::endl;
      EXPECT_C_NEAR(vmat_ref[idx], vmat_calc[idx], eps)
	<< "(" << i << ", " << j << ")" << std::endl;
	*/
    }
  std::cout << std::endl;
  delete[] smat_calc;
  delete[] tmat_calc;
  delete[] zmat_calc;
  delete[] vmat_calc;
}

TEST(MOLINT, overlap) {

  dcomplex old = overlap(C(1.2,0.2), 2, 1, 0,
			 1.0,  0.0,   0.0,
			 C(0.7,-0.5), 1, 1, 2,
			 0.0, 0.3, 0.0);
  dcomplex calc = gto_overlap(2, 1, 0,
			   1.0, 0.0, 0.0,
			   dcomplex(1.2, 0.2),
			   1, 1, 2,
			   0.0, 0.3, 0.0,
			      dcomplex(0.7, -0.5));
  EXPECT_C_EQ(old, calc);

  old = overlap(C(1.2,0.2), 2, 1, 0,
		0.0,  0.0,   0.0,
		C(0.7,-0.5), 2, 1, 0,
		0.0, 0.0, 0.0);
  calc = gto_overlap(2, 1, 0,
		     0.0, 0.0, 0.0,
		     dcomplex(1.2, 0.2),
		     2, 1, 0,
		     0.0, 0.0, 0.0,
		     dcomplex(0.7, -0.5));
  EXPECT_TRUE(abs(calc)>0.0001);  
  EXPECT_C_EQ(old, calc);

}
TEST(MOLINT, kinetic) {

  dcomplex old = kinetic(C(1.2,0.2), 0, 0, 0,
			 0.0,  0.0,   0.0,
			 C(0.7,-0.5), 0, 0, 0,
			 0.0, 0.0, 0.0);
  dcomplex calc = gto_kinetic(0, 0, 0,
			      0.0, 0.0, 0.0,
			      dcomplex(1.2, 0.2),
			      0, 0, 0,
			      0.0, 0.0, 0.0,
			      dcomplex(0.7, -0.5));
  EXPECT_TRUE(abs(calc)>0.0001);  
  EXPECT_C_EQ(old, calc);

  old = kinetic(C(1.2,0.2), 0, 1, 0,
		0.0,  0.0,   0.0,
		C(0.7,-0.5), 0, 1, 0,
		0.0, 0.0, 0.0);
  calc = gto_kinetic(0, 1, 0,
		     0.0, 0.0, 0.0,
		     dcomplex(1.2, 0.2),
		     0, 1, 0,
		     0.0, 0.0, 0.0,
		     dcomplex(0.7, -0.5));
  EXPECT_TRUE(abs(calc)>0.0001);  
  EXPECT_C_EQ(old, calc);

  old = kinetic(dcomplex(1.2,0.2), 2, 0, 3,
		0.2,  0.3,   0.4,
		dcomplex(0.7,-0.5), 1, 1, 2,
		0.0, 0.0, 0.0);
  calc = gto_kinetic(2, 0, 3,
		     0.2, 0.3, 0.4,
		     dcomplex(1.2, 0.2),
		     1, 1, 2,
		     0.0, 0.0, 0.0,
		     dcomplex(0.7, -0.5));
  EXPECT_TRUE(abs(calc)>0.0001);  
  EXPECT_C_EQ(old, calc);

}
TEST(MOLINT, IncGamma) {

  /*
    python code :
    from scipy.special import hyp1f1
    def inc_gamma(m, z):
        return 1.0/(2*m+1) * hyp1f1(m+0.5, m+1.5, -z)
    print inc_gamma(0, 50.0-5.50j)
    print inc_gamma(1, 50.0-5.50j)
    print inc_gamma(2, 50.0-5.50j)
   */

  dcomplex* us = IncompleteGamma(10, dcomplex(0.1, -0.2));
  EXPECT_C_EQ(dcomplex(0.96392503253589557,
		       +0.06262986651323893),		       
	      us[0]);
  EXPECT_C_EQ(dcomplex(0.100986690286, 0.0166268474338),
	      us[4]);
  delete us;

  dcomplex* vs1 = IncompleteGamma(10, dcomplex(21.0, 15.5));
  EXPECT_C_EQ(dcomplex(-3.59002512901*pow(10.0, -6),
		       -0.000190939738972),
	      vs1[2]);
  
  dcomplex* vs = IncompleteGamma(10, dcomplex(50.0, -5.5));

  double step(0.00001);
  EXPECT_C_EQ(dcomplex(0.124767690392,
		       0.00684158939256),
	      vs[0]);
  EXPECT_C_EQ(dcomplex(0.0012253247264,
		       0.00020320161383),
	      vs[1]);
  EXPECT_C_EQ(dcomplex(3.5657718077638774*step,
		       1.0018397403429266*step),
	      vs[2]);
  delete vs;

  dcomplex* v0s = IncompleteGamma(10, 0.0);
  EXPECT_C_EQ(1.0, v0s[0]);
  EXPECT_C_EQ(1.0/3.0, v0s[1]);
  EXPECT_C_EQ(1.0/5.0, v0s[2]);


  
}

TEST(MOLINT, nuclear_attraction) {

  dcomplex zetaA(1.2, -0.3);
  dcomplex zetaB(0.6, 0.0);  
  dcomplex old = nuclear_attraction(0.0,  0.0,  0.0,
				    1.0,
				    0, 0, 0, zetaA,
				    0.0, 0.0, 0.0,
				    1.0,
				    0, 0, 0, zetaB,
				    0.0, 0.0, 0.0);
  dcomplex calc = gto_nuclear_attraction(0, 0, 0, 
					 0.0, 0.0, 0.0,
					 zetaA,
					 0, 0, 0,
					 0.0, 0.0, 0.0,
					 zetaB,
					 0.0, 0.0, 0.0);  
  EXPECT_C_EQ(old, calc);

  old = nuclear_attraction(0.0,  0.0,  0.0,
			   1.0,
			   1, 0, 0, zetaA,
			   0.0, 0.0, 0.0,
			   1.0,
			   0, 0, 0, zetaB,
			   0.0, 0.0, 0.0);

  calc = gto_nuclear_attraction(1, 0, 0, 
				0.0, 0.0, 0.0,
				zetaA,
				0, 0, 0,
				0.0, 0.0, 0.0,
				zetaB,
				0.0, 0.0, 0.0);
  EXPECT_C_EQ(old, calc);  

  old = nuclear_attraction(0.3,  0.0,  0.0,
			   1.0,
			   0, 0, 0, zetaA,
			   0.0, 0.0, 0.0,
			   1.0,
			   0, 0, 0, zetaB,
			   0.0, 0.0, 0.0);

  calc = gto_nuclear_attraction(0, 0, 0, 
				0.3, 0.0, 0.0,
				zetaA,
				0, 0, 0,
				0.0, 0.0, 0.0,
				zetaB,
				0.0, 0.0, 0.0);
  EXPECT_C_EQ(old, calc);  

  old = nuclear_attraction(0.3,  0.0,  0.0,
			   1.0,
			   1, 0, 0, zetaA,
			   0.0, 0.0, 0.0,
			   1.0,
			   0, 0, 0, zetaB,
			   0.0, 0.0, 0.0);

  calc = gto_nuclear_attraction(1, 0, 0, 
				0.3, 0.0, 0.0,
				zetaA,
				0, 0, 0,
				0.0, 0.0, 0.0,
				zetaB,
				0.0, 0.0, 0.0);
  EXPECT_C_EQ(old, calc);  

  old = nuclear_attraction(0.3,  0.0,  0.0,
			   1.0,
			   1, 2, 0, zetaA,
			   0.0, dcomplex(0.1, 0.01), 0.0,
			   1.0,
			   0, 2, 0, zetaB,
			   0.3, 0.0, 0.0);

  calc = gto_nuclear_attraction(1, 2, 0, 
				0.3, 0.0, 0.0,
				zetaA,
				0, 2, 0,
				0.0, dcomplex(0.1, 0.01), 0.0,
				zetaB,
				0.3, 0.0, 0.0);
  EXPECT_C_EQ(old, calc);  

  old = nuclear_attraction(0.0,  0.0,  0.0,
			   1.0,
			   0, 0, 1, zetaA,
			   0.0, 0.0, 0.0,
			   1.0,
			   0, 0, 1, zetaB,
			   0.0, 0.0, 0.0);

  calc = gto_nuclear_attraction(0, 0, 1, 
				0.0, 0.0, 0.0,
				zetaA,
				0, 0, 1,
				0.0, 0.0, 0.0,
				zetaB,
				0.0, 0.0, 0.0);
  EXPECT_C_EQ(old, calc);  

}
TEST(MOLINT, Time) {

  

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
  delete vs;
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

  c3 xyz(0.1, 0.2, 0.3); C zeta(1.1, -0.3);

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

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
