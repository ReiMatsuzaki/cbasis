#include <gtest/gtest.h>
#include "gtest_plus.hpp"

// header only
#include "typedef.hpp"
#include "mult_array.hpp"

// mathematical functions
#include "fact.hpp"
#include "erfc.hpp"
#include "lgamma.hpp"
#include "angmoment.hpp"
#include "mol_func.hpp"
#include "eigen_plus.hpp"
#include "bmatset.hpp"

using namespace std;
using namespace l2func;
using namespace Eigen;

TEST(EigenPlus, Canonical) {

  MatrixXcd A2(2,2);
  A2 <<
    1.0, 0.5,
    0.5, 1.5;

  double eps(0.00001);
  MatrixXcd A3(3,3);
  A3 <<
    1.0,     0.5,     0.5+eps,
    0.5,     1.5,     1.5+eps,
    0.5+eps, 1.5+eps, 1.5+eps;
  MatrixXcd F2(2,2); 
  F2 <<
    1.4, 0.2,
    0.2, 0.4;
  MatrixXcd F3(3,3);
  F3 <<
    1.4,     0.2,     0.2+eps,
    0.2,     0.4,     0.4+eps,
    0.2+eps, 0.4+eps, 0.4+eps;
  
  CM C2, C3, C2_c;
  CV eig2, eig3, eig2_c;
  generalizedComplexEigenSolve(F2, A2, &C2, &eig2);
  generalizedComplexEigenSolve(F3, A3, &C3, &eig3);
  //CEigenSolveCarnonicalNum(F3, A3, 2, &C2_c, &eig2_c);
  CEigenSolveCanonical(F3, A3, 0.00001, &C2_c, &eig2_c);
  
  cout << eig2 << endl;
  cout << endl;
  cout << eig3 << endl;
  cout << endl;
  cout << eig2_c << endl;

  EXPECT_C_NEAR(eig2[0], eig2_c[0], eps);
  EXPECT_C_NEAR(eig2[1], eig2_c[1], eps);

}

TEST(BMat, ReadWrite) {

  BMat bmat1;

  MatrixXcd M10(2, 3); 
  M10 <<
    1.0, 1.1, 1.2,
    2.1, 2.2, 2.3;
  bmat1[make_pair(1, 0)] = M10;

  MatrixXcd M12(3, 1); 
  M12 << 3.1, 3.2, 3.3;
  bmat1[make_pair(1, 2)] = M12;

  string filename("tmp.dat");
  BMatWrite(bmat1, filename);

  BMat bmat2;
  BMatRead(bmat2, filename);

  EXPECT_EQ(bmat1.size(), bmat2.size());
  EXPECT_MATXCD_EQ(bmat1[make_pair(1, 0)],
		   bmat2[make_pair(1, 0)]);
  EXPECT_MATXCD_EQ(bmat1[make_pair(1, 2)],
		   bmat2[make_pair(1, 2)]);
  
}

TEST(BMatSet, SetGet) {

  BMatSet sets(new _BMatSet());

  MatrixXcd s00(2,2); s00 << 1.0, 1.1, 1.2, 1.3;
  sets->SetMatrix("s", 0, 0, s00);

  MatrixXcd s01(2,3); s01 << 2.0, 2.1, 2.2, 2.3, 0.1, 0.1;
  sets->SetMatrix("s", 0, 1, s01);

  MatrixXcd t01(2,3); t01 << 2.6, 2.1, 2.5, 0.3, 0.11, 0.12;
  sets->SetMatrix("t", 0, 1, t01);

  const MatrixXcd s01ref = sets->GetMatrix("s", 0, 1);

  EXPECT_C_EQ(2.0, s01ref(0, 0));
  EXPECT_C_EQ(2.1, s01ref(0, 1));
  EXPECT_C_EQ(2.2, s01ref(0, 2));
  EXPECT_C_EQ(2.3, s01ref(1, 0));
  EXPECT_C_EQ(0.1, s01ref(1, 1));
  EXPECT_C_EQ(0.1, s01ref(1, 2));

}
TEST(BMatSet, Plus) {

  BMatSet sets(new _BMatSet);
  MatrixXcd s00(2,2); s00 << 1.0, 1.1, 1.2, 1.3;
  sets->SetMatrix("s", 0, 0, s00);

  sets->SelfAdd("s", 0, 0, 0, 1, 0.33);
  EXPECT_C_EQ(1.1+0.33, sets->GetValue("s", 0, 0, 0, 1));

}
TEST(BMatSet, Exception) {

  BMatSet sets(new _BMatSet);
  
  MatrixXcd s00 = MatrixXcd::Zero(2, 2);
  
  sets->SetMatrix("s", 0, 0, s00);
  EXPECT_NO_THROW(sets->SetMatrix("s", 2, 1, s00));
  EXPECT_NO_THROW(sets->SetMatrix("s", -1, 1, s00));
  
#ifdef ARG_NO_CHECK
  EXPECT_NO_THROW(sets->GetMatrix("t", 0, 0));
  EXPECT_NO_THROW(sets->GetMatrix("s", -1, 0));
  EXPECT_NO_THROW(sets->GetMatrix("s", 1, 2));
#else
  EXPECT_ANY_THROW(sets->GetMatrix("t", 0, 1));
  EXPECT_ANY_THROW(sets->GetMatrix("s", -1, 0));
  EXPECT_ANY_THROW(sets->GetMatrix("s", 1, 2));
#endif  

  try {
    sets->GetMatrix("s", -1, 0);
  } catch(const runtime_error& e) {
    cout << e.what() << endl;
  }

}
TEST(BMatSet, BMat_swap) {

  BMat s0;
  s0[make_pair(0, 0)] = MatrixXcd::Ones(4, 5);
  s0[make_pair(1, 1)] = MatrixXcd::Zero(4, 5);

  BMat s1;
  swap(s0, s1);
  
  EXPECT_MATXCD_EQ(MatrixXcd::Ones(4, 5), s1[make_pair(0, 0)]);
  EXPECT_MATXCD_EQ(MatrixXcd::Zero(4, 5), s1[make_pair(1, 1)]);

}
TEST(BmatSet, Swap) {

  BMatSet sets(new _BMatSet);
  MatrixXcd s00(2, 2); 
  s00 << 1.0, 1.1, 1.2, 1.3;
  MatrixXcd s01(2, 2); 
  s01 << 2.0, 2.1, 2.2, 2.3;

  sets->SetMatrix("s", 0, 0, s00);
  sets->SetMatrix("s", 0, 1, s01);
  
  BMatSet set2(new _BMatSet);
  swap(set2, sets);

  EXPECT_C_EQ(set2->GetValue("s", 0, 0, 0, 0), 1.0);
  EXPECT_C_EQ(set2->GetValue("s", 0, 0, 0, 1), 1.1);
  EXPECT_C_EQ(set2->GetValue("s", 0, 0, 1, 0), 1.2);
  EXPECT_C_EQ(set2->GetValue("s", 0, 0, 1, 1), 1.3);

  EXPECT_C_EQ(set2->GetValue("s", 0, 1, 0, 0), 2.0);
  EXPECT_C_EQ(set2->GetValue("s", 0, 1, 0, 1), 2.1);
  EXPECT_C_EQ(set2->GetValue("s", 0, 1, 1, 0), 2.2);
  EXPECT_C_EQ(set2->GetValue("s", 0, 1, 1, 1), 2.3);
  
}

TEST(MultArray, MultArray1) {

  MultArray<int, 1> xs(10);
  xs.SetRange(-1, 3);
  for(int i = -1; i <= 3; i++)
    xs(i) = i*2;

  for(int i = -1; i <= 3; i++)
    EXPECT_EQ(i*2, xs(i));

}
TEST(MultArray, MultArray2) {

  MultArray<int, 2> xs(100);
  xs.SetRange(-2, 3, -1, 4);

  for(int i = -2; i <= 3; i++)
    for(int j = -1; j <= 4; j++)
      xs(i, j) = 10*i+j;

  for(int i = -2; i <= 3; i++)
    for(int j = -1; j <= 4; j++)
      EXPECT_EQ(10*i+j, xs(i, j));

}
TEST(MultArray, MultArray2Exception) {

  MultArray<int, 2> xs(100);
  xs.SetRange(-2, 3, -1, 4);

#ifdef ARG_NO_CHECK  
#else
  EXPECT_ANY_THROW(xs(-3, 1));
  EXPECT_ANY_THROW(xs(4, 4));
  EXPECT_ANY_THROW(xs(0, 5));  
  EXPECT_ANY_THROW(xs(-3, 4));
#endif

}
TEST(MultArray, MultArray3) {

  MultArray<int, 3> xs(216);
  xs.SetRange(-2, 3,
	      -1, 4,
	      0, 5);
  EXPECT_TRUE(xs.size() == 216) << xs.size();
  for(int i = -2; i <= 3; i++)
    for(int j = -1; j <= 4; j++)
      for(int k = 0; k <= 5; k++) {
	xs(i, j, k) = 100*i+10*j+k;
      }

  for(int i = -2; i <= 3; i++)
    for(int j = -1; j <= 4; j++)
      for(int k = 0; k <= 5; k++) {
	EXPECT_EQ(100*i+10*j+k, xs(i, j, k));
      }

}
TEST(MultArray, MultArray4) {

  MultArray<int, 4> xs(1000);
  xs.SetRange(-2, 3,
	      -1, 4,
	      0, 5,
	      2, 4);
  for(int i = -2; i <= 3; i++)
    for(int j = -1; j <= 4; j++)
      for(int k = 0; k <= 5; k++)
	for(int l = 2; l <= 4; l++)  {
	  try{ 
	    xs(i, j, k, l) = 1000*i+100*j+10*k+l;
	  } catch(const exception& e) {
	    EXPECT_TRUE(false) << e.what();
	  }
	}


  for(int i = -2; i <= 3; i++)
    for(int j = -1; j <= 4; j++)
      for(int k = 0; k <= 5; k++)
	for(int l = 2; l <= 4; l++) {
	  try{ 
	    EXPECT_EQ(1000*i+100*j+10*k+l, xs(i, j, k, l));
	  } catch(const exception& e) {
	    EXPECT_TRUE(false) << e.what();
	  }
	}

}

TEST(Fact, Factorial) {

  EXPECT_ANY_THROW(Factorial(-1));
  EXPECT_EQ(1, Factorial(0));
  EXPECT_EQ(1, Factorial(1));
  EXPECT_EQ(2, Factorial(2));
  EXPECT_EQ(6, Factorial(3));
  EXPECT_EQ(24, Factorial(4));
  EXPECT_EQ(120, Factorial(5));
  EXPECT_EQ(720, Factorial(6));

  EXPECT_ANY_THROW(DoubleFactorial(-1));
  EXPECT_EQ(1,   DoubleFactorial(0));
  EXPECT_EQ(1,   DoubleFactorial(1));
  EXPECT_EQ(2,   DoubleFactorial(2));
  EXPECT_EQ(3,   DoubleFactorial(3));
  EXPECT_EQ(8,   DoubleFactorial(4));
  EXPECT_EQ(15,  DoubleFactorial(5));
  EXPECT_EQ(48,  DoubleFactorial(6));
  EXPECT_EQ(105, DoubleFactorial(7));
   
}
TEST(Fact, Combination) {
  EXPECT_EQ(10, Combination(5, 2));
  EXPECT_EQ(1, Combination(4, 0));
}

TEST(Erfc, real_erfc) {

  using namespace erfc_mori;

  // this function is forbidden
  //  erfc_add_Eh_q<int>(1, 2);
  
  double y, x, expect;
  ErfcCalcData calc_data;

  x = 1.0;
  expect =0.157299207050285130658779364917390740703933002034;
  //erfc_d(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);
  
  x = 1.0;
  x = 1.0 / 100.0;
  expect = 0.98871658444415038308409047645193078905089904517;
  //  erfc_d(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);

  x = 3.0;
  expect = 0.0000220904969985854413727761295823203798477070873992;
  //  erfc_d(x, y, calc_data);
  // erfc(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);

}
TEST(Erfc, complex_erfc) {

  using namespace erfc_mori;

  dcomplex y;
  ErfcCalcData calc_data;
  double eps = 10.0 * machine_eps();

  dcomplex x(1, -1);
  dcomplex y_expect( -0.31615128169794764488027108024367,
	       +0.190453469237834686284108861969162);

  Erfc(x, y, calc_data);

  EXPECT_DOUBLE_EQ(y.real(), y_expect.real());
  EXPECT_NEAR(y.real(), y_expect.real(), eps);
  EXPECT_DOUBLE_EQ(y.imag(), y_expect.imag());
  EXPECT_NEAR(y.imag(), y_expect.imag(), eps);

  x = dcomplex(0.0157073173118206757532953533099,
	 0.9998766324816605986389071277312);
  y_expect = dcomplex(0.95184545524179913420177473658805,
		-1.64929108965086517748934245403332);

  Erfc(x, y, calc_data);
  EXPECT_TRUE(calc_data.convergence);
  EXPECT_NEAR(     y.real(), y_expect.real(), eps);
  EXPECT_NEAR(     y.imag(), y_expect.imag(), eps);

  x = dcomplex(0.001564344650402308690101053194671668923139,
	 0.009876883405951377261900402476934372607584);
  y_expect = dcomplex(0.9982346553205423153337357292658472915601,
		-0.0111452046101524188315708507537751407281);

  Erfc(x, y, calc_data);
  EXPECT_TRUE(calc_data.convergence);
  EXPECT_NEAR(y.real(), y_expect.real(), eps);
  EXPECT_NEAR(y.imag(), y_expect.imag(), eps);

}

TEST(lgamma, lower_gamma) {
  double eps = pow(10.0, -13.0);
  dcomplex a(1.1, -0.3);

  EXPECT_ANY_THROW(LowerGamma(-1, a));

  EXPECT_C_NEAR(1.0-exp(-a), LowerGamma(1, a), eps);
  EXPECT_C_NEAR(-a*exp(-a) + 1.0-exp(-a),
		LowerGamma(2, a), 
		eps);

}

int lm(int l, int m) {
  return lm_index(l, m);
}
void test_mod_spherical_bessel_nd(int max_n, dcomplex x) {

  double h(abs(x)*0.001);
  dcomplex ih(0, h);
  dcomplex* v00 = new dcomplex[2*max_n+2];
  dcomplex* vp0 = new dcomplex[2*max_n+2];
  dcomplex* vm0 = new dcomplex[2*max_n+2];
  dcomplex* v0p = new dcomplex[2*max_n+2];
  dcomplex* v0m = new dcomplex[2*max_n+2];

  ModSphericalBessel(x,    max_n, v00);
  ModSphericalBessel(x+h,  max_n, vp0);
  ModSphericalBessel(x-h,  max_n, vm0);
  ModSphericalBessel(x+ih, max_n, v0p);
  ModSphericalBessel(x-ih, max_n, v0m);

  dcomplex ii(0, 1);
  for(int n = 0; n < max_n+1; n++) {
    dcomplex nd = (vp0[n] - vm0[n] + ii*v0m[n] - ii*v0p[n]) / (4.0*h);
    dcomplex calc = v00[max_n+1+n];
    EXPECT_C_EQ(0.0, (nd-calc)/nd)
      << "x   :" << x    << endl
      << "nd  :" << nd   << endl
      << "calc:" << calc << endl
      << "n =  " << n;
  }
  
  delete[] v00;
  delete[] vp0;
  delete[] vm0;
  delete[] v0p;
  delete[] v0m;

}
void test_mod_spherical_bessel(dcomplex x) {

  int max_n(3);
  dcomplex* vs = new dcomplex[2*max_n+2];
  ModSphericalBessel(x, max_n, vs);

  dcomplex y0, y1, y2, y3, d0, d1, d2;
  y0 = sinh(x)/x;
  d0 = cosh(x)/x -sinh(x)/(x*x);
  y1 = (x*cosh(x)-sinh(x))/(x*x);
  d1 = (cosh(x)+x*sinh(x)-cosh(x))/(x*x) -2.0*(x*cosh(x)-sinh(x))/(x*x*x);
  y2 = ((x*x+3.0)*sinh(x)-3.0*x*cosh(x))/(x*x*x);
  d2 = ((2.0*x*sinh(x) + (x*x+3.0)*cosh(x) -3.0*cosh(x) -3.0*x*sinh(x)))/(x*x*x) +
    -3.0*((x*x+3.0)*sinh(x)-3.0*x*cosh(x))/(x*x*x*x);
  y3 = ((x*x*x+15.0*x)*cosh(x)-(6.0*x*x+15.0)*sinh(x))/pow(x, 4.0);
  //  d3 = ((3.0*x*x+15.0)*cosh(x) + (x*x*x+15.0*x)*sinh(x) - 12.0*x*sinh(x) -(6.0*x*x+15.0)*cosh(x)) / pow(x, 4.0) - 4.0 * ((x*x*x+15.0*x)*cosh(x)-(6.0*x*x+15.0)*sinh(x)) / pow(x, 5.0);
  EXPECT_C_EQ(y0, vs[0])         << x;
  EXPECT_C_EQ(d0, vs[max_n+1+0]) << x;
  EXPECT_C_EQ(y1, vs[1])         << x;
  EXPECT_C_EQ(d1, vs[max_n+1+1]) << x;
  EXPECT_C_EQ(y2, vs[2])         << x;
  EXPECT_C_EQ(d2, vs[max_n+1+2]) << x;
  EXPECT_C_EQ(y3, vs[3])         << x;
  /*
  EXPECT_C_EQ(0.0, (d3-vs[max_n+1+3])/d3)
    << "expect:" << d3            << endl 
    << "actual:" << vs[max_n+1+3] << endl
    << "x = " << x;
    */

  delete[] vs;

}
void test_mod_spherical_bessel_rec(dcomplex x) {

  int max_n(5);
  dcomplex* vs = new dcomplex[2*max_n+2];
  ModSphericalBessel(x, max_n, vs);  

  for(int n = 2; n < max_n; n++) {
    EXPECT_C_EQ(1.0, (vs[n-1]-vs[n+1])*x/ ((2*n+1.0)*vs[n]))
      << "l : " << (vs[n-1]-vs[n+1])*x << endl
      << "r : " << (2*n+1.0)*vs[n]     << endl
      << "x = " << x                   << endl
      << "n = " << n;
    //    EXPECT_C_EQ( (2*n+1.0)*vs[n], vs[n-1]*x-vs[n+1]*x) << x;
  }  

  delete[] vs;

}
void test_mod_spherical_bessel_0(){

  int max_n(5);
  dcomplex* vs = new dcomplex[2*max_n+2];
  ModSphericalBessel(0.0, max_n, vs);

  for(int n = 0; n <= max_n; n++) {
    dcomplex in0 = n==0 ? 1.0 : 0.0;
    EXPECT_C_EQ(in0, vs[n]) << n;
  }

  delete[] vs;
  
}
TEST(Angmoment, mod_spherical_bessel) {
  
  dcomplex x;

  x =5.54904;
  test_mod_spherical_bessel(x);

  x = dcomplex(1.1, 0.0);
  test_mod_spherical_bessel(x);
  test_mod_spherical_bessel_rec(x);
  // test_mod_spherical_bessel_nd(5, x);

  x = dcomplex(1.1, 0.3);
  test_mod_spherical_bessel(x);
  test_mod_spherical_bessel_rec(x);
  // test_mod_spherical_bessel_nd(5, x);

  x = dcomplex(0.1, 0.002);
  test_mod_spherical_bessel(x);  
  test_mod_spherical_bessel_rec(x);
  // test_mod_spherical_bessel_nd(5, x);

  x = dcomplex(10.0, 2.5);
  test_mod_spherical_bessel(x);
  test_mod_spherical_bessel_rec(x);
  test_mod_spherical_bessel_nd(5, x);

  x = dcomplex(0.022, 0.001);
  test_mod_spherical_bessel(x);
  test_mod_spherical_bessel_rec(x);

  // -- small x --
  test_mod_spherical_bessel_0();

  x = dcomplex(0.019, 0.001);
  test_mod_spherical_bessel_rec(x);

  x = dcomplex(0.01, 0.01);
  test_mod_spherical_bessel_rec(x);

  x = dcomplex(0.0012, 0.0001);
  test_mod_spherical_bessel_rec(x);

  x = dcomplex(0.00012, 0.0);
  test_mod_spherical_bessel_rec(x);
  
}
TEST(Angmoment, mod_spherical_bessel_n0) {

  int max_n(0);
  dcomplex* vs = new dcomplex[2*max_n+2];
  dcomplex x;

  // large x
  x = 1.1;
  ModSphericalBessel(x, max_n, vs);
  dcomplex y0 = sinh(x)/x;
  dcomplex d0 = cosh(x)/x -sinh(x)/(x*x);
  EXPECT_C_EQ(y0, vs[0]);
  EXPECT_C_EQ(d0, vs[1]);

  
  // small x
  x = 0.001;
  ModSphericalBessel(x, max_n, vs);
  y0 = sinh(x)/x;
  d0 = cosh(x)/x -sinh(x)/(x*x);
  EXPECT_C_EQ(y0, vs[0]);
  EXPECT_C_EQ(d0, vs[1]);


  delete[] vs;
  
}
void test_associated_legendre(int lmax, dcomplex x) {

  dcomplex* vs = new dcomplex[num_lm_pair(lmax)];
  AssociatedLegendre(x, lmax, vs);  
  EXPECT_C_EQ(1.0, vs[lm(0, +0)]) << 0 << 0 << x;

  if(lmax > 0) {
    EXPECT_C_EQ(-sqrt(1.0-x*x), vs[lm(1, +1)]) << 1 << 1 << x;
    EXPECT_C_EQ(x,              vs[lm(1, +0)]) << 1 << 0 << x;;  
  }

  if(lmax > 1) {
    EXPECT_C_EQ(0.5*(3.0*x*x-1.0),  vs[lm(2, 0)]) << 2 << 0 << x;
    EXPECT_C_EQ(-3.0*x*sqrt(1.0-x*x), vs[lm(2, 1)]) << 2 << 1 << x;
    EXPECT_C_EQ(3.0*(1.0-x*x),      vs[lm(2,+2)]) << 2 << 2 << x;
  }

  for(int L = 0; L <= lmax; L++) 
    for(int M = 1; M <= L; M++) 
      EXPECT_C_EQ(vs[lm(L, -M)], pow(-1.0, M)*Factorial(L-M)*1.0/Factorial(L+M)*vs[lm(L, M)]) << L << M << x;

  delete[] vs;  
}
TEST(Angmoment, associated_legendre) {
  
  test_associated_legendre(4, 1.0);
  test_associated_legendre(4, 0.999);
  test_associated_legendre(4, -1.0);
  test_associated_legendre(4, 0.344);

  test_associated_legendre(0, 0.33);
  test_associated_legendre(1, 0.21);
}
void test_real_spherical_harm(dcomplex t, dcomplex p) {

  int max_l(10);
  dcomplex* vs = new dcomplex[num_lm_pair(max_l)];
  RealSphericalHarmonics(t, p, max_l, vs);

  ostringstream oss;
  oss << "theta = " << t << ",  phi = " << p;
  string msg = oss.str();
  
  dcomplex a = 1.0/sqrt(4.0*M_PI);
  EXPECT_C_EQ(a,                         vs[lm(0, 0)])<< msg;
  EXPECT_C_EQ(a*sqrt(3.0)*cos(t),        vs[lm(1, 0)])<< msg;
  EXPECT_C_EQ(a*sqrt(3.0)*sin(t)*cos(p), vs[lm(1,+1)])<< msg;
  EXPECT_C_EQ(a*sqrt(3.0)*sin(t)*sin(p), vs[lm(1,-1)])<< msg;
  
  EXPECT_C_EQ(a*sqrt(5.0/4.0)*(3.0*cos(t)*cos(t)-1.0),   vs[lm(2, 0)])<< msg;
  EXPECT_C_EQ(a*sqrt(15.0)*sin(t)*cos(t)*cos(p),         vs[lm(2,+1)])<< msg;
  EXPECT_C_EQ(a*sqrt(15.0)*sin(t)*cos(t)*sin(p),         vs[lm(2,-1)])<< msg;
  EXPECT_C_EQ(a*sqrt(15.0/4.0)*sin(t)*sin(t)*cos(2.0*p), vs[lm(2,+2)])<< msg;
  EXPECT_C_EQ(a*sqrt(15.0/4.0)*sin(t)*sin(t)*sin(2.0*p), vs[lm(2,-2)])<< msg;

  delete[] vs;
}
TEST(Angmoment, real_spherical_harm) {

  dcomplex t = dcomplex(1.1, 0.2);
  dcomplex p = dcomplex(1.1, 0.1);
  test_real_spherical_harm(t, p); 
  test_real_spherical_harm(0.0, 0.0); 

}
TEST(Angmoment, cg_coef) {
  int j1 = 2;
  int m1 = 1;
  EXPECT_DOUBLE_EQ(-1.0/sqrt(2*j1+1),
		   cg_coef(j1,j1,m1,-m1, 0,0));
  EXPECT_NEAR(pow(Factorial(2*j1),2)/(6 * sqrt(Factorial(8)*1.0)),
	      cg_coef(j1,j1,m1,-m1,2*j1,0),
	      0.0000000001);
  std::cout << cg_coef(1,1,1,0,0,0) << std::endl;
}
TEST(Angmoment, lm_pair) {

  for(int L = 0; L < 4; L++)
    for(int M = -L; M <= L; M++)
      EXPECT_TRUE(is_lm_pair(L, M));

  EXPECT_FALSE(is_lm_pair(0, -1));
  EXPECT_FALSE(is_lm_pair(2, -3));
  EXPECT_FALSE(is_lm_pair(2, 4));

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
TEST(Angmoment, gto_00_r) {

  double pi(M_PI);  
  int num_r(5);
  dcomplex* rs = new dcomplex[num_r];
  dcomplex* vs = new dcomplex[num_r];
  
  for(int i = 0; i < num_r; i++) {
    rs[i] = 0.1 + 0.2*i;
  }
  dcomplex zeta(1.0, 0.0);
  int L(0);
  int M(0);
  dcomplex* work = new dcomplex[num_lm_pair(L) + 2*L + 2];
  gto_00_r(0.0, 0.0, 0.0, L, M, zeta, rs, num_r, work, vs);
  for(int i = 0; i < num_r; i++) {
    dcomplex r(rs[i]);
    dcomplex calc(sqrt(4.0*pi)*exp(-zeta*r*r));    
    EXPECT_C_EQ(vs[i], calc) << i;    
  }
  delete[] vs;
  delete[] rs;
  delete[] work;
}

TEST(MolFunc, IncGamma) {

  /*
    python code :
    from scipy.special import hyp1f1
    def inc_gamma(m, z):
        return 1.0/(2*m+1) * hyp1f1(m+0.5, m+1.5, -z)
    print inc_gamma(0, 50.0-5.50j)
    print inc_gamma(1, 50.0-5.50j)
    print inc_gamma(2, 50.0-5.50j)
   */

  dcomplex* us = new dcomplex[11];
  IncompleteGamma(10, dcomplex(0.1, -0.2), us);
  EXPECT_C_EQ(dcomplex(0.96392503253589557,
		       +0.06262986651323893),		       
	      us[0]);
  EXPECT_C_EQ(dcomplex(0.100986690286, 0.0166268474338),
	      us[4]);

  IncompleteGamma(10, dcomplex(21.0, 15.5), us);
  EXPECT_C_EQ(dcomplex(-3.59002512901*pow(10.0, -6),
		       -0.000190939738972),
	      us[2]);
  
  IncompleteGamma(10, dcomplex(50.0, -5.5), us);

  double step(0.00001);
  EXPECT_C_EQ(dcomplex(0.124767690392,
		       0.00684158939256),
	      us[0]);  
  EXPECT_C_EQ(dcomplex(0.0012253247264,
		       0.00020320161383),
	      us[1]);
  EXPECT_C_EQ(dcomplex(3.5657718077638774*step,
		       1.0018397403429266*step),
	      us[2]);

  IncompleteGamma(10, 0.0, us);
  EXPECT_C_EQ(1.0, us[0]);
  EXPECT_C_EQ(1.0/3.0, us[1]);
  EXPECT_C_EQ(1.0/5.0, us[2]);

  delete us;

}
TEST(MolFunc, ExpIncGamma) {

  /*
    see ref_molfunc.py
    python code :
    from scipy.special import hyp1f1
    def inc_gamma(m, z):
        return 1.0/(2*m+1) * hyp1f1(1; m+1.5; z)
    print inc_gamma(0, 50.0-5.50j)
    print inc_gamma(1, 50.0-5.50j)
    print inc_gamma(2, 50.0-5.50j)
   */

  dcomplex* us = new dcomplex[11];
  dcomplex ref;
  ExpIncompleteGamma(10, dcomplex(40.0, 40.0), us);
  ref = dcomplex(0.00624843963457, -0.0063295853572);
  EXPECT_C_EQ(ref, us[0]);
  ref = dcomplex(0.00625050716077, -0.00617138734381);
  EXPECT_C_EQ(ref, us[1]);
  ref = dcomplex(0.00624851650343, -0.0060170894781);
  EXPECT_C_EQ(ref, us[2]);

  ExpIncompleteGamma(10, dcomplex(5.0, -4.3), us);
  ref = dcomplex(0.0547155856638, 0.056197386053);
  EXPECT_C_EQ(ref, us[0]);
  ref = dcomplex(0.0571173928686, 0.0435012192617);
  EXPECT_C_EQ(ref, us[1]);
  ref = dcomplex(0.0540860523735, 0.0334636392627);
  EXPECT_C_EQ(ref, us[2]);

  delete us;

}

TEST(Eigen, GenEig) {

  int n(3);
  Eigen::MatrixXcd H(n, n);
  Eigen::MatrixXcd S(n, n);
  H <<
    dcomplex(1.0, 0.1), 0.2, dcomplex(0.1, -0.1),
    0.2,                0.3, 0.4,
    dcomplex(0.2, -0.1),0.4, dcomplex(1.0, -0.2);
  S <<
    1.0, 0.2, 0.2,
    0.2, 0.1, 0.1,
    0.2, 0.1, 1.0;
  SymGenComplexEigenSolver solver(H, S);

  for(int i = 0; i < n; i++) {

    VectorXcd Ci = solver.eigenvectors().col(i);
    VectorXcd lexp = H*Ci;
    VectorXcd rexp = S*Ci*solver.eigenvalues()(i);
    EXPECT_C_EQ(0.0, (lexp-rexp).norm()) << i;

  }

}
TEST(Eigen, GenEig2) {

  int n(3);
  Eigen::MatrixXcd H(n, n);
  Eigen::MatrixXcd S(n, n);
  H <<
    dcomplex(1.0, 0.1), 0.2, dcomplex(0.1, -0.1),
    0.2,                0.3, 0.4,
    dcomplex(0.2, -0.1),0.4, dcomplex(1.0, -0.2);
  S <<
    1.0, 0.2, 0.2,
    0.2, 0.1, 0.1,
    0.2, 0.1, 1.0;
  MatrixXcd C;
  VectorXcd eig;

  generalizedComplexEigenSolve(H, S, &C, &eig);

  for(int i = 0; i < n; i++) {

    VectorXcd Ci = C.col(i);
    VectorXcd lexp = H*Ci;
    VectorXcd rexp = S*Ci*eig(i);
    EXPECT_C_EQ(0.0, (lexp-rexp).norm()) << i;

  }

}
TEST(Eigen, sort) {

  int num(10);
  VectorXcd xs = VectorXcd::Random(num);
  MatrixXcd cs = MatrixXcd::Random(num, num);
  
  SortEigs(xs, cs);

  for(int i = 0; i < num-1; i++) 
    EXPECT_TRUE(xs(i).real() < xs(i+1).real());

}


int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}

