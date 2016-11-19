#include <iostream>
#include <Eigen/Core>
#include <gtest/gtest.h>
#include "gtest_plus.hpp"
#include "eigen_plus.hpp"
#include "fact.hpp"
#include "timestamp.hpp"

using namespace std;
using namespace cbasis;
using namespace Eigen;

TEST(EigenPlus, Canonical) {
  PrintTimeStamp("EigenPlus", NULL);
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
  CEigenSolveCanonicalNum(F3, A3, 2, &C2_c, &eig2_c);
  
  CV resid = F3 * C2_c.col(0) - A3 * C2_c.col(0) * eig2_c[0];
  EXPECT_C_NEAR(0.0, resid(0), eps);
  EXPECT_C_NEAR(0.0, resid(1), eps);
  EXPECT_C_NEAR(0.0, resid(2), eps);
  
  resid = F3 * C2_c.col(1) - A3 * C2_c.col(1) * eig2_c[1];
  EXPECT_C_NEAR(0.0, resid(0), eps);
  EXPECT_C_NEAR(0.0, resid(1), eps);
  EXPECT_C_NEAR(0.0, resid(2), eps);

}
TEST(EigenPlus, GenEig) {

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
TEST(EigenPlus, GenEig2) {

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
TEST(EigenPlus, sort) {

  int num(10);
  VectorXcd xs = VectorXcd::Random(num);
  MatrixXcd cs = MatrixXcd::Random(num, num);
  
  SortEigs(xs, cs, TakeReal);

  for(int i = 0; i < num-1; i++) 
    EXPECT_TRUE(xs(i).real() < xs(i+1).real());

}
TEST(EigenPlus, sort_reverse) {

  int num(10);
  VectorXcd xs = VectorXcd::Random(num);
  MatrixXcd cs = MatrixXcd::Random(num, num);
  
  SortEigs(xs, cs, TakeReal, true);

  for(int i = 0; i < num-1; i++) 
    EXPECT_TRUE(xs(i).real() > xs(i+1).real());

  
}

TEST(Fact, iabs) {

  EXPECT_EQ(0, iabs(0));
  EXPECT_EQ(1, iabs(1));
  EXPECT_EQ(1, iabs(-1));
  EXPECT_EQ(4, iabs(4));
  EXPECT_EQ(4, iabs(-4));
  
}
TEST(Fact, ipow) {
  
  EXPECT_EQ(-1, ipow(-1, 1));
  EXPECT_EQ(-1, ipow(-1, -3));
  EXPECT_EQ(1, ipow(-1, -6));
  EXPECT_EQ(1, ipow(-1, 0));

  EXPECT_EQ(1, ipow(3, 0));
  EXPECT_EQ(2, ipow(2, 1));
  EXPECT_EQ(9, ipow(3, 2));
  EXPECT_EQ(81, ipow(3, 4));
  
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

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}


