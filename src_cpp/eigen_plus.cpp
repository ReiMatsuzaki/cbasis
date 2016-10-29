#include <sstream>
#include <stdexcept>
#include <iostream>

#include <Eigen/Eigenvalues>
#include <Eigen/Core>

#include "typedef.hpp"
#include "macros.hpp"

#include "eigen_plus.hpp"

using namespace std;
using namespace Eigen;

MatrixXcd
m13cd(dcomplex m00, dcomplex m01, dcomplex m02) {
  MatrixXcd M(1, 3);
  M << m00, m01, m02;
  return M;
}

MatrixXcd
m33cd(dcomplex m00, dcomplex m01, dcomplex m02,
      dcomplex m10, dcomplex m11, dcomplex m12,
      dcomplex m20, dcomplex m21, dcomplex m22) {
  MatrixXcd M(3, 3);
  M << m00, m01, m02, m10, m11, m12, m20, m21, m22;
  return M;
}

VectorXi v1i(int i) {
  VectorXi v(1); v << i;
  return v;
}
VectorXi v3i(int i,int j,int k) {
  VectorXi v(3); v<<i,j,k;
  return v;
}
VectorXcd v1cd(dcomplex v) {
  VectorXcd vec(1); vec << v;
  return vec;
}
VectorXcd v3cd(dcomplex v0, dcomplex v1, dcomplex v2) {
  VectorXcd vec(3); vec << v0, v1, v2;
  return vec;
}

double TakeReal(dcomplex x) {
  return x.real();
}
double TakeAbs(dcomplex x) {
  return abs(x);
}
std::complex<double> cnorm(const CV& v) {
  Eigen::ArrayXcd a = v.array();
  return sqrt((a*a).sum());}
void complex_normalize(CV& v) {
  Eigen::ArrayXcd u = v.array();
  std::complex<double> cnorm = sqrt((u*u).sum());
  v = v / cnorm;
}
void col_cnormalize(CM& c) {

  int n = c.cols();
  for(int j = 0; j < n; j++) {

    std::complex<double> cn = cnorm(c.col(j));

    for(int i = 0; i < n; i++)
      c(i, j) /= cn;

  }
}
void matrix_inv_sqrt(const CM& s, CM* s_inv_sqrt) {
  // compute eigen value problem of S.
  /// S D = D V;
  Eigen::ComplexEigenSolver<CM> es;
  es.compute(s, true);
  CV v = es.eigenvalues();
  CM c = es.eigenvectors();

  // ensure ||c_i|| = 1;
  col_cnormalize(c);

  // lambda_ij = delta_ij / sqrt(v_i);
  CV tmp = v.array().inverse().sqrt();
  CM lambda_mat = tmp.asDiagonal();

  // S^(-1/2) = D diag{1/sqrt(v_i)} D^T
  CM c_tr = c;
  c_tr.transposeInPlace();
  *s_inv_sqrt = c * lambda_mat * c_tr;
}
void SortEigs(Eigen::VectorXcd& eigs, Eigen::MatrixXcd& eigvecs,
	      double (*to_real)(dcomplex), bool reverse) {
  
  int n = eigs.size();
  if(n != eigvecs.cols()) {
    string msg; SUB_LOCATION(msg);
    msg += ": size mismatch";
    throw runtime_error(msg);
  }

  for(int i = 1; i < n; i++) {
    for(int j = i-1; j > -1; j--) {
      if(not reverse && to_real(eigs[j]) > to_real(eigs[j+1])) {
	dcomplex tmp = eigs[j];
	eigs[j] = eigs[j+1];
	eigs[j+1] = tmp;
	eigvecs.col(j).swap(eigvecs.col(j+1));
      }
      if(reverse && to_real(eigs[j]) < to_real(eigs[j+1])) {
	dcomplex tmp = eigs[j];
	eigs[j] = eigs[j+1];
	eigs[j+1] = tmp;
	eigvecs.col(j).swap(eigvecs.col(j+1));
      }
    }
  }
    
}
void generalizedComplexEigenSolve(const CM& f, const CM& s, CM* c, CV* eig){

  if(f.rows() != s.rows() || f.rows() == 0) {
    string msg; SUB_LOCATION(msg);
    stringstream oss;
    oss << msg << ": invalid matrix size for f and s." << endl
	<< "s = " << s.rows() << s.cols() << endl
	<< "f = " << f.rows() << f.cols() << endl;
    throw runtime_error(oss.str());
  }

  // s2inv means S^(-1/2)
  CM s2inv;
  matrix_inv_sqrt(s, &s2inv);
  
  // fp means F' = S^(-1/2)FS^(-1/2)
  CM fp = s2inv * f * s2inv;

  // solve F'C' = C diag{e_i}
  Eigen::ComplexEigenSolver<CM> es;
  es.compute(fp, true);
  *c = s2inv *  es.eigenvectors();
  *eig = es.eigenvalues();

  SortEigs(*eig, *c, TakeReal);
    
}
void CanonicalMatrix(const CM& S, double eps, CM* res) {

  int num(S.rows());
  if(num != S.cols()) {
    string msg; SUB_LOCATION(msg);
    msg += ": size mismatch."; 
    throw runtime_error(msg);
  }

  ComplexEigenSolver<CM> es;
  es.compute(S, true);
  CV s;
  CM U;
  s.swap(es.eigenvalues());
  U.swap(es.eigenvectors());
  SortEigs(s, U, TakeAbs, true);
  
  int num_non0(0);
  for(int j = 0; j < num; j++)
    if(abs(s[j]) > eps) {
      num_non0++;
    }
  
  MatrixXcd X(num, num_non0);
  int j_idx(0);
  for(int j = 0; j < num; j++) {
    if(abs(s[j]) > eps) {
      for(int i = 0; i < num; i++) {
	X(i, j_idx) = U(i, j) / sqrt(s(j));
      }
      j_idx++;
    }
  }
  
  *res = CM::Zero(1, 1);
  res->swap(X);

}
void CanonicalMatrixNum(const CM& S, int num0, CM* res) {

  int num_all(S.rows());
  
  if(num_all != S.cols()) {
    string msg; SUB_LOCATION(msg);
    msg += ": size mismatch."; 
    throw runtime_error(msg);
  }
  
  if(num_all < num0) {
    string msg; SUB_LOCATION(msg);
    msg += "num0 must be lesser than num_all";
    throw runtime_error(msg);
  }

  ComplexEigenSolver<CM> es;
  es.compute(S, true);
  CV s;
  CM U;
  s.swap(es.eigenvalues());
  U.swap(es.eigenvectors());
  SortEigs(s, U, TakeAbs, true);

  MatrixXcd X = MatrixXcd::Zero(num_all, num0);
  for(int j = 0; j < num0; j++) {
    for(int i = 0; i < num_all; i++ ) {
      X(i, j) = U(i, j) / sqrt(s(j));
    }
  }
  
  *res = CM::Zero(1, 1);
  res->swap(X);
  
}
void CEigenSolveCanonical(const CM& F, const CM& S, double eps, CM* c, CV* eig) {

  MatrixXcd X;
  CanonicalMatrix(S, eps, &X);
  CM Fp = X.transpose() * F * X;

  ComplexEigenSolver<CM> es;
  es.compute(Fp, true);
  *c = X * es.eigenvectors();
  eig->swap(es.eigenvalues());
  SortEigs(*eig, *c, TakeReal, false);

}
void CEigenSolveCanonicalNum(const CM& F, const CM& S, int num0,
			      CM* c, CV* eig) {

  MatrixXcd X;
  CanonicalMatrixNum(S, num0, &X);
  CM Fp = X.transpose() * F * X;

  ComplexEigenSolver<CM> es;
  es.compute(Fp, true);
  *c = X * es.eigenvectors();
  eig->swap(es.eigenvalues());
  SortEigs(*eig, *c, TakeReal, false);
  
}

SymGenComplexEigenSolver::SymGenComplexEigenSolver() {}
SymGenComplexEigenSolver::SymGenComplexEigenSolver(const CM& A, const CM& B) {

  this->compute(A, B);

}
void SymGenComplexEigenSolver::compute(const CM& A, const CM& B) {

  generalizedComplexEigenSolve(A, B, &this->eigenvectors_, &this->eigenvalues_);

}
CV SymGenComplexEigenSolver::eigenvalues() const {
  return eigenvalues_;
}
CM SymGenComplexEigenSolver::eigenvectors() const {
  return eigenvectors_;
}
