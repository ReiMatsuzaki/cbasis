#include <sstream>
#include <iostream>

#include <Eigen/Eigenvalues>
#include <Eigen/Core>

#include "typedef.hpp"
#include "macros.hpp"

#include "eigen_plus.hpp"

using namespace std;
using namespace Eigen;

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
void SortEigs(VectorXcd& eigs, MatrixXcd& eigvecs) {
  
  int n = eigs.size();
  if(n != eigvecs.cols()) {
    string msg; SUB_LOCATION(msg);
    msg += ": size mismatch";
    throw runtime_error(msg);
  }

  for(int i = 1; i < n; i++) {
    for(int j = i-1; j > -1; j--) {
      if(eigs[j].real() > eigs[j+1].real()) {

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

  SortEigs(*eig, *c);
    
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
  SortEigs(s, U);
  
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

  int num1(S.rows());
  
  if(num1 != S.cols()) {
    string msg; SUB_LOCATION(msg);
    msg += ": size mismatch."; 
    throw runtime_error(msg);
  }
  
  if(num1 < num0) {
    string msg; SUB_LOCATION(msg);
    msg += "num0 must be lesser than num1";
    throw runtime_error(msg);
  }

  ComplexEigenSolver<CM> es;
  es.compute(S, true);
  CV s;
  CM U;
  s.swap(es.eigenvalues());
  U.swap(es.eigenvectors());
  SortEigs(s, U);
  
  MatrixXcd X(num1, num0);
  for(int j = 0; j < num0; j++) {
    for(int i = 0; i < num1; i++ ) {
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
  SortEigs(*eig, *c);
}
void CEigenSolveCarnonicalNum(const CM& F, const CM& S, int num0,
			      CM* c, CV* eig) {

  cout << 1 << endl;
  MatrixXcd X;
  CanonicalMatrixNum(S, num0, &X);
  cout << 1.5 << endl;
  CM Fp = X.transpose() * F * X;

  cout << 2 << endl;
  ComplexEigenSolver<CM> es;
  es.compute(Fp, true);
  cout << 3 << endl;
  *c = X * es.eigenvectors();
  eig->swap(es.eigenvalues());
  cout << 4 << endl;
  SortEigs(*eig, *c);
  
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
