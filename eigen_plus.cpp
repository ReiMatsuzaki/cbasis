#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include "eigen_plus.hpp"
#include "macros.hpp"
using namespace std;
using namespace Eigen;

typedef complex<double> dcomplex;

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
  if(n != eigvecs.cols() || n != eigvecs.rows()) {
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
