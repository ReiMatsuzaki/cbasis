#ifndef EIGNN_PLUS_H
#define EIGNN_PLUS_H


#include <complex>
#include <Eigen/Core>

namespace {
  typedef Eigen::MatrixXcd CM;
  typedef Eigen::VectorXcd CV;
}


std::complex<double> cnorm(const CV& v);
void complex_normalize(CV& v);
void col_cnormalize(CM& c);
void matrix_inv_sqrt(const CM& s, CM* s_inv_sqrt);
void SortEigs(Eigen::VectorXcd& eigs, Eigen::MatrixXcd& eigvecs);
void generalizedComplexEigenSolve(const CM& f, const CM& s, CM* c, CV* eig);

class SymGenComplexEigenSolver {
private:
  CM eigenvectors_;
  CV eigenvalues_;
  
public:
  SymGenComplexEigenSolver();
  SymGenComplexEigenSolver(const CM& _A, const CM& _B);
  void compute(const CM& _A, const CM& B);
  CV eigenvalues() const;
  CM eigenvectors() const;
};



#endif
