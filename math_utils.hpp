#ifndef MATH_UTILS_TEMPLATE_H
#define MATH_UTILS_TEMPLATE_H

#include <complex>
#include <sstream>
#include <vector>
#include "macros.hpp"

namespace l2func {

  void set3(int* xs);

  template<class F>
  class MultArray3 {
  private:
    F* data_;
    int num_;
    int n0_[3];
    int n1_[3];
  public:
    MultArray3(F* data, int nx0, int nx1, int ny0, int ny1, int nz0, int nz1) {
      n0_[0] = nx0; n0_[1] = ny0; n0_[2] = nz0;
      n1_[0] = nx1; n1_[1] = ny1; n1_[2] = nz1;
      num_ = (nx1-nx0+1)*(ny1-ny0+1)*(nz1-nz0+1);
      data_ = data;
    }
    int idx(int nx, int ny, int nz) {
      return ((nx - n0_[0]) +
	      (ny - n0_[1]) * (n1_[0] - n0_[0] + 1) + 
	      (nz - n0_[2]) * (n1_[0] - n0_[0] + 1) * (n1_[1] - n0_[1] + 1)); 
    }
    void check_index(int nx, int ny, int nz) {
      int index = this->idx(nx, ny, nz);
      if(index < 0   || num_-1 < index ||
	 nx < n0_[0] || n1_[0] < nx ||
	 ny < n0_[1] || n1_[1] < ny ||
	 nz < n0_[2] || n1_[2] < nz) {
	std::string msg;
	std::stringstream ss;
	SUB_LOCATION(msg);
	ss << "index: (" << nx << ", " << ny << ", " << nz << ") "
	   << index << std::endl;
	msg += ss.str();
	throw std::out_of_range(msg);
      }
    }
    void set(int nx, int ny, int nz, F v) {
      data_[idx(nx, ny, nz)] = v;
    }
    void set_safe(int nx, int ny, int nz, F v) {
      check_index(nx, ny, nz);
      set(nx, ny, nz, v);
    }
    F get(int nx, int ny, int nz) {
      return data_[idx(nx, ny, nz)];
    }
    F get_safe(int nx, int ny, int nz) {
      check_index(nx, ny, nz);
      return this->get(nx, ny, nz);
    }
  };

  typedef std::complex<double> dcomplex;

  int Factorial(int n);
  double DFactorial(int n);
  int DoubleFactorial(int n);
  double DDoubleFactorial(int n);

  template<class F>
  class array3 {
  public:
    typedef F Field;

    // ---- member field ----
  private:
    F xs_[3];

  public:
    // ---- constructors ----
    array3(F x, F y, F z) {
      this->xs_[0] = x;
      this->xs_[1] = y;
      this->xs_[2] = z;
    }
    template<class F2>
    array3(const array3<F2>& o) {
      for(int i = 0; i < 3; i++) 
	this->xs_[i] = o[i];
    }

    // ---- accessors ----
    F operator [](int i) const { return xs_[i]; }
    F x() const { return xs_[0]; }
    F y() const { return xs_[1]; }
    F z() const { return xs_[2]; }
    void set_x(F x) { xs_[0] = x; }
    void set_y(F y) { xs_[1] = y; }
    void set_z(F z) { xs_[2] = z; }
  };

  template<class F> F ConjugateIfPossible(F x);

}

#endif
