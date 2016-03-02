#ifndef MATH_UTILS_TEMPLATE_H
#define MATH_UTILS_TEMPLATE_H

#include <complex>

namespace l2func {

/*
  template<class F, int n>
  class array {
    // ---- member field ----
  private:
    F xs_[n];

  public:
    // ---- constructors ----
    array(F* xs) {
      for(int i = 0; i < n; i++) 
	this->xs_[i] = xs[i];
    }
    template<class F2>
    array(const array<F2, n>& o) {
      for(int i = 0; i < n; i++) 
	this->xs_[i] = o[i];
    }

    // ---- accessors ----
    F operator [](int i) const { return xs_[i]; }  
  };
*/

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

  /*
  template<class F, int N>
  class MultiArray;

  template<class F>
  class MultiArray<F, 3> {
  private:
    int num_[3];
    F* data_;
  public:
    MultiArray(int num[3]) {
      int cumprod(1);
      for(int i = 0; i < 3; i++) {
	num_[i] = num[i];
	cumprod *= num[i];
      }
      data_ = new F[cumprod];
    }
    ~MultiArray() {
      delete data_;
    }
    void set(int i0, int i1, int i2, F val) {
      data_[i0 + *num_[1]]
    }
  };
  */

  template<class F> F ConjugateIfPossible(F x);

}

#endif
