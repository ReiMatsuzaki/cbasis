#include <boost/format.hpp>
#include <unsupported/Eigen/FFT>
#include <iostream>
#include <complex>

using namespace boost;
using namespace Eigen;
using namespace std;

typedef complex<double> dcomplex;

dcomplex gau_t(double a, double t0, double w0, double t) {
  /*
    f(t) = exp[-a(t-t0)^2 + iw0 t]
   */
  return exp(dcomplex(-a*pow(t-t0, 2), w0*t));
}
dcomplex fft_gau_w(double a, double t0, double w0, double w) {
  /*
    f(t) = exp[-a(t-t0)^2 + iw0 t]
    Int_{-oo}^{+oo} f(t) e^{-iwt} dt
   */
  double x = (w-w0) / (4.0*a);
  dcomplex expo(x*(-w+w0), x*(-4.0*a*t0));
  return exp(expo) * sqrt(M_PI/a);
}
dcomplex inv_gau_w(double a, double t0, double w0, double w) {
  /*
    f(t) = exp[-a(t-t0)^2 + iw0 t]
    Int_{-oo}^{+oo} f(t) e^{+iwt} dt
   */  
  return fft_gau_w(a, t0, w0, -w);
}
int main () {

  double a0 = 0.77;
  double t0 = 3.0;
  double w0 = 0.8;
  double T = 8.0;
  int num = 200;
  double dt = T / (num);
  vector<dcomplex> t_vec;
  cout << "t,f(t)\n";
  for(int n = 0; n < num; n++) {
    double t = n * dt;
    dcomplex y = gau_t(a0, t0, w0, t);
    t_vec.push_back(y);
    cout << format("%f,%f,%f\n") % t % y.real() % y.imag();
  }
  
  vector<dcomplex> w_vec;
  vector<dcomplex> w_inv_vec;
  FFT<double> fft;
  fft.fwd(w_vec, t_vec);
  fft.inv(w_inv_vec, t_vec);

  cout << "w,fft\n";
  for(int n = 0; n < 20; n++) {
    double w = 2.0*M_PI / T * n;
    dcomplex y = w_vec[n] * dt;
    dcomplex z = fft_gau_w(a0, t0, w0, w);
    cout << format("%f,%f\n") % w % abs(y-z);
  }

  cout << "w,inv_fft\n";
  for(int n = 0; n < 20; n++) {
    double w = 2.0*M_PI / T * n;
    dcomplex calc = w_inv_vec[n] * T;
    dcomplex ref = inv_gau_w(a0, t0, w0, w);
    cout << format("%f,%f\n") % w % abs(calc-ref);
    //    cout << format("%f,%f,%f,%f,%f\n")
    //      % w % calc.real() % calc.imag()
    //      %     ref.real() %  ref.imag();

  }
}
