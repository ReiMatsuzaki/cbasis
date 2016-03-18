#ifndef ANGMOMENT_HPP
#define ANGMOMENT_HPP

#include <exception>
#include "typedef.hpp"

namespace l2func {

  // ==== Exception class ====

  class ExceptionBadYlm : public std::exception {
  public:
    std::string msg_;
    ExceptionBadYlm(int L, int M, std::string msg);
    ~ExceptionBadYlm() throw();
    virtual const char* what() const throw();
  };

  // ==== Utilities ====
  double cg_coef(int j1, int j2, int m1, int m2, int j3, int m3);
  int lm_index(int L, int M);
  int num_lm_pair(int max_l);
  class Array5Dim {
  private:
    double* xs_;
    int n0_;
    int n1_;
    int n2_;
    int n3_;
    int n4_;    
  public:
    Array5Dim(int n0, int n1, int n2, int n3, int n4) {
      n0_ = n0;
      n1_ = n1;
      n2_ = n2;
      n3_ = n3;
      n4_ = n4;
      xs_ = new double[n0*n1*n2*n3*n4];
    }
    ~Array5Dim() {
      delete xs_;
    }
    double get(int i0, int i1, int i2, int i3, int i4) {
      return xs_[i0 + i1*n0_ + i2*n0_*n1_ + i3*n0_*n1_*n2_ + i4*n0_*n1_*n2_*n3_];
    }
    void set(int i0, int i1, int i2, int i3, int i4, double x) {
      xs_[i0 + i1*n0_ + i2*n0_*n1_ + i3*n0_*n1_*n2_ + i4*n0_*n1_*n2_*n3_] = x;
    }
  };

  // ==== Spherical GTO single expansion ====
  double GTOExpansionCoef(int l, int m, int lp, int lppp,
			  int Jp, int Mp, int Jpp, int Mpp);
  void ModSphericalBessel(dcomplex x, int max_n, dcomplex* res);
  void AssociatedLegendre(dcomplex x, int max_l, dcomplex* res);
  void RealSphericalHarmonics(dcomplex theta, dcomplex phi, int max_l, dcomplex* res);  
  //  void SphericalHarm(int J, int M, dcomplex theta, dcomplex phi);  
  
  // <Y(Jpp,Mpp) | GTO(zeta; (x,y,z); (l, m)>
  void gto_00_r(dcomplex x, dcomplex y, dcomplex z,		    
		     int Jpp, int Mpp,
		     dcomplex zeta,
		     dcomplex* r, int num_r, dcomplex* work, dcomplex* res);
  void gto_lm_r_center(int l, int m, int Jpp, int Mpp,
		       dcomplex* r, int num_r, dcomplex zeta, dcomplex* res);
  void gto_lm_r(int l, int m,       
		dcomplex x, dcomplex y, dcomplex z,      
		int Jpp, int Mpp,   
		dcomplex zeta,
		int lppp_max,
		dcomplex* r, int num_r,
		dcomplex* work,
		dcomplex* res);

  /*
  class GTOCenterExpansion {
  private:
    int lmax_;
    double* coef;
  public:
    void Init(int lmax) {
      lmax_ = lmax;
    }
  };
  */
  
}

#endif
