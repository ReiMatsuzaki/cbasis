#ifndef DELTA_TEMPLATE_H
#define DELTA_TEMPLATE_H

#include "func.hpp"
#include <complex>

/**
   Represent Dirac's Delta function.
   δ(r-r_0)
   whose width is very small value h.
*/

namespace l2func {

  template<class F>
  class DiracDelta {
    
  public:
    // ---- type ----
    typedef F Field;

  private:
    // ---- Field Member ----
    double r0_; // location of delta function
    double h_;  // width (most calculation are done as h=0)

  public:
    // ---- Constructors ----
    DiracDelta();
    DiracDelta(double _r0);
    DiracDelta(const DiracDelta<F>& o);

    // ---- Accessor ----
    double r0() const { return this->r0_; }
    void set_r0(double r0) { this->r0_ = r0; }

    // ---- Method ----
    F at(F x) const;
  };

  // ==== External ====
  template<class F>
  std::ostream& operator << (std::ostream& os, const DiracDelta<F>& a);

  // ==== Traits ====
  struct delta_tag :public func_tag {};
  template<class F> struct func_traits<DiracDelta<F> > { 
    typedef delta_tag func_tag; };
  template<class F> struct is_l2func<DiracDelta<F> > {};

  typedef DiracDelta<double> RDelta;
  typedef DiracDelta<std::complex<double> > CDelta;

}

#endif
