#ifndef DELTA_TEMPLATE_H
#define DELTA_TEMPLATE_H

#include "prim.hpp"

namespace l2func {

  template<class F>
  class DiracDelta {
    /**
       Represent Dirac's Delta function.
       δ(r-r_0)
       whose width is very small value h.
     */
    
  public:
    // ---- typedef ----
    typedef F Field;

  private:
    // ---- Field Member ----
    double r0_; // location of delta function
    double h_;  // width (most calculation are done as h=0)

  public:
    // ---- Big 3 ----
    DiracDelta();
    DiracDelta(double _r0);
    DiracDelta(const DiracDelta& o);

    // ---- Accessor ----
    double r0() const { return this->r0_; }
    void set_r0(double r0) { this->r0_ = r0; }
    F at(F x) const;

    // ---- C.C. ----
    //    void SetComplexConjugate(const DiracDelta& o);
    //    DiracDelta ComplexConjugate() const;
  };

  // ==== TpyeDef ====
  typedef DiracDelta<CD> CDelta;

  // ==== CIP ====
  template<class F> F CIP(const DiracDelta<F>& a, const CSTO& b);
  template<class F> F CIP(const DiracDelta<F>& a, const CGTO& b);
  template<class F> F CIP(const CSTO& b, const DiracDelta<F>& a);
  template<class F> F CIP(const CGTO& b, const DiracDelta<F>& a);

  // ==== Prim ====
  template<> struct IsPrimitive<DiracDelta<CD> > {};
}

#endif
