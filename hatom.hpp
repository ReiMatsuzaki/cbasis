#ifndef HATOM_HPP
#define HATOM_HPP

#include <stdexcept>
#include "prim.hpp"
#include "lcomb.hpp"
#include "op.hpp"

namespace {
  using std::string;
}


namespace l2func {

  template<class F=double>
  class HLikeAtom {
  private:
    // ---------- Field ------------
    int n_;
    F z_;
    int l_;

    // ---------- type -------------
    typedef LinearComb<ExpBasis<F,1> > STOs;
    
  public:
    // -------- constructor --------
    HLikeAtom(int _n, F _z, int _l);

    // -------- Getter -------------
    int n() const { return n_; }
    F z() const   { return z_; }
    int l() const { return l_; }

    // -------- operator -----------
    template<class Prim> Op<Prim> Hamiltonian() const {

      Op<Prim> op;
      op.Add(-0.5, OpDDr2<Prim>());
      if(l_ != 0) {
	F c(l_ * (l_ + 1) * 0.5);
	op.Add(c, OpRM<Prim>(-2));
      }
      op.Add(-z_, OpRM<Prim>(-1));
      return op;
    }
    template<class Prim> Op<Prim> HMinusEnergy(F ene) const {
      Op<Prim> op = this->Hamiltonian<Prim>();
      op.AddOther(OpCst<Prim>(-ene));
      return op;
    }

    // -------- state vector -------
    STOs EigenState() const;
    STOs DipoleInitLength(int l1) const;
    STOs DipoleInitVelocity(int l1) const;
    double EigenEnergy() const;
  };
}
#endif



