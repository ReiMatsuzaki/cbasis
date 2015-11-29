#ifndef HATOM_HPP
#define HATOM_HPP

#include <stdexcept>

#include "exp_func.hpp"
#include "lin_func.hpp"

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
    int l_;
    F z_;

  public:
    // ---------- type -------------
    typedef ExpFunc<F, 1> STO;
    typedef LinFunc<ExpFunc<F, 1> > STOs;
    typedef OpAdd<OpScalarProd<F, OpD2>,
		  OpAdd<OpScalarProd<F, OpRm>,
			OpScalarProd<F, OpRm> > > HamiltonianOp;
    typedef OpAdd<HamiltonianOp, OpRm> HMinusEnergyOp;

  public:
    // -------- constructor --------
    HLikeAtom();
    HLikeAtom(int _n, int _l): n_(_n), l_(_l), z_(1) {}
    HLikeAtom(int _n, int _l, F _z): n_(_n), l_(_l), z_(_z) {}

    // -------- Getter -------------
    int n() const { return n_; }
    int l() const { return l_; }
    F z() const   { return z_; }

    // -------- operator -----------
    HamiltonianOp Hamiltonian() const {

      return AddOp(ProdOp(       -F(1)/F(2),    OpD2()),
		   AddOp(ProdOp(  -z_,          OpRm(-1)),
			 ProdOp( F(l_*l_-l_)/F(2), OpRm(-2))));

    }
    template<class OpT> OpT HMinusEnergy(F ene) const {

      return AddOp(Hamiltonian(), 
		   ProdOp(-ene, OpRm(0)));

    }

    // -------- state vector -------
    STOs EigenState() const;
    STOs DipoleInitLength(int l1) const;
    STOs DipoleInitVelocity(int l1) const;
    F EigenEnergy() const;
  };
}
#endif



