#ifndef HATOM_HPP
#define HATOM_HPP

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
    typedef LinFunc<STO> STOs;
    typedef OpAdd<OpScalarProd<F, OpD2>,
		  OpAdd<OpScalarProd<F, OpRm>,
			OpScalarProd<F, OpRm> > > HamiltonianOp;
    typedef OpAdd<HamiltonianOp, 
		  OpScalarProd<F, OpRm> > HMinusEnergyOp;

  public:
    // -------- constructor --------
    HLikeAtom();
    HLikeAtom(int _n, int _l);
    HLikeAtom(int _n, int _l, F _z);

    // -------- Getter -------------
    int n() const { return n_; }
    int l() const { return l_; }
    F z() const   { return z_; }

    // -------- operator -----------
    STOs EigenState() const;
    
    LinFunc<STOs> DipoleInitLength(int l1) const;

    LinFunc<STOs> DipoleInitVelocity(int l1) const;

    F EigenEnergy() const;

    HamiltonianOp Hamiltonian() const;

    HMinusEnergyOp HMinusEnergy(F ene) const;

  };
}
#endif



