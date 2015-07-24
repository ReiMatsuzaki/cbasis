#include "hatom.hpp"
#include <boost/lexical_cast.hpp>
#include "macros.hpp"
#include "prim.hpp"
#include "lcomb.hpp"
#include "op.hpp"


namespace l2func {

  template<class F>
  HLikeAtom<F>::HLikeAtom() :
    n_(1), z_(1.0), l_(0) {}

  template<class F>
  HLikeAtom<F>::HLikeAtom(int _n, F _z, int _l) :
  n_(_n), z_(_z), l_(_l) {}
  
  // ----------- eigen state ------------
  template<class F>
  LinearComb<ExpBasis<F,1> >
  HLikeAtom<F>::EigenState() const {

    STOs lc;

    if( n_ == 1 && l_ ==0) {
      
      lc += ExpBasis<F,1>(2.0, 1, 1.0);
      
    } else if ( n_ == 2 && l_ == 0) {
      
      lc += ExpBasis<F,1>(1.0/sqrt(2.0), 1, 0.5);
      lc += ExpBasis<F,1>(-1.0/(2.0*sqrt(2.0)), 2, 0.5);
      
    } else if ( n_ == 2 && l_ == 1) {
      
      lc += ExpBasis<F,1>(1.0/(2.0 * sqrt(6.0)), 2, 0.5);

    } else if (n_ == 3 && l_ == 0) {

      lc += ExpBasis<F,1>(4.0 / (81.0 * sqrt(30.0)), 3, 1.0/3.0);
      
    } else if (n_ == 3 && l_ == 1) {
      double c = 8.0 / (27.0 * sqrt(6.0));
      double z = 1.0 / 3.0;
      lc += ExpBasis<F,1>(c,      2, z);
      lc += ExpBasis<F,1>(-c/6.0, 3, z);

    } else if (n_ == 3 && l_ == 2) {

      double c = (1.0 / 81.0) * sqrt(8.0/15.0);
      lc += ExpBasis<F, 1>(c, 3, 1/3.0);
      
    } else {

      string msg;
      SUB_LOCATION(msg);
      msg += "\ninputted n and l is not supported for HLikeAtom::EigenState\n";
      msg += "n: ";
      msg += boost::lexical_cast<string>(n_);
      msg += "\nl: ";
      msg += boost::lexical_cast<string>(l_);
      throw std::invalid_argument(msg);
      
    }

    return lc;
  }
  template<class F>
  LinearComb<ExpBasis<F,1> >
  HLikeAtom<F>::DipoleInitLength(int l1) const {

    if(l_ != l1 + 1 && l_ != l1 - 1) {

      string msg;
      msg = "|l0 - l1| != 1 in HAtomDipoleInitL";
      throw std::invalid_argument(msg);

    }

    typedef ExpBasis<F, 1> Prim;
    typedef LinearComb<Prim> LC;

    LC psi_n = this->EigenState();
    LC res = OpRM<Prim>(1)(psi_n);
    return res;
    
  }
  template<class F>
  LinearComb<ExpBasis<F,1> >
  HLikeAtom<F>::DipoleInitVelocity(int l1) const {

    if(l_ != l1 + 1 && l_ != l1 - 1) {
      string msg;
      msg = "|l0 - l1| != 1 in HAtomDipoleInitV";
      throw std::invalid_argument(msg);
    }

    typedef ExpBasis<F,1> Prim;
    typedef LinearComb<Prim> LC;
    
    LC psi_n = this->EigenState();
    LC a = OpDDr<Prim>()(psi_n);

    if( l_ > 0) {
      double coef;
      LC b;
      if(l_ < l1) 
	coef = - (l_ + 1);
      else 
        coef = l_;
      
      Op<Prim> op;
      op.Add(coef, OpRM<Prim>(-1));
      b = op(psi_n);
      a += b;
    }
    return a;    	
  }
  template<class F>
  double HLikeAtom<F>::EigenEnergy() const {
    return -1.0 / (2.0 * n_ * n_); 
  }

  // ----------- explicit instance ---------
  typedef std::complex<double> CD;
  template class HLikeAtom<double>;
  template class HLikeAtom<CD>;

  /*
  template Op<RSTO> HLikeAtom<double>::Hamiltonian<RSTO>();
  template Op<CSTO> HLikeAtom<CD>::Hamiltonian<CSTO>();
  template Op<RGTO> HLikeAtom<double>::Hamiltonian<RGTO>();
  template Op<CGTO> HLikeAtom<CD>::Hamiltonian<CGTO>();

  template Op<RSTO> HLikeAtom<double>::HMinusEnergy<RSTO>(double e);
  template Op<CSTO> HLikeAtom<CD>::HMinusEnergy<CSTO>(CD);
  template Op<RGTO> HLikeAtom<double>::HMinusEnergy<RGTO>(double);
  template Op<CGTO> HLikeAtom<CD>::HMinusEnergy<CGTO>(CD);
  */
  
}
