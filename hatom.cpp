#include "hatom.hpp"
#include <boost/lexical_cast.hpp>
#include "macros.hpp"
#include "prim.hpp"
#include "lcomb.hpp"
#include "op.hpp"
#include "op_func.hpp"


namespace l2func {
  
  // ----------- eigen state ------------

  template<class F>
  typename HLikeAtom<F>::STOs HLikeAtom<F>::EigenState() const {

    STOs lc;

    if( n_ == 1 && l_ ==0) {
      
      lc.Add(F(1), STO(F(2), 1, F(1)));
      
    } else if ( n_ == 2 && l_ == 0) {

      lc.Add(F(1), STO(F(1)/sqrt(F(2)), 1, F(1)/F(2)));
      lc.Add(F(1), STO(-F(1)/(F(2)*sqrt(F(2))), 2, F(1)/F(2)));
      
    } else if ( n_ == 2 && l_ == 1) {

      lc.Add(F(1), STO(F(1) / (F(2)*sqrt(F(6))), 2, F(1)/F(2)));

    } else if (n_ == 3 && l_ == 0) {

      lc.Add(F(1), STO(F(4)/(F(81)*sqrt(F(30))), 3, F(1)/F(3)));
      
    } else if (n_ == 3 && l_ == 1) {
      F c = F(8) / (F(27) * sqrt(F(6)));
      F z = F(1) / F(3);
      lc.Add(F(1), STO(c, 2, z));
      lc.Add(F(1), STO(-c/F(6), 3, z));

    } else if (n_ == 3 && l_ == 2) {

      F c = F(1) / F(81) * sqrt(F(8)/F(15));
      lc.Add(F(1), STO(c, 3, F(1)/F(3)));
      
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
  typename HLikeAtom<F>::STOs HLikeAtom<F>::DipoleInitLength(int l1) const {

    if(l_ != l1 + 1 && l_ != l1 - 1) {

      string msg;
      msg = "|l0 - l1| != 1 in HAtomDipoleInitL";
      throw std::invalid_argument(msg);

    }

    STOs psi_n = this->EigenState();
    STOs res = Expand(OP(OpRm(1), psi_n));
    return res;
    
  }
  template<class F>
  typename HLikeAtom<F>::STOs HLikeAtom<F>::DipoleInitVelocity(int l1) const {

    if(l_ != l1 + 1 && l_ != l1 - 1) {
      string msg;
      msg = "|l0 - l1| != 1 in HAtomDipoleInitV";
      throw std::invalid_argument(msg);
    }
    
    STOs psi_n = this->EigenState();
    STOs res = Expand(OP(OpD1(), psi_n));

    if( l_ > 0) {
      F coef = F(l_ < l1 ? - (l_ + 1) : l_);
      STOs b = Expand(OP(OpRm(-1), psi_n));
      b.SetScalarProd(coef);
      res.Add(b);
    }

    return res;    	
  }
  template<class F>
  F HLikeAtom<F>::EigenEnergy() const {
    return -F(1) / (F(2) * n_ * n_); 
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
