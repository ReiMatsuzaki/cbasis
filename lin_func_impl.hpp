#ifndef LIN_FUNC_IMPL_TEMPLATE_H
#define LIN_FUNC_IMPL_TEMPLATE_H

#include "lin_func.hpp"

namespace l2func {

  template<class FuncT>
  void LinFunc<FuncT>::Add(const LinFunc<FuncT>& o) {

    std::copy(o.begin(), o.end(), std::back_inserter(this->coef_func_list_));

  }

  template<class FuncT>
  typename LinFunc<FuncT>::Field LinFunc<FuncT>::at(typename FuncT::Field x) const {

    typedef typename LinFunc<FuncT>::Field F;
    F acc(0);
    for(typename LinFunc<FuncT>::const_iterator it = this->begin();
	it != this->end(); it++) {

      acc += it->first * it->second.at(x);

    }
    return acc;
  }

  template<class FuncT>
  void LinFunc<FuncT>::SetComplexConjugate() {

    for(typename LinFunc<FuncT>::iterator it = this->begin();
	it != this->end(); ++it) {

      it->first = ConjugateIfPossible(it->first);
      it->second.SetComplexConjugate();

    }
  }

  template<class FuncT>
  void LinFunc<FuncT>::SetScalarProd(F c) {

    for(typename LinFunc<FuncT>::iterator it = this->begin();
	it != this->end(); ++it) {

      it->first *= c;

    }
    
  }

  template<class FuncT>
  void LinFunc<FuncT>::SetRmProd(int n) {

    for(typename LinFunc<FuncT>::iterator it = this->begin();
	it != this->end(); ++it) {

      it->second.SetRmProd(n);

    }
  }


  // ==== Externals ====
  template<class FuncT>
  std::ostream& operator << (std::ostream& os, const LinFunc<FuncT>& a) {

    for(typename LinFunc<FuncT>::const_iterator it = a.begin(); 
	it != a.end(); ++it) {
      
      os << it->first << "(" << it->second << ") + ";

    }

    return os;

  }

  template<class FuncT>
  LinFunc<FuncT> Expand(const LinFunc<LinFunc<FuncT> >& a) {

    LinFunc<FuncT> res;
    for(typename LinFunc<LinFunc<FuncT> >::const_iterator it_a = a.begin();
	it_a != a.end(); ++it_a) {
      for(typename LinFunc<LinFunc<FuncT> >::const_iterator it = it_a->begin();
	  it != it_a->end(); ++it) {
	res.Add(it_a->first * it->first, it->second);
      }
    }
    return res;
  }
}

#endif
