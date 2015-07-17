#include "lcomb.hpp"

namespace l2func {

  template<class Prim>
  LinearComb<Prim>::LinearComb() { IsPrimitive<Prim>();}
  template<class Prim>
  LinearComb<Prim>::LinearComb(int n):
    cf_list_(n) {}
  template<class Prim>
  LinearComb<Prim>::LinearComb(const Prim& prim):
    cf_list_(1) {
    cf_list_[0] = make_pair(1.0, prim);
  }

  template<class Prim>
  typename Prim::Field LinearComb<Prim>::at_x(Field x) const {

    Field acc(0);

    typedef typename LinearComb<Prim>::const_iterator IT;
    for(IT it = this->begin(), end_it = this->end(); it != end_it; ++it) 
      acc += it->first * it->second.at_x(x);

    return acc;    
  }

  template<class Prim>
  void LinearComb<Prim>::AddCoefPrim
  (Field c, const Prim& o) {
    cf_list_.push_back(make_pair(c, o));
  }
  template<class Prim>
  void LinearComb<Prim>::AddOther
  (const LinearComb<Prim>& o) {

    for(const_iterator it = o.begin(),
	  end = o.end(); it != end; ++it) {
      Field c = it->first;
      Prim  f = it->second;
      cf_list_.push_back(make_pair(c, f));
    }

  }

  template<class Prim>
  void LinearComb<Prim>::Add(Field c, const Prim& o) {
    this->AddCoefPrim(c, o);
  }
  template<class Prim>
  void LinearComb<Prim>::Add(const LinearComb<Prim>& o) {
    this->AddOther(o);
  }
 
  template<class Prim>
  void LinearComb<Prim>::operator += (pair<Field, Prim> cf) {
    cf_list_.push_back(cf);
  }

  template<class Prim>
  void LinearComb<Prim>::operator += (const Prim& f) {
    cf_list_.push_back(make_pair(Field(1), f));
  }
  template<class Prim>
  void LinearComb<Prim>::operator += (const LinearComb<Prim>& o) {
    for(const_iterator it = o.begin(), it_end = o.end();
	it != it_end; ++it ) {
      cf_list_.push_back(*it);
    }
  }

  // explicit instance
  template class LinearComb<RSTO>;
  template class LinearComb<CSTO>;
  template class LinearComb<RGTO>;
  template class LinearComb<CGTO>;
}
