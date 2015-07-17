#include <boost/bind.hpp>
#include "op.hpp"

namespace {
  using boost::bind;
}

namespace l2func {

  template<class Prim> Op<Prim>::Op() {}

  template<class Prim> int Op<Prim>::size() const {
    return coef_op_list_.size();
  } 
  template<class Prim> typename Op<Prim>::const_iterator 
  Op<Prim>::begin() const {
    return coef_op_list_.begin();
  }
  template<class Prim> typename Op<Prim>::const_iterator 
  Op<Prim>::end() const {
    return coef_op_list_.end();
  }
  template<class Prim> typename Op<Prim>::iterator Op<Prim>::begin() {
    return coef_op_list_.begin();
  }
  template<class Prim> typename Op<Prim>::iterator Op<Prim>::end() {
    return coef_op_list_.end();
  }

  template<class Prim> void Op<Prim>::AddFunc(Func op) {
    Field one(1);
    this->AddCoefFunc(one, op);    
  }
  template<class Prim> void Op<Prim>::AddCoefFunc(Field c, Func f) {
    coef_op_list_.push_back(make_pair(c, f));
  }
  template<class Prim> void Op<Prim>::AddOther(const Op<Prim>& o) {
    
    for(const_iterator it = o.begin(), end = o.end(); it != end; ++it) {

      coef_op_list_.push_back(*it);

    }

  }
  template<class Prim> void Op<Prim>::Add(Func op) {
    this->AddFunc(op);
  }
  template<class Prim> void Op<Prim>::Add(Field c, Func op) {
    this->AddCoefFunc(c, op);
  }
  template<class Prim> void Op<Prim>::Add(const Op<Prim>& o) {
    this->AddOther(o);
  }
  template<class Prim> void Op<Prim>::ScalarProduct(Field c) {
    
    for(iterator it = this->begin(), end = this->end();
	  it != end; ++it) {

      it->first *= c;

    }

  }

  template<class Prim> LinearComb<Prim> 
  Op<Prim>::OperatePrim(const Prim& a) const {
    LC acc;
    for(const_iterator it = coef_op_list_.begin(),
	  end = coef_op_list_.end(); it != end; ++it) {
      Field c = it->first;
      Func op = it->second;
      LC tmp = op(a);
      tmp.ScalarProduct(c);
      acc += tmp;
    }
    return acc;
  }
  template<class Prim> LinearComb<Prim> 
  Op<Prim>::OperateLC(const LC& a) const {

      LC res;

      for( typename LC::const_iterator it = a.begin(),
	     end = a.end(); it != end; ++it) {
	typename LC::Field c = it->first;
	Prim      f = it->second;
	LC op_f = this->OperatePrim(f);
	op_f.ScalarProduct(c);
	res += op_f;
      }
      
      return res;
    }
  template<class Prim> LinearComb<Prim> 
  Op<Prim>::operator() (const Prim& a) const  {
    return OperatePrim(a);
  }
  template<class Prim> LinearComb<Prim> 
  Op<Prim>::operator() (const LinearComb<Prim>& a) const {
    return OperateLC(a);
  }

  template<class Prim> Op<Prim> OpRM(int m)  {
    Op<Prim> acc;
    acc.Add(bind(OperateRm<Prim>, m, _1));
    return acc;
  }
  template<class Prim> Op<Prim> OpCst(typename Prim::Field c) {
    Op<Prim> acc;
    acc.Add(bind(OperateCst<Prim>, c, _1));
    return acc;
  }
  template<class Prim> Op<Prim> OpDDr() {
    Op<Prim> acc;
    acc.Add(bind(OperateDDr<Prim>, _1));
    return acc;
  }
  template<class Prim> Op<Prim> OpDDr2() {
    Op<Prim> acc;
    acc.Add(bind(OperateDDr2<Prim>, _1));
    return acc;
  }

  // --------- explicit instance ------------
  template class Op<RSTO>;
  template class Op<CSTO>;
  template class Op<RGTO>;
  template class Op<CGTO>;

  template Op<RSTO> OpRM(int);
  template Op<CSTO> OpRM(int);
  template Op<RGTO> OpRM(int);
  template Op<CGTO> OpRM(int);

  template Op<RSTO> OpCst(double);
  template Op<CSTO> OpCst(std::complex<double>);
  template Op<RGTO> OpCst(double);
  template Op<CGTO> OpCst(std::complex<double>);

  template Op<RSTO> OpDDr();
  template Op<CSTO> OpDDr();
  template Op<RGTO> OpDDr();
  template Op<CGTO> OpDDr();

  template Op<RSTO> OpDDr2();
  template Op<CSTO> OpDDr2();
  template Op<RGTO> OpDDr2();
  template Op<CGTO> OpDDr2();

}


