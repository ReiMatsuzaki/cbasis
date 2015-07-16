#include "op.hpp"

namespace l2func {

  template<class Prim> Op<Prim>::Op() {}
  template<class Prim> void Op<Prim>::AddOp(Func op) {
    Field one(1);
    this->AddCoefOp(one, op);    
  }
  template<class Prim> void Op<Prim>::AddCoefOp(Field c, Func f) {
    coef_op_list_.push_back(make_pair(c, f));
  }
  template<class Prim> void Op<Prim>::Add(Func op) {
    this->AddOp(op);
  }
  template<class Prim> void Op<Prim>::Add(Field c, Func op) {
    this->AddCoefOp(c, op);
  }
  template<class Prim> int Op<Prim>::size() const {
    return coef_op_list_.size();
  } 
  template<class Prim> LinearComb<Prim> Op<Prim>::OperatePrim(const Prim& a) const {
    LC acc;
    for(cIT it = coef_op_list_.begin(),
	  end = coef_op_list_.end(); it != end; ++it) {
      Field c = it->first;
      Func op = it->second;
      LC tmp = op(a);
      tmp.ScalarProduct(c);
      acc += tmp;
    }
    return acc;
  }
  template<class Prim> LinearComb<Prim> Op<Prim>::OperateLC(const LC& a) const {

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
  template<class Prim> LinearComb<Prim> Op<Prim>::operator() (const Prim& a) const  {
    return OperatePrim(a);
  }
  template<class Prim> LinearComb<Prim> Op<Prim>::operator() (const LinearComb<Prim>& a) const {
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


}


