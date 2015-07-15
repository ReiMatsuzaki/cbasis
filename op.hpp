#ifndef OP_HPP
#define OP_HPP

#include <vector>
#include <boost/function.hpp>
#include "lcomb.hpp"

namespace {
  using std::vector;
  using std::pair;
  using std::make_pair;
}

namespace l2func {
  
  template<class Prim>
  class LinearOp {
  public:
    // -------------- typedef -------------
    typedef LinearComb<Prim> LC;
    typedef function<LC (const Prim&)> OP;
    typedef typename Prim::Field Field;
    typedef vector<pair<Field, OP> > F_OPs;
    typedef typename F_OPs::const_iterator cIT;
  private:
    F_OPs coef_op_list_;
  public:
    // ------- Constructors -----------
    LinearOp() { IsPrimitive<Prim>(); }
    // ------- setter -----------------
    void Add(Field c, OP op) {
      coef_op_list_.push_back(make_pair(c, op));
    }
    // ------- getter -----------------
    int size() const { 
      return coef_op_list_.size(); 
    }
    // -------- operate ---------------
    LC OperatePrim(const Prim& a) const {
      LC acc;
      for(cIT it = coef_op_list_.begin(),
	    end = coef_op_list_.end(); it != end; ++it) {
	Field c = it->first;
	OP op = it->second;
	LC tmp = op(a);
	tmp.ScalarProduct(c);
	acc += tmp;
      }
      return acc;
    }
    LC OperateLC(const LC& a) const {

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
    LC operator () (const Prim& a) const {
      return OperatePrim(a);
    }
   
    LC operator () (const LC& a) const {
      return OperateLC(a);
    }
  };
  
  template<class Prim>
  LinearOp<Prim> LinearOpRM(int m) {
    LinearOp<Prim> acc;
    typename Prim::Field one(1);
    acc.Add(one, OpRm<Prim>(m));
    return acc;
  }

}



/*
#include <vector>
#include <boost/functio.hpp>
#include "lcomb.hpp"

namespace {
  using std::vector;
  using std::pair;
  using std::make_pair;
  
  using boost::function;
}

namespace l2func {

  template<class Prim>
  class LinearOp {

  public:
    // ---------- typedef ----------------
    typedef function<LinearComb<Prim>(const Prim&)> OP;
    typedef Prim::Field Field;
    typedef vector<pair<Field, OP> >::const_iterator const_iterator;
    
  private:
    // ---------- Member Field ----------
    vector<pair<Field, OP> > coef_op_list_;

  public:
    // ------- Constructor --------------
    LinearOp() { IsPrimitive<Prim>(); }
    LinearOp(OP op) {
      coef_op_list_.resize(1);
      coef_op_list_[0] = make_pair(Field(1), op);
    }
    
    // ------- Getter -------------------
    int size() const { return coef_op_list_.size(); }
    const_iterator begin() const { return coef_op_list_.begin(); }
    const_iterator end()   const { return coef_op_list_.end(); }
    
    // ------- Setter ------------------
    void operator += (const pair<Field, Prim>& coef_op) {
      coef_op_list_.push_back(coef_op);
    }
    void operator += (const LinearOp& o) {
      
      for(const_iterator it = o.begin(); it != o.end(); ++it) 
	coef_op_list_.push_back(*it);
      
    }
    pair<Field, Prim>& operator [] (int i) { return coef_op_list_[i]; }

    // ------- Operation ---------------
    LinearComb<Prim> operator() (const Prim& a) {
      LinearComb<Prim> res;
      for(const_iterator it = coef_op_list_.begin(),
	    it_end = coef_op_list_.end(); it != it_end; ++it) {
	Field c = it->first;
	OP    f = it->second;
	res += c * f(a);
      }
      return res;
  };
}
*/

#endif
