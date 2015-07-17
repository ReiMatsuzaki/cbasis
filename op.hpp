#ifndef OP_HPP
#define OP_HPP

#include <vector>
#include <boost/function.hpp>
#include "lcomb.hpp"

namespace {
  using std::vector;
  using std::pair;
  using boost::function;
}

namespace l2func {
  
  template<class Prim>
  class Op {
  public:
    // -------------- typedef -------------
    typedef LinearComb<Prim> LC;
    typedef function<LC (const Prim&)> Func;
    typedef typename Prim::Field Field;
    typedef vector<pair<Field, Func> > F_Funcs;
    typedef typename F_Funcs::const_iterator const_iterator;
    typedef typename F_Funcs::iterator iterator;
  private:
    F_Funcs coef_op_list_;
  public:
    // ------- Constructors -----------
    Op();
    // ------- getter -----------------
    int size() const;
    const_iterator begin() const;
    const_iterator end() const;
    iterator begin();
    iterator end();
    // ------- setter -----------------
    void AddFunc(Func op);
    void AddCoefFunc(Field c, Func op);
    void AddOther(const Op<Prim>&);
    void Add(Func op);
    void Add(Field c, Func op);
    void Add(const Op<Prim>&);
    void ScalarProduct(Field c);
    // -------- operate ---------------
    LC OperatePrim(const Prim& a) const;
    LC OperateLC(const LC& a) const;
    LC operator () (const Prim& a) const;
    LC operator () (const LC& a) const;
  };

  // ------- Factory Functions ------------
  template<class Prim> Op<Prim> OpRM(int m);
  template<class Prim> Op<Prim> OpCst(typename Prim::Field c);
  template<class Prim> Op<Prim> OpDDr();
  template<class Prim> Op<Prim> OpDDr2();

}

#endif
