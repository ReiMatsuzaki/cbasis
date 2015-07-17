#ifndef LCOMB_HPP
#define LCOMB_HPP

#include <vector>
#include <iostream>
#include "prim.hpp"

namespace{
  using std::vector;
  using std::pair;
  using std::make_pair;
}

namespace l2func {

  // ============== Linear Combination ==================
  template<class Prim>
  class LinearComb {

  public:
    // ----------- typedef -------------------
    typedef typename Prim::Field Field;
    typedef vector<pair<Field, Prim> > VFP;
    typedef typename VFP::const_iterator const_iterator;
    
  private:
    // ----------- Field Member ---------------
    VFP cf_list_; // c: coefficient, f:function
    
  public:
    // ---------- Constructor ---------------
    LinearComb();
    LinearComb(int n);
    LinearComb(const Prim& prim);
    LinearComb<Prim> Clone() const;
    
    // ---------- Getter ------------------
    int size() const { return cf_list_.size(); }
    const_iterator begin() const { return cf_list_.begin(); }
    const_iterator end() const { return cf_list_.end(); }
    pair<Field, Prim>& operator [] (int i) {
      return cf_list_[i]; }
    const Prim& prim_i (int i) const {
      return cf_list_[i].second; 
    }
    const Prim prim_i_copied(int i) const {
      return cf_list_[i].second;
    }
    Field coef_i (int i) const { 
      return cf_list_[i].first; }
    Field at(Field x) const;

    // ---------- Setter -----------------
    void AddCoefPrim(Field, const Prim&);
    void AddPair(pair<Field, Prim> cf);
    void AddOther(const LinearComb<Prim>&);
    void Add(Field, const Prim&);
    void Add(pair<Field, Prim>);
    void Add(const LinearComb<Prim>&);
    void operator += (pair<Field, Prim> cf);
    void operator += (const Prim& f);
    void operator += (const LinearComb<Prim>& o);
    void ScalarProduct(Field c);

  };

  // =============== Functions ==================
  template<class Prim>
  pair<typename Prim::Field, Prim> operator *
  (typename Prim::Field c, const Prim& f) {
   
    return make_pair(c, f);
  }

  // --------------- Complex Inner Product ------
  template<class Prim1, class Prim2>
  typename Prim1::Field CIP(const LinearComb<Prim1>& a,
			    const LinearComb<Prim2>& b){
    typedef LinearComb<Prim1> LC1;
    typedef LinearComb<Prim2> LC2;
    typedef typename LC1::const_iterator IT1;
    typedef typename LC2::const_iterator IT2;

    typename Prim1::Field acc(0);
    for(IT1 it_a = a.begin(), end_a = a.end(); 
	it_a != end_a; ++it_a)
      for(IT2 it_b = b.begin(), end_b = b.end(); 
	  it_b != end_b; ++it_b) 
	acc += it_a->first * it_b->first * 
	  CIP(it_a->second, it_b->second);

    return acc;
  }
  template<class Prim, class L2Func>
  typename Prim::Field CIP(const LinearComb<Prim>& lc1,
			   const L2Func& o) {

    //    IsPrimitive<Prim>();

    typename Prim::Field acc(0);
    typedef typename LinearComb<Prim>::const_iterator IT;
    for(IT it = lc1.begin(), it_end = lc1.end(); it != it_end; ++it) {
      acc += it->first * CIP(it->second, o);
    }
    return acc;
  }
  template<class Prim, class L2Func>
  typename Prim::Field CIP( const L2Func& o,
			    const LinearComb<Prim>& lc1) {
    return CIP(lc1, o);
  }

  // --------- Operator ----------------
  template<class Prim>
  LinearComb<Prim> OperateDDr (const Prim& f);
  template<class Prim>
  LinearComb<Prim> OperateDDr2 (const Prim& f);

  // --------- derivative basis -------------
  // assuming "a" is normalized basis set.
  template<class F>
  LinearComb<ExpBasis<F,1> > D1Normalized(const ExpBasis<F, 1>& a);
  template<class F>
  LinearComb<ExpBasis<F,2> > D1Normalized(const ExpBasis<F, 2>& a);
  template<class F>
  LinearComb<ExpBasis<F,1> > D2Normalized(const ExpBasis<F, 1>& a);
  template<class F>
  LinearComb<ExpBasis<F,2> > D2Normalized(const ExpBasis<F, 2>& a);
  
}

#endif
