#ifndef LIN_FUNC_TEMPLATE_H
#define LIN_FUNC_TEMPLATE_H

#include <vector>
#include "linspace.hpp"
#include "func.hpp"


/**
   Finite linear combination of functions 
*/

namespace l2func {

  template<class FuncT>
  class LinFunc {    
  public:
    // ---- type ----
    typedef typename FuncT::Field Field;
    typedef std::vector<std::pair<Field, FuncT> > FFs;
    typedef typename FFs::const_iterator const_iterator;
    typedef typename FFs::iterator iterator;

  private:
    typedef Field F;
    
  private:
    // ---- Field ----
    FFs coef_func_list_;
    
  public:
    // ---- constructors ----
    LinFunc() { is_l2func<FuncT>(); }

    // ---- accessors ----
    int size() const { return coef_func_list_.size(); }
    const_iterator begin() const { return coef_func_list_.begin(); }
    const_iterator end() const { return coef_func_list_.end(); }
    iterator begin() { return coef_func_list_.begin(); }
    iterator end() { return coef_func_list_.end(); }

    // ---- methods ----
    Field at(F x) const;
    void Add(F c, FuncT f) { coef_func_list_.push_back(std::make_pair(c, f));}
    void Add(const LinFunc<FuncT>& o);
    void SetComplexConjugate();
    //    void SetDerivParam();
    void SetScalarProd(F c);
    void SetRmProd(int n);

  };  

  // ==== Traits ====
  struct linfunc_tag :public func_tag {};
  template<class T> struct func_traits<LinFunc<T> > {
    typedef linfunc_tag func_tag;
  };
  template<class T> struct is_l2func<LinFunc<T> > {};
  template<class T> struct is_compound<LinFunc<T> > : public boost::true_type {};
  template<class T> struct is_fundamental<LinFunc<T> > : public boost::false_type {};

  // ==== Externals ====
  //  template<class FuncT>
  //  std::ostream& operator << (std::ostream& os, const LinFunc<FuncT>& a);

  template<class FuncT>
  LinFunc<FuncT> Expand(const LinFunc<LinFunc<FuncT> >& a) {

    LinFunc<FuncT> res;
    for(typename LinFunc<LinFunc<FuncT> >::const_iterator it_a = a.begin();
	it_a != a.end(); ++it_a) {
      for(typename LinFunc<FuncT>::const_iterator it = it_a->second.begin();
	  it != it_a->second.end(); ++it) {
	res.Add(it_a->first * it->first, it->second);
      }
    }
    return res;
  }

  //  template<class FuncT>
  //  LinFunc<FuncT> Expand(const LinFunc<LinFunc<FuncT> >& a);

}
#endif
