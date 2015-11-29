#ifndef FUNC_TEMPLATE_H
#define FUNC_TEMPLATE_H

/**
   General L2-function definitions
 */

namespace l2func {

  // ==== General Func ====
  template<typename L2FuncT> struct func_traits;
  template<typename L2FuncT> struct is_l2func;

  // base class for l2func tag 
  struct func_tag {};

  
  // ==== Add ====
  template<class FuncA, class FuncB>
  struct FuncAdd {
    typedef typename FuncA::Field Field;
    FuncA funcA;
    FuncB funcB;
    FuncAdd(FuncA a, FuncB b): funcA(a), funcB(b) {}
  };
  template<class FuncA, class FuncB>
  FuncAdd<FuncA, FuncB> AddFunc(const FuncA& a, const FuncB& b) { 
    return FuncAdd<FuncA, FuncB>(a, b);
  }
  struct func_add_tag : public func_tag{};
  template<class FuncA, class FuncB> struct func_traits<FuncAdd<FuncA, FuncB> > {
    typedef func_add_tag func_tag;
  };
  template<class FuncA, class FuncB> struct is_l2func<FuncAdd<FuncA, FuncB> > {};


  // ==== Scalar Product ====
  template<class FuncT>
  struct FuncProd {
    typedef typename FuncT::Field Field;
    Field c;
    FuncT f;
    FuncProd(Field _c, FuncT _f): c(_c), f(_f) {}
  };
  template<class FuncT>
  FuncProd<FuncT> ProdFunc(typename FuncT::Field c, FuncT f) {
    return FuncProd<FuncT>(c, f);
  }
  struct func_prod_tag : public func_tag {};
  template<class FuncT> struct func_traits<FuncProd<FuncT> > {
    typedef func_prod_tag func_tag;
  };
  template<class FuncT> struct is_l2func<FuncProd<FuncT> > {};
    
  // used to represent calculation with normalization
  //  enum ENormalized { Normalized };

}
#endif
