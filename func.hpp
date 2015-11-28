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

  // used to represent calculation with normalization
  //  enum ENormalized { Normalized };

}
#endif
