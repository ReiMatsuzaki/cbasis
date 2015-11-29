#include "lin_func.hpp"
#include "math_utils.hpp"
#include "exp_func.hpp"
#include "cut_exp.hpp"

#include "lin_func_impl.hpp"

namespace l2func {

  // ==== Explicit Decralation ====
  template class  LinFunc<RSTO>;
  template class  LinFunc<CSTO>;
  template class  LinFunc<RGTO>;
  template class  LinFunc<CGTO>;

  template class  LinFunc<CutRSTO>;
  template class  LinFunc<CutCSTO>;
  template class  LinFunc<CutRGTO>;
  template class  LinFunc<CutCGTO>;
}
