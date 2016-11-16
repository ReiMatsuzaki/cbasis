#ifndef INT_EXP_H
#define INT_EXP_H

#include "../utils/typedef.hpp"

namespace cbasis {

  dcomplex STOInt_Rplus(int n, dcomplex a);
  dcomplex GTOInt_Rplus(int n, dcomplex a);
  dcomplex STO_GTOInt_Rplus(int n, dcomplex a, dcomplex b);

  dcomplex GTOInt_R(int n, dcomplex a);
  dcomplex STO_GTOInt_R(int n, dcomplex a, dcomplex b);
   
}

#endif
