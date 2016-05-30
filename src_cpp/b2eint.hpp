#ifndef B2EINT_H
#define B2EINT_H

//#include "symmolint.hpp"
#include <vector>
#include "typedef.hpp"

namespace l2func {

  /**
    Interface for store of two electron integrals.
   */
  class IB2EInt {
  public:
    virtual ~IB2EInt();

    virtual bool Get(int *ib, int *jb, int *kb, int *lb,
		     int *i, int *j, int *k, int *l, int *type, dcomplex *val) = 0;
    virtual bool Set(int ib, int jb, int kb, int lb,
	     int i, int j, int k, int l, dcomplex val) = 0;
    virtual void Reset() = 0;
  };

  class B2EIntMem :public IB2EInt {
  private:
    int capacity_; // capacity of each array.
    int size_;     // size of data
    int idx_;      // used for Get function.
    std::vector<int> ibs, jbs, kbs, lbs;
    std::vector<int> is, js, ks, ls;
    std::vector<int> ts;
    std::vector<dcomplex> vs;
  public:
    B2EIntMem(int num);
    ~B2EIntMem();
    bool Get(int *ib, int *jb, int *kb, int *lb,
	     int *i, int *j, int *k, int *l, int *type, dcomplex *val);
    bool Set(int ib, int jb, int kb, int lb,
	     int i, int j, int k, int l, dcomplex val);
    void Reset();
  };

}
#endif
