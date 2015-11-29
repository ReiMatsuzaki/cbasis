#ifndef OP_TEMPLATE_H
#define OP_TEMPLATE_H

#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/sequence.hpp>

namespace l2func {

  // ==== General Operator ====
  // should be move other source file.
  template<typename OpT> struct op_traits;
  template<typename OpT> struct is_op;
  struct op_tag {};

  // ==== Scalar Product ====
  template<class F, class OpT>
  struct OpScalarProd {
    F c;
    OpT op;
    OpScalarProd(F _c, OpT _op): c(_c), op(_op) {}
  };
  struct scalar_prod_tag : public op_tag {};
  template<class F, class OpT> struct is_op<OpScalarProd<F, OpT> > {};
  template<class F, class OpT> struct op_traits<OpScalarProd<F, OpT> > {
    typedef scalar_prod_tag op_tag;
  };
  template<class F, class OpT>
  std::ostream& operator << (std::ostream& os, const OpScalarProd<F, OpT>& a) {
    os << a.c << " " << a.op;
    return os;
  }
  template<class F, class OpT>
  OpScalarProd<F, OpT> ProdOp(F c, OpT f) { return OpScalarProd<F,OpT>(c,f);}

  // ==== Add ====
  template<class OpA, class OpB>
  struct OpAdd {
    OpA opA;
    OpB opB;
    OpAdd(OpA _opA, OpB _opB) : opA(_opA), opB(_opB) {}
  };
  struct op_add_tag : public op_tag {};
  template<class OpA, class OpB> struct is_op<OpAdd<OpA, OpB> > {};
  template<class OpA, class OpB> struct op_traits<OpAdd<OpA, OpB> > {
    typedef op_add_tag op_tag;
  };
  template<class OpA, class OpB>
  std::ostream& operator << (std::ostream& os, const OpAdd<OpA, OpB>& a) {
    os << a.opA << " + " << a.opB;
    return os;
  }
  template<class OpA, class OpB>
  OpAdd<OpA, OpB> AddOp(OpA a, OpB b) { return OpAdd<OpA, OpB>(a, b); }

  template<int N, class OpT>
  struct NTermOp {
    typedef OpAdd<OpT, typename NTermOp<N-1, OpT>::type> type;
  };
  template<class OpT>
  struct NTermOp<1, OpT> {
    typedef OpT type;
  };

  template<int N, class OpT>
  struct NLinOp {
    typedef OpAdd<OpScalarProd<typename OpT::Field, OpT>,
		  typename NLinOp<N-1, OpT>::type> type;
  };
  template<class OpT>
  struct NLinOp<1, OpT> {
    typedef OpScalarProd<typename OpT::Field, OpT> type;
  };
  
  // ==== D1 operator ====
  class OpD1 {
  public:
    OpD1() {}
  };
  struct d1_tag : public op_tag {};
  template<> struct is_op<OpD1> {};
  template<> struct op_traits<OpD1> {
    typedef d1_tag op_tag;
  };
  std::ostream& operator << (std::ostream& os, const OpD1&) {
    os << (char*)"OpD1";
    return os;
  }

  // ==== D2 operator ====
  class OpD2 {
  public:
    OpD2() {}
  };
  struct d2_tag : public op_tag {};
  template<> struct is_op<OpD2> {};
  template<> struct op_traits<OpD2> {
  typedef d2_tag op_tag;
};
  std::ostream& operator << (std::ostream& os, const OpD2&) {
  os << (char*)"OpD2";
  return os;
}
  
  // ==== Rm operator ====
  class OpRm {
  private:
    int m_;
  public:
    OpRm(int _m):m_(_m) {}
    int m() const {return m_; }
  };
  struct rm_tag : public op_tag {};
  template<> struct is_op<OpRm> {};
  template<> struct op_traits<OpRm> {
  typedef rm_tag op_tag;
};
  std::ostream& operator << (std::ostream& os, const OpRm& a) {
  os << "OpRm(" << a.m() << ")";
  return os;
}

  // ==== Linear Combination ====
  struct printer {
    template<typename F, typename OpT>
    void operator()(std::pair<F, OpT>& a) const {
      std::cout << a.first << a.second<< " + ";
    }
  };

  template<class F, class T0, class T1, class T2> class OpLin {
  private:
    // ---- type ----
    typedef boost::fusion::vector<std::pair<F, T0>, 
				  std::pair<F, T1>, 
				  std::pair<F, T2> > DataT;

    std::vector<F> cs_;
    DataT ops_;
  public:
    OpLin(const std::vector<F>& cs, T0 t0, T1 t1, T2 t2):
      ops_(std::make_pair(cs[0], t0), 
	   std::make_pair(cs[1], t1), 
	   std::make_pair(cs[2], t2)) {}
  
    void print() {
      //    std::cout << boost::fusion::at_c<0>(this->ops_) << std::endl;;
      boost::fusion::for_each(this->ops_, printer());
    }
    const DataT& CoefOp() const { return ops_; }
  };
  struct linop_tag : public op_tag {};
  template<class F, class T1, class T2, class T3> 
  struct is_op<OpLin<F, T1, T2, T3> > {};
  template<class F, class T1, class T2, class T3> 
  struct op_traits<OpLin<F, T1, T2, T3> > {
    typedef linop_tag op_tag;
  };

/*
struct streamer {
  template<typename F, typename OpT>
  std::ostream& operator()(std::ostream& os, const std::pair<F, OpT>& a) const {
    os << a.first << a.second << "+";
    return os;
  }
};

template<class F, class T1, class T2, class T3> 
std::ostream& operator << (std::ostream& os, const OpLin<F, T1, T2, T3>& a) {
  //  a.print();
  return boost::fusion::fold(a.CoefOp(), os, streamer());
}
*/
}
#endif
