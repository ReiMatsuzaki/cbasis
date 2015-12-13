#ifndef LINSPACE_TEMPLATE_H
#define LINSPACE_TEMPLATE_H

#include <complex>
#include <boost/array.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/or.hpp>
#include <boost/static_assert.hpp>

namespace l2func {

  // ==== L2 function space ====  
  // ---- base class ----
  template<class _Field, class _Coord>
  struct L2funcComponent {
    typedef _Field Field;
    typedef _Coord Coord;
  };

  // ---- func and op ----
  template<class _Field, class _Coord>
  struct Func : public L2funcComponent<_Field, _Coord> {};

  template<class _Field, class _Coord>
  struct Op : public L2funcComponent<_Field, _Coord> {};

  // ---- meta function ----
  template<class A>
  struct is_lin_separable: boost::mpl::false_ {};

  template<class A, class B>
  struct in_same_space : boost::mpl::and_<
    boost::is_same<typename A::Field, typename B::Field>,
    boost::is_same<typename A::Coord, typename B::Coord> > {};

  // ==== Linear space operation ====
  // ---- zero ----
  template<class _Field, class _Coord> 
  class ZeroFunc : public Func<_Field, _Coord> {};

  // ---- one ----
  struct One {};
  template<class Field, class Coord>
  class OneOp : public Op<Field, Coord>, One {};

  // ---- func_add_func ----
  template<class A, class B>
  struct FuncAdd : public Func<typename A::Field, typename A::Coord> {
    BOOST_STATIC_ASSERT((in_same_space<A, B>::value));
    A left;
    B right;
    FuncAdd(A a, B b): left(a), right(b) {}
  };

  template<class A, class B>
  struct is_lin_separable<FuncAdd<A, B> > : boost::mpl::true_ {};

  template<class A, class B>
  FuncAdd<A, B> func_add_func(const A& a, const B& b) {
    return FuncAdd<A, B>(a, b);
  }
  template<class B>
  B func_add_func(const ZeroFunc<typename B::Field, typename B::Coord>&, const B& b) {
    return b; }
  template<class A>
  A func_add_func(const A& a, const ZeroFunc<typename A::Field, typename A::Coord>&) {
    return a; }
  

  // ---- scalar_mult_func ----
  template<class A>
  struct ScalarFuncMult : public Func<typename A::Field, typename A::Coord> {
    typename A::Field scalar;
    A func;
    ScalarFuncMult(typename A::Field c, A a): scalar(c), func(a) {}
  };

  template<class A>
  struct is_lin_separable<ScalarFuncMult<A> > : boost::mpl::true_ {};

  template<class A>
  ScalarFuncMult<A> scalar_mult_func(typename A::Field c, const A& a) {
    return ScalarFuncMult<A>(c, a);
  }
  
  // ---- op_add_op ----
  template<class A, class B>
  struct OpAdd :public Op<typename A::Field, typename A::Coord> {
    BOOST_STATIC_ASSERT((in_same_space<A, B>::value));
    A left;
    B right;
    OpAdd(A a, B b): left(a), right(b) {}
  };
  
  template<class A, class B>
  struct is_lin_separable<OpAdd<A, B> > : boost::mpl::true_ {};

  template<class A, class B>
  OpAdd<A, B> op_add_op(const A& a, const B& b) { 
    return OpAdd<A, B>(a, b);
  }

  // ---- scalar_mult_op ----
  template<class A>
  struct ScalarOpMult :public Op<typename A::Field, typename A::Coord> {
    typename A::Field scalar;
    A op;
    ScalarOpMult(typename A::Field c, A a): scalar(c), op(a) {}    
  };

  template<class A>
  struct is_lin_separable<ScalarOpMult<A> > : boost::mpl::true_{};

  template<class A>
  ScalarOpMult<A> scalar_mult_op(typename A::Field c, const A& a) {
    return ScalarOpMult<A>(c, a);
  }

  // ==== inner product ====
  // ---- primitive ----
  template<class A, class OpT, class B>
  typename A::Field cip_impl_prim(const A&, const OpT&, const B&) {
    //    BOOST_STATIC_ASSERT((in_same_space<A, B>::value));
    //    BOOST_STATIC_ASSERT((in_same_space<A, OpT>::value));
    //    BOOST_STATIC_ASSERT((boost::mpl::false_::value));
    return typename A::Field(777);
  }

  // ---- linearly separate ----
  template<class A, class OpT, class B>
  typename A::Field 
  cip_impl_lin(const A& a, const OpT& o, const B& b, 
	       typename boost::disable_if<is_lin_separable<A> >::type* =0,
	       typename boost::disable_if<is_lin_separable<OpT> >::type* =0,
	       typename boost::disable_if<is_lin_separable<B> >::type* =0) {
    return cip_impl_prim(a, o, b);
  }

  template<class A,class OpT, class B, class C>
  typename A::Field 
  cip_impl_lin(const A& a, const OpT& o, const FuncAdd<B, C>& bc,
	       typename boost::disable_if<is_lin_separable<A> >::type* =0,
	       typename boost::disable_if<is_lin_separable<OpT> >::type* =0) {
    return cip_impl_lin(a, o, bc.left) + cip_impl_lin(a, o, bc.right);
  }

  template<class A, class OpT, class B>
  typename A::Field 
  cip_impl_lin(const A& a, const OpT& o, const ScalarFuncMult<B>& b,
	       typename boost::disable_if<is_lin_separable<A> >::type* =0,
	       typename boost::disable_if<is_lin_separable<OpT> >::type* =0) {
    return b.scalar * cip_impl_lin(a, o, b.func);
  }

  template<class A, class OpT, class OpS, class B>
  typename A::Field
  cip_impl_lin(const A& a, const OpAdd<OpT, OpS>& o, const B& b,
	       typename boost::disable_if<is_lin_separable<A> >::type* =0) {
    return cip_impl_lin(a, o.left, b) + cip_impl_lin(a, o.right, b);
  }
  template<class A, class OpT, class B>
  typename A::Field
  cip_impl_lin(const A& a, const ScalarOpMult<OpT>& o, const B& b,
	       typename boost::disable_if<is_lin_separable<A> >::type* =0) {
    return o.scalar * cip_impl_lin(a, o.op, b);
  }

  template<class A, class B, class OpT, class C>
  typename A::Field 
  cip_impl_lin(const FuncAdd<A, B>& ab, const OpT& o, const C& c) {
    return cip_impl_lin(ab.left, o, c) + cip_impl_lin(ab.right, o, c);
  }

  template<class A, class OpT, class B>
  typename A::Field 
  cip_impl_lin(const ScalarFuncMult<A>& a, const OpT& o, const B& b) {
    return a.scalar * cip_impl_lin(a.func, o, b);
  }


  // ---- general interface ----
  template<class A, class B>
  typename A::Field cip(const A& a, const B& b) {
    BOOST_STATIC_ASSERT((in_same_space<A, B>::value));
    // return cip_impl_lin(a, OneOp<typename A::Field, typename A::Coord>(), b);
    return cip_impl_lin(a, One(), b);
  }

  template<class A, class O, class B>
  typename A::Field cip(const A& a, const O& o, const B& b) {
    BOOST_STATIC_ASSERT((in_same_space<A, B>::value));
    // return cip_impl_lin(a, OneOp<typename A::Field, typename A::Coord>(), b);
    return cip_impl_lin(a, o, b);
  }


  // ==== For practice ====
  class RealBasisOn1D : public Func<double, double> { 
  public:
    double param; 
    RealBasisOn1D(double p) : param(p) {}
  };
  class RealBasisOn3D : public Func<double, boost::array<double, 3> > {};
  class ComplexBasisOn1D : public Func< std::complex<double>, double> {};

  template<>double cip_impl_prim<RealBasisOn1D, One, RealBasisOn1D>
  (const RealBasisOn1D& a, const One&, const RealBasisOn1D& b) {
    return a.param * b.param;
  }

  template<class T> struct is_fundamental;
  template<class T> struct is_compound;

}

#endif
