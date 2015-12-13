#ifndef LINSPACE_TEMPLATE_H
#define LINSPACE_TEMPLATE_H

#include <complex>
#include <boost/array.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/or.hpp>
#include <boost/static_assert.hpp>

namespace l2func {

  // ==== Basic component in L2function ====  
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
  // ---- add ----
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

  // ---- scalar mult ----
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


  // ==== inner product ====
  // ---- primitive ----
  template<class A, class B>
  typename A::Field cip_impl_prim(const A&, const B&) {
    BOOST_STATIC_ASSERT((in_same_space<A, B>::value));
  }

  // ---- linearly separate ----
  template<class A, class B>
  typename A::Field 
  cip_impl_lin(const A& a, const B& b, 
	       typename boost::disable_if<is_lin_separable<A> >::type* =0,
	       typename boost::disable_if<is_lin_separable<B> >::type* =0) {
    return cip_impl_prim(a, b);
  }

  template<class A, class B, class C>
  typename A::Field 
  cip_impl_lin(const A& a, const FuncAdd<B, C>& bc,
	       typename boost::disable_if<is_lin_separable<A> >::type* =0) {
    return cip_impl_lin(a, bc.left) + cip_impl_lin(a, bc.right);
  }

  template<class A, class B>
  typename A::Field 
  cip_impl_lin(const A& a, const ScalarFuncMult<B>& b,
	       typename boost::disable_if<is_lin_separable<A> >::type* =0) {
    return b.scalar * cip_impl_lin(a, b.func);
  }

  template<class A, class B, class C>
  typename A::Field 
  cip_impl_lin(const FuncAdd<A, B>& ab, const C& c) {
    return cip_impl_lin(ab.left, c) + cip_impl_lin(ab.right, c);
  }

  template<class A, class B>
  typename A::Field 
  cip_impl_lin(const ScalarFuncMult<A>& a, const B& b) {
    return a.scalar * cip_impl_lin(a.func, b);
  }


  // ---- general interface ----
  template<class A, class B>
  typename A::Field cip(const A& a, const B& b) {
    BOOST_STATIC_ASSERT((in_same_space<A, B>::value));
    return cip_impl_lin(a, b);
  }


  // ==== For practice ====
  class RealBasisOn1D : public Func<double, double> { 
  public:
    double param; 
    RealBasisOn1D(double p) : param(p) {}
  };
  class RealBasisOn3D : public Func<double, boost::array<double, 3> > {};
  class ComplexBasisOn1D : public Func<std::complex<double>, double> {};

  template<>
  double cip_impl_prim<RealBasisOn1D, RealBasisOn1D>(const RealBasisOn1D& a,
						     const RealBasisOn1D& b) {
    return a.param * b.param;
  }

  template<class T> struct is_fundamental;
  template<class T> struct is_compound;

}

#endif
