#include <iostream>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <Eigen/Core>

#include "erfc.hpp"
#include "r1_lc.hpp"
#include "r1basis.hpp"

using namespace Eigen;
using namespace std;
using namespace cbasis;
using namespace erfc_mori;
using namespace boost::python;

dcomplex TDot(const VectorXcd& a, const VectorXcd& b) {
  return (a.array() * b.array()).sum();
}
void BindMath() {
  def("tdot", &TDot);
  def("erfc", &erfc<dcomplex>);
  def("erfcx", &erfcx<dcomplex>);
}
void BindR1LC() {
  
  typedef _LC_EXPs<1> _LC_STOs;
  register_ptr_to_python<LC_STOs>();
  def("LC_STOs", &Create_LC_STOs);
  class_<_LC_STOs>("_LC_STOs", init<>())
    .def("size", &_LC_EXPs<1>::size)
    .def("c", (dcomplex(_LC_EXPs<1>::*)(int)const)&_LC_EXPs<1>::c)
    .def("n", (int(_LC_EXPs<1>::*)(int)const)&_LC_EXPs<1>::n)
    .def("z", (dcomplex(_LC_EXPs<1>::*)(int)const)&_LC_EXPs<1>::z)
    .def("add",  &_LC_EXPs<1>::Add, return_self<>())
    .def("str",  &_LC_EXPs<1>::str)
    .def("at_r", &_LC_EXPs<1>::AtR)
    .def("clone", &_LC_EXPs<1>::Clone)
    .def("conj", &_LC_EXPs<1>::Conj)
    .def("__str__", &_LC_EXPs<1>::str);

  
  typedef _LC_EXPs<2> _LC_GTOs;
  register_ptr_to_python<LC_GTOs>();
  def("LC_GTOs", &Create_LC_GTOs);
  class_<_LC_GTOs>("_LC_GTOs", init<>())
    .def("size", &_LC_EXPs<2>::size)
    .def("c", (dcomplex(_LC_EXPs<2>::*)(int)const)&_LC_EXPs<2>::c)
    .def("n", (int(_LC_EXPs<2>::*)(int)const)&_LC_EXPs<2>::n)
    .def("z", (dcomplex(_LC_EXPs<2>::*)(int)const)&_LC_EXPs<2>::z)
    .def("add",  &_LC_EXPs<2>::Add, return_self<>())
    .def("str",  &_LC_EXPs<2>::str)
    .def("at_r", &_LC_EXPs<2>::AtR)
    .def("clone", &_LC_EXPs<2>::Clone)
    .def("conj", &_LC_EXPs<2>::Conj)
    .def("__str__", &_LC_EXPs<2>::str);

}
void BindR1Basis() {

  def("sto_int", &STOInt);
  def("gto_int", &GTOInt);
  def("sto_gto_int", &STO_GTOInt);

  /*
  def("int_lc", &EXPIntLC<1,1>);
  def("int_lc", &EXPIntLC<1,2>);
  def("int_lc", &EXPIntLC<2,1>);
  def("int_lc", &EXPIntLC<2,2>);
  */
  //  def("int_lc", (dcomplex (*)(LC_STOs, int, LC_STOs))EXPIntLC);
  //  def("int_lc", (dcomplex (*)(LC_GTOs, int, LC_GTOs))EXPIntLC);

  def("calc_vec", &CalcVec<1,1>);
  def("calc_vec", &CalcVec<1,2>);
  def("calc_vec", &CalcVec<2,1>);
  def("calc_vec", &CalcVec<2,2>);
  
  def("calc_rm_mat", &CalcRmMat<1,1>);
  def("calc_rm_mat", &CalcRmMat<1,2>);
  def("calc_rm_mat", &CalcRmMat<2,1>);
  def("calc_rm_mat", &CalcRmMat<2,2>);
  
  def("calc_d2_mat", &CalcD2Mat<1,1>);
  def("calc_d2_mat", &CalcD2Mat<1,2>);
  def("calc_d2_mat", &CalcD2Mat<2,1>);
  def("calc_d2_mat", &CalcD2Mat<2,2>);
  
  typedef _EXPs<1> _STOs;
  typedef _EXPs<2> _GTOs;

  register_ptr_to_python<STOs>();
  def("STOs", &Create_STOs);
  class_<_STOs>("_STOs", init<>())
    .def("size",      &_EXPs<1>::size)
    .def("basis",     &_EXPs<1>::basis)
    .def("is_prim",       &_EXPs<1>::IsPrim)
    .def("is_prim_all",   &_EXPs<1>::IsPrimAll)
    .def("is_normal",     &_EXPs<1>::IsNormal)
    .def("is_normal_all", &_EXPs<1>::IsNormalAll)
    .def("has_coef",      &_EXPs<1>::HasCoef)
    .def("has_coef_all",  &_EXPs<1>::HasCoefAll)
    .def("exp_power", &_EXPs<1>::exp_power)
    .def("add",       &_EXPs<1>::AddPrim, return_self<>())
    .def("add",       &_EXPs<1>::AddPrims, return_self<>())
    .def("add",       &_EXPs<1>::AddLC, return_self<>())
    .def("add_not_normal", &_EXPs<1>::AddNotNormalPrim, return_self<>())
    .def("add_not_normal", &_EXPs<1>::AddNotNormalLC, return_self<>())    
    .def("setup",       &_EXPs<1>::SetUp, return_self<>())
    .def("str",         &_EXPs<1>::str)
    .def("at_r",        &_EXPs<1>::AtR)
    .def("at_r",        &_EXPs<1>::AtR_One)
    .def("clone",       &_EXPs<1>::Clone)
    .def("conj",        &_EXPs<1>::Conj)
    .def("calc_rm_mat", &_EXPs<1>::CalcRmMat)
    .def("calc_d2_mat", &_EXPs<1>::CalcD2Mat)
    .def("calc_vec",    &_EXPs<1>::CalcVecSTO)
    .def("calc_vec",    &_EXPs<1>::CalcVecGTO)
    .def("__str__",     &_EXPs<1>::str);

  register_ptr_to_python<GTOs>();
  def("GTOs", &Create_GTOs);
  class_<_GTOs>("_GTOs", init<>())
    .def("size",      &_EXPs<2>::size)
    .def("basis",     &_EXPs<2>::basis)
    .def("is_prim",       &_EXPs<2>::IsPrim)
    .def("is_prim_all",   &_EXPs<2>::IsPrimAll)
    .def("is_normal",     &_EXPs<2>::IsNormal)
    .def("is_normal_all", &_EXPs<2>::IsNormalAll)
    .def("has_coef",      &_EXPs<2>::HasCoef)
    .def("has_coef_all",  &_EXPs<2>::HasCoefAll)    
    .def("exp_power", &_EXPs<2>::exp_power)
    .def("add",       &_EXPs<2>::AddPrim, return_self<>())
    .def("add",       &_EXPs<2>::AddPrims, return_self<>())
    .def("add",       &_EXPs<2>::AddLC, return_self<>())
    .def("add_not_normal", &_EXPs<2>::AddNotNormalPrim, return_self<>())
    .def("add_not_normal", &_EXPs<2>::AddNotNormalLC, return_self<>())        
    .def("setup",     &_EXPs<2>::SetUp, return_self<>())
    .def("str",   &_EXPs<2>::str)
    .def("at_r",  &_EXPs<2>::AtR)
    .def("at_r",  &_EXPs<2>::AtR_One)
    .def("clone", &_EXPs<2>::Clone)
    .def("conj",  &_EXPs<2>::Conj)
    .def("calc_rm_mat", &_EXPs<2>::CalcRmMat)
    .def("calc_d2_mat", &_EXPs<2>::CalcD2Mat)
    .def("calc_vec",    &_EXPs<2>::CalcVecSTO)
    .def("calc_vec",    &_EXPs<2>::CalcVecGTO)
    .def("__str__",     &_EXPs<2>::str);
}

BOOST_PYTHON_MODULE(r1basis_bind) {

  Py_Initialize();
  BindMath();
  BindR1LC();
  BindR1Basis();
  
}
