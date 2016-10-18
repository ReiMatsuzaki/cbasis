#include <iostream>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <Eigen/Core>

#include "r1_lc.hpp"
#include "r1basis.hpp"

using namespace Eigen;
using namespace std;
using namespace cbasis;
using namespace boost::python;

int add(int a, int b) {
  return a+b;
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
    .def("conj", &_LC_EXPs<1>::Conj);
  
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
    .def("conj", &_LC_EXPs<2>::Conj);
}

void BindR1Basis() {

  def("sto_int", &STOInt);
  def("gto_int", &GTOInt);
  def("int_lc", (dcomplex (*)(LC_STOs, int, LC_STOs))EXPIntLC);
  def("int_lc", (dcomplex (*)(LC_GTOs, int, LC_GTOs))EXPIntLC);

  typedef _EXPs<1> _STOs;
  typedef _EXPs<2> _GTOs;

  register_ptr_to_python<STOs>();
  def("STOs", &Create_STOs);
  class_<_STOs>("_STOs", init<>())
    .def("size",      &_EXPs<1>::size)
    .def("basis",     &_EXPs<1>::basis)
    .def("only_prim", &_EXPs<1>::OnlyPrim)
    .def("add",       &_EXPs<1>::AddPrim, return_self<>())
    .def("add",       &_EXPs<1>::AddPrims, return_self<>())
    .def("add_lc",    &_EXPs<1>::AddLC, return_self<>())
    .def("setup",     &_EXPs<1>::SetUp, return_self<>())
    .def("str",   &_EXPs<1>::str)
    .def("at_r",  &_EXPs<1>::AtR)
    .def("clone", &_EXPs<1>::Clone)
    .def("conj",  &_EXPs<1>::Conj)
    .def("calc_rm_mat", &_EXPs<1>::CalcRmMat)
    .def("calc_d2_mat", &_EXPs<1>::CalcD2Mat)
    .def("calc_vec",    &_EXPs<1>::CalcVecSTO)
    .def("calc_vec",    &_EXPs<1>::CalcVecGTO);

  register_ptr_to_python<GTOs>();
  def("GTOs", &Create_GTOs);
  class_<_GTOs>("_STOs", init<>())
    .def("size",      &_EXPs<2>::size)
    .def("basis",     &_EXPs<2>::basis)
    .def("only_prim", &_EXPs<2>::OnlyPrim)
    .def("add",       &_EXPs<2>::AddPrim, return_self<>())
    .def("add",       &_EXPs<2>::AddPrims, return_self<>())
    .def("add_lc",    &_EXPs<2>::AddLC, return_self<>())
    .def("setup",     &_EXPs<2>::SetUp, return_self<>())
    .def("str",   &_EXPs<2>::str)
    .def("at_r",  &_EXPs<2>::AtR)
    .def("clone", &_EXPs<2>::Clone)
    .def("conj",  &_EXPs<2>::Conj)
    .def("calc_rm_mat", &_EXPs<2>::CalcRmMat)
    .def("calc_d2_mat", &_EXPs<2>::CalcD2Mat)
    .def("calc_vec",    &_EXPs<2>::CalcVecSTO)
    .def("calc_vec",    &_EXPs<2>::CalcVecGTO);

  /*
  register_ptr_to_python<GTOs>();
  def("_GTOs", &Create_GTOs);
  class_<_GTOs>("_GTOs", init<>())
    .def("size",      &_GTOs::size)
    .def("basis",     &_GTOs::basis)
    .def("only_prim", &_GTOs::OnlyPrim)
    .def("add",       &_GTOs::AddPrim, return_self<>())
    .def("add",       &_GTOs::AddPrims, return_self<>())
    .def("add_lc",    &_GTOs::AddLC, return_self<>())
    .def("setup",     &_GTOs::SetUp, return_self<>())
    .def("str",  &_GTOs::str)
    .def("at_r", &_GTOs::AtR)
    .def("clone", &_GTOs::Clone)
    .def("conj", &_GTOs::Conj)
    .def("calc_rm_mat", &_GTOs::CalcRmMat)
    .def("calc_d2_mat", &_GTOs::CalcD2Mat);
  */
}

BOOST_PYTHON_MODULE(r1basis_bind) {

  Py_Initialize();
  BindR1LC();
  BindR1Basis();
  def("add", &add);
  
}
