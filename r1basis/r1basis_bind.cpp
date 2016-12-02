#include <iostream>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <Eigen/Core>

#include "../math/erfc.hpp"
#include "../math/int_exp.hpp"
#include "../utils/eigen_plus.hpp"
#include "../src_cpp/angmoment.hpp"
#include "r1_lc.hpp"
#include "r1basis.hpp"
#include "opt_green.hpp"

using namespace Eigen;
using namespace std;
using namespace cbasis;
using namespace erfc_mori;
using namespace boost::python;

void BindMath() {
  def("tdot", &TDot);
  def("erfc", &erfc<dcomplex>);
  def("erfcx", &erfcx<dcomplex>);
  def("cg_coef", &cg_coef);
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
    .def("set_c", &_LC_EXPs<1>::set_c)
    .def("set_n", &_LC_EXPs<1>::set_n)
    .def("set_z", &_LC_EXPs<1>::set_z)        
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
    .def("set_c", &_LC_EXPs<2>::set_c)
    .def("set_n", &_LC_EXPs<2>::set_n)
    .def("set_z", &_LC_EXPs<2>::set_z)            
    .def("add",  &_LC_EXPs<2>::Add, return_self<>())
    .def("str",  &_LC_EXPs<2>::str)
    .def("at_r", &_LC_EXPs<2>::AtR)
    .def("clone", &_LC_EXPs<2>::Clone)
    .def("conj", &_LC_EXPs<2>::Conj)
    .def("__str__", &_LC_EXPs<2>::str);

}
template<int m>
void BindEXPs() {

  //  string name = (m==1 ? "_STOs" : "_GTOs");
  //const char *_name = name.c_str();
  const char *_name = (m==1 ? "_STOs" : "_GTOs");
  class_<_EXPs<m> >(_name, init<>())
    .def("size",  &_EXPs<m>::size)
    .def("basis", &_EXPs<m>::basis)
    .def("is_prim",       &_EXPs<m>::IsPrim)
    .def("is_prim_all",   &_EXPs<m>::IsPrimAll)
    .def("is_normal",     &_EXPs<m>::IsNormal)
    .def("is_normal_all", &_EXPs<m>::IsNormalAll)
    .def("has_coef",      &_EXPs<m>::HasCoef)
    .def("has_coef_all",  &_EXPs<m>::HasCoefAll)
    .def("exp_power", &_EXPs<m>::exp_power)
    .def("add",       &_EXPs<m>::AddPrim, return_self<>())
    .def("add",       &_EXPs<m>::AddPrims, return_self<>())
    .def("add",       &_EXPs<m>::AddLC, return_self<>())
    .def("add_not_normal", &_EXPs<m>::AddNotNormalPrim, return_self<>())
    .def("add_not_normal", &_EXPs<m>::AddNotNormalLC, return_self<>())
    .def("replace",     &_EXPs<m>::ReplaceLC, return_self<>())
    .def("setup",       &_EXPs<m>::SetUp, return_self<>())
    .def("str",         &_EXPs<m>::str)
    .def("at_r",        &_EXPs<m>::AtR)
    .def("at_r",        &_EXPs<m>::AtR_One)
    .def("clone",       &_EXPs<m>::Clone)
    .def("conj",        &_EXPs<m>::Conj)
    .def("calc_rm_mat", (MatrixXcd(_EXPs<m>::*)(int)const)&_EXPs<m>::CalcRmMat)
    .def("calc_d2_mat", (MatrixXcd(_EXPs<m>::*)()const)&_EXPs<m>::CalcD2Mat)
    .def("calc_vec",    (VectorXcd(_EXPs<m>::*)(LC_STOs)const)&_EXPs<m>::CalcVecSTO)
    .def("calc_vec",    (VectorXcd(_EXPs<m>::*)(LC_GTOs)const)&_EXPs<m>::CalcVecGTO)
    .def("__str__",     &_EXPs<m>::str);
    /*
    .def("calc_rm_mat", (void(_EXPs<m>::*)(int, MatrixXcd&))&_EXPs<m>::CalcRmMat)
    .def("calc_d2_mat", (void(_EXPs<m>::*)(MatrixXcd&))&_EXPs<m>::CalcD2Mat)
    .def("calc_vec",    (void(_EXPs<m>::*)(LC_STOs,VectorXcd&))&_EXPs<m>::CalcVec)
    .def("calc_vec",    (void(_EXPs<m>::*)(LC_GTOs,VectorXcd&))&_EXPs<m>::CalcVec)
    
    */
}
void BindR1Basis() {

  def("sto_int", &STOInt_Rplus);
  def("gto_int", &GTOInt_Rplus);
  def("sto_gto_int", &STO_GTOInt_Rplus);

  def("calc_vec", (VectorXcd(*)(STOs, LC_STOs))&CalcVec<1,1>);
  def("calc_vec", (VectorXcd(*)(STOs, LC_GTOs))&CalcVec<1,2>);
  def("calc_vec", (VectorXcd(*)(GTOs, LC_STOs))&CalcVec<2,1>);
  def("calc_vec", (VectorXcd(*)(GTOs, LC_GTOs))&CalcVec<2,2>);
  
  def("calc_rm_mat", (MatrixXcd(*)(STOs, int, STOs))&CalcRmMat<1,1>);
  def("calc_rm_mat", (MatrixXcd(*)(STOs, int, GTOs))&CalcRmMat<1,2>);
  def("calc_rm_mat", (MatrixXcd(*)(GTOs, int, STOs))&CalcRmMat<2,1>);
  def("calc_rm_mat", (MatrixXcd(*)(GTOs, int, GTOs))&CalcRmMat<2,2>);
  
  def("calc_d2_mat", (MatrixXcd(*)(STOs, STOs))&CalcD2Mat<1,1>);
  def("calc_d2_mat", (MatrixXcd(*)(STOs, GTOs))&CalcD2Mat<1,2>);
  def("calc_d2_mat", (MatrixXcd(*)(GTOs, STOs))&CalcD2Mat<2,1>);
  def("calc_d2_mat", (MatrixXcd(*)(GTOs, GTOs))&CalcD2Mat<2,2>);

  typedef _EXPs<1> _STOs;
  register_ptr_to_python<STOs>();
  def("STOs", &Create_STOs);
  BindEXPs<1>();

  typedef _EXPs<2> _GTOs;
  register_ptr_to_python<GTOs>();
  def("GTOs", &Create_GTOs);
  BindEXPs<2>();

}
void BindOptGreen() {
  class_<OptGreen<1,1,1> >("OptGreen_SSS", init<LC_STOs,STOs,LC_STOs,int,double>())
    .def("L00", &OptGreen<1,1,1>::L00, return_internal_reference<>())
    .def("calc_S0_L00_R0", &OptGreen<1,1,1>::Calc_S0_L00_R0);
}
BOOST_PYTHON_MODULE(r1basis_bind) {

  Py_Initialize();
  BindMath();
  BindR1LC();
  BindR1Basis();
  BindOptGreen();
}
