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
    .def("c", &_LC_EXPs<1>::c)
    .def("n", &_LC_EXPs<1>::n)
    .def("z", &_LC_EXPs<1>::z)
    .def("add",  &_LC_EXPs<1>::Add, return_self<>())
    .def("str",  &_LC_EXPs<1>::str)
    .def("at_r", &_LC_EXPs<1>::AtR)
    .def("conj", &_LC_EXPs<1>::Conj);
  
  typedef _LC_EXPs<2> _LC_GTOs;
  register_ptr_to_python<LC_GTOs>();
  def("LC_GTOs", &Create_LC_GTOs);
  class_<_LC_GTOs>("_LC_GTOs", init<>())
    .def("size", &_LC_EXPs<2>::size)
    .def("c", &_LC_EXPs<2>::c)
    .def("n", &_LC_EXPs<2>::n)
    .def("z", &_LC_EXPs<2>::z)
    .def("add",  &_LC_EXPs<2>::Add, return_self<>())
    .def("str",  &_LC_EXPs<2>::str)
    .def("at_r", &_LC_EXPs<2>::AtR)
    .def("conj", &_LC_EXPs<2>::Conj);


  /* 
     
     
    .def("d_at_r", &_LC_STOs::DAtR)
    .def("d2_at_r", &_LC_STOs::D2AtR)
  */

  /*
    
  class_<_LC_STOs>("_LC_STOs", init<>())
    .def("size", &_LC_STOs::size)
    .def("at_r", &_LC_STOs::AtR)
    .def("d_at_r", &_LC_STOs::DAtR)
    .def("d2_at_r", &_LC_STOs::D2AtR)
    .def("add",     &_LC_STOs::Add)
    .def("conj",     &_LC_STOs::Conj);

  register_ptr_to_python<LC_GTOs>();
  class_<_LC_GTOs>("_LC_GTOs", init<>())
    .def("size", &_LC_GTOs::size)
    .def("at_r", &_LC_GTOs::AtR)
    .def("d_at_r", &_LC_GTOs::DAtR)
    .def("d2_at_r", &_LC_GTOs::D2AtR)
    .def("add",     &_LC_GTOs::Add)
    .def("conj",     &_LC_GTOs::Conj);
  */
  
}

void BindR1Basis() {
  
  register_ptr_to_python<GTOs>();
  def("GTOs", &Create_GTOs);
  class_<_GTOs>("_GTOs", init<>())
    .def("size",      &_GTOs::size)
    .def("basis",     &_GTOs::basis)
    .def("only_prim", &_GTOs::OnlyPrim)
    .def("add",       &_GTOs::AddPrim, return_self<>())
    .def("add",       &_GTOs::AddPrims, return_self<>())
    .def("add_lc",    &_GTOs::AddLC, return_self<>())
    .def("setup",     &_GTOs::SetUp, return_self<>());
    //    .def("str",  &_GTOs::str)
    //    .def("at_r", &_GTOs::AtR)
    //    .def("conj", &_GTOs::Conj);
  

}

BOOST_PYTHON_MODULE(r1basis_bind) {

  Py_Initialize();
  BindR1LC();
  BindR1Basis();
  def("add", &add);
  
}
