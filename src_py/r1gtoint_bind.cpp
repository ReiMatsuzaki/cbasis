#include <iostream>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <Eigen/Core>

#include "../src_cpp/r1gtoint.hpp"
#include "../src_cpp/opt_alpha.hpp"

using namespace Eigen;
using namespace std;
using namespace l2func;
using namespace boost::python;

void PrintMatrixXi(const MatrixXi& x) {
  cout << x<< endl;
}
void PrintVectorXi(const VectorXi& x) {
  cout << x<< endl;
}

tuple R1GTOs_AtR(R1GTOs* gtos, const VectorXcd& cs, const VectorXcd& rs) {
  VectorXcd* res = gtos->AtR(cs, rs);
  return make_tuple(res);
}
void print_c(const VectorXcd& xs){
  cout << xs << endl;
}
void print_i(const VectorXi& xs) {
  cout << xs << endl;
}
tuple OptAlphaShift_py(const R1STOs driv,
		       const VectorXi opt_idx_i,
		       dcomplex ene,
		       double h,
		       int max_iter,
		       double eps,
		       double eps_canonical,
		       R1GTOs& gtos,
		       dcomplex z_shift) {

  bool convq;
  dcomplex alpha;
  dcomplex z_shift_res(z_shift);
  OptAlphaShift(driv, opt_idx_i, ene, h, max_iter, eps, eps_canonical,
		gtos, z_shift_res,
		&convq, &alpha);
  return make_tuple(convq, alpha, z_shift_res);
}

BOOST_PYTHON_MODULE(r1gtoint_bind) {

  Py_Initialize();

  def("print_matrixxi", PrintMatrixXi);
  def("print_vectorxi", PrintVectorXi);
  
  class_<R1STO>("R1STO", init<dcomplex, int, dcomplex>())
    .def_readwrite("c", &R1STO::c)
    .def_readwrite("n", &R1STO::n)
    .def_readwrite("z", &R1STO::z)
    .def(self_ns::str(self))
    .def(self_ns::repr(self));
  class_<R1GTO>("R1GTO", init<dcomplex, int, dcomplex>())
    .def_readwrite("c", &R1GTO::c)
    .def_readwrite("n", &R1GTO::n)
    .def_readwrite("z", &R1GTO::z)
    .def(self_ns::str(self))
    .def(self_ns::repr(self));

  class_<vector<R1STO> >("vector_R1STO")
    .def("__getitem__", (R1STO const&(vector<R1STO>::*)(vector<R1STO>::size_type) const)&vector<R1STO>::at, return_value_policy<copy_const_reference>())
    .def("append",  &vector<R1STO>::push_back)
    .def("__len__", &vector<R1STO>::size);

  class_<VecMap>("VecMap")
    .def(map_indexing_suite<VecMap>());
  class_<MatMap>("MatMap")
    .def(map_indexing_suite<MatMap>());

  class_<R1GTOs>("R1GTOs", init<int>())
    .def("add", (void(R1GTOs::*)(dcomplex,int,dcomplex))&R1GTOs::Add)
    .def("add", (void(R1GTOs::*)(int,dcomplex))&R1GTOs::Add)
    .def("add", (void(R1GTOs::*)(int,const VectorXcd&))&R1GTOs::Add)
    .def("normalize", &R1GTOs::Normalize)
    .def("basis",(R1GTO&(R1GTOs::*)(int))&R1GTOs::basis,
	 return_internal_reference<>())
    .def("mat", &R1GTOs::mat, return_internal_reference<>())
    .def("vec", &R1GTOs::vec, return_internal_reference<>())
    .def("calc_mat", &R1GTOs::CalcMat)
    .def("calc_vec", (void(R1GTOs::*)(const R1GTOs&))&R1GTOs::CalcVec)
    .def("calc_vec", (void(R1GTOs::*)(const R1STOs&))&R1GTOs::CalcVec)
    .def(self_ns::str(self))
    .def(self_ns::repr(self))
    .def("at_r_cpp",
       	 (VectorXcd*(R1GTOs::*)
	  (const VectorXcd&, const VectorXcd&))&R1GTOs::AtR,
	 return_value_policy<manage_new_object>());


  class_<R1STOs>("R1STOs")
    .def("add",  (void(R1STOs::*)(dcomplex,int,dcomplex))&R1STOs::Add)
    .def("add",  (void(R1STOs::*)(int,dcomplex))&R1STOs::Add)
    .def("basis",  &R1STOs::basis,
	 return_internal_reference<>());

  def("calc_alpha", (dcomplex (*)(const R1STOs&, R1GTOs&, dcomplex, double))CalcAlpha);
  def("opt_alpha_shift", OptAlphaShift_py);
  def("solve_alpha", SolveAlpha);
  
  //def("print_c", print_c);
  //def("print_i", print_i);

}
