 #include <iostream>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <Eigen/Core>

#include "../src_cpp/r1gtoint.hpp"
#include "../src_cpp/op_driv.hpp"
#include "../src_cpp/opt_alpha.hpp"

using namespace Eigen;
using namespace std;
using namespace cbasis;
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
/*
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
*/
tuple R1GTOs_MatSTV_1(R1GTOs* gtos, int L) {
  MatrixXcd s,t,v;
  gtos->CalcMatSTV(L, s, t, v);
  return make_tuple(s, t, v);
}
tuple R1GTOs_MatSTV_2(R1GTOs* gi, R1GTOs* gj, int L) {
  MatrixXcd s,t,v;
  gi->CalcMatSTV(*gj, L, s, t, v);
  return make_tuple(s, t, v);
}
MatrixXcd* R1GTOs_MatSTO_1(R1GTOs* gs, const R1STOs& sto) {
  MatrixXcd* mat = new MatrixXcd();
  gs->CalcMatSTO(sto, *mat);
  return mat;
}
MatrixXcd* R1GTOs_MatSTO_2(R1GTOs* gi, R1GTOs* gj, const R1STOs& sto) {
  MatrixXcd* mat = new MatrixXcd();
  gi->CalcMatSTO(*gj, sto, *mat);
  return mat;
}
VectorXcd* R1GTOs_VecSTO(R1GTOs* gs, const R1STOs& pot) {
  VectorXcd* vec = new VectorXcd();
  gs->CalcVec(pot, *vec);
  return vec;
}
VectorXcd* SolveAlpha_Py(IDriv* driv, IOp* op, const R1GTOs& gs) {
  VectorXcd* ptr = new VectorXcd();
  SolveAlpha(driv, op, gs, *ptr);
  return ptr;
}
void BindBasic() {
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

}
void BindOpt() {

  // ==== Driven term ====
  class_<IDriv, boost::noncopyable>("IDriv", no_init);
  class_<DrivSTO, bases<IDriv> >("DrivSTO", init<const R1STOs&>());


  // ==== Operator ====
  class_<IOp, boost::noncopyable>("IOp", no_init);
  class_<OpCoulomb, bases<IOp> >("OpCoulomb", init<int, dcomplex>());
  class_<OpCoulombShort, bases<IOp> >("OpCoulombShort", init<int, dcomplex, const R1STOs&>());


  // ==== Optimization Result ====
  class_<OptResult>("OptRes", no_init)
    .def_readwrite("conv_q", &OptResult::conv_q)
    .def_readwrite("zs",     &OptResult::zs)
    .def_readwrite("val",    &OptResult::val)
    .def_readwrite("dz",     &OptResult::dz)
    .def_readwrite("grad",   &OptResult::grad)
    .def_readwrite("hess",   &OptResult::hess)
    .def(self_ns::str(self))
    .def(self_ns::repr(self));
    
  // ==== Optimizer ====
  class_<IOptimizer, boost::noncopyable>("IOptimizer", no_init)
    .def("optimize",
	 (OptResult*(IOptimizer::*)(dcomplex))&IOptimizer::Optimize,
	 return_value_policy<manage_new_object>())
    .def("optimize",
	 (OptResult*(IOptimizer::*)(const VectorXcd&))&IOptimizer::Optimize,
	 return_value_policy<manage_new_object>());
  
  class_<OptNewton, bases<IOptimizer> >("OptNewton",
					init<int, double, IOptTarget*, int>());

  // ==== Opt target ====
  class_<IOptTarget, boost::noncopyable>("IOptTarget", no_init);
  class_<OptAlpha, bases<IOptTarget> >("OptAlpha", init<IDriv*,IOp*,R1GTOs&>());
  class_<OptAlphaShift, bases<IOptTarget> >("OptAlphaShift", init<IDriv*,IOp*,R1GTOs&,
 			const VectorXi&>());
  class_<OptAlphaPartial, bases<IOptTarget> >("OptAlphaPartial", init<IDriv*, IOp*, R1GTOs&, const VectorXi&>());

  def("calc_alpha", CalcAlpha);
  def("solve_alpha", SolveAlpha);
  def("solve_alpha", SolveAlpha_Py, return_value_policy<manage_new_object>());
  

}
BOOST_PYTHON_MODULE(r1gtoint_bind) {

  Py_Initialize();

  BindBasic();
  BindOpt();
 
  class_<R1GTOs>("R1GTOs")
    .def("add", (void(R1GTOs::*)(int,dcomplex))&R1GTOs::Add)
    .def("add", (void(R1GTOs::*)(int,const VectorXcd&))&R1GTOs::Add)
    .def("add", (void(R1GTOs::*)(int,const VectorXcd&, const MatrixXcd&))&R1GTOs::Add)
    .def("set", (void(R1GTOs::*)(const VectorXcd&))&R1GTOs::Set)
    .def("set", (void(R1GTOs::*)(int, const VectorXcd&))&R1GTOs::Set)
    .def("set_conj", &R1GTOs::SetConj)
    .def("set_one_deriv", &R1GTOs::SetOneDeriv)
    .def("set_two_deriv", &R1GTOs::SetTwoDeriv)
    .def("z_prim", &R1GTOs::z_prim)
    .def("n_prim", &R1GTOs::n_prim)
    .def("calc_mat_stv", R1GTOs_MatSTV_1)
    .def("calc_mat_stv", R1GTOs_MatSTV_2)
    .def("calc_mat_stv", (void(R1GTOs::*)(int, MatrixXcd&, MatrixXcd&, MatrixXcd&) const)&R1GTOs::CalcMatSTV)
    .def("calc_mat_stv", (void(R1GTOs::*)(const R1GTOs&, int, MatrixXcd&, MatrixXcd&, MatrixXcd&) const)&R1GTOs::CalcMatSTV)
    .def("calc_mat_sto", R1GTOs_MatSTO_1, return_value_policy<manage_new_object>())
    .def("calc_mat_sto", R1GTOs_MatSTO_2, return_value_policy<manage_new_object>())
    .def("calc_mat_sto", (void(R1GTOs::*)(const R1STOs&, MatrixXcd&) const)&R1GTOs::CalcMatSTO)
    .def("calc_mat_sto", (void(R1GTOs::*)(const R1GTOs&, const R1STOs&, MatrixXcd&) const)&R1GTOs::CalcMatSTO)
    .def("calc_vec_sto", R1GTOs_VecSTO,   return_value_policy<manage_new_object>())
    .def("calc_vec_sto", (void(R1GTOs::*)(const R1STOs&, VectorXcd&) const)&R1GTOs::CalcVec)
    .def("setup", &R1GTOs::SetUp)
    .def("normalize", &R1GTOs::Normalize)
    .def("at_r_cpp",
       	 (VectorXcd*(R1GTOs::*)
	  (const VectorXcd&, const VectorXcd&))&R1GTOs::AtR,
	 return_value_policy<manage_new_object>())
    .def("deriv_at_r_cpp",
       	 (VectorXcd*(R1GTOs::*)
	  (const VectorXcd&, const VectorXcd&) const)&R1GTOs::DerivAtR,
	 return_value_policy<manage_new_object>())
    .def("deriv_2_at_r_cpp",
	 (VectorXcd*(R1GTOs::*)
	  (const VectorXcd&, const VectorXcd&) const)&R1GTOs::Deriv2AtR,
	 return_value_policy<manage_new_object>())
    .def(self_ns::str(self))
    .def(self_ns::repr(self));

  class_<R1STOs>("R1STOs")
    .def("add",  (void(R1STOs::*)(dcomplex,int,dcomplex))&R1STOs::Add)
    .def("add",  (void(R1STOs::*)(int,dcomplex))&R1STOs::Add)
    .def("basis",  &R1STOs::basis,
	 return_internal_reference<>())
    .def("at_r_cpp", (VectorXcd*(R1STOs::*)(const VectorXcd&) const)&R1STOs::AtR,
	 return_value_policy<manage_new_object>())
    .def(self_ns::str(self))
    .def(self_ns::repr(self));

  /*
  def("calc_alpha", (dcomplex (*)(const R1STOs&, R1GTOs&, dcomplex, double))CalcAlpha);
  def("calc_alpha", (dcomplex (*)(const R1STOs&, R1GTOs&, dcomplex, const R1STOs&))CalcAlpha);
  def("opt_alpha_shift", OptAlphaShift_py);
  def("solve_alpha", (VectorXcd (*)(const R1STOs&, R1GTOs&, dcomplex, double))SolveAlpha);
  def("solve_alpha", (VectorXcd (*)(const R1STOs&, R1GTOs&, dcomplex, const R1STOs&))SolveAlpha);
  */
  
  //def("print_c", print_c);
  //def("print_i", print_i);

}
