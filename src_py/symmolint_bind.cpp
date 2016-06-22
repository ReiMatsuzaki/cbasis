#include <iostream>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/numpy.hpp>
#include <Eigen/Core>

#include "../src_cpp/macros.hpp"

#include "../src_cpp/eigen_plus.hpp"
#include "../src_cpp/b2eint.hpp"
#include "../src_cpp/bmatset.hpp"
#include "../src_cpp/symmolint.hpp"
#include "../src_cpp/mo.hpp"

namespace {
  using namespace boost::python;
  using namespace l2func;
  using namespace std;
  using namespace Eigen;
  typedef std::complex<double> F;

  namespace bp = boost::python;
  namespace np = boost::numpy;
}

// ==== Boost/Python utilities ====
// http://nonbiri-tereka.hatenablog.com/entry/2014/02/02/104717
template<class T>
boost::python::list PyListFromCppVector(vector<T>& cpp_vector) {
  boost::python::list py_list;
  typename vector<T>::const_iterator it;
  for(it = cpp_vector.begin(); it != cpp_vector.end(); ++it) {
    py_list.append(*it);
  }
  return py_list;
}

MatrixXcd& BMat_get(BMat* bmat, tuple args) {
  int i = extract<int>(args[0]);
  int j = extract<int>(args[1]);
  return (*bmat)[make_pair(i, j)];
}
BMatSet* SymGTOs_CalcMat(SymGTOs* gtos) {
  BMatSet* res = new BMatSet(gtos->sym_group->order());
  gtos->CalcMat(res);
  return res;
}
BMatSet* SymGTOs_CalcMatOther(SymGTOs* a, SymGTOs* b, bool calc_coulombq) {
  BMatSet* res = new BMatSet(a->sym_group->order());
  a->CalcMatOther(*b, calc_coulombq, res);
  return res;
}
IB2EInt* SymGTOs_CalcERI(SymGTOs* gtos, int method) {
  int num(gtos->size_basis());
  IB2EInt* eri = new B2EIntMem(pow(num, 4));
  gtos->CalcERI(eri, method);
  return eri;
}
IB2EInt* py_CalcERI(SymGTOs* gi, SymGTOs* gj, SymGTOs* gk, SymGTOs* gl) {
  IB2EInt *eri = new B2EIntMem();
  CalcERI(*gi, *gj, *gk, *gl, eri);
  return eri;
}
BMat* py_CalcSEHamiltonian(MO mo, IB2EInt* eri, Irrep I0, int i0) {

  BMat *res = new BMat();
  CalcSEHamiltonian(mo, eri, I0, i0, res);
  return res;

}
tuple SymGTOs_AtR_Ylm(SymGTOs* gtos,
		      int L, int M,  int irrep,
		      const VectorXcd& cs_ibasis,
		      VectorXcd rs) {
  
  VectorXcd* ys = new VectorXcd();
  VectorXcd* ds = new VectorXcd();
  gtos->AtR_Ylm(L, M, irrep, cs_ibasis, rs, ys, ds);
  return make_tuple(ys, ds);
  
}
VectorXcd* SymGTOs_CorrectSign(SymGTOs* gtos, int L, int M,  int irrep,
			       const VectorXcd& cs_ibasis) {
  VectorXcd* cs = new VectorXcd(cs_ibasis); // call copy constructor
  gtos->CorrectSign(L, M, irrep, *cs);
  return cs;
}
const MatrixXcd& BMatSets_getitem(BMatSet* self, tuple name_i_j) {
  string name = extract<string>(name_i_j[0]);
  Irrep i = extract<int>(name_i_j[1]);
  Irrep j = extract<int>(name_i_j[2]);
  return self->GetMatrix(name, i, j);  
}
tuple generalizedComplexEigenSolve_py(const CM& F, const CM& S) {

  MatrixXcd C;
  VectorXcd eig;
  generalizedComplexEigenSolve(F, S, &C, &eig);
  return make_tuple(eig, C);

  /*
  int n(F.cols());
  if((F.rows() != n) ||
     (S.cols() != n) ||
     (S.rows() != n)) {
    string msg; SUB_LOCATION(msg); 
    msg += ": matrix F and S must be same size square matrix. ";
    throw runtime_error(msg);
  }
  */
}
MatrixXcd* CanonicalMatrix_py(const MatrixXcd& S, double eps) {

  MatrixXcd* X = new MatrixXcd();
  CanonicalMatrix(S, eps, X);
  return X;
}

MO CalcRHF1(SymGTOs& gtos, int nele, int max_iter, double eps, int debug_lvl) {
  bool is_conv;
  MO mo = CalcRHF(gtos, nele, max_iter, eps, &is_conv, debug_lvl);
  return mo;
}
MO CalcRHF2(pSymmetryGroup sym, BMatSet& mat_set, IB2EInt* eri, int nele, 
	    int max_iter, double eps, int debug_lvl) {
  bool is_conv;
  MO mo = CalcRHF(sym, mat_set, eri, nele, max_iter, eps, &is_conv, debug_lvl);
  return mo;
}
boost::python::list MO_get_num_occ(MO mo) {
  boost::python::list num_occ = PyListFromCppVector(mo->num_occ_irrep);
  return num_occ;
}

BOOST_PYTHON_MODULE(symmolint_bind) {

  Py_Initialize();
  np::initialize();
  
  def("canonical_matrix", CanonicalMatrix_py,
      return_value_policy<manage_new_object>());
  def("ceig", generalizedComplexEigenSolve_py);

  register_ptr_to_python<pSymmetryGroup>();
  class_<SymmetryGroup>("SymmetryGroup", no_init)
    .add_property("order", &SymmetryGroup::order)
    .add_property("name", &SymmetryGroup::name)
    .def("get_irrep", &SymmetryGroup::GetIrrep)
    .def("get_irrep_name", &SymmetryGroup::GetIrrepName)
    .def("__str__", &SymmetryGroup::str)
    .def("__repr__", &SymmetryGroup::str);
  def("C1", SymmetryGroup::C1);
  def("Cs", SymmetryGroup::Cs);
  def("C2h", SymmetryGroup::C2h);
  def("D2h", SymmetryGroup::D2h);
  def("C4", SymmetryGroup::C4);

  class_<vector<int> >("vector_i")
    .def(vector_indexing_suite<vector<int> >());

  class_<BVec>("BVec", no_init)
    .def(map_indexing_suite<BVec>());

  class_<BMat>("BMat", no_init)
    .def("__getitem__", BMat_get, return_internal_reference<>());
  
  class_<BMatSet>("BMatSet", init<int>())
    .def("get_matrix",  &BMatSet::GetMatrix, return_internal_reference<>())
    .def("__getitem__", BMatSets_getitem, return_internal_reference<>())
    .def("set_matrix", &BMatSet::SetMatrix);

  class_<IB2EInt, boost::noncopyable>("IB2EInt", no_init);
  class_<B2EIntMem, bases<IB2EInt> >("B2EIntMem", init<>())
    .def("get", &IB2EInt::Get)
    .def("set", &IB2EInt::Set)
    .def("reset", &IB2EInt::Reset)
    .def("at", &IB2EInt::At)
    .def("exist", &IB2EInt::Exist)
    .def("read",  &IB2EInt::Read)
    .def("write", &IB2EInt::Write)
    .def("size", &IB2EInt::size)
    .def("capacity", &IB2EInt::capacity);
  
  class_<Reduction>("Reduction", init<int, MatrixXcd>())
    .add_property("irrep",        &Reduction::irrep)
    .def("coef_iat_ipn", &Reduction::get_coef_iat_ipn,
	 return_internal_reference<>())
    .def("__str__", &Reduction::str)
    .def("__repr__", &Reduction::str);

  class_<SubSymGTOs>("SubSymGTOs", init<pSymmetryGroup>())
    .def("xyz", &SubSymGTOs::AddXyz, return_self<>())
    .def("ns",  &SubSymGTOs::AddNs,  return_self<>())
    .def("rds", &SubSymGTOs::AddRds,  return_self<>())
    .def("zeta",&SubSymGTOs::AddZeta, return_self<>())
    .def("get_zeta", &SubSymGTOs::zeta)
    .def("x", &SubSymGTOs::x)
    .def("y", &SubSymGTOs::y)
    .def("z", &SubSymGTOs::z)
    .def("nx", &SubSymGTOs::nx)
    .def("ny", &SubSymGTOs::ny)
    .def("nz", &SubSymGTOs::nz)
    .def("__str__", &SubSymGTOs::str)
    .def("__repr__", &SubSymGTOs::str);
  def("sub_s", Sub_s);
  def("sub_pz", Sub_pz);
  def("sub_two_sgto", Sub_TwoSGTO);

  class_<SymGTOs>("SymGTOs", init<pSymmetryGroup>())
    .def("sub", &SymGTOs::AddSub, return_self<>())
    .def("atom", &SymGTOs::AddAtom, return_self<>())
    .def("setup", &SymGTOs::SetUp, return_self<>())
    .def("set_cc", &SymGTOs::SetComplexConj)
    .def("calc_mat", SymGTOs_CalcMat,
	 return_value_policy<manage_new_object>())    
    .def("calc_mat", SymGTOs_CalcMatOther,
	 return_value_policy<manage_new_object>())
    .def("calc_eri", SymGTOs_CalcERI,
	 return_value_policy<manage_new_object>())
    .def("at_r_ylm_cpp", SymGTOs_AtR_Ylm)
    .def("correct_sign", SymGTOs_CorrectSign,
	 return_value_policy<manage_new_object>())
    //    .def("at_r_ylm", SymGTOs_AtR_Ylm)    
    .def("__str__", &SymGTOs::str)
    .def("__repr__", &SymGTOs::str);
  def("calc_ERI", py_CalcERI, return_value_policy<manage_new_object>());

  register_ptr_to_python<MO>();
  class_<_MO>("MO", no_init)
    .def_readonly("H", &_MO::H)
    .def_readonly("S", &_MO::S)
    .def_readonly("C", &_MO::C)
    .def_readonly("eigs", &_MO::eigs)
    .def("num_occ", MO_get_num_occ)
    .def_readonly("energy", &_MO::energy);
  
  def("calc_RHF", CalcRHF1);
  def("calc_RHF", CalcRHF2);
  def("calc_SEHamiltonian", py_CalcSEHamiltonian,
      return_value_policy<manage_new_object>());
  def("calc_alpha", CalcAlpha);
  def("pi_total_crosssection", PITotalCrossSection);
}

