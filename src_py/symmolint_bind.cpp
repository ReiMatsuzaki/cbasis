#include <iostream>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
// #include <boost/numpy.hpp>
#include <Eigen/Core>

#include "../src_cpp/macros.hpp"
#include "../src_cpp/eigen_plus.hpp"
#include "../src_cpp/angmoment.hpp"
#include "../src_cpp/b2eint.hpp"
#include "../src_cpp/bmatset.hpp"
#include "../src_cpp/one_int.hpp"
#include "../src_cpp/two_int.hpp"
#include "../src_cpp/symmolint.hpp"
#include "../src_cpp/mo.hpp"


namespace {
  using namespace boost::python;
  using namespace cbasis;
  using namespace std;
  using namespace Eigen;
  typedef std::complex<double> F;

  namespace bp = boost::python;
  // namespace np = boost::numpy;
}

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
void BMat_set(BMat* bmat, int i, int j, MatrixXcd& mat) {
  (*bmat)[make_pair(i, j)] = mat;
}
void BMatSetZero(BMat& C, BMat& mat) {

  for(int I = 0; I < (int)C.size(); I++) {
    pair<Irrep, Irrep> II(I, I);
    MatrixXcd& CII = C[II];
    mat[II] = MatrixXcd::Zero(CII.rows(), CII.cols());
  }
  
}
void BMatSetZero(SymGTOs g0, SymGTOs g1, BMat& bmat) {

  //
  // g0 and g1 are assumed to have same symmetry
  //

  if(! g0->sym_group->IsSame(g1->sym_group)) {
    string msg; SUB_LOCATION(msg);
    msg += ": symmetry group of g0 and g1 are different";
    throw runtime_error(msg);
  }
  
  for(int irrep = 0; irrep < g0->sym_group->order(); irrep++) {
    
    int n0(g0->size_basis_isym(irrep));
    int n1(g1->size_basis_isym(irrep));
    bmat[make_pair(irrep, irrep)] = MatrixXcd::Zero(n0, n1);
  }

}

void AddAngmoemnt() {

  def("cg_coef", cg_coef);

}

boost::python::list MO_get_num_occ(MO mo) {
  boost::python::list num_occ = PyListFromCppVector(mo->num_occ_irrep);
  return num_occ;
}
MO CalcRHF1(SymGTOs gtos, int nele, int max_iter, double eps, int debug_lvl) {
  bool is_conv;
  MO mo = CalcRHF(gtos, nele, max_iter, eps, &is_conv, debug_lvl);
  return mo;
}
MO CalcRHF2(pSymmetryGroup sym, BMatSet mat_set, B2EInt eri, int nele, 
	    int max_iter, double eps, int debug_lvl) {
  bool is_conv;
  MO mo = CalcRHF(sym, mat_set, eri, nele, max_iter, eps, &is_conv, debug_lvl);
  return mo;
}
BMat* py_CalcSEHamiltonian(MO mo, B2EInt eri, Irrep I0, int i0) {

  BMat *res = new BMat();
  CalcSEHamiltonian(mo, eri, I0, i0, res);
  return res;

}
BMat* py_JK(B2EInt eri, BMat& C, Irrep I0, int i0, dcomplex c_J, dcomplex c_K) {

  BMat *mat = new BMat();
  BMatSetZero(C, *mat);
  
  AddJK(eri, C, I0, i0, c_J, c_K, *mat);

  return mat;
}
BMat* py_J(B2EInt eri, VectorXcd& ca, Irrep irrep_a, 
	   SymGTOs g0, SymGTOs g1, dcomplex coef_J) {

  BMat *mat = new BMat();
  BMatSetZero(g0, g1, *mat);
  
  AddJ(eri, ca, irrep_a, coef_J, *mat);

  return mat;

}
BMat* py_K(B2EInt eri, VectorXcd& ca, Irrep ir_a,
	   SymGTOs g0, SymGTOs g1, dcomplex coef_K) {

  BMat *mat = new BMat();
  BMatSetZero(g0, g1, *mat);
  
  AddK(eri, ca, ir_a, coef_K, *mat);

  return mat;

}

void AddMO() {

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

  def("calc_JK", py_JK, return_value_policy<manage_new_object>());

  def("calc_J", py_J, return_value_policy<manage_new_object>());
  def("calc_K", py_K, return_value_policy<manage_new_object>());

  def("calc_alpha", CalcAlpha);
  def("pi_total_crosssection", PITotalCrossSection);

}

void AddSymmetryGroup() {

  register_ptr_to_python<pSymmetryGroup>();
  class_<SymmetryGroup>("SymmetryGroup", no_init)
    .add_property("order", &SymmetryGroup::order)
    .add_property("name", &SymmetryGroup::name)
    .add_property("irrep_s", &SymmetryGroup::irrep_s)
    .add_property("irrep_x", &SymmetryGroup::irrep_x)
    .add_property("irrep_y", &SymmetryGroup::irrep_y)
    .add_property("irrep_z", &SymmetryGroup::irrep_z)
    .def("get_irrep", &SymmetryGroup::GetIrrep)
    .def("get_irrep_name", &SymmetryGroup::GetIrrepName)
    .def("__str__", &SymmetryGroup::str)
    .def("__repr__", &SymmetryGroup::str);
  def("C1", SymmetryGroup::C1);
  def("Cs", SymmetryGroup::Cs);
  def("C2h", SymmetryGroup::C2h);
  def("D2h", SymmetryGroup::D2h);
  def("C4", SymmetryGroup::C4);

}

MatrixXcd* CanonicalMatrix_py(const MatrixXcd& S, double eps) {

  MatrixXcd* X = new MatrixXcd();
  CanonicalMatrix(S, eps, X);
  return X;
}

tuple CEigenSolveCanonicalNum_py(const CM& F, const CM& S, int num0) {

  MatrixXcd C;
  VectorXcd eig;
  CEigenSolveCanonicalNum(F, S, num0, &C, &eig);
  return make_tuple(eig, C);

}

tuple generalizedComplexEigenSolve_py(const CM& F, const CM& S) {

  MatrixXcd C;
  VectorXcd eig;
  generalizedComplexEigenSolve(F, S, &C, &eig);
  return make_tuple(eig, C);
}

void AddLinearAlgebra() {

  def("canonical_matrix", CanonicalMatrix_py,
      return_value_policy<manage_new_object>());
  def("ceig", generalizedComplexEigenSolve_py);
  def("ceig_canonical", CEigenSolveCanonicalNum_py);

  class_<vector<int> >("vector_i")
    .def(vector_indexing_suite<vector<int> >());

  class_<BVec>("BVec", no_init)
    .def(map_indexing_suite<BVec>());

  class_<BMat>("BMat", init<>())
    .def("write", BMatWrite)
    .def("read", BMatRead)
    .def("set_matrix", BMat_set)
    .def("__getitem__", BMat_get, return_internal_reference<>());

  register_ptr_to_python<BMatSet>();
  class_<_BMatSet>("BMatSet", init<int>())
    .def("get_matrix",  &_BMatSet::GetMatrix, return_internal_reference<>())
    .def("__getitem__", &_BMatSet::GetBlockMatrix, return_internal_reference<>())
    .def("set_matrix", &_BMatSet::SetMatrix);

}

void AddB2EInt() {

  register_ptr_to_python<B2EInt>();
  class_<IB2EInt, boost::noncopyable>("IB2EInt", no_init);
  class_<B2EIntMem, bases<IB2EInt> >("B2EIntMem", init<>())
    .def("get", &IB2EInt::Get)
    .def("set", &IB2EInt::Set)
    .def("reset", &IB2EInt::Reset)
    .def("at", &IB2EInt::At)
    .def("exist", &IB2EInt::Exist)
    .def("write", &IB2EInt::Write)
    .def("size", &IB2EInt::size)
    .def("capacity", &IB2EInt::capacity);
  def("ERI_read", ERIRead);
  
}

tuple SymGTOs_AtR_Ylm(SymGTOs gtos,
		      int L, int M,  int irrep,
		      const VectorXcd& cs_ibasis,
		      VectorXcd rs) {
  
  VectorXcd* ys = new VectorXcd();
  VectorXcd* ds = new VectorXcd();
  gtos->AtR_Ylm(L, M, irrep, cs_ibasis, rs, ys, ds);
  return make_tuple(ys, ds);
  
}
VectorXcd* SymGTOs_CorrectSign(SymGTOs gtos, int L, int M,  int irrep,
			       const VectorXcd& cs_ibasis) {
  VectorXcd* cs = new VectorXcd(cs_ibasis); // call copy constructor
  gtos->CorrectSign(L, M, irrep, *cs);
  return cs;
}
void AddSymMolInt() {

  class_<Reduction>("Reduction", init<int, MatrixXcd>())
    .add_property("irrep",        &Reduction::irrep)
    .def("coef_iat_ipn", &Reduction::get_coef_iat_ipn,
	 return_internal_reference<>())
    .def("__str__", &Reduction::str)
    .def("__repr__", &Reduction::str);

  class_<SubSymGTOs>("SubSymGTOs", init<>())
    .def("sym", &SubSymGTOs::SetSym, return_self<>())
    .def("xyz", (void(SubSymGTOs::*)(dcomplex, dcomplex, dcomplex))
	 &SubSymGTOs::AddXyz, return_self<>())
    .def("xyz", (void(SubSymGTOs::*)(Vector3cd))&SubSymGTOs::AddXyz, return_self<>())
    .def("ns",  (void(SubSymGTOs::*)(Vector3i))
	   &SubSymGTOs::AddNs,  return_self<>())
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
  def("sub_solid_sh", Sub_SolidSH_M);
  def("sub_solid_sh", Sub_SolidSH_Ms);


  register_ptr_to_python<SymGTOs>();
  class_<_SymGTOs>("SymGTOs", init<>())
    .def("create", CreateSymGTOs).staticmethod("create")
    .add_property("sym_group", &_SymGTOs::sym_group)
    .def("sym", &_SymGTOs::SetSym, return_self<>())
    .def("atom", &_SymGTOs::AddAtom, return_self<>())
    .def("sub", &_SymGTOs::AddSub, return_self<>())
    .def("set_atoms", &_SymGTOs::SetAtoms, return_self<>())
    .def("clone", &_SymGTOs::Clone)
    .def("conj", &_SymGTOs::Conj)
    .def("setup", &_SymGTOs::SetUp, return_self<>())
    .def("at_r_ylm_cpp", SymGTOs_AtR_Ylm)
    .def("correct_sign", SymGTOs_CorrectSign,
	 return_value_policy<manage_new_object>())
    //    .def("at_r_ylm", _SymGTOs_AtR_Ylm)    
    .def("__str__", &_SymGTOs::str)
    .def("__repr__", &_SymGTOs::str);
}

void AddOneTwoInt() {

  class_<ERIMethod>("ERI_method", init<>())
    .def("use_symmetry",&ERIMethod::set_symmetry,
	 return_self<>())
    .def("coef_R_memo", &ERIMethod::set_coef_R_memo,
	 return_self<>());

  def("calc_mat_complex", CalcMat_Complex);
  def("calc_mat_hermite", CalcMat_Hermite);
  def("calc_mat", CalcMat);
  
  def("calc_ERI_complex", CalcERI_Complex);
  def("calc_ERI_hermite", CalcERI_Hermite);
  def("calc_ERI", CalcERI);

}


BOOST_PYTHON_MODULE(symmolint_bind) {

  Py_Initialize();
  // np::initialize();
  AddAngmoemnt();
  AddMO();
  AddSymmetryGroup();
  AddLinearAlgebra();
  AddB2EInt();
  AddSymMolInt();
  AddOneTwoInt();

}

