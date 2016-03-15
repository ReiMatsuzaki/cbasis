#include <iostream>
#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include <Eigen/Core>
#include "../eigen_plus.hpp"
#include "../symmolint.hpp"

namespace {
  using namespace boost::python;
  using namespace l2func;
  using namespace std;
  using namespace Eigen;
  typedef std::complex<double> F;

  namespace bp = boost::python;
  namespace np = boost::numpy;
}

// -- to be removed
void VectorEigen2NumPy(VectorXcd& a, np::ndarray* b) {

  int num = a.size();
  dcomplex* ptr;
  Map<VectorXcd>(ptr, num) = a;
  *b = np::from_data(ptr,
		     np::dtype::get_builtin<dcomplex>(),
		     bp::make_tuple(num),
		     bp::make_tuple(sizeof(dcomplex)),
		     bp::object());
}
// -- to be removed
void VectorNumPy2Eigen(np::ndarray a, VectorXcd* b) {

  if(a.get_nd() != 1) {
    string msg; SUB_LOCATION(msg);
    msg += "dimension of a must be 1";
    throw runtime_error(msg);
  }

  if(a.get_dtype() != np::dtype::get_builtin<dcomplex>()) {
    string msg; SUB_LOCATION(msg);
    msg += "only complex type is supported";
    throw runtime_error(msg);
  }

  int num = a.shape(0);

  dcomplex* ptr = reinterpret_cast<dcomplex*>(a.get_data());
  *b = Map<VectorXcd>(ptr, num);

}
// -- to be removed
void MatrixEigen2NumPy(MatrixXcd& a, np::ndarray* b) {
  int numi = a.rows();
  int numj = a.cols();

  dcomplex* ptr;
  Map<MatrixXcd>(ptr, numi, numj) = a;

  *b = np::from_data(ptr,
		     np::dtype::get_builtin<dcomplex>(),
		     bp::make_tuple(numi, numj),
		     bp::make_tuple(sizeof(dcomplex), sizeof(dcomplex)*numi),
		     bp::object());
  
}
// -- to be removed
void MatrixNumPy2Eigen(np::ndarray& a, MatrixXcd* b) {

  if(a.get_nd() != 2){
    string msg; SUB_LOCATION(msg);
    msg += "dimension of a must be 2";
    throw runtime_error(msg);
  }

  if(a.get_dtype() != np::dtype::get_builtin<dcomplex>()) {
    string msg; SUB_LOCATION(msg);
    msg += "only complex type is supported";
    throw runtime_error(msg);
  }
  
  int numi = a.shape(0);
  int numj = a.shape(1);
  
  dcomplex* ptr = reinterpret_cast<dcomplex*>(a.get_data());
  typedef Matrix<dcomplex, Dynamic, Dynamic, RowMajor> Mat;
  *b = Map<Mat>(ptr, numi, numj);

}
// -- to be removed
void VectorNumPy2Eigen_3(np::ndarray a, Vector3cd* b) {

  if(a.shape(0) != 3) {
    string msg; SUB_LOCATION(msg);
    msg += ": a must have 3 elements";
    throw runtime_error(msg);    
  }  

  VectorXcd c;
  VectorNumPy2Eigen(a, &c);
  *b = Vector3cd(c(0), c(1), c(2));

}
// -- to be removed
void PrintAsEigen(np::ndarray& a) {

  if(a.get_nd() == 1) {

    VectorXcd b; 
    VectorNumPy2Eigen(a, &b);
    cout << b << endl;

  } else if(a.get_nd() == 2) {

    MatrixXcd b;
    MatrixNumPy2Eigen(a, &b);
    cout << b << endl;
      
  } else {
    string msg; SUB_LOCATION(msg);
    msg += "only 1 or 2 dimension of a is supported";
    throw runtime_error(msg);

  }

}
// -- to be removed
ReductionSets ReductionSetsInit(int irrep, np::ndarray& cs_iat_ipn) {

  MatrixXcd coef;
  try {
    MatrixNumPy2Eigen(cs_iat_ipn, &coef);
  } catch(const runtime_error& e) {
    string msg; SUB_LOCATION(msg);
    msg += "Error on conversion: \n";
    msg += e.what();
    throw runtime_error(msg);
  }

  ReductionSets rds(irrep, coef);
  return rds;

}
// -- to be removed
SubSymGTOs Sub_s_py(Irrep irrep, np::ndarray& xyz, np::ndarray& zs) {

  Vector3cd xyz_in; VectorNumPy2Eigen_3(xyz, &xyz_in);
  VectorXcd zs_in; VectorNumPy2Eigen(zs, &zs_in);

  return Sub_s(irrep, xyz_in, zs_in);
}
// -- to be removed
SubSymGTOs Sub_pz_py(Irrep irrep, np::ndarray& xyz, np::ndarray& zs) {

  Vector3cd xyz_in; VectorNumPy2Eigen_3(xyz, &xyz_in);
  VectorXcd zs_in; VectorNumPy2Eigen(zs, &zs_in);

  return Sub_pz(irrep, xyz_in, zs_in);

}
// -- to be removed
SubSymGTOs Sub_TwoSGTO_py(SymmetryGroup sym, Irrep irrep,
			  np::ndarray xyz, np::ndarray zs) {
  Vector3cd xyz_in; VectorNumPy2Eigen_3(xyz, &xyz_in);
  VectorXcd zs_in;  VectorNumPy2Eigen(zs, &zs_in);

  return Sub_TwoSGTO(sym, irrep, xyz_in, zs_in);
}
// -- to be removed
np::dict SymGTOs_CalcMat(SymGTOs* gtos) {

  BMatMap* res0;
  gtos->CalcMat(res0);

  bp::dict res1;
  for(BMatMap::iterator it = res1.begin(); it != res1.end(); ++it) {
    
    string name = it->first;
  }
  return res1;

}

BOOST_PYTHON_MODULE(symmolint_bind) {

  Py_Initialize();
  np::initialize();
  
  def("print_as_eigen", PrintAsEigen);

  class_<SymmetryGroup>("SymmetryGroup", init<int, string>())
    .add_property("order", &SymmetryGroup::order)
    .add_property("name", &SymmetryGroup::name)
    .def("__str__", &SymmetryGroup::str)
    .def("__repr__", &SymmetryGroup::str);
  def("C1", SymmetryGroup_C1);
  def("Cs", SymmetryGroup_Cs);
  def("Cs_Ap", Cs_Ap);
  def("Cs_App", Cs_App);

  class_<ReductionSets>("ReductionSets", init<int, MatrixXcd>())
    .add_property("irrep",        &ReductionSets::irrep)
    .add_property("coef_iat_ipn", &ReductionSets::get_coef_iat_ipn)
    .def("__str__", &ReductionSets::str)
    .def("__repr__", &ReductionSets::str);

  class_<SubSymGTOs>("SubSymGTOs", init<MatrixXcd, MatrixXi,
		     vector<ReductionSets>, VectorXcd>())
    .add_property("xyz_iat", &SubSymGTOs::get_xyz_iat)
    .add_property("ns_ipn", &SubSymGTOs::get_ns_ipn)
    .add_property("rds", &SubSymGTOs::get_rds)
    .add_property("zeta_iz", &SubSymGTOs::get_zeta_iz)
    .def("__str__", &SubSymGTOs::str)
    .def("__repr__", &SubSymGTOs::str);
  def("sub_s", Sub_s);
  def("sub_pz", Sub_pz);
  def("sub_two_sgto", Sub_TwoSGTO);

  class_<SymGTOs>("SymGTOs", init<SymmetryGroup>())
    .def("add_sub", &SymGTOs::AddSub)
    .def("add_atoms", &SymGTOs::AddAtoms)
    .def("setup", &SymGTOs::SetUp)
    .def("__str__", &SymGTOs::str)
    .def("__repr__", &SymGTOs::str)
  
}

