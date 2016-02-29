#include <boost/python.hpp>
//#include <boost/python/numeric.hpp>
//#include <numpy/ndarray.h>
//#include <Python.h>
#include <boost/numpy.hpp>
#include "../linspace.hpp"
#include "../exp_func.hpp"
#include "../cut_exp.hpp"
#include "../delta.hpp"
#include "../cip_exp.hpp"
#include "../op.hpp"

#include "../gto3dset.hpp"

namespace {
  using namespace boost::python;
  using namespace l2func;
  typedef std::complex<double> F;

  namespace bp = boost::python;
  namespace np = boost::numpy;
  using namespace std;

}

using namespace std;
using namespace boost;
using namespace boost::python;

/*
np::ndarray array_transform(dcmplx* xs, int n) {  
  return np::from_data(xs,
		       np::dtype::get_builtin<dcmplx>(),
		       bp::make_tuple(n),
		       bp::make_tuple(sizeof(dcmplx)),
		       bp::object());

}
		       */
/*
double trace(bp::numeric::array& y, int n) {

  double cumsum(0);
  for(int i = 0; i < n; i++) {
    cumsum += bp::extract<double>(y[i]);
  }
  return cumsum;
}
*/

void np_array_set3(int num, np::ndarray& ys) {

  for(int i = 0; i < num; i++)
    ys[bp::make_tuple(i)] = dcmplx(3.0, 2.0);

}

np::ndarray CalcSMat(const SphericalGTOSet& a, const SphericalGTOSet& b) {
  // dcmplx* vs = a.SMatWithOther(b);
  // return array_transform(vs, a.size() * b.size());

  /*
  dcmplx* xs = new dcmplx[5];
  for(int i = 0; i < 5; i++)
    xs[i] = 1.1*i;
  return array_transform(xs, 5);
  */
/*
  double* xs = new double[5];
  for(int i = 0; i < 5; i++)
    xs[i] = 1.1*i;
  np::ndarray res = np::from_data(xs,
				  np::dtype::get_builtin<double>(),
				  bp::make_tuple(5),
				  bp::make_tuple(sizeof(double)),
				  bp::object());
*/
  dcmplx* vs = a.SMatWithOther(b);
  np::ndarray res = np::from_data(vs,
				  np::dtype::get_builtin<dcmplx>(),
				  bp::make_tuple(a.size() * b.size()),
				  bp::make_tuple(sizeof(dcmplx)),
				  bp::object());
  return res;
}
/*
np::ndarray CalcTMat(const SphericalGTOSet& a, const SphericalGTOSet& b) {
  dcmplx* vs = a.TMatWithOther(b);
  return array_transform(vs, a.size() * b.size());
}
np::ndarray CalcVMat(const SphericalGTOSet& a, const SphericalGTOSet& b,
		     dcmplx q, dcmplx x, dcmplx y, dcmplx z) {
  dcmplx* vs = a.VMatWithOther(b, q, x, y, z);
  return array_transform(vs, a.size() * b.size());
}
*/

BOOST_PYTHON_MODULE(l2func_bind) {

  Py_Initialize();
  np::initialize();
  //  import_array();
  // bp::numeric::array::set_module_and_type("numpy", "ndarray");
  // def("trace", trace);

  def("np_array_set3", np_array_set3);

  class_<array3<F> >("tuple_c3", init<F, F, F>())
    .add_property("x", &array3<F>::x, &array3<F>::set_x)
    .add_property("y", &array3<F>::y, &array3<F>::set_y)
    .add_property("z", &array3<F>::z, &array3<F>::set_z);

  class_<array3<int> >("tuple_i3", init<int, int, int>())
    .add_property("x", &array3<int>::x, &array3<int>::set_x)
    .add_property("y", &array3<int>::y, &array3<int>::set_y)
    .add_property("z", &array3<int>::z, &array3<int>::set_z);
  
  typedef ExpFunc<F,1> STO;
  class_<STO>("STO", init<F, int, F>())
    .add_property("c", &STO::c, &STO::set_c)
    .add_property("n", &STO::n, &STO::set_n)
    .add_property("z", &STO::z, &STO::set_z)
    .def("__str__", &STO::str)
    .def("at", &STO::at);

  typedef ExpFunc<F,2> GTO;
  class_<GTO>("GTO", init<F, int, F>())
    .add_property("c", &GTO::c, &GTO::set_c)
    .add_property("n", &GTO::n, &GTO::set_n)
    .add_property("z", &GTO::z, &GTO::set_z)
    .def("__str__", &GTO::str)
    .def("at", &GTO::at);

  typedef CutExpFunc<F, 1> CutSTO;
  class_<CutSTO>("CutSTO", init<F, int, F, double>())
    .add_property("c", &CutSTO::c, &CutSTO::set_c)
    .add_property("n", &CutSTO::n, &CutSTO::set_n)
    .add_property("z", &CutSTO::z, &CutSTO::set_z)
    .add_property("r0", &CutSTO::r0, &CutSTO::set_r0)
    .def("__str__", &CutSTO::str)
    .def("at", &CutSTO::at);

  /*
  typedef CutExpFunc<F, 2> CutGTO;
  class_<CutGTO>("CutGTO", init<F, int, F, double>())
    .add_property("c", &CutGTO::c, &CutGTO::set_c)
    .add_property("n", &CutGTO::n, &CutGTO::set_n)
    .add_property("z", &CutGTO::z, &CutGTO::set_z)
    .add_property("r0", &CutGTO::r0, &CutGTO::set_r0)
    .def("__str__", &CutGTO::str)
    .def("at", &CutGTO::at);
  */

  def("cip_ss", CIP<STO, STO>);
  def("cip_sg", CIP<STO, GTO>);
  def("cip_gs", CIP<GTO, STO>);
  def("cip_gg", CIP<GTO, GTO>);
  def("cip_cut_ss", CIP<CutSTO, CutSTO>);
  //  def("cip_cut_gg", CIP<CutGTO, CutGTO>);
  

  class_<CRm>("Rm", init<int>())
    .add_property("m", &CRm::m);
  class_<CD1>("D1", init<>());
  class_<CD2>("D2", init<>());
  class_<Cut<F, double> >("Cut", init<double>())
    .add_property("r0", &Cut<F, double>::r0);

  def("cip_s_d2_s", CIP<STO, CD2, STO>);
  def("cip_s_d2_g", CIP<STO, CD2, GTO>);
  def("cip_g_d2_s", CIP<GTO, CD2, STO>);
  def("cip_g_d2_g", CIP<GTO, CD2, GTO>);
  def("cip_cut_s_d2_s", CIP<CutSTO, CD2, CutSTO>);

  def("cip_s_d1_s", CIP<STO, CD1, STO>);
  def("cip_s_d1_g", CIP<STO, CD1, GTO>);
  def("cip_g_d1_s", CIP<GTO, CD1, STO>);
  def("cip_g_d1_g", CIP<GTO, CD1, GTO>);
  def("cip_cut_s_d1_s", CIP<CutSTO, CD1, CutSTO>);

  def("cip_s_rm_s", CIP<STO, CRm, STO>);
  def("cip_s_rm_g", CIP<STO, CRm, GTO>);
  def("cip_g_rm_s", CIP<GTO, CRm, STO>);
  def("cip_g_rm_g", CIP<GTO, CRm, GTO>);
  def("cip_cut_s_rm_s", CIP<CutSTO, CRm, CutSTO>);


  // ==== 3D ====
  class_<SphericalGTOSet>("SphericalGTOSet", init<>())
    .def("add_one_basis", &SphericalGTOSet::AddOneBasis)
    .def("add_basis",     &SphericalGTOSet::AddBasis);
  // def("overlap", overlap);
  //  def("kinetic", kinetic);
  // def("overlap_1D", overlap_1D);
  // def("nuclear_attraction", nuclear_attraction0);

  def("calc_s_mat", CalcSMat);
  //def("calc_t_mat", CalcTMat);
  //def("calc_v_mat", CalcVMat);

  /*
  typedef SphericalGTO<F,F> SphericalGTO3d;
  class_<SphericalGTO3d>("SphericalGTO", init<int, int, c3, F>())
    .add_property("l", &SphericalGTO3d::L)
    .add_property("m", &SphericalGTO3d::M)
    .add_property("zeta", &SphericalGTO3d::zeta)
    .add_property("xyz", &SphericalGTO3d::xyz);
    */

  // class_<OpKE<F, c3> >("KE", init<int>());
  // class_<OpNA<F, c3> >("NA", init<F, c3>())
  //   .add_property("q",   &OpNA<F, c3>::q)
  //   .add_property("xyz", &OpNA<F, c3>::xyz);
  //  class_<OpXyz<F, c3> >("Xyz", init<int, int, int>());
    //    .add_property("n", &OpNA<F, c3>::n)
    //    .add_property("m", &OpNA<F, c3>::m)
    //    .add_property("l", &OpNA<F, c3>::l);  

}

