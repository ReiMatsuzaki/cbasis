#include <boost/python.hpp>
#include "../exp_func.hpp"
#include "../cip.hpp"
#include "../cip_impl.hpp"
#include "../op.hpp"

namespace {
  using namespace boost::python;
  using namespace l2func;
  typedef std::complex<double> F;
}

BOOST_PYTHON_MODULE(l2func_bind) {
  
  typedef ExpFunc<F,1> STO;
  class_<STO>("STO", init<F, int, F>())
    .add_property("c", &STO::c, &STO::set_c)
    .add_property("n", &STO::n, &STO::set_n)
    .add_property("z", &STO::z, &STO::set_z)
    .def("at", &STO::at);

  typedef ExpFunc<F,2> GTO;
  class_<GTO>("GTO", init<F, int, F>())
    .add_property("c", &GTO::c, &GTO::set_c)
    .add_property("n", &GTO::n, &GTO::set_n)
    .add_property("z", &GTO::z, &GTO::set_z)
    .def("at", &GTO::at);

  def("cip_ss", CIP<STO, STO>);
  def("cip_sg", CIP<STO, GTO>);
  def("cip_gs", CIP<GTO, STO>);
  def("cip_gg", CIP<GTO, GTO>);

  class_<OpRm>("Rm", init<int>())
    .add_property("m", &OpRm::m);
  class_<OpD1>("D1", init<>());
  class_<OpD2>("D2", init<>());

  def("cip_s_d2_s", CIP<STO, OpD2, STO>);
  def("cip_s_d2_g", CIP<STO, OpD2, GTO>);
  def("cip_g_d2_s", CIP<GTO, OpD2, STO>);
  def("cip_g_d2_g", CIP<GTO, OpD2, GTO>);

  def("cip_s_d1_s", CIP<STO, OpD1, STO>);
  def("cip_s_d1_g", CIP<STO, OpD1, GTO>);
  def("cip_g_d1_s", CIP<GTO, OpD1, STO>);
  def("cip_g_d1_g", CIP<GTO, OpD1, GTO>);

  def("cip_s_rm_s", CIP<STO, OpRm, STO>);
  def("cip_s_rm_g", CIP<STO, OpRm, GTO>);
  def("cip_g_rm_s", CIP<GTO, OpRm, STO>);
  def("cip_g_rm_g", CIP<GTO, OpRm, GTO>);

}

