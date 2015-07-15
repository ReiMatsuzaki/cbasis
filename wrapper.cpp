#include <boost/python.hpp>
#include "fact.hpp"
#include "prim.hpp"
#include "lcomb.hpp"

/**
 *  This is wrapper for python bindings 
 */

namespace {
  using namespace boost::python;
  using namespace l2func;
  typedef std::complex<double> F;
  typedef ExpBasis<F, 1> STO;
  typedef ExpBasis<F, 2> GTO;
  typedef LinearComb<STO> STOs;
  typedef LinearComb<GTO> GTOs;
}

template<class Prim>
Prim NormalizedBasis(int n, typename Prim::Field z) {
  return Prim(n, z, Normalized);
}

template<class P1, class P2>
F sym_ip(const P1& a, const P2 b) {
  return CIP(a, b);
}


BOOST_PYTHON_MODULE(l2func_bind) {

  def("factorial", fact::Factorial);

  class_<STO>("STO", init<int, F>())
    .def("c", &STO::c)
    .def("n", &STO::n)
    .def("z", &STO::z)
    .def("set_z", &STO::set_z);
  def("normalized_sto", NormalizedBasis<STO>);

  class_<GTO>("GTO", init<int, F>())
    .def("c", &GTO::c)
    .def("n", &GTO::n)
    .def("z", &GTO::z)
    .def("set_z", &GTO::set_z);
  def("normalized_gto", NormalizedBasis<GTO>);

  def("sym_ip_ss", sym_ip<STO,STO>);
  def("sym_ip_sg", sym_ip<STO,GTO>);
  def("sym_ip_gs", sym_ip<GTO,STO>);
  def("sym_ip_gg", sym_ip<GTO,GTO>);

  class_<STOs>("STOs", init<>())
    .def("size", &STOs::size)
    .def("coef_i", &STOs::coef_i)
    .def("add", &STOs::Add)
    .def("prim_i", &STOs::prim_i_copied);

}
