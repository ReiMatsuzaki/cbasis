#include <boost/python.hpp>
#include "fact.hpp"
#include "prim.hpp"
#include "lcomb.hpp"
#include "op.hpp"
#include "hatom.hpp"

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
LinearComb<Prim> NormalizedBasis(int n, typename Prim::Field z) {
  LinearComb<Prim> res(Prim(n, z, Normalized));
  return res;
}

template<class P1, class P2>
F sym_ip(const P1& a, const P2 b) {
  return CIP(a, b);
}

template<class Prim>
LinearComb<Prim> DerivNormalized(int n, typename Prim::Field z) {

  Prim f(n, z, Normalized);
  return D1Normalized(f);

}
template<class Prim>
LinearComb<Prim> DoubleDerivNormalized
(int n, typename Prim::Field z) {
  Prim f(n, z, Normalized);
  return D2Normalized(f);
}


BOOST_PYTHON_MODULE(l2func_bind) {

  def("factorial", fact::Factorial);

  class_<STO>("STO", init<F, int, F>())
    .def("c", &STO::c)
    .def("n", &STO::n)
    .def("z", &STO::z)
    .def("set_z", &STO::set_z)
    .def("at_x", &STO::at);
  

  class_<GTO>("GTO", init<F, int, F>())
    .def("c", &GTO::c)
    .def("n", &GTO::n)
    .def("z", &GTO::z)
    .def("set_z", &GTO::set_z)
    .def("at_x", &GTO::at);

  class_<STOs>("STOs", init<>())
    .def("size", &STOs::size)
    .def("coef_i", &STOs::coef_i)
    .def("add", &STOs::AddOther)
    .def("at_x", &STOs::at)
    .def("add_one", &STOs::AddCoefPrim)
    .def("prim_i", &STOs::prim_i_copied);
  

  class_<GTOs>("GTOs", init<>())
    .def("size", &GTOs::size)
    .def("coef_i", &GTOs::coef_i)
    .def("add", &GTOs::AddOther)
    .def("at_x", &GTOs::at)
    .def("add_one", &GTOs::AddCoefPrim)
    .def("prim_i", &GTOs::prim_i_copied);

  def("scalar_prod_sto", ScalarProductForLC<STO>);
  def("scalar_prod_gto", ScalarProductForLC<GTO>);
  def("normalized_sto", NormalizedBasis<STO>);
  def("normalized_gto", NormalizedBasis<GTO>);

  def("sym_ip_ss", sym_ip<STOs, STOs>);
  def("sym_ip_sg", sym_ip<STOs, GTOs>);
  def("sym_ip_gs", sym_ip<GTOs, STOs>);
  def("sym_ip_gg", sym_ip<GTOs, GTOs>);

  def("d_normalized_sto", DerivNormalized<STO>);
  def("dd_normalized_sto", DoubleDerivNormalized<STO>);
  def("d_normalized_gto", DerivNormalized<GTO>);
  def("dd_normalized_gto", DoubleDerivNormalized<GTO>);

  class_<Op<STO> >("Op_sto", init<>())
    .def("add", &Op<STO>::AddOther)
    .def("operate", &Op<STO>::OperateLC);

  class_<Op<GTO> >("Op_gto", init<>())
    .def("add", &Op<GTO>::AddCoefOp)
    .def("operate", &Op<GTO>::OperateLC);

  def("scalar_prod_op_sto", ScalarProductOp<STO>);
  def("scalar_prod_op_gto", ScalarProductOp<GTO>);

  def("op_rm_sto", OpRM<STO>);
  def("op_cst_sto", OpCst<STO>);
  def("op_ddr_sto", OpDDr<STO>);
  def("op_ddr2_sto", OpDDr2<STO>);

  def("op_rm_gto", OpRM<GTO>);
  def("op_cst_gto", OpCst<GTO>);
  def("op_ddr_gto", OpDDr<GTO>);
  def("op_ddr2_gto", OpDDr2<GTO>);

  class_<HLikeAtom<F> >("HLikeAtom", init<int, F, int>())
    .def("hamiltonian_sto", &HLikeAtom<F>::Hamiltonian<STO>)
    .def("hamiltonian_gto", &HLikeAtom<F>::Hamiltonian<GTO>)
    .def("h_minus_energy_sto", &HLikeAtom<F>::HMinusEnergy<STO>)
    .def("h_minus_energy_gto", &HLikeAtom<F>::HMinusEnergy<GTO>)
    .def("n", &HLikeAtom<F>::n)
    .def("z", &HLikeAtom<F>::z)
    .def("l", &HLikeAtom<F>::l)
    .def("eigen_state", &HLikeAtom<F>::EigenState)
    .def("dipole_init_length", &HLikeAtom<F>::DipoleInitLength)
    .def("dipole_init_velocity", &HLikeAtom<F>::DipoleInitVelocity)
    .def("eigen_energy", &HLikeAtom<F>::EigenEnergy);
  
}
