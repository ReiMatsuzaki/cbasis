#ifndef MO_H
#define MO_H

#include <vector>
#include <Eigen/Core>
#include <boost/shared_ptr.hpp>
#include "typedef.hpp"
#include "symgroup.hpp"
#include "symmolint.hpp"

namespace l2func {

  class _MO {
  public:
    _MO();
    pSymmetryGroup sym;
    BMat H, S, J, K, F, C, P;
    BVec eigs;
    std::vector<int> num_occ_irrep;
    std::vector<Irrep> irrep_list;
    dcomplex energy;
  };
  typedef boost::shared_ptr<_MO> MO;

  // ---- RHF ----
  vector<int> CalcOccNum(const BVec& eigs, int num_sym, int num_orb);
  MO CalcOneEle(pSymmetryGroup sym, BMatSet& mat_set, int debug_lvl = 0);
  MO CalcRHF(SymGTOs& gtos, int nele, int max_iter, double eps, bool *is_conv,
	     int debug_lvl = 0);
  MO CalcRHF(pSymmetryGroup sym, BMatSet& mat_set, IB2EInt* eri, int nele, 
	     int max_iter, double eps, bool *is_conv, int debug_lvl=0);
  void CalcSEHamiltonian(MO mo, IB2EInt* eri, Irrep I0, int i0, BMat* hmat);
  dcomplex CalcAlpha(MO mo, BMatSet& mat_set, Irrep I0, int i0, BMat& h_stex, double w);
  double PITotalCrossSection(dcomplex alpha, double w, int num_occ_ele);  
}

#endif

