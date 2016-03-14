#ifndef SYMMOLINT_H
#define SYMMOLINT_H

#include "math_utils.hpp"
#include <Eigen/Core>
#include <string>
#include <map>

namespace l2func {

  // ==== Each integrate ====

  // ==== Type def ====
  typedef int Irrep;
  typedef std::map<std::pair<Irrep, Irrep>, Eigen::MatrixXcd> BMat;
  typedef std::map<std::string, BMat> BMatMap;

  void swap(BMat& a, BMat& b);
  void swap(BMatMap& a, BMatMap& b);

  // ==== Irrep Group ====
  // ---- Class ----
  class SymmetryGroup {
  private:
    int order_;
  public:
    SymmetryGroup(int order);
    int order()  const { return order_; }
    void CheckIrrep(Irrep a) const;
    // Irrep GetIrrep(int n) const;
    bool IncludeScalar_2(Irrep a, Irrep b) const;
    // bool IncludeScalar_3(const Irrep& a, const Irrep& b, const Irrep& c) const;
    bool IncludeZ_2(Irrep a, Irrep b) const;
  };
  SymmetryGroup SymmetryGroup_Cs();
  SymmetryGroup SymmetryGroup_C1();
   
  // ==== AO Reduction Sets ====
  struct ReductionSets {
    Irrep sym;
    Eigen::MatrixXcd coef_iat_ipn;

    // ---- for calculation  ----
    Eigen::VectorXcd coef_iz;
    int offset;

    ReductionSets(int _sym, Eigen::MatrixXcd _coef):
      sym(_sym), coef_iat_ipn(_coef), offset(0) {}
    void set_zs_size(int num_zs) {
      coef_iz = Eigen::VectorXcd::Zero(num_zs);
    }
  };

  // ==== Sub sets of SymGTOs ====
  struct SubSymGTOs {
    // ---- data ----
    Eigen::MatrixXcd xyz_iat;
    Eigen::MatrixXi  ns_ipn;
    std::vector<ReductionSets> rds;
    Eigen::VectorXcd zeta_iz;

    // ---- for calculation  ----
    int maxn;

    int size_at() const { return xyz_iat.cols();}
    int size_pn() const { return ns_ipn.cols(); }
    int size_cont() const { return rds.size(); }
    int size_zeta() const { return zeta_iz.rows(); }
    SubSymGTOs(Eigen::MatrixXcd xyz, Eigen::MatrixXi ns,
	       std::vector<ReductionSets> cs, Eigen::VectorXcd zs);
  };

  // ---- Helper ----
  SubSymGTOs Sub_s(Irrep sym, Eigen::Vector3cd xyz, Eigen::VectorXcd zs);
  SubSymGTOs Sub_pz(Irrep sym, Eigen::Vector3cd xyz, Eigen::VectorXcd zs);

  // ==== SymGTOs ====
  class SymGTOs {
  public:
    SymmetryGroup sym_group;
    std::vector<SubSymGTOs> subs;
    Eigen::MatrixXcd xyzq_iat;
    bool setupq;
  public:
    // ---- Constructors ----
    SymGTOs(SymmetryGroup _sym_group);

    // ---- Accessors ----
    int size_atom() const;
    int size_basis_isym(Irrep isym) const;

    // ---- Add information ----
    void SetAtoms(Eigen::MatrixXcd _xyzq_iat);
    void AddSub(SubSymGTOs);

    void SetUp();
  private:
    void SetOffset();
    void Normalize();
  public:
    // ---- Calculation ----
    void loop();
    void CalcMat(BMatMap* res);
    void STVMat(BMatMap* res);
    void ZMat(BMatMap* res);

    void AtR_Ylm_add_center(int L, int M, const Eigen::VectorXcd& rs,
			    const Eigen::MatrixXcd& cs_irrep_ibasis,  
			    Eigen::VectorXcd* vs,
			    std::vector<SubSymGTOs>::const_iterator it);
    void AtR_Ylm(int L, int M, const Eigen::VectorXcd& rs,
		 const Eigen::MatrixXcd& cs_irrep_ibasis, Eigen::VectorXcd* vs );
  };
}

#endif
