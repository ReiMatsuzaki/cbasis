#ifndef SYMMOLINT_H
#define SYMMOLINT_H

#include <string>
#include <vector>
#include <Eigen/Core>
#include "typedef.hpp"
#include "bmatset.hpp"

namespace l2func {

  // ==== Type def ====
  typedef int Irrep;
  
  // ==== Symmetry Group ====
  // ---- Class ----
  class SymmetryGroup {
  private:
    int order_;
    std::string name_;
  public:
    SymmetryGroup(int order, std::string name);
    int order()  const { return order_; }
    std::string name() const { return name_; }
    std::string str() const;
    void Display() const;
    void CheckIrrep(Irrep a) const;
    // Irrep GetIrrep(int n) const;
    bool IncludeScalar_2(Irrep a, Irrep b) const;
    // bool IncludeScalar_3(const Irrep& a, const Irrep& b, const Irrep& c) const;
    bool IncludeZ_2(Irrep a, Irrep b) const;
  };
  SymmetryGroup SymmetryGroup_Cs();
  Irrep Cs_Ap(); 
  Irrep Cs_App();
  SymmetryGroup SymmetryGroup_C1();

  // ==== AO Reduction Sets ====
  struct ReductionSets {
    Irrep sym;
    Eigen::MatrixXcd coef_iat_ipn;

    // ---- for calculation  ----
    Eigen::VectorXcd coef_iz;
    int offset;

    // ReductionSets() {}
    ReductionSets(int _sym, Eigen::MatrixXcd _coef):
      sym(_sym), coef_iat_ipn(_coef), offset(0) {}
    std::string str() const;
    void Display() const;
    Irrep irrep() const { return sym; }
    const Eigen::MatrixXcd& get_coef_iat_ipn() const { return coef_iat_ipn; }
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
    bool setupq;
    int maxn;

    // ---- new ----
    SubSymGTOs();
    void SetUp();
    std::string str() const;
    void AddXyz(Eigen::Vector3cd xyz);
    void AddNs(Eigen::Vector3i ns);
    void AddZeta(const Eigen::VectorXcd& zs);
    void AddRds(const ReductionSets& rds);
    
    // ---- old ----
    SubSymGTOs(Eigen::MatrixXcd xyz, Eigen::MatrixXi ns,
	       std::vector<ReductionSets> cs, Eigen::VectorXcd zs);
    int size_at() const { return xyz_iat.cols();}
    int size_pn() const { return ns_ipn.cols(); }
    int size_cont() const { return rds.size(); }
    int size_zeta() const { return zeta_iz.rows(); }
    const Eigen::MatrixXcd& get_xyz_iat() { return xyz_iat; }
    const Eigen::MatrixXi&  get_ns_ipn() const { return  ns_ipn; }
    const ReductionSets&  get_rds(int i) const { return  rds[i]; }
    const Eigen::VectorXcd& get_zeta_iz() const { return  zeta_iz; }    
    void Display() const;
  };

  // ---- Helper ----
  SubSymGTOs Sub_s(Irrep sym, Eigen::Vector3cd xyz, Eigen::VectorXcd zs);
  SubSymGTOs Sub_pz(Irrep sym, Eigen::Vector3cd xyz, Eigen::VectorXcd zs);
  SubSymGTOs Sub_TwoSGTO(SymmetryGroup sym, Irrep irrep,
			 Eigen::Vector3cd xyz, Eigen::VectorXcd zs);

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
    std::string str() const;

    // ---- Add information ----
    void SetAtoms(Eigen::MatrixXcd _xyzq_iat);
    void AddAtom(Eigen::Vector3cd _xyz, dcomplex q);
    void AddSub(SubSymGTOs);
    void SetUp();
  private:
    void SetOffset();
    void Normalize();
  public:
    // ---- Calculation ----
    // -- to be removed
    void loop();
    void CalcMat(BMatSet* res);
    // -- to be removed
    void STVMat(BMatMap* res);
    // -- to be removed
    void ZMat(BMatMap* res);
    void AtR_Ylm(int L, int M,  int irrep,
		 const Eigen::VectorXcd& cs_ibasis,
		 const Eigen::VectorXcd& rs, Eigen::VectorXcd* vs );

  };
}

#endif
