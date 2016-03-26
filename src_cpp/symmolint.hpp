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
  struct Reduction {
    Irrep irrep;
    Eigen::MatrixXcd coef_iat_ipn;

    // ---- for calculation  ----
    Eigen::VectorXcd coef_iz;
    int offset;

    // ReductionSets() {}
    Reduction(int _irrep, Eigen::MatrixXcd _coef):
      irrep(_irrep), coef_iat_ipn(_coef), offset(0) {}
    std::string str() const;
    void Display() const;
    
    const Eigen::MatrixXcd& get_coef_iat_ipn() const { return coef_iat_ipn; }
    inline int size_at() const { return coef_iat_ipn.rows(); }
    inline int size_pn() const { return coef_iat_ipn.cols(); }
    void set_zs_size(int num_zs) {
      coef_iz = Eigen::VectorXcd::Zero(num_zs);
    }
  };

  // ==== Sub sets of SymGTOs ====
  struct SubSymGTOs {
    // ---- type ----
    typedef std::vector<Reduction>::const_iterator cRdsIt;
    typedef std::vector<Reduction>::iterator RdsIt;

    // ---- Calculation data ----
    Eigen::MatrixXcd xyz_iat;
    Eigen::MatrixXi  ns_ipn;
    Eigen::VectorXcd zeta_iz;
    std::vector<Reduction> rds;

    // ---- for calculation  ----
    bool setupq;
    int maxn;

    // ---- Constructors ----
    SubSymGTOs();

    // ---- Accessors ----
    std::string str() const;
    inline int nx(int ipn) const { return ns_ipn(0, ipn); }
    inline int ny(int ipn) const { return ns_ipn(1, ipn); }
    inline int nz(int ipn) const { return ns_ipn(2, ipn); }
    inline dcomplex x(int iat) const { return xyz_iat(0, iat); }
    inline dcomplex y(int iat) const { return xyz_iat(1, iat); }
    inline dcomplex z(int iat) const { return xyz_iat(2, iat); }
    inline dcomplex zeta(int iz) const { return zeta_iz[iz]; }
    inline cRdsIt begin_rds() const { return rds.begin(); }
    inline cRdsIt end_rds() const { return rds.end(); }
    inline RdsIt begin_rds() { return rds.begin(); }
    inline RdsIt end_rds() { return rds.end(); }
    void AddXyz(Eigen::Vector3cd xyz);
    void AddNs(Eigen::Vector3i ns);
    void AddZeta(const Eigen::VectorXcd& zs);
    void AddRds(const Reduction& rds);
    inline int size_at() const { return xyz_iat.cols();}
    inline int size_pn() const { return ns_ipn.cols(); }
    inline int size_zeta() const { return zeta_iz.rows(); }

    // ---- SetUp ----
    // -- calculate inner information and check values.
    void SetUp();
    
    
    // ---- old ----
    SubSymGTOs(Eigen::MatrixXcd xyz, Eigen::MatrixXi ns,
	       std::vector<Reduction> cs, Eigen::VectorXcd zs);
    int size_cont() const { return rds.size(); }
    const Eigen::MatrixXcd& get_xyz_iat() { return xyz_iat; }
    const Eigen::MatrixXi&  get_ns_ipn() const { return  ns_ipn; }
    const Reduction&  get_rds(int i) const { return  rds[i]; }
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
    inline dcomplex x_at(int i) const { return xyzq_iat(0, i); }
    inline dcomplex y_at(int i) const { return xyzq_iat(1, i); }
    inline dcomplex z_at(int i) const { return xyzq_iat(2, i); }
    inline dcomplex q_at(int i) const { return xyzq_iat(3, i); }

    // ---- Add information ----
    void SetAtoms(Eigen::MatrixXcd _xyzq_iat);
    void AddAtom(Eigen::Vector3cd _xyz, dcomplex q);
    void AddSub(SubSymGTOs);

    // ---- SetUp ----
    void SetUp();
  private:
    void SetOffset();
    void Normalize();

  public:
    // ---- Calculation ----
    // -- to be removed --
    void loop();
    // -- not uesd now --
    int max_n() const;
    // -- matrix calculation --
    void CalcMatOther(SymGTOs& o, bool calc_coulomb, BMatSet*);
    void CalcMat(BMatSet* res);
    // -- Radial wave function --
    void AtR_Ylm(int L, int M,  int irrep,
		 const Eigen::VectorXcd& cs_ibasis,
		 const Eigen::VectorXcd& rs,
		 Eigen::VectorXcd* vs,
		 Eigen::VectorXcd* dvs);
    // -- Correction of wave function sign --
    void CorrectSign(int L, int M, int irrep, Eigen::VectorXcd& cs);

  };
}

#endif
