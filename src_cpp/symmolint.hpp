#ifndef SYMMOLINT_H
#define SYMMOLINT_H

#include <string>
#include <vector>
#include <Eigen/Core>
#include "typedef.hpp"
#include "bmatset.hpp"
#include "symgroup.hpp"
#include "b2eint.hpp"

namespace l2func {

  // ==== AO Reduction ====
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
    pSymmetryGroup sym_group;
    std::vector<dcomplex> x_iat;
    std::vector<dcomplex> y_iat;
    std::vector<dcomplex> z_iat;
    std::vector<int> nx_ipn;
    std::vector<int> ny_ipn;
    std::vector<int> nz_ipn;
    Eigen::VectorXcd zeta_iz;
    std::vector<Reduction> rds;

    // ---- for calculation  ----
    Eigen::MatrixXi ip_iat_ipn;
    Eigen::MatrixXi ip_jg_kp;       // GTO[j] = G[i][GTO[k]]
    Eigen::MatrixXi sign_ip_jg_kp; // sign of above relation
    
    bool setupq;
    int maxn;
    //    int maxnx;

    // ---- Constructors ----
    SubSymGTOs(pSymmetryGroup);

    // ---- Accessors ----
    std::string str() const;
    inline int nx(int ipn) const { return nx_ipn[ipn]; }
    inline int ny(int ipn) const { return ny_ipn[ipn]; }
    inline int nz(int ipn) const { return nz_ipn[ipn]; }
    inline dcomplex x(int iat) const { return x_iat[iat]; }
    inline dcomplex y(int iat) const { return y_iat[iat]; }
    inline dcomplex z(int iat) const { return z_iat[iat]; }
    inline dcomplex zeta(int iz) const { return zeta_iz[iz]; }
    inline cRdsIt begin_rds() const { return rds.begin(); }
    inline cRdsIt end_rds() const { return rds.end(); }
    inline RdsIt begin_rds() { return rds.begin(); }
    inline RdsIt end_rds() { return rds.end(); }
    void AddXyz(Eigen::Vector3cd xyz);
    void AddNs(Eigen::Vector3i ns);
    void AddZeta(const Eigen::VectorXcd& zs);
    void AddRds(const Reduction& rds);   
    inline int size_at() const { return x_iat.size();}
    inline int size_pn() const { return nx_ipn.size(); }
    inline int size_prim() const { return this->size_at() * this->size_pn(); }
    inline int size_zeta() const { return zeta_iz.rows(); }

    // ---- SetUp ----
    // -- calculate inner information and check values.
    void SetUp();

    // ---- old ----
    //    SubSymGTOs(Eigen::MatrixXcd xyz, Eigen::MatrixXi ns,
    //	       std::vector<Reduction> cs, Eigen::VectorXcd zs);
    int size_cont() const { return rds.size(); }
    //    const Eigen::MatrixXcd& get_xyz_iat() { return xyz_iat; }
    //    const Eigen::MatrixXi&  get_ns_ipn() const { return  ns_ipn; }
    const Reduction&  get_rds(int i) const { return  rds[i]; }
    const Eigen::VectorXcd& get_zeta_iz() const { return  zeta_iz; }    
    void Display() const;
  };

  // ---- Helper ----
  SubSymGTOs Sub_s(pSymmetryGroup sym, Irrep irrep,
		   Eigen::Vector3cd xyz, Eigen::VectorXcd zs);
  SubSymGTOs Sub_pz(pSymmetryGroup sym, Irrep irrep,
		    Eigen::Vector3cd xyz, Eigen::VectorXcd zs);
  SubSymGTOs Sub_TwoSGTO(pSymmetryGroup sym, Irrep irrep,
			 Eigen::Vector3cd xyz, Eigen::VectorXcd zs);
  SubSymGTOs Sub_mono(pSymmetryGroup sym, Irrep irrep,
		      Eigen::Vector3cd xyz, Eigen::Vector3i ns, Eigen::VectorXcd zs);

  // ==== SymGTOs ====
  class SymGTOs {
  public:
    pSymmetryGroup sym_group;
    std::vector<SubSymGTOs> subs;
    Eigen::MatrixXcd xyzq_iat;
    bool setupq;
  public:
    // ---- Constructors ----
    SymGTOs(pSymmetryGroup _sym_group);

    // ---- Accessors ----    
    int size_atom() const;
    int size_basis() const;
    int size_basis_isym(Irrep isym) const;
    std::string str() const;    
    inline dcomplex x_at(int i) const { return xyzq_iat(0, i); }
    inline dcomplex y_at(int i) const { return xyzq_iat(1, i); }
    inline dcomplex z_at(int i) const { return xyzq_iat(2, i); }
    inline dcomplex q_at(int i) const { return xyzq_iat(3, i); }
    int max_num_prim() const;

    // ---- Add information ----
    void SetAtoms(Eigen::MatrixXcd _xyzq_iat);
    void AddAtom(Eigen::Vector3cd _xyz, dcomplex q);
    void AddSub(SubSymGTOs);
    void SetComplexConj(SymGTOs&);

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
    void CalcERI(IB2EInt* eri, int method=0);
    // -- Radial wave function --
    void AtR_Ylm(int L, int M,  int irrep,
		 const Eigen::VectorXcd& cs_ibasis,
		 const Eigen::VectorXcd& rs,
		 Eigen::VectorXcd* vs,
		 Eigen::VectorXcd* dvs);
    // -- Correction of wave function sign --
    void CorrectSign(int L, int M, int irrep, Eigen::VectorXcd& cs);
  };
  
  void CalcERI(SymGTOs& gi, SymGTOs& gj,SymGTOs& gk,SymGTOs& gl, IB2EInt* eri);
}

#endif
