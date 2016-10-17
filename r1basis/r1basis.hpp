#ifndef R1GTOINT_H
#define R1GTOINT_H

#include <vector>
#include <map>
#include <Eigen/Core>

#include "../src_cpp/typedef.hpp"
#include "r1_lc.hpp"

// TODO
// CalcVecの引数をR1GTOsから、std::vector<R1GTO>に変更。
// それぞれのメッソドをAccessor等に分けて書く。


namespace cbasis {


  dcomplex GTOInt(int n, dcomplex a);
  dcomplex GTOIntLC(LC_GTOs a, int m, LC_GTOs b);

  class _GTOs;
  typedef boost::shared_ptr<_GTOs> GTOs;
  GTOs Create_GTOs();

  class _GTOs {
  private:
    // ---- Member field ----
    std::vector<LC_GTOs> basis_;
    bool setupq_;
    
  public:
    // ---- Constructors ----
    _GTOs();
    
    // ---- Getter ----
    int size()  const { return basis_.size(); }
    LC_GTOs basis(int i) const { return basis_[i]; }
    bool OnlyPrim() const;
    Eigen::VectorXcd AtR(const Eigen::VectorXcd&, 
			 const Eigen::VectorXcd&) const;
    Eigen::VectorXcd DAtR(const Eigen::VectorXcd&,
			  const Eigen::VectorXcd&) const;
    std::string str() const;
    
    // ---- Setter ----
    _GTOs* AddPrim(int n, dcomplex z);
    _GTOs* AddPrims(int n, Eigen::VectorXcd zs);
    _GTOs* AddLC(LC_GTOs lc);
    _GTOs* SetUp();
    
    // ---- Create ----
    GTOs Conj() const;
    GTOs Clone() const;
    
    // ---- Calculate ----
    Eigen::MatrixXcd CalcRmMat(int m) const;
    Eigen::MatrixXcd CalcD2Mat()      const;
    Eigen::VectorXcd CalcVecSTO(LC_STOs stos) const;
    Eigen::VectorXcd CalcVecGTO(LC_GTOs gtos) const;

  };
  

  dcomplex NPrimeGTO(dcomplex nterm, int n, dcomplex z);
  dcomplex NDoublePrimeGTO(dcomplex nterm, int n, dcomplex z);

  struct R1GTO {
    R1GTO(dcomplex _c, int _n, dcomplex _z);
    dcomplex c;
    int n;
    dcomplex z;
    bool operator=(const R1GTO& o) const { return this==&o; }
    bool operator!=(const R1GTO& o) const { return this!=&o; }
  };
  struct R1STO {
    R1STO(dcomplex _c, int _n, dcomplex _z);
    dcomplex c;
    int n;
    dcomplex z;  
    bool operator=( const R1STO& o) const { return this==&o; }
    bool operator!=(const R1STO& o) const { return this!=&o; }  
    bool operator=( const R1STO& o) { return this==&o; }
    bool operator!=(const R1STO& o) { return this!=&o; }  
  };

  std::ostream& operator<<(std::ostream& out, const R1GTO& basis);
  std::ostream& operator<<(std::ostream& out, const R1STO& basis);
  
  //  typedef std::map<std::string, Eigen::MatrixXcd> MatMap;
  typedef std::vector<R1GTO>::const_iterator CIt;
  //  typedef std::map<std::string, Eigen::VectorXcd> VecMap;
  class R1GTOs;
  class R1STOs;  

  class MatVecMap {
  private:
    std::map<std::string, Eigen::MatrixXcd> mat_;
    std::map<std::string, Eigen::VectorXcd> vec_;
  public:
    const Eigen::MatrixXcd& mat(std::string lbl) const;
    Eigen::MatrixXcd& mat(std::string lbl);
    bool exist_mat(std::string lbl) const { return mat_.find(lbl) != mat_.end(); }
    void InitMatIfNecessary(std::string lbl, int ni, int nj);
    const Eigen::VectorXcd& vec(std::string lbl) const;
    Eigen::VectorXcd& vec(std::string lbl);
    bool exist_vec(std::string lbl) const { return vec_.find(lbl) != vec_.end(); }
    void InitVecIfNecessary(std::string lbl, int n);
  };

  void CalcGTOInt(int maxn, dcomplex a, dcomplex* res);

  
  class R1GTOs {
    // setup_q_==false => coef, offset, num_basis_, num_prim is not setted.
  public:
    struct Prim {
      int n;
      dcomplex z;
      Prim(int _n, dcomplex _z) : n(_n), z(_z) {}
    };
    struct Contraction {
      std::vector<Prim> prim;
      Eigen::MatrixXcd  coef;
      int offset;
      int size_prim() const { return coef.cols(); }
      int size_basis() const { return coef.rows(); }
    };
  public:
    //    bool coef_set_q_;        // coefcient is setup or not
    bool setup_q_; 
    std::string coef_type_;  // coefficient type (Nothing, normalized, derivativ)
    std::vector<Contraction> conts_;
    int num_basis_;
    int num_prim_;

    typedef std::vector<Prim>::iterator ItPrim;
    typedef std::vector<Prim>::const_iterator cItPrim;
    typedef std::vector<Contraction>::iterator ItCont;
    typedef std::vector<Contraction>::const_iterator cItCont;
  public:

    // ---- Constructors ----
    R1GTOs();
    void swap(R1GTOs& o);

    // ---- Utils ----
    int max_n() const;
    bool IsSameStructure(const R1GTOs& o) const;

    // ---- Getter -----
    //    int L() const { return L_; }
    int calc_size_basis() const;
    int calc_size_prim() const;
    int size_basis() const { return num_basis_; }
    int size_prim() const { return num_prim_; }
    const Prim& prim(int i) const;    
    Prim& prim(int i);
    dcomplex z_prim(int i) { return this->prim(i).z; }
    int      n_prim(int i) { return this->prim(i).n; }
    bool setup_q() const { return this->setup_q_; }
    std::string coef_type() const { return coef_type_; }

    // ---- Setter ----
    void Add(int n, dcomplex zeta);
    void Add(int n, const Eigen::VectorXcd& zs);    
    void Add(int n, const Eigen::VectorXcd& zs, const Eigen::MatrixXcd& coef);
    void Set(const Eigen::VectorXcd& zs);
    void Set(int n, const Eigen::VectorXcd& zs);

    // ---- basis convert ----
    void SetConj(const R1GTOs& o);
    void SetOneDeriv(const R1GTOs& o);
    void SetTwoDeriv(const R1GTOs& o);


    // ---- Matrix/Vector ----
    void CreateMat(const R1GTOs& o, Eigen::MatrixXcd& m) const;
    void CreateVec(Eigen::VectorXcd& m) const;
    void CalcMatSTV(const R1GTOs& o, int L,  MatVecMap& mat_vec,
		    std::string s_lbl, std::string t_lbl, std::string v_lbl) const;
    void CalcMatSTV(const R1GTOs& o, int L, Eigen::MatrixXcd& S, 
		    Eigen::MatrixXcd& T, Eigen::MatrixXcd& V) const;
    void CalcMatSTV(int L, MatVecMap& mat_vec,
		    std::string s_lbl, std::string t_lbl, std::string v_lbl) const;
    void CalcMatSTV(int L, Eigen::MatrixXcd& S, 
		    Eigen::MatrixXcd& T, Eigen::MatrixXcd& V) const;
    void CalcMatSTO(const R1GTOs& o, const R1STOs& v, MatVecMap& res, std::string) const;
    void CalcMatSTO(const R1STOs& v, MatVecMap& res, std::string) const;
    void CalcMatSTO(const R1GTOs& o, const R1STOs& v, Eigen::MatrixXcd&) const;
    void CalcMatSTO(const R1STOs& v, Eigen::MatrixXcd&) const;
    void CalcVec(const R1STOs& o, MatVecMap& mat_vec, std::string label="m") const;
    void CalcVec(const R1STOs& v, Eigen::VectorXcd& m) const;

    // ---- Other ----
    void SetUp();
    void Normalize();
    void AtR(const Eigen::VectorXcd&,
	     const Eigen::VectorXcd&, Eigen::VectorXcd*);
    void AtR(const Eigen::VectorXcd&, Eigen::VectorXcd*);
    Eigen::VectorXcd* AtR(const Eigen::VectorXcd&,
			  const Eigen::VectorXcd&);
    Eigen::VectorXcd* AtR(const Eigen::VectorXcd&);
    void DerivAtR(const Eigen::VectorXcd&,
		  const Eigen::VectorXcd&, Eigen::VectorXcd*) const;
    void DerivAtR(const Eigen::VectorXcd&, Eigen::VectorXcd*) const;
    Eigen::VectorXcd* DerivAtR(const Eigen::VectorXcd&,
			       const Eigen::VectorXcd&) const;    
    Eigen::VectorXcd* DerivAtR(const Eigen::VectorXcd&) const;
    void Deriv2AtR(const Eigen::VectorXcd&,
		  const Eigen::VectorXcd&, Eigen::VectorXcd*) const;
    void Deriv2AtR(const Eigen::VectorXcd&, Eigen::VectorXcd*) const;

    Eigen::VectorXcd* Deriv2AtR(const Eigen::VectorXcd&,
			       const Eigen::VectorXcd&) const;    
    Eigen::VectorXcd* Deriv2AtR(const Eigen::VectorXcd&) const;
			   
  };
  class R1STOs {
  public:
    bool normalized_q_;
    std::vector<R1STO> stos_;
  public:
    R1STOs(): normalized_q_(false) {}
    int size_basis() const { return stos_.size();}
    const R1STO& basis(int i) const { return stos_[i]; }
    bool normalized_q() const { return normalized_q_; }
    void Add(dcomplex c, int n, dcomplex zeta);
    void Add(int n, dcomplex zeta);
    void AtR(const Eigen::VectorXcd&, Eigen::VectorXcd*) const;
    Eigen::VectorXcd* AtR(const Eigen::VectorXcd&)       const;

  };

  std::ostream& operator<<(std::ostream& out, const R1GTOs& basis);
  std::ostream& operator<<(std::ostream& out, const R1STOs& basis);
  
}

#endif
