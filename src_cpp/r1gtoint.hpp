#ifndef R1GTOINT_H
#define R1GTOINT_H

#include <vector>
#include <map>
#include <Eigen/Core>
#include "typedef.hpp"

// TODO
// CalcVecの引数をR1GTOsから、std::vector<R1GTO>に変更。
// それぞれのメッソドをAccessor等に分けて書く。


namespace l2func {

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
    typedef std::vector<Prim>::iterator ItPrim;
    typedef std::vector<Prim>::const_iterator cItPrim;
    typedef std::vector<Contraction>::iterator ItCont;
    typedef std::vector<Contraction>::const_iterator cItCont;
  public:
    bool coef_set_q_;         // coefcient is setup or not
    std::string coef_type_;        // coefficient type (Nothing, normalized, derivativ)
    std::vector<Contraction> conts_;
  public:

    // ---- Constructors ----
    R1GTOs();
    void swap(R1GTOs& o);

    // ---- Utils ----
    int max_n() const;
    bool IsSameStructure(const R1GTOs& o) const;

    // ---- Getter -----
    //    int L() const { return L_; }
    int size_basis() const;
    int size_prim() const;
    const Prim& prim(int i) const;    
    Prim& prim(int i);
    dcomplex z_prim(int i) { return this->prim(i).z; }
    int      n_prim(int i) { return this->prim(i).n; }
    bool coef_set_q() const { return coef_set_q_; }

    // ---- Setter ----
    void Add(int n, dcomplex zeta);
    void Add(int n, const Eigen::VectorXcd& zs);    
    void Add(int n, const Eigen::VectorXcd& zs, const Eigen::MatrixXcd& coef);
    void Set(int n, const Eigen::VectorXcd& zs);

    // ---- basis convert ----
    void SetConj(const R1GTOs& o);

    // ---- Calculation ----
    void Normalize();
    void CalcMatSTV(const R1GTOs& o, int L,  MatVecMap& mat_vec,
		    std::string s_lbl, std::string t_lbl, std::string v_lbl) const;
    void CalcMatSTV(int L, MatVecMap& mat_vec,
		    std::string s_lbl, std::string t_lbl, std::string v_lbl) const;
    void CalcMatSTO(const R1GTOs& o, const R1STOs& v, MatVecMap& res, std::string) const;
    void CalcMatSTO(const R1STOs& v, MatVecMap& res, std::string) const;
    void CalcVec(const R1STOs& o, MatVecMap& mat_vec, std::string label="m") const;

    void AtR(const Eigen::VectorXcd&,
	     const Eigen::VectorXcd&, Eigen::VectorXcd*);
    Eigen::VectorXcd* AtR(const Eigen::VectorXcd&,
			  const Eigen::VectorXcd&);
    
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
  };

  std::ostream& operator<<(std::ostream& out, const R1GTOs& basis);
  std::ostream& operator<<(std::ostream& out, const R1STOs& basis);
  
}

#endif
