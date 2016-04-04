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

  typedef std::map<std::string, Eigen::MatrixXcd> MatMap;
  typedef std::vector<R1GTO>::const_iterator CIt;
  typedef std::map<std::string, Eigen::VectorXcd> VecMap;
  class R1GTOs;
  class R1STOs;

  void CalcGTOInt(int maxn, dcomplex a, dcomplex* res);

  class R1GTOs {
  public:
    bool normalized_q_;
    std::vector<R1GTO> gtos_;
    int L_;
    std::vector<dcomplex> buf_;
  public:
    R1GTOs(int _L);
    int L() const { return L_; }
    int size_basis() const { return gtos_.size();}
    const R1GTO& basis(int i) const { return gtos_[i]; }
    R1GTO& basis(int i) { return gtos_[i]; }
    bool normalized_q() const { return normalized_q_; }
    void Add(dcomplex c, int n, dcomplex zeta);
    void Add(int n, dcomplex zeta);
    void Add(int n, const Eigen::VectorXcd& zs);    
    void Set(int n, const Eigen::VectorXcd& zs);
    int max_n() const;
    void Normalize();
    void Reserve();
    void Reserve(int n);
    void CalcMat(MatMap* res);
    MatMap* CalcMatNew();
    void CalcVec(R1GTOs& o, VecMap* res);    
    void CalcVecSTO(const R1STOs&, VecMap* res);
    VecMap* CalcVecSTONew(const R1STOs&);
    void AtR(const Eigen::VectorXcd&,
	     const Eigen::VectorXcd&, Eigen::VectorXcd*);
    Eigen::VectorXcd* AtR(const Eigen::VectorXcd&,
			  const Eigen::VectorXcd&);
    void swap(R1GTOs& o);
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
