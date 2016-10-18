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
  
  template<int m>
  class _EXPs {
  public:
    // ---- typedef ----
    typedef typename _LC_EXPs<m>::LC_EXPs LC_EXPs;
    typedef boost::shared_ptr<_EXPs<m> > EXPs;
    
  private:
    // ---- Member field ----
    std::vector<LC_EXPs> basis_;
    bool setupq_;
    
  public:
    // ---- Constructors ----
    _EXPs();
    
    // ---- Getter ----
    int size()  const { return basis_.size(); }
    LC_EXPs basis(int i) const { return basis_[i]; }
    bool OnlyPrim() const;
    Eigen::VectorXcd AtR(const Eigen::VectorXcd&, 
			 const Eigen::VectorXcd&) const;
    Eigen::VectorXcd DAtR(const Eigen::VectorXcd&,
			  const Eigen::VectorXcd&) const;
    std::string str() const;
    
    // ---- Setter ----
    _EXPs<m>* AddPrim(int n, dcomplex z);
    _EXPs<m>* AddPrims(int n, Eigen::VectorXcd zs);
    _EXPs<m>* AddLC(LC_EXPs lc);
    _EXPs<m>* SetUp();
    
    // ---- Create ----
    EXPs Conj() const;
    EXPs Clone() const;
    
    // ---- Calculate ----
    Eigen::MatrixXcd CalcRmMat(int M) const;
    Eigen::MatrixXcd CalcD2Mat()      const;
    Eigen::VectorXcd CalcVecSTO(LC_STOs stos) const;
    Eigen::VectorXcd CalcVecGTO(LC_GTOs gtos) const;

  };    

  typedef boost::shared_ptr<_EXPs<1> > STOs;
  typedef boost::shared_ptr<_EXPs<2> > GTOs;
  template<int m>
  boost::shared_ptr<_EXPs<m> > Create_EXPs();
  STOs Create_STOs();
  GTOs Create_GTOs();
  dcomplex STOInt(int n, dcomplex a);
  dcomplex GTOInt(int n, dcomplex a);
  dcomplex EXPIntLC(LC_STOs a, int m, LC_STOs b);
  dcomplex EXPIntLC(LC_GTOs a, int m, LC_GTOs b);
 
}

#endif
