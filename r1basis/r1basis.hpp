#ifndef R1GTOINT_H
#define R1GTOINT_H

#include <vector>
#include <map>
#include <Eigen/Core>

#include "../src_cpp/typedef.hpp"
#include "r1_lc.hpp"

#define COEF_NO 0	   // coefficient is nothing
#define COEF_NOT_NORMAL 1 // not normalized
#define COEF_NORMAL 2     // normalized


//
// represent set of linear combination of STO/GTO
// you can use normalized basis set or not normalized basis set
//

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
    std::vector<int>     coef_type_;
    bool setupq_;
    
  public:
    // ---- Constructors ----
    _EXPs();
    
    // ---- Getter ----
    int size()  const { return basis_.size(); }
    LC_EXPs basis(int i) const { return basis_[i]; }
    bool OnlyPrim() const;
    bool IsNormal() const;
    bool HasCoef() const;
    int exp_power() const { return m; }
    Eigen::VectorXcd AtR(const Eigen::VectorXcd&, 
			 const Eigen::VectorXcd&) const;
    Eigen::VectorXcd DAtR(const Eigen::VectorXcd&,
			  const Eigen::VectorXcd&) const;
    std::string str() const;
    
    // ---- Setter ----
    _EXPs<m>* AddPrim(int n, dcomplex z);
    _EXPs<m>* AddPrims(int n, Eigen::VectorXcd zs);
    _EXPs<m>* AddLC(LC_EXPs lc);
    _EXPs<m>* AddNotNormalPrim(dcomplex c, int n, dcomplex z);    
    _EXPs<m>* AddNotNormalLC(LC_EXPs lc);
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

  template<int m1, int m2>
  Eigen::MatrixXcd CalcRmMat(_EXPs<m1>& a,
			     int M,
			     _EXPs<m2>& b);
  template<int m1, int m2>
  Eigen::MatrixXcd CalcD2Mat(_EXPs<m1>& a,
			     _EXPs<m2>& b);  
  /*
  template<int m1, int m2>
  Eigen::MatrixXcd CalcRMMat(boost::shared_ptr<_EXPs<m1> > a,
			     int M,
			     boost::shared_ptr<_EXPs<m2> >b);
  template<int m1, int m2>
  Eigen::MatrixXcd CalcD2Mat(boost::shared_ptr<_EXPs<m1> > a,
			     boost::shared_ptr<_EXPs<m2> >b);
  */
  /*    
  Eigen::MatrixXcd CalcRMMatSS(STOs a, int M, STO b);
  Eigen::MatrixXcd CalcRMMatSG(STOs a, int M, GTO b);
  Eigen::MatrixXcd CalcRMMatGS(GTOs a, int M, STO b);
  Eigen::MatrixXcd CalcRMMatGG(GTOs a, int M, GTO b);

  Eigen::MatrixXcd CalcD2MatSS(STOs a, STO b);
  Eigen::MatrixXcd CalcD2MatSG(STOs a, GTO b);
  Eigen::MatrixXcd CalcD2MatGS(GTOs a, STO b);
  Eigen::MatrixXcd CalcD2MatGG(GTOs a, GTO b);

  Eigen::VectorXcd CalcVecSS(STOs a, LC_STOs b);
  Eigen::VectorXcd CalcVecSG(STOs a, LC_GTOs b);
  Eigen::VectorXcd CalcVecGS(GTOs a, LC_STOs b);
  Eigen::VectorXcd CalcVecGG(GTOs a, LC_GTOs b);  
  */
 
}

#endif
