#ifndef R1GTOINT_H
#define R1GTOINT_H

#include <vector>
#include <map>
#include <Eigen/Core>
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>

#include "../utils/typedef.hpp"
#include "r1_lc.hpp"

#define COEF_NO 0	   // coefficient is nothing
#define COEF_NOT_NORMAL 1 // not normalized
#define COEF_NORMAL 2     // normalized

// represent set of linear combination of STO/GTO
// you can use normalized basis set or not normalized basis set
//

namespace cbasis {

  template<int m1, int m2>
  dcomplex EXPInt(int n, dcomplex a, dcomplex b);

  template<int m>
  dcomplex NormalizationTermContinue(int n, dcomplex a);
  
  template<int m>
  class _EXPs : public boost::enable_shared_from_this<_EXPs<m> > {
  public:
    // ---- typedef ----
    typedef typename _LC_EXPs<m>::LC_EXPs LC_EXPs;
    typedef boost::shared_ptr<_EXPs<m> > EXPs;
    
  private:
    // ---- Member field ----
    std::vector<LC_EXPs> basis_;
    std::vector<int>     coef_type_;
  
  public:
    // ---- Constructors ----
    _EXPs();
  
    // ---- Getter ----
    int size()  const { return basis_.size(); }
    LC_EXPs basis(int i) const { return basis_[i]; }
    bool IsPrim(int i) const;
    bool IsPrimAll() const;
    bool IsNormal(int i) const;
    bool IsNormalAll() const;
    bool HasCoef(int i) const;
    bool HasCoefAll() const;
    int exp_power() const { return m; }
    Eigen::VectorXcd AtR(const Eigen::VectorXcd& rs, 
			 const Eigen::VectorXcd& cs) const;
    Eigen::VectorXcd DAtR(const Eigen::VectorXcd& rs,
			  const Eigen::VectorXcd& cs) const;
    dcomplex AtR_One(dcomplex r, const Eigen::VectorXcd& cs) const;
    std::string str() const;
    EXPs self();
    
    // ---- Setter ----
    _EXPs<m>* AddPrim(int n, dcomplex z);
    _EXPs<m>* AddPrims(int n, Eigen::VectorXcd zs);
    _EXPs<m>* AddLC(LC_EXPs lc);
    _EXPs<m>* AddNotNormalPrim(dcomplex c, int n, dcomplex z);    
    _EXPs<m>* AddNotNormalLC(LC_EXPs lc);
    _EXPs<m>* SetUp();
    _EXPs<m>* ReplaceLC(int i, LC_EXPs lc);
    
    _EXPs<m>* SetStructure(int num, int n);
    _EXPs<m>* Clear() {this->basis_.clear(); this->coef_type_.clear(); return this;}
    // ---- Create ----
    EXPs Conj() const;
    EXPs Clone() const;
    void DerivOneZeta(EXPs other) const;
    void DerivTwoZeta(EXPs other) const;
  
    // ---- Calculate ----
    Eigen::MatrixXcd CalcRmMat(int M) const;
    Eigen::MatrixXcd CalcD2Mat()      const;
    Eigen::VectorXcd CalcVecSTO(LC_STOs stos) const;
    Eigen::VectorXcd CalcVecGTO(LC_GTOs gtos) const;

    void InitVec(Eigen::VectorXcd&);
    void InitMat(Eigen::MatrixXcd&);
    void CalcVec(LC_STOs stos, Eigen::VectorXcd&);
    void CalcVec(LC_GTOs gtos, Eigen::VectorXcd&);    
    void CalcRmMat(int M        , Eigen::MatrixXcd&);
    void CalcD2Mat(               Eigen::MatrixXcd&);

    
  };

  typedef boost::shared_ptr<_EXPs<1> > STOs;
  typedef boost::shared_ptr<_EXPs<2> > GTOs;
  template<int m>
  boost::shared_ptr<_EXPs<m> > Create_EXPs();
  STOs Create_STOs();
  GTOs Create_GTOs();
  
  template<int m1, int m2>
  Eigen::VectorXcd CalcVec(typename _EXPs<m1>::EXPs a,
			   typename _EXPs<m2>::LC_EXPs b);

  template<int m1, int m2>
  Eigen::MatrixXcd CalcRmMat(typename _EXPs<m1>::EXPs a,
			     int M,
			     typename _EXPs<m2>::EXPs b);
			     
  template<int m1, int m2>
  Eigen::MatrixXcd CalcD2Mat(typename _EXPs<m1>::EXPs a,
			     typename _EXPs<m2>::EXPs b);

  // ==== Calculation of matrix/vector (speed) ====
  template<int m1>
  void InitVec(typename _EXPs<m1>::EXPs a, Eigen::VectorXcd& m);
  template<int m1, int m2>
  void CalcVec(typename _EXPs<m1>::EXPs a,
	       typename _EXPs<m2>::LC_EXPs b,
	       Eigen::VectorXcd& vec);

  template<int m1, int m2>
  void initMat(typename _EXPs<m1>::EXPs a,
	       typename _EXPs<m2>::EXPs b, Eigen::MatrixXcd& mat);
  template<int m1, int m2>
  void CalcRmMat(typename _EXPs<m1>::EXPs a,
		 int M,
		 typename _EXPs<m2>::EXPs b,
		 Eigen::MatrixXcd& mat);
			     
  template<int m1, int m2>
  void CalcD2Mat(typename _EXPs<m1>::EXPs a,
		 typename _EXPs<m2>::EXPs b,
		 Eigen::MatrixXcd& mat);

  template<int m>
  void CalcSTVMat(typename _EXPs<m>::EXPs us, std::vector<dcomplex>& buf,
		  Eigen::MatrixXcd *S, Eigen::MatrixXcd *T, Eigen::MatrixXcd *V);  
}

#endif
