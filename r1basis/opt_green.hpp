#ifndef OPT_GREEN_H
#define OPT_GREEN_H

#include <vector>
#include <Eigen/Core>
#include "../utils/typedef.hpp"
#include "r1basis.hpp"

namespace cbasis {
  template<int m_s, int m, int m_r>
  class OptGreen {
  public:
    typedef typename _EXPs<m_s>::LC_EXPs LC_EXPs_S;
    typedef typename _EXPs<m>::EXPs      EXPs;
    typedef typename _EXPs<m_r>::LC_EXPs LC_EXPs_R;
    typedef Eigen::MatrixXcd Mat;
  private:

    LC_EXPs_S driv_S_;
    EXPs      basis_;
    LC_EXPs_R driv_R_;
    int L_;
    double Z_;
    double W_;    
    
    Eigen::MatrixXcd S, D2, R1, R2;
    //    Eigen::MatrixXcd S00, S01, S02, S11;
    std::vector<dcomplex> buf;
    Mat L00_, L01_, L02_, L11_;
    Mat S0_, S1_, S2_, R0_, R1_, R2_;
  public:

    // ==== Constructors ====
    OptGreen(LC_EXPs_S _driv_S, EXPs _basis, LC_EXPs_R _driv_R,
	     int _L, double _Z, double _W);

    // ==== Accessors ====
    int L() const { return L_; }
    void setL(int _L) { L_=_L; }
    double Z() const { return Z_;}
    void setZ(double _Z) { Z_=_Z; }
    double W() const { return W_;}
    void setW(double _W) { W_=_W; }
    const Mat& L00() const { return L00_; }

    // ==== Calculations ====
    // -- compute L00 and S0, R0
    void Calc_S0_L00_R0();
    void CalcV(dcomplex *val);
    void CalcVGH(dcomplex *val, Eigen::VectorXcd *grad, Eigen::MatrixXcd *hess);
  };	   
	   
}
#endif
