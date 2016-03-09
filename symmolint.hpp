#ifndef SYMMOLINT_H
#define SYMMOLINT_H

#include "math_utils.hpp"
#include <Eigen/Core>
#include <string>

namespace l2func {

  /*
  struct Symmetry {
    std::string name;    
    Symmetry(std::string);
  };
  */

  typedef int Symmetry;

  struct Contraction {
    int sym;
    Eigen::MatrixXcd coef_iat_ipn;
    Contraction(int _sym, Eigen::MatrixXcd _coef): sym(_sym), coef_iat_ipn(_coef){}
  };

  struct SubSymGTOs {
    Eigen::MatrixXcd xyz_iat;
    Eigen::MatrixXi  ns_ipn;
    std::vector<Contraction> contraction_icont;
    Eigen::VectorXcd zeta_iz;
    Eigen::MatrixXcd coef_iz_icont;

    int size_at() const { return xyz_iat.cols();}
    int size_pn() const { return ns_ipn.cols(); }
    int size_cont() const { return contraction_icont.size(); }
    int size_zeta() const { return zeta_iz.rows(); }
    //    int size_prim() const { return size_pn() * size_zeta() * size_at(); }
    SubSymGTOs(Eigen::MatrixXcd xyz, Eigen::MatrixXi ns,
	       std::vector<Contraction> cs, Eigen::VectorXcd zs) :
      xyz_iat(xyz), ns_ipn(ns), contraction_icont(cs), zeta_iz(zs) {

      coef_iz_icont = Eigen::MatrixXcd::Zero(this->size_zeta(), this->size_cont());

    }
  };

  class SymGTOs {
  public:
    std::vector<Symmetry> symmetries;
    std::vector<SubSymGTOs> sub_sym_gtos;
  public:
    void AddSymmetry(Symmetry sym);
    void AddSub(Eigen::MatrixXcd xyz, Eigen::MatrixXi ns,
	       std::vector<Contraction> cs, Eigen::VectorXcd zs) {
      sub_sym_gtos.push_back(SubSymGTOs(xyz, ns, cs, zs));
    }
    void loop();
  };

  /*
  struct Shell {
    dcomplex zeta;
    Eigen::MatrixXcd xyz_iat;
    Eigen::MatrixXi ns_iprim;
    std::vector<Contraction> contractions;
    int offset;

    int size_at_prim() const { return xyz_iat.size() * ns_iprim.size(); }
    int size_at() const { return xyz_iat.size(); }
    int size_prim() const { return ns_iprim.size(); }
    int size_contractions() const { return contractions.size();}
    Shell(dcomplex _zeta, Eigen::MatrixXcd _xyz, Eigen::MatrixXi _ns) :
      zeta(_zeta), xyz_iat(_xyz), ns_iprim(_ns) {}
    void AddContraction(const Contraction& cont) {
      contractions.push_back(cont);
    }
  };

  class SymGTOs {
  private:
    std::vector<Symmetry> symmetries_;
    std::vector<Shell> shells_;
  public:
    const std::vector<Symmetry> symmetries() const { return symmetries_; }
    const std::vector<Shell> shells() const { return shells_; }
    void AddSymmetry(Symmetry sym);
    void AddShell(const Shell& shell);
    void loop();
  };
  */

}

#endif
