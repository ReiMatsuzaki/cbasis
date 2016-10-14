#ifndef R1_LC_HPP
#define R1_LC_HPP

#include <vector>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include "../src_cpp/typedef.hpp"

/**
   Linear combination of functions
 */

namespace cbasis {

  template<int M>
  class _LC_EXPs {
  private:
    // ---- Member field ----
    std::vector<dcomplex> cs;
    std::vector<int>      ns;
    std::vector<dcomplex> zs;
    
    // ---- typedef ----
    typedef boost::shared_ptr<_LC_EXPs<M> > LC_EXPs;
    
  public:
    // ---- Constructors ----

    // ---- Utils ----
    int max_n() const;
    int power() const;

    // ---- Getter ----
    int size() const;
    Eigen::VectorXcd AtR(const Eigen::VectorXcd&) const;
    Eigen::VectorXcd DAtR(const Eigen::VectorXcd&) const;
    //Eigen::VectorXcd D2AtR(const Eigen::VectorXcd&) const;

    // ---- Setter ----
    _LC_EXPs<M>* Add(dcomplex c, int n, dcomplex z);

    // ---- Generate other ----
    LC_EXPs Conj() const;
    
    // ---- IO ----
    std::string str() const;
    
  };

  typedef boost::shared_ptr<_LC_EXPs<1> > LC_STOs;
  typedef boost::shared_ptr<_LC_EXPs<2> > LC_GTOs;
  template<int M>
  boost::shared_ptr<_LC_EXPs<M> > Create_LC_EXPs();
  LC_STOs Create_LC_STOs();
  LC_GTOs Create_LC_GTOs();

  class _LC_GTOs {
  private:
    // ---- Member field ----
    std::vector<dcomplex> cs;
    std::vector<int>      ns;
    std::vector<dcomplex> zs;
    
  public:
    // ---- Constructors ----

    // ---- Utils ----
    int max_n() const;
    // bool IsSameStructure(const LC_STOs& o) const;

    // ---- Getter ----
    int size() const;
    void AtR(const Eigen::VectorXcd&, Eigen::VectorXcd&) const;
    void DAtR(const Eigen::VectorXcd&, Eigen::VectorXcd&) const;
    void D2AtR(const Eigen::VectorXcd&, Eigen::VectorXcd&) const;

    // ---- Setter ----
    _LC_GTOs* Add(dcomplex c, int n, dcomplex z);

    // ---- Generate other ----
    LC_GTOs Conj() const;
    
  };
  
  
}


#endif
