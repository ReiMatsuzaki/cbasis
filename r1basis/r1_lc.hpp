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

  class _LC_STOs;
  typedef boost::shared_ptr<_LC_STOs> LC_STOs;
  LC_STOs create_LC_STOs();

  class _LC_STOs {
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
    void AddOne(dcomplex c, int n, dcomplex z);
    _LC_STOs* Add(dcomplex c, int n, dcomplex z);

    // ---- Generate other ----
    LC_STOs Conj() const;
    
    // ---- IO ----
    std::string str() const;
    
  };

  class _LC_GTOs;
  typedef boost::shared_ptr<_LC_GTOs> LC_GTOs;
  LC_GTOs create_LC_GTOs();

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
