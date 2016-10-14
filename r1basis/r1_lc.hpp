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

    // ---- Setter ----
    _LC_EXPs<M>* Add(dcomplex c, int n, dcomplex z);

    // ---- Generate other ----
    LC_EXPs Conj() const;
    
    // ---- IO ----
    std::string str() const;
    
  };

  typedef boost::shared_ptr<_LC_EXPs<1> > LC_STOs;
  typedef boost::shared_ptr<_LC_EXPs<2> > LC_GTOs;
  LC_STOs Create_LC_STOs();
  LC_GTOs Create_LC_GTOs();
  
}


#endif
