#ifndef BMATSET_H
#define  BMATSET_H

#include <map>
#include <Eigen/Core>
#include "typedef.hpp"

/**
   BMatSet is Block Matrix set.
   if ARG_NO_CHECK is defined, given argument for map and vector is not checked.
*/

namespace l2func {

  typedef std::map<std::pair<int, int>, Eigen::MatrixXcd> BMat;
  typedef std::map<std::string, BMat> BMatMap;

  class BMatSet {
  private:
    int     block_num_;
    BMatMap mat_map_;
  public:
    BMatSet();
    BMatSet(int _block_num);
    void SetMatrix(std::string name, int i, int j, Eigen::MatrixXcd& a);
    const Eigen::MatrixXcd& GetMatrix(std::string name, int i, int j);
    void SelfAdd(std::string name, int ib, int jb, int i, int j, dcomplex v);
    dcomplex GetValue(std::string name, int ib, int jb, int i, int j);
    void swap(BMatSet& o);
    std::string str() const;
  };  
  
  void swap(BMat& a, BMat& b);
  void swap(BMatSet& a, BMatSet& b);

}

#endif
