#ifndef BMATSET_H
#define  BMATSET_H

#include <map>
#include <Eigen/Core>
#include <boost/shared_ptr.hpp>
#include "../utils/typedef.hpp"

/**
   BMatSet is Block Matrix set.
   if ARG_NO_CHECK is defined, given argument for map and vector is not checked.
*/

namespace cbasis {

  typedef std::map<int, Eigen::VectorXcd> BVec;
  typedef std::map<std::pair<int, int>, Eigen::MatrixXcd> BMat;
  typedef std::map<std::string, BMat> BMatMap;

  void BMatRead(BMat& bmat, std::string fn);
  void BMatWrite(BMat& bmat, std::string fn);

  class _BMatSet {
  private:
    int     block_num_;
    BMatMap mat_map_;
  public:
    _BMatSet();
    _BMatSet(int _block_num);
    int block_num() const { return block_num_; }
    void SetMatrix(std::string name, int i, int j, Eigen::MatrixXcd& a);
    bool Exist(std::string, int i, int j);
    const Eigen::MatrixXcd& GetMatrix(std::string name, int i, int j);
    const BMat& GetBlockMatrix(std::string name);
    void SelfAdd(std::string name, int ib, int jb, int i, int j, dcomplex v);
    dcomplex GetValue(std::string name, int ib, int jb, int i, int j);
    void swap(_BMatSet& o);
    std::string str() const;
  };

  typedef boost::shared_ptr<_BMatSet> BMatSet;

  void swap(BMat& a, BMat& b);
}

#endif
