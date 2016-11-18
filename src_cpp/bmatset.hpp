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

  class BVec {
  public:
    typedef std::map<int, Eigen::VectorXcd>::iterator iterator;
    typedef std::map<int, Eigen::VectorXcd>::const_iterator const_iterator;
  private:
    std::string name_;
    std::map<int, Eigen::VectorXcd> map_;
  public:
    BVec() {}
    BVec(std::string _name) :name_(_name) {}
    iterator begin() { return map_.begin(); }
    const_iterator begin() const { return map_.begin(); }
    iterator end() { return map_.end(); }
    const_iterator end() const { return map_.end(); }
    void set_name(std::string _name) { name_ = _name; }
    std::string get_name() const { return name_; }
    bool has_irrep(int irrep) { return map_.find(irrep) != map_.end(); }
    Eigen::VectorXcd& operator()(int irrep) {return map_[irrep];}
    const Eigen::VectorXcd& operator()(int irrep) const {return map_.find(irrep)->second;}
    Eigen::VectorXcd& operator[] (int irrep) { return this->operator()(irrep); }
    const Eigen::VectorXcd& operator[](int irrep) const {return this->operator()(irrep); }
    
  };
  std::ostream& operator << (std::ostream& os, const BVec& a);
  
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
