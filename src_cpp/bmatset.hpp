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

  // ==== Block Vector ====
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
    bool has_block(int irrep) { return map_.find(irrep) != map_.end(); }
    int size() const { return map_.size(); }
    Eigen::VectorXcd& at(int irrep) {return map_[irrep];}
    const Eigen::VectorXcd& at(int irrep) const {return map_.find(irrep)->second;}    
    Eigen::VectorXcd& operator()(int irrep) {return at(irrep); }
    const Eigen::VectorXcd& operator()(int irrep) const {return at(irrep); }
    Eigen::VectorXcd& operator[] (int irrep) { return at(irrep); }
    const Eigen::VectorXcd& operator[](int irrep) const {return at(irrep); }
    void Write(std::string filename) const;
    void Read(std::string filename);
    
  };
  std::ostream& operator << (std::ostream& os, const BVec& a);

  // ==== Block Matrix ====
  class BMat {
  public:
    typedef std::pair<int, int> Key;
    typedef Eigen::MatrixXcd Value;
    typedef std::map<Key, Value> Map;
    typedef Map::iterator iterator;
    typedef Map::const_iterator const_iterator;
    
  private:
    std::string name_;
    Map map_;
  public:
    BMat() {}
    BMat(std::string _name): name_(_name) {}
    iterator begin() { return map_.begin(); }
    const_iterator begin() const { return map_.begin(); }
    iterator end() { return map_.end(); }
    const_iterator end() const { return map_.end(); }
    void set_name(std::string _name) { name_ = _name; }
    std::string get_name() const { return name_; }
    bool has_block(Key ijrrep) {
      return map_.find(ijrrep) != map_.end();
    }    
    bool has_block(int irrep, int jrrep) {
      return this->has_block(std::make_pair(irrep, jrrep));
    }
    int size() const { return map_.size(); }
    Value& operator()(int irrep, int jrrep) {
      return map_[std::make_pair(irrep, jrrep)];
    }
    const Value& operator()(int irrep, int jrrep) const {
      return map_.find(std::make_pair(irrep, jrrep))->second;
    }
    Value& operator[] (Key k) { return map_[k]; }
    const Value& operator[](Key k) const {return map_.find(k)->second; }
    
    void Write(std::string filename) const;
    void Read(std::string filename);
  };
  std::ostream& operator << (std::ostream& os, const BMat& a);

  // ==== Old ====
  void BMatRead(BMat::Map& bmat, std::string fn);
  void BMatWrite(BMat::Map& bmat, std::string fn);
  typedef std::map<std::string, BMat::Map> BMatMap;
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
    const BMat::Map& GetBlockMatrix(std::string name);
    void SelfAdd(std::string name, int ib, int jb, int i, int j, dcomplex v);
    dcomplex GetValue(std::string name, int ib, int jb, int i, int j);
    void swap(_BMatSet& o);
    std::string str() const;
  };
  typedef boost::shared_ptr<_BMatSet> BMatSet;
}

#endif
