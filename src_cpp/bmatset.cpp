#include <string>
#include <ostream>
#include <fstream>
#include <iostream>
#include "macros.hpp"
#include "bmatset.hpp"

using namespace std;
using namespace Eigen;

namespace l2func {


  // ==== BMat ====
  // typedef map<pair<int, int>, Matrix> BMat 
  void BMatRead(BMat& bmat, string fn) {

    ifstream f(fn, ios::in|ios::binary);
    
    if(!f) {
      string msg; SUB_LOCATION(msg);
      msg += ": file not found";
      throw runtime_error(msg);
    }

    int num;
    f.read((char*)&num, sizeof(int));
    for(int i = 0; i < num; i++) {

      int irrep, jrrep, n, m;
      f.read((char*)&irrep, sizeof(int));
      f.read((char*)&jrrep, sizeof(int));
      f.read((char*)&n, sizeof(int));
      f.read((char*)&m, sizeof(int));
      MatrixXcd M(n, m);
      for(int i = 0; i < n; i++) {
	for(int j = 0; j < m; j++) {
	  dcomplex v;
	  f.read((char*)&v, sizeof(dcomplex));
	  M(i, j) = v;
	}
      }
      bmat[make_pair(irrep, jrrep)] = MatrixXcd::Zero(1, 1);
      bmat[make_pair(irrep, jrrep)].swap(M);
    }

  }
  void BMatWrite(BMat& bmat, string fn) {

    ofstream f;
    f.open(fn, ios::out|ios::binary|ios::trunc);
    
    if(!f) {
      string msg; SUB_LOCATION(msg); msg+=": file not found";
      throw runtime_error(msg);
    }

    int num = bmat.size();
    f.write((char*)&num, sizeof(int));

    for(BMat::iterator it = bmat.begin(); it != bmat.end(); ++it) {
      int irrep = it->first.first;
      int jrrep = it->first.second;
      f.write((char*)&irrep, sizeof(int));
      f.write((char*)&jrrep, sizeof(int));
      MatrixXcd& M = it->second;
      int n = M.rows();
      int m = M.cols();
      f.write((char*)&n, sizeof(int));
      f.write((char*)&m, sizeof(int));
      for(int i = 0; i < n; i++) {
	for(int j = 0; j < m; j++) {
	  f.write((char*)&M.coeff(i, j), sizeof(dcomplex));
	}
      }
    }
  }

  // ==== BlockMatrixSets ====
  BMatSet::BMatSet(): block_num_(1) {}
  BMatSet::BMatSet(int _block_num): block_num_(_block_num) {}
  void BMatSet::SetMatrix(string name, int i, int j, MatrixXcd& mat) {
#ifndef ARG_NO_CHECK
#endif

    MatrixXcd tmp;
    mat_map_[name][make_pair(i, j)] = tmp;
    mat_map_[name][make_pair(i, j)].swap(mat);

  }
  bool BMatSet::Exist(string name, int i, int j) {

    if(i < 0 || block_num_ <= i ||
       j < 0 || block_num_ <= j) {
      string msg; SUB_LOCATION(msg);
      stringstream ss;
      ss << msg;
      ss << ": Error for int i or j" << endl;
      ss << "i: " << i << endl;
      ss << "j: " << j << endl;
      ss << "block_num: " << block_num_ << endl;
      throw runtime_error(ss.str());
    }

    if(mat_map_.find(name) == mat_map_.end()) 
      return false;

    if(mat_map_[name].find(make_pair(i, j)) == mat_map_[name].end()) 
      return false;

    MatrixXcd& mat = mat_map_[name][make_pair(i, j)];
    if(mat.rows() == 0 || mat.cols() == 0)
      return false;
    
    return true;
  }
  const MatrixXcd& BMatSet::GetMatrix(string name, int i, int j) {
    return mat_map_[name][make_pair(i, j)];
  }
  const BMat& BMatSet::GetBlockMatrix(std::string name) {
    return mat_map_[name];
  }
  void BMatSet::SelfAdd(string name, int i, int j, int a, int b, dcomplex v) {

#ifndef ARG_NO_CHECK
    if(i < 0 || block_num_ <= i ||
       j < 0 || block_num_ <= j) {
      string msg; SUB_LOCATION(msg);
      msg += "Error for int i or j.";
      throw runtime_error(msg);
    }

    if(mat_map_.find(name) == mat_map_.end()) {
      string msg; SUB_LOCATION(msg);
      msg += " : key not found.\n" ;
      msg += "name : " + name;
      throw runtime_error(msg);
    }

    if(mat_map_[name].find(make_pair(i, j)) == mat_map_[name].end()) {
      string msg; SUB_LOCATION(msg);
      msg += " : matrix (i,j) is not found.\n";
      throw runtime_error(msg);
    }
#endif

    mat_map_[name][make_pair(i, j)](a, b) += v;
  }
  dcomplex BMatSet::GetValue(string name, int i, int j, int a, int b) {

#ifndef ARG_NO_CHECK
    if(mat_map_.find(name) == mat_map_.end()) {
      string msg; SUB_LOCATION(msg);
      msg += ": key not found.\n" ;
      msg += "name : " + name;
      throw runtime_error(msg);
    }

    if(i < 0 || block_num_ <= i ||
       j < 0 || block_num_ <= j) {
      string msg; SUB_LOCATION(msg);
      stringstream oss;
      oss << msg << endl << "Error for int i or j" << endl;
      oss << "i = " << i << endl;
      oss << "j = " << j << endl;
      oss << "block_num = " << block_num_ << endl;
      throw runtime_error(oss.str());
    }

    if(mat_map_[name].find(make_pair(i, j)) == mat_map_[name].end()) {
      string msg; SUB_LOCATION(msg);
      msg += ": matrix (i,j) is not found.\n";
      throw runtime_error(msg);
    }
#endif

    return mat_map_[name][make_pair(i, j)](a, b);
  } 
  void BMatSet::swap(BMatSet& o) {
    
    typedef BMatMap::iterator It;
    BMatMap tmp_o;
    for(It it_o = o.mat_map_.begin(); it_o != o.mat_map_.end();) {
      tmp_o[it_o->first] = BMat();
      ::swap(tmp_o[it_o->first], it_o->second);
      o.mat_map_.erase(it_o++);
    }
    
    for(It it_this = this->mat_map_.begin(); it_this != this->mat_map_.end();) {
      o.mat_map_[it_this->first] = BMat();
      o.mat_map_[it_this->first].swap(it_this->second);
      this->mat_map_.erase(it_this++);
    }

    for(It it_tmp = tmp_o.begin(); it_tmp != tmp_o.end();) {
      this->mat_map_[it_tmp->first] = BMat();
      this->mat_map_[it_tmp->first].swap(it_tmp->second);
      ++it_tmp;
    }

    int tmp = o.block_num_;
    o.block_num_ = this->block_num_;
    this->block_num_ = tmp;

  }
  string BMatSet::str() const {
    ostringstream oss;
    for(BMatMap::const_iterator it = mat_map_.begin(); it != mat_map_.end(); ++it) {
      oss << it->first << endl;
      for(BMat::const_iterator it_mat = it->second.begin(); it_mat != it->second.end();
	  ++it_mat) {
	oss << it_mat->first.first << ", " << it_mat->first.second<< endl
	    << it_mat->second << endl;
      }
    }
    return oss.str();
  }
	      
  // ==== External ====
  void swap(BMat& a, BMat& b) {

    typedef map<pair<int, int>, MatrixXcd>::iterator It;
    map<pair<int, int>, MatrixXcd> tmp_a;
    for(It it_a = a.begin(); it_a != a.end();) {
      tmp_a[it_a->first] = MatrixXcd();
      tmp_a[it_a->first].swap(it_a->second);
      a.erase(it_a++);
    }

    for(It it_b = b.begin(); it_b != b.end();) {
      a[it_b->first] = MatrixXcd();
      a[it_b->first].swap(it_b->second);
      b.erase(it_b++);
    }

    for(It it_tmp = tmp_a.begin(); it_tmp != tmp_a.end();) {
      b[it_tmp->first] = MatrixXcd();
      b[it_tmp->first].swap(it_tmp->second);
      it_tmp++;
    }
  }
  void swap(BMatSet& a, BMatSet& b) {
    a.swap(b);
  }

}
