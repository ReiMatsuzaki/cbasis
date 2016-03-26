#ifndef MULT_ARRAY_H
#define MULT_ARRAY_H

// header library
#include <string>
#include <iostream>
#include <ostream>
#include "macros.hpp"

namespace l2func {

  template<class F, int N>
  class MultArray {};

  template<class F>
  class MultArray<F, 2> {
  private:
    static const int N = 2;
    F* data_;
    int data_num_;  // size of data_ (Capacity)
    int num_;       // Prod_i(n1_[i] - n0_[i] + 1)
    int n0_[N];
    int n1_[N];
  public:
    MultArray(int _num0) {
      data_ = new F[_num0];
      data_num_ = _num0;
      num_ = _num0;
    }
    ~MultArray() {
      delete[] data_;
    }
    void SetRange(int nx0, int nx1, int ny0, int ny1) {
      n0_[0] = nx0; n0_[1] = ny0;
      n1_[0] = nx1; n1_[1] = ny1;
      num_ = (nx1-nx0+1)*(ny1-ny0+1);
      if(data_num_ < num_) {
	delete[] data_;
	data_num_ = num_;
	data_ = new F[num_];
      }
    }
    int idx(int nx, int ny) {
      return ((nx - n0_[0]) * (n1_[1] - n0_[1] + 1) + 
	      (ny - n0_[1]));
    }
    void check_index(int nx, int ny) {
      int index = this->idx(nx, ny);
      if(index < 0   || num_-1 < index ||
	 nx < n0_[0] || n1_[0] < nx ||
	 ny < n0_[1] || n1_[1] < ny) {

	std::string msg;
	std::stringstream ss;
	SUB_LOCATION(msg);
	ss << ": index: (" << nx << ", " << ny << ") "
	   << index << std::endl;
	msg += ss.str();
	throw std::runtime_error(msg);
      }
    }
    void set(int nx, int ny, F v) {
      data_[idx(nx, ny)] = v;
    }
    void set_safe(int nx, int ny, F v) {
      try{
	check_index(nx, ny);
      } catch(const std::exception& e) {
	throw e;
      }
      set(nx, ny, v);
    }
    int size() const { return num_; }
    F& get(int nx, int ny) {
      return data_[idx(nx, ny)];
    }
    F& get_safe(int nx, int ny) {
      try{
	check_index(nx, ny);
      } catch(const std::runtime_error& e) {
	throw e;
      }
      return this->get(nx, ny);
    }
    F& operator()(int nx, int ny) {

#ifndef ARG_NO_CHECK
      try {
	this->check_index(nx, ny);
      } catch(const std::runtime_error& e) {
	std::string msg; SUB_LOCATION(msg);
	msg += ": ";
	msg += e.what();
	throw std::runtime_error(msg);
      }
#endif
      return get_safe(nx, ny);
    }
  };

  template<class F>
  class MultArray<F, 3> {
  private:
    static const int N = 3;
    F* data_;
    int data_num_;
    int num_;
    int n0_[N];
    int n1_[N];
  public:
    MultArray(int _num0) {
      data_ = new F[_num0];
      data_num_ = _num0;
      num_ = _num0;
    }
    ~MultArray() {
      delete[] data_;
    }
    void SetRange(int nx0, int nx1, int ny0, int ny1, int nz0, int nz1) {
      n0_[0] = nx0; n0_[1] = ny0; n0_[2] = nz0; 
      n1_[0] = nx1; n1_[1] = ny1; n1_[2] = nz1; 
      num_ = (nx1-nx0+1)*(ny1-ny0+1)*(nz1-nz0+1);
      if(data_num_ < num_) {
	delete[] data_;
	data_num_ = num_;
	data_ = new F[num_];
      }
    }
    int idx(int nx, int ny, int nz) {
      return ((nx - n0_[0]) * (n1_[2] - n0_[2] + 1) * (n1_[1] - n0_[1] + 1) + 
	      (ny - n0_[1]) * (n1_[2] - n0_[2] + 1) + 
	      (nz - n0_[2]));
    }
    void check_index(int nx, int ny, int nz) {
      int index = this->idx(nx, ny, nz);
      if(index < 0   || num_-1 < index ||
	 nx < n0_[0] || n1_[0] < nx ||
	 ny < n0_[1] || n1_[1] < ny ||
	 nz < n0_[2] || n1_[2] < nz) {
	std::string msg;
	std::stringstream ss;
	SUB_LOCATION(msg);
	ss << "index: (" << nx << ", " << ny << ", " << nz << ") "
	   << index << std::endl;
	msg += ss.str();
	throw std::runtime_error(msg);
      }
    }
    void set(int nx, int ny, int nz, F v) {
      data_[idx(nx, ny, nz)] = v;
    }
    void set_safe(int nx, int ny, int nz, F v) {
      check_index(nx, ny, nz);
      set(nx, ny, nz, v);
    }
    int size() const { return num_; }
    F& get(int nx, int ny, int nz) {
      return data_[idx(nx, ny, nz)];
    }
    F& get_safe(int nx, int ny, int nz) {
      check_index(nx, ny, nz);
      return this->get(nx, ny, nz);
    }
    F& operator()(int nx, int ny, int nz) {

#ifndef ARG_NO_CHECK
      try {
	this->check_index(nx, ny, nz);
      } catch(const std::runtime_error& e) {
	std::string msg; SUB_LOCATION(msg);
	msg += ": ";
	msg += e.what();
	throw std::runtime_error(msg);
      }
#endif

      return get_safe(nx, ny, nz);
    }
  };

  template<class F>
  class MultArray<F, 4> {
  private:
    static const int N = 4;
    F* data_;
    int data_num_;
    int num_;
    int n0_[N];
    int n1_[N];
  public:
    MultArray(int _num0) {
      data_ = new F[_num0];
      data_num_ = _num0;
      num_ = _num0;
    }
    ~MultArray() {
      delete[] data_;
    }
    void SetRange(int nx0, int nx1, int ny0, int ny1,
		  int nz0, int nz1, int nw0, int nw1) {
      n0_[0] = nx0; n0_[1] = ny0; n0_[2] = nz0; n0_[3] = nw0; 
      n1_[0] = nx1; n1_[1] = ny1; n1_[2] = nz1; n1_[3] = nw1; 
      num_ = (nx1-nx0+1)*(ny1-ny0+1)*(nz1-nz0+1)*(nw1-nw0+1);
      if(data_num_ < num_) {
	delete[] data_;
	data_num_ = num_;
	data_ = new F[num_];
      }
    }
    int idx(int nx, int ny, int nz, int nw) {
      int num1 = n1_[1] - n0_[1] + 1;
      int num2 = n1_[2] - n0_[2] + 1;
      int num3 = n1_[3] - n0_[3] + 1;
      return ((nx - n0_[0]) * num3 * num2 * num1 + 
	      (ny - n0_[1]) * num3 * num2 + 
	      (nz - n0_[2]) * num3 +
	      (nw - n0_[3]));
    }
    void check_index(int nx, int ny, int nz, int nw) {
      int index = this->idx(nx, ny, nz, nw);
      if(index < 0   || num_-1 < index ||
	 nx < n0_[0] || n1_[0] < nx ||
	 ny < n0_[1] || n1_[1] < ny ||
	 nz < n0_[2] || n1_[2] < nz ||
	 nw < n0_[3] || n1_[3] < nw) {
	std::string msg;
	std::stringstream ss;
	SUB_LOCATION(msg);
	ss << ": index: (" << nx << ", " << ny << ", " << nz << ", " << nw << ") "
	   << index << std::endl;
	for(int i = 0; i < N; i++)
	  ss << n0_[i] << ", " << n1_[i] << " ";
	msg += ss.str();
	throw std::out_of_range(msg);
      }
    }
    void set(int nx, int ny, int nz, int nw, F v) {
      data_[idx(nx, ny, nz, nw)] = v;
    }
    void set_safe(int nx, int ny, int nz, int nw, F v) {
      check_index(nx, ny, nz, nw);
      set(nx, ny, nz, nw, v);
    }
    F& get(int nx, int ny, int nz, int nw) {
      return data_[idx(nx, ny, nz, nw)];
    }
    F& get_safe(int nx, int ny, int nz, int nw) {
      check_index(nx, ny, nz, nw);
      return this->get(nx, ny, nz, nw);
    }
    F& operator()(int nx, int ny, int nz, int nw) {

#ifndef ARG_NO_CHECK
      try {
	this->check_index(nx, ny, nz, nw);
      } catch(const std::runtime_error& e) {
	std::string msg; SUB_LOCATION(msg);
	msg += ": ";
	msg += e.what();
	throw std::runtime_error(msg);
      }
#endif

      return get_safe(nx, ny, nz, nw);
    }
  };

}
#endif
