#ifndef MULT_ARRAY_H
#define MULT_ARRAY_H

#include <string>
#include <ostream>
#include <iostream>
#include <stdexcept>
#include "macros.hpp"

/**
   MultArray gives multi indexed array. This file is header library for inline optimization.
*/

namespace l2func {

  template<class F, int N>
  class MultArray {};

  template<class F>
  class MultArray<F, 1> {
  private:
    static const int N = 1;
    F* data_;
    int data_num_;
    int num_;
    int n0_[N];
    int n1_[N];
  public:
    MultArray(int _num) {
      data_ = new F[_num];
      data_num_ = _num;
      num_ = _num;
    }
    ~MultArray() {
      delete[] data_;
    }
    int size() const { return num_; }
    void SetRange(int n0, int n1) {
      n0_[0] = n0; n1_[0] = n1;
      int num_ = n1-n0+1;
      if(data_num_ < num_) {
	delete[] data_;
	data_num_ = num_;
	data_ = new F[num_];
      }
    }
    F& operator()(int nx) {
      int index = nx - n0_[0];

#ifndef ARG_NO_CHECK

      if(index < 0   || num_-1 < index ||
	 nx < n0_[0] || n1_[0] < nx) {
	std::string msg; SUB_LOCATION(msg);
	std::ostringstream oss;
	oss << std::endl << msg;
	oss << "index: (" << nx << ") "
	   << index << std::endl;
	throw std::runtime_error(oss.str());
      }
#endif
    return data_[index];
    }
  };

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
    int size() const { return num_; }
    void SetRange(int nx0, int nx1, int ny0, int ny1) {
      n0_[0] = nx0; n0_[1] = ny0;
      n1_[0] = nx1; n1_[1] = ny1;
      int num_ = (nx1-nx0+1)*(ny1-ny0+1);
      if(data_num_ < num_) {
	delete[] data_;
	data_num_ = num_;
	data_ = new F[num_];
      }
    }
    F& operator()(int nx, int ny) {
      int index = ((nx - n0_[0]) * (n1_[1] - n0_[1] + 1) + 
		   (ny - n0_[1]));
#ifndef ARG_NO_CHECK

      if(index < 0   || num_-1 < index ||
	 nx < n0_[0] || n1_[0] < nx ||
	 ny < n0_[1] || n1_[1] < ny) { 
	std::string msg;
	std::stringstream ss;
	SUB_LOCATION(msg);
	ss << "index: (" << nx << ", " << ny << ") "
	   << index << std::endl;
	msg += ss.str();
	throw std::runtime_error(msg);
      }
#endif
      return data_[index];
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
      this->num_ = (nx1-nx0+1)*(ny1-ny0+1)*(nz1-nz0+1);
      if(data_num_ < num_) {
	delete[] data_;
	data_num_ = num_;
	data_ = new F[num_];
      }
    }
    int size() const { return num_; }
    F& operator()(int nx, int ny, int nz) {
      int index = ((nx - n0_[0]) * (n1_[2] - n0_[2] + 1) * (n1_[1] - n0_[1] + 1) + 
		   (ny - n0_[1]) * (n1_[2] - n0_[2] + 1) + 
		   (nz - n0_[2]));
#ifndef ARG_NO_CHECK

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
#endif

      return data_[index];
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
    void SetValue(F val) {
      for(int i = 0; i < data_num_; i++)
	data_[i] = val;
    }
    F& operator()(int nx, int ny, int nz, int nw) {

      int num1 = n1_[1] - n0_[1] + 1;
      int num2 = n1_[2] - n0_[2] + 1;
      int num3 = n1_[3] - n0_[3] + 1;
      int index = ((nx - n0_[0]) * num3 * num2 * num1 + 
		   (ny - n0_[1]) * num3 * num2 + 
		   (nz - n0_[2]) * num3 +
		   (nw - n0_[3]));

#ifndef ARG_NO_CHECK
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
#endif
      return data_[index];
    }
  };

}
/*
namespace future {
  
  template<class F, int N>
  class MultArray {

    // ---- Typedef ----
  public:
    typedef MultArray<F, N-1> Element;

    // ---- Field member ----
    Element* data;
    int step;
    int n0;
    int n1;

    // ---- Constructors ----
    MultArray(int num): elements(num), n0(0), n1(num-1) {}
    static Create
    ~MultArray() {}
    int size() const { return elements.size(); }
    Element& SetRange(int _n0, int _n1) { n0 = _n0; n1 = _n1; }
    void CheckArg(int index, int n) {
      int num(this->size());
      if(index < 0   || num-1 < index) {
	std::string msg;
	std::stringstream ss;
	SUB_LOCATION(msg);
	ss << "index: " << n << ", "
	   << index << std::endl;
	msg += ss.str();
	throw std::runtime_error(msg);
      }
    }
    Element& operator[](int n) {
      int index = (n-n0);
#ifndef ARG_NO_CHECK
      CheckArg(index, n);
#endif
      return elements[index];
    }
    const Element& operator[](int n) const {
      int index = n-n0;
#ifndef ARG_NO_CHECK
      CheckArg(index, n);
#endif
      return elements[index];      
    }
  };
  
}
*/
#endif
