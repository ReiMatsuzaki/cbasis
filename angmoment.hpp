#ifndef ANGMOMENT_HPP
#define ANGMOMENT_HPP

#include <sstream>
#include <exception>
#include <boost/exception/all.hpp>
#include <boost/throw_exception.hpp>


namespace l2func {

  // ==== Exception class ====
  class ExceptionBadYlm : public std::exception, public boost::exception {
  public:
    std::string msg_;
    ExceptionBadYlm(int L, int M);
    ~ExceptionBadYlm() throw();
    virtual const char* what() const throw();
  };

}

#endif
