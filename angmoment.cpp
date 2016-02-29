#include "angmoment.hpp"

namespace l2func {

  ExceptionBadYlm::ExceptionBadYlm(int L, int M) :std::exception() {
      std::stringstream ss;
      ss << "Unphysical (L, M) pair. (L, M) = (" << L << ", " << M << ")";
      msg_ = ss.str();
  }
  ExceptionBadYlm::~ExceptionBadYlm() throw() {}
  const char* ExceptionBadYlm::what() const throw() {
    return msg_.c_str();
  }

}
