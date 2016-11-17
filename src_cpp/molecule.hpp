#ifndef MOLECULE_H
#define MOLECULE_H

#include <iostream>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include "../utils/typedef.hpp"

namespace cbasis {

  class _Molecule {
  public:
    std::vector<std::pair<Eigen::Vector3cd, dcomplex> > at_q_list;
  public:
    Eigen::Vector3cd& At(int i) { return at_q_list[i].first; }
    const Eigen::Vector3cd& At(int i) const { return at_q_list[i].first; }
    dcomplex x(int i) const { return at_q_list[i].first[0]; }
    dcomplex y(int i) const { return at_q_list[i].first[1]; }
    dcomplex z(int i) const { return at_q_list[i].first[2]; }
    dcomplex& q(int i) { return at_q_list[i].second; }
    const dcomplex& q(int i) const { return at_q_list[i].second; }
    int size() const { return at_q_list.size(); }
    void Add(Eigen::Vector3cd at, dcomplex q);
    std::string str() const;
  };
  typedef boost::shared_ptr<_Molecule> Molecule;
}

#endif
