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
    struct Atom {
      std::string name;
      Eigen::Vector3cd xyz;
      dcomplex q;
      const std::string& get_name() const { return name; }
      const Eigen::Vector3cd& get_xyz() const { return xyz; }
      dcomplex get_q() const { return q; }
    };
    std::vector<Atom> atoms;
  public:
    Eigen::Vector3cd& At(int i) { return atoms[i].xyz; }
    const Eigen::Vector3cd& At(int i) const { return atoms[i].xyz; }
    const std::string& name(int i) const { return this->atoms[i].get_name(); }
    dcomplex x(int i) const { return this->At(i)[0]; }
    dcomplex y(int i) const { return this->At(i)[1]; }
    dcomplex z(int i) const { return this->At(i)[2]; }
    dcomplex& q(int i) { return this->atoms[i].q; }
    const dcomplex& q(int i) const { return this->atoms[i].q; }
    int size() const { return atoms.size(); }
    void Add(Eigen::Vector3cd at, dcomplex q);
    void Add(std::string name, Eigen::Vector3cd at, dcomplex q);
    std::string str() const;
  };
  typedef boost::shared_ptr<_Molecule> Molecule;
}

#endif
