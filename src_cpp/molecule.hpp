#ifndef MOLECULE_H
#define MOLECULE_H

#include <iostream>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include "../utils/typedef.hpp"

namespace cbasis {

  struct _Atom {
    std::string name;
    dcomplex q;
    std::vector<Eigen::Vector3cd> xyz_list;
  };
  typedef boost::shared_ptr<_Atom> Atom;

  class _Molecule {
  public:
    typedef std::map<std::string, Atom> Map;
    Map atoms_map_;

    // -- data for calculation --
    std::vector<std::string> names;
    std::vector<dcomplex> qs;
    std::vector<Eigen::Vector3cd> xyzs;
  public:
    Eigen::Vector3cd& At(int i) { return xyzs[i]; }
    const Eigen::Vector3cd& At(int i) const { return xyzs[i]; }
    const std::string& name(int i) const { return names[i]; }
    dcomplex x(int i) const { return this->At(i)[0]; }
    dcomplex y(int i) const { return this->At(i)[1]; }
    dcomplex z(int i) const { return this->At(i)[2]; }
    dcomplex& q(int i) { return this->qs[i]; }
    const dcomplex& q(int i) const { return this->qs[i]; }
    int size() const { return names.size(); }
    void Add(std::string name, dcomplex q, std::vector<Eigen::Vector3cd>& xyz_list);
    void Add(std::string name, dcomplex q, Eigen::Vector3cd xyz);
    void Add(std::string name, dcomplex q, Eigen::MatrixXcd xyz_list);
    Atom atom(std::string name) { return atom_map_[name]; }
    
    std::string str() const;
  };
  typedef boost::shared_ptr<_Molecule> Molecule;
}

#endif
