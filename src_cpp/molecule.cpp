#include <stdexcept>
#include "molecule.hpp"
#include "../utils/macros.hpp"

using namespace std;
using namespace Eigen;

namespace cbasis {
  void _Molecule::Add(string name, dcomplex q, vector<Vector3cd>& xyz_list) {

    if(atoms_map_.find(name) != atoms_map_.end()) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + " : atom " + name + " already exist." ;
      throw runtime_error(msg);
    }
    Atom x(new _Atom);
    x->name = name;
    x->q = q;
    x->xyz_list = xyz_list;
    atoms_map_[name] = x;

    typedef vector<Vector3cd>::iterator It;
    for(It it = xyz_list.begin(); it != xyz_list.end(); ++it) {
      this->names.push_back(name);
      this->qs.push_back(q);
      this->xyzs.push_back(*it);
    }
  }
  void _Molecule::Add(string name, dcomplex q, Vector3cd xyz) {

    vector<Vector3cd> xyz_list;
    xyz_list.push_back(xyz);
    this->Add(name, q, xyz_list);

  }  
  void _Molecule::Add(std::string name, dcomplex q, Eigen::MatrixXcd xyz_mat) {

    if(xyz_mat.cols() != 3) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + ": cols must be 3";
    }

    vector<Vector3cd> xyz_list;
    for(int i = 0; i < xyz_mat.rows(); i++)
      xyz_list.push_back(Vector3cd(xyz_mat(i, 0), xyz_mat(i, 1), xyz_mat(i, 2)));
    this->Add(name, q, xyz_list);
    
  }
  string _Molecule::str() const {
    ostringstream oss;
    oss << "==== Molecule ====" << endl;

    typedef Map::const_iterator It;
    for(It it = atoms_map_.begin(); it != atoms_map_.end(); ++it) {
      int i = distance(atoms_map_.begin(), it);
      oss << "atom" << i << ": {"
	  << it->first << ", " << it->second->q << ", [";
      const vector<Vector3cd>& xyz_list = it->second->xyz_list;
      for(vector<Vector3cd>::const_iterator jt = xyz_list.begin();
	  jt != xyz_list.end(); ++jt) {
	const Vector3cd& xyz = *jt;
	oss << "[" << xyz[0] << xyz[1] << xyz[2] << "] ";
      }
      oss << "]}" ;
    }
    return oss.str();
  }
}
