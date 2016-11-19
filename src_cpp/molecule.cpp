#include "molecule.hpp"

using namespace std;
using namespace Eigen;

namespace cbasis {
  void _Molecule::Add(string name, Vector3cd at, dcomplex q) {
    Atom atom = {name, at, q};
    atoms.push_back(atom);
  }
  string _Molecule::str() const {
    ostringstream oss;
    oss << "==== Molecule ====" << endl;
    typedef vector<Atom>::const_iterator It;
    for(It it = atoms.begin(); it != atoms.end(); ++it ) {
      int i = distance(atoms.begin(), it);
      oss << "atom" << i << ": {"
	  << this->name(i) << ", ("
	  << this->x(i) << ", "
	  << this->y(i) << ", "
	  << this->z(i) << "), "
	  << this->q(i) << "}" << endl;
    }
    oss << "==================" << endl;
    return oss.str();
  }
}
