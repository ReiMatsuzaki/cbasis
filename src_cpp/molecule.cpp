#include "molecule.hpp"

using namespace std;
using namespace Eigen;

namespace cbasis {
  void _Molecule::Add(Vector3cd at, dcomplex q) {
    at_q_list.push_back(make_pair(at, q));
  }
  string _Molecule::str() const {
    ostringstream oss;
    oss << "==== Molecule ====" << endl;
    typedef std::vector<pair<Vector3cd, dcomplex> >::const_iterator It;
    for(It it = at_q_list.begin(); it != at_q_list.end(); ++it ) {
      oss << it->first << " : " << it->second << endl;
    }
    return oss.str();
  }
}
