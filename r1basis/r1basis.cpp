#include "r1basis.hpp"
#include "r1_lc.hpp"


using namespace std;
using namespace Eigen;

namespace cbasis {

  // ==== member field ====
  _GTOs::_GTOs() {
    this->setupq_ =false;
  }
  bool _GTOs::OnlyPrim() const {

    bool acc(true);

    for(vector<LC_GTOs>::const_iterator
	  it = this->basis_.begin();
	it != this->basis_.end();
	++it) {
      acc &= (*it)->size() == 1;
    }

    return acc;
    
  }
  _GTOs* _GTOs::AddPrim(int n, dcomplex z) {
    this->setupq_ = false;
    LC_GTOs g = Create_LC_GTOs();
    g->Add(1.0, n, z);
    this->basis_.push_back(g);
    return this;
  }
  _GTOs* _GTOs::AddPrims(int n, Eigen::VectorXcd zs) {
    this->setupq_ = false;
    for(int i = 0; i < zs.size(); i++)
      this->AddPrim(n, zs[i]);
    return this;
  }
  _GTOs* _GTOs::AddLC(LC_GTOs lc) {
    this->setupq_ =false;
    this->basis_.push_back(lc);
    return this;
  }
  _GTOs* _GTOs::SetUp() {
    this->setupq_ = true;
    return this;
  }

  // ==== create =====
  GTOs Create_GTOs() {

    GTOs ptr(new _GTOs());
    return ptr;
    
  }



}
