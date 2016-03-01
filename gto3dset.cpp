#include "gto3dset.hpp"

namespace l2func {

  // typedef std::complex<double> dcmplx;
  // typedef SphericalGTO<dcmplx, array3<dcmplx> > Basis;

  SphericalGTOSet::SphericalGTOSet() {}
  SphericalGTOSet::~SphericalGTOSet() {
    for(iterator it = basis_list_.begin(); 
	it != basis_list_.end(); ++it) {
      delete(*it);
    }
  }
  int SphericalGTOSet::size() const { return basis_list_.size(); }
  const SGTO& SphericalGTOSet::basis(int i) const {
    return *basis_list_[i];
  }
  void SphericalGTOSet::AddOneBasis(int L, int M, dcomplex x, dcomplex y, dcomplex z, dcomplex zeta) {
    SGTO* ptr;
    try {
      ptr = new SGTO(L, M, c3(x, y, z), zeta);
    } catch(const ExceptionBadYlm& e) {
      BOOST_THROW_EXCEPTION(e);
    }
    this->basis_list_.push_back(ptr);
  }
  void SphericalGTOSet::AddBasis(int L, dcomplex x, dcomplex y, dcomplex z, dcomplex zeta) {

    for (int M = -L; M <= L; M++) {
      this->AddOneBasis(L, M, x, y, z, zeta);
    }

  }

  template<class TOp>
  dcomplex* CalcMat(const SphericalGTOSet& a, const TOp& op, const SphericalGTOSet& b) {
    int numi = a.size();
    int numj = b.size();
    dcomplex* vs = new dcomplex[numi * numj];

    if(&a == &b) {
      for(int i = 0; i < numi; i++) {
	vs[i*numj + i] = CIP(a.basis(i), op, a.basis(i));
	for(int j = 0; j < i; j++) {
	  dcomplex v = CIP(a.basis(i), op, a.basis(j));
	  vs[i*numj + j] = v;
	  vs[j*numj + i] = v;
	}
      }
    } else {
      for(int i = 0; i < numi; i++) 
	for(int j = 0; j < numj; j++) 
	  vs[i*numj + j] = CIP(a.basis(i), op, b.basis(j));
    }
    return vs;
  }
  dcomplex* SphericalGTOSet::SMat(const SphericalGTOSet& o) const {
    return CalcMat(*this, OpXyz<dcomplex, c3>(0, 0, 0), o);
  }
  dcomplex* SphericalGTOSet::TMat(const SphericalGTOSet& o) const {
    return CalcMat(*this, OpKE<dcomplex, c3>(), o);
  }
  
  dcomplex* SphericalGTOSet::VMat(dcomplex q, dcomplex x, dcomplex y, dcomplex z,
				const SphericalGTOSet& o) const {
    OpNA<dcomplex, c3> op(q, c3(x, y, z));
    return CalcMat(*this, op, o);
  }
  dcomplex* SphericalGTOSet::XyzMat(int nx, int ny, int nz,
				  const SphericalGTOSet& o) const {
    OpXyz<dcomplex, c3> op(nx, ny, nz);
    return CalcMat(*this, op, o);
  }

}
