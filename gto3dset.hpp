#ifndef GTO3DSET_H
#define GTO3DSET_H

#include <vector>
//#include <boost/python.hpp>
//#include <boost/numpy.hpp>
#include "gto3d.hpp"

namespace l2func {

  typedef std::complex<double> dcmplx;
  typedef SphericalGTO<dcmplx, dcmplx> SGTO;

  class SphericalGTOSet {
  private:
    
    std::vector<SGTO*> basis_list_;
  public:
    typedef vector<SGTO*>::const_iterator const_iterator;
    typedef vector<SGTO*>::iterator iterator;

  public:
    SphericalGTOSet();
    ~SphericalGTOSet();
    int size() const;
    const SGTO& basis(int i) const;
    void AddOneBasis(int L, int M, dcmplx x, dcmplx y, dcmplx z, dcmplx zeta);
    void AddBasis(int L, dcmplx x, dcmplx y, dcmplx z, dcmplx zeta);
    dcmplx* SMat(const SphericalGTOSet& o) const;
    dcmplx* TMat(const SphericalGTOSet& o) const;
    dcmplx* VMat(dcmplx q, dcmplx x, dcmplx y, dcmplx z, const SphericalGTOSet& o) const;
    dcmplx* XyzMat(int nx, int ny, int nz, const SphericalGTOSet& o) const;
		   
    // np::ndarray SMatWtihOther(const SphericalGTOSet& o);
  };

}

#endif
