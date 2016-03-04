#ifndef MOLINT_HPP
#define MOLINT_HPP

#include <vector>
#include <boost/array.hpp>
#include "math_utils.hpp"

namespace l2func {

  using std::vector;

  dcomplex* IncompleteGamma(int max_m, dcomplex z);

  dcomplex coef_d(dcomplex zetap,
		  dcomplex wPk, dcomplex wAk, dcomplex wBk,
		  int nAk, int nBk, int Nk);
  void calc_d_coef(int num_ni, int num_nj, int num_n,
		   dcomplex zetaP, dcomplex wPx, dcomplex xi, dcomplex xj,
		   dcomplex** dsx);

  dcomplex gto_overlap(int nAx, int nAy, int nAz, 
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB);

  dcomplex gto_kinetic(int nAx, int nAy, int nAz,
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB);

  dcomplex gto_moment_z(int nAx, int nAy, int nAz,
			dcomplex wAx, dcomplex wAy, dcomplex wAz,
			dcomplex zetaA,
			int nBx, int nBy, int nBz, 
			dcomplex wBx, dcomplex wBy, dcomplex wBz,
			dcomplex zetaB);

  dcomplex gto_nuclear_attraction(int nAx, int nAy, int nAz, 
				  dcomplex wAx, dcomplex wAy, dcomplex wAz,
				  dcomplex zetaA,
				  int nBx, int nBy, int nBz, 
				  dcomplex wBx, dcomplex wBy, dcomplex wBz,
				  dcomplex zetaB,
				  dcomplex wCx, dcomplex wCy, dcomplex wCz);

  /*
  class GTOShell {
    typedef boost::array<int, 3> i3;
  public:
    dcomplex zeta;
    dcomplex x;
    dcomplex y;
    dcomplex z;
    vector<i3> nml;
    vector<dcomplex> 
  };
  */
  class GTOs {
  private:
    typedef boost::array<dcomplex, 3> dc3;
    typedef boost::array<int, 3>      i3;
  public:
    // ish   : index of shell
    // iprim : index of primitive GTO in ish shell.
    // icont : index of contracted basis in ish shell.
    //
    vector<dcomplex>    zeta_ish; 
    vector<dcomplex>    x_ish;
    vector<dcomplex>    y_ish;
    vector<dcomplex>    z_ish;
    vector<vector<int> > nx_ish_iprim;
    vector<vector<int> > ny_ish_iprim;
    vector<vector<int> > nz_ish_iprim;
    
    vector<int> offset_ish; // offset_ish[ish] + ibasis gives global index    
    vector<vector<vector<dcomplex> > > coef_ish_icont_iprim;

  public:
    GTOs();
    void Add(dcomplex _zeta,
	     dcomplex x, dcomplex y, dcomplex z,
	     vector<int> _nx, vector<int> _ny, vector<int> _nz, 
	     vector<vector<dcomplex> > _coef);
    int size_basis() const;
    int size_sh()  const;
    int size_basis_ish(int ish) const;
    int size_prim_ish(int ish) const;
    void AddSphericalGTO(int L, dcomplex x, dcomplex y, dcomplex z, dcomplex _zeta);
    void Normalize();
    dcomplex overlap(int ish, int iprim, int jsh, int jprim) const;
    dcomplex kinetic(int ish, int iprim, int jsh, int jprim) const;
    dcomplex* SMat() const;
    dcomplex* SMat2() const;
    dcomplex* TMat() const;
    void CalcMat(dcomplex** s, dcomplex** t, dcomplex** dz, dcomplex** v);
    // void VMat(dcomplex** v);
  };
  class MolePot {
  public:
    vector<dcomplex> cq;
    vector<dcomplex> cx;
    vector<dcomplex> cy;
    vector<dcomplex> cz;
    void Add(dcomplex _cq, dcomplex _cx, dcomplex _cy, dcomplex _cz) {
      cq.push_back(_cq);
      cx.push_back(_cx);
      cy.push_back(_cy);
      cz.push_back(_cz);
    }
  };
  /*
  class CartGTOs {
  public:
    vector<dcomplex> zeta;
    vector<int>      nx;
    vector<int>      ny;
    vector<int>      nz;
    vector<dcomplex> wx;
    vector<dcomplex> wy;
    vector<dcomplex> wz;

  public:
    CartGTOs();
    ~CartGTOs();
    void AddGTO(dcomplex _zeta, int _nx, int _ny, int _nz,
		dcomplex _wx, dcomplex _wy, dcomplex _wz) {
      zeta.push_back(_zeta);
      nx.push_back(_nx);
      ny.push_back(_ny);
      nz.push_back(_nz);
      wx.push_back(_wx);
      wy.push_back(_wy);
      wz.push_back(_wz);
    }
    int size() const {
      return zeta.size();
    }    
  };
  */
  class MatrixSet {
  public:
    int numa;
    int numb;
    dcomplex* S;
    dcomplex* V;
    dcomplex* T;
    dcomplex* Dz;
  };

  /*
  MatrixSet CalcMat(const CartGTOs& a, const MolePot& v, const CartGTOs& b);
*/
}


#endif
