#ifndef MOLINT_HPP
#define MOLINT_HPP

#include <vector>
#include <boost/array.hpp>
#include "math_utils.hpp"
#include <map>

namespace l2func {

  using std::vector;

  void IncompleteGamma(int max_m, dcomplex z, dcomplex* res);

  dcomplex coef_d(dcomplex zetap,
		  dcomplex wPk, dcomplex wAk, dcomplex wBk,
		  int nAk, int nBk, int Nk);

  MultArray3<dcomplex> calc_d_coef(int max_ni, int max_nj, int max_n,
				   dcomplex zetaP, dcomplex wPx,
				   dcomplex xi, dcomplex xj, dcomplex* buf);

  void calc_d_coef(int max_ni, int max_nj, int max_n,
		   dcomplex zetaP, dcomplex wPx,
		   dcomplex xi, dcomplex xj,
		   dcomplex* buffer, dcomplex* res);
  
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

  class MatrixSet {
  public:
    int nbasis_i_;
    int nbasis_j_;
    std::map<std::string, dcomplex*> mat_map_;
  public:
    MatrixSet(int nbi, int nbj);
    ~MatrixSet();
    int size_basis_i() const;
    int size_basis_j() const;
    void set(std::string label, dcomplex* ptr);    
    dcomplex* get(std::string label);
    std::map<std::string, dcomplex*>& get_map() { return mat_map_; }
    //np::ndarray get_py(std::string label);
  };  
  class GTOs {
  private:
    //    typedef boost::array<dcomplex, 3> dc3;
    //    typedef boost::array<int, 3>      i3;
    void Normalize();
    void Add(dcomplex _zeta,
	     dcomplex x, dcomplex y, dcomplex z,
	     vector<int> _nx, vector<int> _ny, vector<int> _nz, 
	     vector<vector<dcomplex> > _coef);

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
    vector<int> l_ish;
    vector<vector<int> > m_ish_iprim;
    
    vector<int> offset_ish; // offset_ish[ish] + ibasis gives global index    
    vector<dcomplex>                   coef_ylm_ish;
    vector<vector<vector<dcomplex> > > coef_ish_icont_iprim;

    vector<dcomplex> x_iat;
    vector<dcomplex> y_iat;
    vector<dcomplex> z_iat;
    vector<dcomplex> q_iat;

  public:
    GTOs();
    int size_basis() const {
      int cumsum(0);
      for(int ish = 0; ish < this->size_sh(); ish++) {
	cumsum += this->size_basis_ish(ish);      
      }
      return cumsum;
    }
    int size_sh()  const {
      return zeta_ish.size();    
    }
    int size_basis_ish(int ish) const {
      return coef_ish_icont_iprim[ish].size();    
    }
    int size_prim_ish(int ish) const {
      return nx_ish_iprim[ish].size();
    }
    int size_prim() const {
      int cumsum(0);
      for(int ish = 0; ish < this->size_sh(); ish++) {
	cumsum += this->size_prim_ish(ish);
      }
      return cumsum;
    }
    int size_atom() const {
      return x_iat.size();
    }

    void AddSphericalGTOs(int L, dcomplex x, dcomplex y, dcomplex z, dcomplex _zeta);
    void AddOneSphericalGTO(int L, int M,
			    dcomplex x, dcomplex y, dcomplex z, dcomplex _zeta);
    void AddAtom(dcomplex q, dcomplex x, dcomplex y, dcomplex z);
    //void Normalize();
    //    dcomplex overlap(int ish, int iprim, int jsh, int jprim) const;
    //    dcomplex kinetic(int ish, int iprim, int jsh, int jprim) const;
    //    dcomplex* SMat() const;
    //    dcomplex* SMat2() const;
    //    dcomplex* TMat() const;
    MatrixSet Calc();
    void CalcMat(dcomplex** s, dcomplex** t, dcomplex** dz, dcomplex** v);
    MatrixSet CalcZMatOther(const GTOs&);
    void Show() const;
    void AtR_Ylm(int l, int m, dcomplex* rs, int num_r, dcomplex* cs, dcomplex* res);
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

  /*
  MatrixSet CalcMat(const CartGTOs& a, const MolePot& v, const CartGTOs& b);
*/
}


#endif
