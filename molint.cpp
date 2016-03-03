#include <iostream>
#include "molint.hpp"
#include "macros.hpp"

namespace l2func {

  // ==== Utility ====
  dcomplex dist2(dcomplex dx, dcomplex dy, dcomplex dz) {
    return dx*dx+dy*dy+dz*dz;
  }

  double div(int a, int b) {
    return a*1.0/b;
  }

  // ==== Incomplete Gamma ====
  // K.Ishida J.Comput.Chem. 25, (2004), 739
  // F1 and F2 algorithms 
  // Warning:
  // Eq. (24)(25)(26) are incorrect.

  dcomplex* IncompleteGamma_F1(int max_m, dcomplex z) {

    double x = real(z);
    double y = imag(z);
    double eps(pow(10.0, -10.0));

    if(x < -eps || y < -eps) {
      std::string msg;
      SUB_LOCATION(msg);
      msg += "Re[z] and Im[z] must positive in F1 algorithms";
      throw std::runtime_error(msg);
    }

    double pi(M_PI);
    double z2 = x*x+y*y;
    double az = abs(z);
    double c = sqrt(pi*(az+x)/(8.0*z2));
    double s = sqrt(pi*(az-x)/(8.0*z2));

    int NF_max(50);
    int NF(0);
    int NF_0(10);
    double* anR = new double[NF_max];
    double* anI = new double[NF_max];
    
    double phiR = -x/(2.0*z2)*exp(-x)*cos(y) + y/(2.0*z2)*exp(-x)*sin(y);
    double phiI = +y/(2.0*z2)*exp(-x)*cos(y) + x/(2.0*z2)*exp(-x)*sin(y);
    anR[0] = phiR;
    anI[0] = phiI;
    double sum_anR(anR[0]);
    double sum_anI(anI[0]);
    double delta(pow(10.0, -15.0));
    for(int n = 1; n <= NF_max; n++) {
      anR[n] = -(2.0*n-1) * (x/(2.0*z2)*anR[n-1] + y/(2.0*z2)*anI[n-1]);
      anI[n] = +(2.0*n-1) * (y/(2.0*z2)*anR[n-1] - x/(2.0*z2)*anI[n-1]);      
      sum_anR += anR[n];
      sum_anI += anI[n];
      if(n > NF_0)
	if(c * delta > std::abs(anR[n]) && s*delta > std::abs(anI[n])) {
	  NF = n;
	  break;
	}
    }

    double* fmR = new double[max_m+1];
    double* fmI = new double[max_m+1];
    fmR[0] = +c + sum_anR;
    fmI[0] = -s + sum_anI;
    for(int m = 1; m <= max_m; m++) {
      fmR[m] = (2*m-1.0) * (x/(2.0*z2)*fmR[m-1] + y/(2.0*z2)*fmI[m-1]);
      fmI[m] = (2*m-1.0) * (x/(2.0*z2)*fmI[m-1] - y/(2.0*z2)*fmR[m-1]);
    }

    double* bmR = new double[max_m+1];
    double* bmI = new double[max_m+1];
    bmR[0] = 0.0;
    bmI[0] = 0.0;
    //    bmR[1] = phiR;
    //    bmI[1] = phiI;
    for(int m = 1; m <= max_m; m++) {
      bmR[m] = phiR + (2*m-1.0) * (x/(2.0*z2)*bmR[m-1] + y/(2.0*z2)*bmI[m-1]);
      bmI[m] = phiI + (2*m-1.0) * (x/(2.0*z2)*bmI[m-1] - y/(2.0*z2)*bmR[m-1]);
    }
    //    std::cout << "bmR[2]: " << bmR[2] << std::endl;    

    dcomplex* Fm = new dcomplex[max_m+1];
    for(int m = 0; m <= max_m; m++) 
      Fm[m] = dcomplex(fmR[m]+bmR[m], fmI[m]+bmI[m]);

    delete anR; delete anI; delete fmR; delete fmI; delete bmR; delete bmI;

    return Fm;

  }

  dcomplex* IncompleteGamma_F2(int max_m, dcomplex z) {

    double x = real(z);
    double y = imag(z);
    double eps(pow(10.0, -10.0));
    
    if(x < -eps) {
      std::string msg;
      SUB_LOCATION(msg);
      msg += "Re[z] must positive in F2 algorithms";
      throw std::runtime_error(msg);
    }

    if(y < -eps) {
      dcomplex* conj_res = IncompleteGamma_F1(max_m, dcomplex(x, -y));
      for(int m = 0; m <= max_m; m++)
	conj_res[m] = conj(conj_res[m]);
      return conj_res;
    }

    //    double pi(M_PI);
    //    double z2 = x*x+y*y;

    int NR(47);

    dcomplex* Bn = new dcomplex[NR+1];
    Bn[0] = 1.0;
    Bn[1] = 1.0 + 0.5*z;
    for(int n = 2; n <= NR; n++) 
      Bn[n] = Bn[n-1] + z*z/(4*(2*n-1)*(2*n-3)*1.0)*Bn[n-2];

    dcomplex* Fm = new dcomplex[max_m+1];
    for(int m = 0; m <= max_m; m++) {

      dcomplex* An = new dcomplex[NR+1];
      An[0] = 1.0;
      double t1 = div(2*m+1, 2*m+3);
      double t2 = div(2*m+1, (2*m+3)*(2*m+5));
      double t3 = div((2*m+1)*((2*m+1)*(2*m+1)+44), 60*(2*m+3)*(2*m+5)*(2*m+7));
      An[1] = Bn[1] - t1*z;
      An[2] = Bn[2] - t1*z - t2*z*z;
      An[3] = Bn[3] - t1*z - t2*z*z - t3*z*z*z;
      for(int n = 4; n <= NR; n++) {
	double F1 = div(2*n-2*m-5, 2*(2*n-3)*(2*n+2*m+1));
	double F2 = div(1, 4*(2*n-1)*(2*n-3));
	double F3 = -F1/(4*(2*n-3)*(2*n-5));
	double E = -F1;
	An[n] = (1.0+F1*z)*An[n-1] + (E + F2*z)*z*An[n-2] + F3*z*z*z*An[n-3];
      }
      Fm[m] = div(1, 2*m+1) * An[NR]/Bn[NR];
    }

    return Fm;
    
  }

  dcomplex* IncompleteGamma(int max_m, dcomplex z) {
    
    double x = real(z);
    double y = imag(z);
    double eps(pow(10.0, -10.0));

    if(x < -eps) {
      std::string msg;
      SUB_LOCATION(msg);
      msg += "negative Re[z] is not implemented yet.";
      throw std::runtime_error(msg);
    } 
    
    dcomplex* res;
    if(x > -eps && y > -eps) {
      if(x < 21.0 && x+y < 37.0) 
	res = IncompleteGamma_F2(max_m, z);
      else
	res = IncompleteGamma_F1(max_m, z);
    } else {
      res = IncompleteGamma(max_m, dcomplex(x, -y));
      for(int m = 0; m <= max_m; m++)
	res[m] = conj(res[m]);
    }

    return res;

  }

  dcomplex coef_d(dcomplex zetap,
		  dcomplex wPk, dcomplex wAk, dcomplex wBk,
		  int nAk, int nBk, int Nk) {
    
    //    std::cout << nAk << nBk << Nk <<std::endl;
		  
    if(nAk == 0 && nBk == 0 && Nk == 0)
      return 1.0;

    if(Nk < 0 || Nk > nAk + nBk) 
      return 0.0;

    if(nAk > 0) 
      return 1.0/(2.0*zetap) * coef_d(zetap, wPk, wAk, wBk, nAk-1, nBk, Nk-1) +
	(wPk - wAk)          * coef_d(zetap, wPk, wAk, wBk, nAk-1, nBk, Nk)   +
	(Nk + 1.0)           * coef_d(zetap, wPk, wAk, wBk, nAk-1, nBk, Nk+1);
    else 
      return 1.0/(2.0*zetap) * coef_d(zetap, wPk, wAk, wBk, nAk, nBk-1, Nk-1) +
	(wPk - wBk)          * coef_d(zetap, wPk, wAk, wBk, nAk, nBk-1, Nk)   +
	(Nk + 1.0)           * coef_d(zetap, wPk, wAk, wBk, nAk, nBk-1, Nk+1);
    
  }

  dcomplex gto_overlap(int nAx, int nAy, int nAz, 
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB) {

    dcomplex zetaP = zetaA + zetaB;
    dcomplex wPx    = (zetaA * wAx + zetaB * wBx) / zetaP;
    dcomplex wPy    = (zetaA * wAy + zetaB * wBy) / zetaP;
    dcomplex wPz    = (zetaA * wAz + zetaB * wBz) / zetaP;
		       
    dcomplex eAB = exp(-zetaA*zetaB/zetaP*dist2(wAx-wBx, wAy-wBy, wAz-wBz));

    return (eAB * pow(M_PI/zetaP, 1.5) *
	    coef_d(zetaP, wPx, wAx, wBx, nAx, nBx, 0) *
	    coef_d(zetaP, wPy, wAy, wBy, nAy, nBy, 0) *
	    coef_d(zetaP, wPz, wAz, wBz, nAz, nBz, 0));

  }

  dcomplex gto_kinetic(int nAx, int nAy, int nAz,
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB) {
    dcomplex res(0.0);
    dcomplex s000 = gto_overlap(nAx, nAy, nAz, wAx, wAy, wAz, zetaA,
				nBx, nBy, nBz, wBx, wBy, wBz, zetaB);
    dcomplex sp00 = gto_overlap(nAx,   nAy, nAz, wAx, wAy, wAz, zetaA,
				nBx+2, nBy, nBz, wBx, wBy, wBz, zetaB);
    dcomplex s0p0 = gto_overlap(nAx, nAy, nAz, wAx, wAy, wAz, zetaA,
				nBx, nBy+2, nBz, wBx, wBy, wBz, zetaB);
    dcomplex s00p = gto_overlap(nAx, nAy, nAz, wAx, wAy, wAz, zetaA,
				nBx, nBy, nBz+2, wBx, wBy, wBz, zetaB);
    res += -2.0*zetaB*(2.0*nBx+2*nBy+2*nBz+3.0) * s000;
    res += 4.0*zetaB*zetaB*(sp00 + s0p0 + s00p);
    if(nBx > 1)
      res += 1.0*nBx*(nBx-1)* gto_overlap(nAx,   nAy, nAz, wAx, wAy, wAz, zetaA,
					  nBx-2, nBy, nBz, wBx, wBy, wBz, zetaB);
    
    if(nBy > 1)
      res += 1.0*nBy*(nBy-1)* gto_overlap(nAx, nAy,   nAz, wAx, wAy, wAz, zetaA,
					  nBx, nBy-2, nBz, wBx, wBy, wBz, zetaB);     

    if(nBz > 1)
      res += 1.0*nBz*(nBz-1)* gto_overlap(nAx, nAy, nAz,   wAx, wAy, wAz, zetaA,
					  nBx, nBy, nBz-2, wBx, wBy, wBz, zetaB);
    

    return -0.5*res;
  }

  dcomplex gto_moment_z(int nAx, int nAy, int nAz,
			dcomplex wAx, dcomplex wAy, dcomplex wAz,
			dcomplex zetaA,
			int nBx, int nBy, int nBz, 
			dcomplex wBx, dcomplex wBy, dcomplex wBz,
			dcomplex zetaB) {
    
    return (gto_overlap(nAx, nAy, nAz+1, wAx, wAy, wAz, zetaA,
			nBx, nBy, nBz,   wBx, wBy, wBz, zetaB)
	    + wAz *
	    gto_overlap(nAx, nAy, nAz, wAx, wAy, wAz, zetaA,
			nBx, nBy, nBz, wBx, wBy, wBz, zetaB));
    
  }
  
  dcomplex coef_R(dcomplex zetaP,
		  dcomplex wPx, dcomplex wPy, dcomplex wPz,
		  dcomplex cx,  dcomplex cy,  dcomplex cz,
		  int mx, int my, int mz, int j, dcomplex* Fjs) {
    //    std::cout << mx << my << mz << j << std::endl;

    if(mx == 0 && my == 0 && mz == 0) {
      //      dcomplex arg = zetaP * dist2(wPx-cx, wPy-cy, wPz-cz);
      //      dcomplex* Fjs = IncompleteGamma(j, arg);
      dcomplex Fj = Fjs[j];
      return pow(-2.0*zetaP, j) * Fj;
    }

    if(mx > 0) {
      dcomplex res(0.0);
      if(mx > 1) 
	res += (mx-1.0) * coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx-2, my, mz, j+1, Fjs);
      res += (wPx-cx) *   coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx-1, my, mz, j+1, Fjs);
      return res;
    }
    if(my > 0) {
      dcomplex res(0.0);
      if(my > 1) 
	res += (my-1.0) * coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx, my-2, mz, j+1, Fjs);
      res += (wPy-cy) *   coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx, my-1, mz, j+1, Fjs);
      return res;
    }
    if(mz > 0) {
      dcomplex res(0.0);
      if(mz > 1) 
	res += (mz-1.0) * coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx, my, mz-2, j+1, Fjs);
      res += (wPz-cz) *   coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx, my, mz-1, j+1, Fjs);
      return res;
    }
    std::string msg;
    SUB_LOCATION(msg);
    msg += " one of mx, my, mz is negative integer.";
    throw std::runtime_error(msg);    
    return 0.0;    
  }
  dcomplex gto_nuclear_attraction(int nAx, int nAy, int nAz, 
				  dcomplex wAx, dcomplex wAy, dcomplex wAz,
				  dcomplex zetaA,
				  int nBx, int nBy, int nBz, 
				  dcomplex wBx, dcomplex wBy, dcomplex wBz,
				  dcomplex zetaB,
				  dcomplex wCx, dcomplex wCy, dcomplex wCz) {

    dcomplex zetaP = zetaA + zetaB;
    dcomplex wPx    = (zetaA * wAx + zetaB * wBx) / zetaP;
    dcomplex wPy    = (zetaA * wAy + zetaB * wBy) / zetaP;
    dcomplex wPz    = (zetaA * wAz + zetaB * wBz) / zetaP;
		       
    dcomplex eAB = exp(-zetaA*zetaB/zetaP*dist2(wAx-wBx, wAy-wBy, wAz-wBz));

    dcomplex* dxs = new dcomplex[nAx+nBx+1];
    for(int nx = 0; nx <= nAx+nBx; nx++)
      dxs[nx] = coef_d(zetaP, wPx, wAx, wBx, nAx, nBx, nx);

    dcomplex* dys = new dcomplex[nAy+nBy+1];
    for(int ny = 0; ny <= nAy+nBy; ny++)
      dys[ny] = coef_d(zetaP, wPy, wAy, wBy, nAy, nBy, ny);

    dcomplex* dzs = new dcomplex[nAz+nBz+1];
    for(int nz = 0; nz <= nAz+nBz; nz++) 
      dzs[nz] = coef_d(zetaP, wPz, wAz, wBz, nAz, nBz, nz);
    
    dcomplex* Fjs = IncompleteGamma(nAx+nAy+nAz+nBx+nBy+nBz,
				    zetaP * dist2(wPx-wCx, wPy-wCy, wPz-wCz));

    dcomplex cumsum(0);    
    for(int nx = 0; nx <= nAx+nBx; nx++)
      for(int ny = 0; ny <= nAy+nBy; ny++)
	for(int nz = 0; nz <= nAz+nBz; nz++) {
	  cumsum += (-2.0*M_PI/zetaP *
		     dxs[nx] * dys[ny] * dzs[nz] * 
		     coef_R(zetaP, wPx, wPy, wPz, wCx, wCy, wCz,
			    nx, ny, nz, 0, Fjs));
	}
    delete dxs;
    delete dys;
    delete dzs;

    /*
    dcomplex cumsum(0);    
    for(int nx = 0; nx <= nAx+nBx; nx++)
      for(int ny = 0; ny <= nAy+nBy; ny++)
	for(int nz = 0; nz <= nAz+nBz; nz++) {
	  cumsum += (-2.0*M_PI/zetaP *
		     coef_d(zetaP, wPx, wAx, wBx, nAx, nBx, nx) * 
		     coef_d(zetaP, wPy, wAy, wBy, nAy, nBy, ny) * 
		     coef_d(zetaP, wPz, wAz, wBz, nAz, nBz, nz) * 
		     coef_R(zetaP, wPx, wPy, wPz, wCx, wCy, wCz,
			    nx, ny, nz, 0));
	}
	*/
    return eAB * cumsum;
  //  dcomplex coef_d(int max_mA, int max_mB, dcomplex* ds, dcomplex zetap,
  }
  GTOs::GTOs() {
    offset_ish.push_back(0);
  }
  void GTOs::Add(dcomplex _zeta,
		 dcomplex x, dcomplex y, dcomplex z,
		 vector<int> _nx, vector<int> _ny, vector<int> _nz, 
		 vector<vector<dcomplex> > _coef) {
    // new shell index
    int ish = this->size_sh();

    // offset for the new shell
    int offset0 = offset_ish[ish];

    int basis_size0 = _coef.size();
    offset_ish.push_back(offset0 + basis_size0);

    zeta_ish.push_back(_zeta);
    x_ish.push_back(x);
    y_ish.push_back(y);
    z_ish.push_back(z);
    nx_ish_iprim.push_back(_nx);
    ny_ish_iprim.push_back(_ny);
    nz_ish_iprim.push_back(_nz);
    coef_ish_icont_iprim.push_back(_coef);
    
  }
  void GTOs::AddSphericalGTO(int L, dcomplex x, dcomplex y, dcomplex z,
			     dcomplex _zeta) {

    if(L == 0) {
      vector<int> n0; n0.push_back(0);
      vector<dcomplex> c1; c1.push_back(1.0);      
      vector<vector<dcomplex> > c; c.push_back(c1);      
      this->Add(_zeta, x, y, z, n0, n0, n0, c);     
    }
    else if(L == 1) {
      vector<int> nx(3), ny(3), nz(3);
      nx[0] = 1; ny[0] = 0; nz[0] = 0;
      nx[1] = 0; ny[1] = 1; nz[1] = 0;
      nx[2] = 0; ny[2] = 0; nz[2] = 1;
      vector<vector<dcomplex> > cs(3, vector<dcomplex>(3));
      cs[0][0] = 0.0; cs[0][1] = 1.0; cs[0][2] = 0.0; // M =-1
      cs[1][0] = 0.0; cs[1][1] = 0.0; cs[1][2] = 1.0; // M = 0
      cs[2][0] = 1.0; cs[2][1] = 0.0; cs[2][2] = 0.0; // M = 1

      this->Add(_zeta, x, y, z, nx, ny, nz, cs);     

    }
    else if(L == 2) {
      vector<int> nx(6), ny(6), nz(6);
      nx[0] = 2; ny[0] = 0; nz[0] = 0; 
      nx[1] = 0; ny[1] = 2; nz[1] = 0; 
      nx[2] = 0; ny[2] = 0; nz[2] = 2; 
      nx[3] = 1; ny[3] = 1; nz[3] = 0; 
      nx[4] = 1; ny[4] = 0; nz[4] = 1;
      nx[5] = 0; ny[5] = 1; nz[5] = 1;
      
      vector<vector<dcomplex> > cs(5, vector<dcomplex>(6));
      // M = -2
      cs[0][0] = 0.0; cs[0][1] = 0.0; cs[0][2] = 0.0;
      cs[0][3] = 1.0; cs[0][4] = 0.0; cs[0][5] = 0.0; 

      // M = -1
      cs[1][0] = 0.0; cs[1][1] = 0.0; cs[1][2] = 0.0;
      cs[1][3] = 0.0; cs[1][4] = 0.0; cs[1][5] = 1.0; 

      // M = 0
      cs[2][0] = -1.0; cs[2][1] = -1.0; cs[2][2] = 2.0;
      cs[2][3] = +0.0; cs[2][4] = +0.0; cs[2][5] = 0.0; 

      // M = 1
      cs[3][0] = 0.0; cs[3][1] = 0.0; cs[3][2] = 0.0;
      cs[3][3] = 0.0; cs[3][4] = 1.0; cs[3][5] = 0.0; 

      // M = 2
      cs[4][0] = 1.0; cs[4][1] =-1.0; cs[4][2] = 0.0;
      cs[4][3] = 0.0; cs[4][4] = 0.0; cs[4][5] = 0.0; 

      this->Add(_zeta, x, y, z, nx, ny, nz, cs);     
    } else {
      std::string msg;
      SUB_LOCATION(msg);      
      msg += "Not implemented yet for this L";
      throw std::runtime_error(msg);      
    }

  }

  int GTOs::size_basis() const {
    int cumsum(0);
    for(int ish = 0; ish < this->size_sh(); ish++) {
      cumsum += this->size_basis_ish(ish);      
    }
    return cumsum;
  }
  int GTOs::size_sh()  const {
    return zeta_ish.size();    
  }
  int GTOs::size_basis_ish(int ish) const {
    return coef_ish_icont_iprim[ish].size();    
  }
  int GTOs::size_prim_ish(int ish) const {
    return nx_ish_iprim[ish].size();
  }
  void GTOs::Normalize() {

    for(int ish = 0; ish < this->size_sh(); ish++) {
      for(int ibasis = 0; ibasis < this->size_basis_ish(ish); ibasis++) {
	dcomplex cumsum(0.0);
	for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++) {
	  for(int jprim = 0; jprim < this->size_prim_ish(ish); jprim++) {
	    dcomplex s = gto_overlap(nx_ish_iprim[ish][iprim],
				     ny_ish_iprim[ish][iprim],
				     nz_ish_iprim[ish][iprim],
				     x_ish[ish], y_ish[ish], z_ish[ish],
				     zeta_ish[ish],
				     nx_ish_iprim[ish][jprim],
				     ny_ish_iprim[ish][jprim],
				     nz_ish_iprim[ish][jprim],
				     x_ish[ish], y_ish[ish], z_ish[ish],
				     zeta_ish[ish]);
	    s *= coef_ish_icont_iprim[ish][ibasis][iprim];
	    s *= coef_ish_icont_iprim[ish][ibasis][jprim];
	    cumsum += s;
	  }
	}
	dcomplex scale = 1.0/sqrt(cumsum);	
	for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++) 
	  coef_ish_icont_iprim[ish][ibasis][iprim] *= scale;	
      }
    }    
  }
  dcomplex* GTOs::SMat() const {

    int nb = this->size_basis();    
    dcomplex* S = new dcomplex[nb*nb];

    int num_sh = this->size_sh();    
    int idx(0);
    for(int ish = 0; ish < num_sh; ish++) {
      for(int jsh = 0; jsh < num_sh; jsh++) {
	int n_basis_ish = this->size_basis_ish(ish);
	int n_basis_jsh = this->size_basis_ish(jsh);
	int n_prim_i = this->size_prim_ish(ish);
	int n_prim_j = this->size_prim_ish(jsh);
	for(int icont = 0; icont < n_basis_ish; icont++)
	  for(int jcont = 0; jcont < n_basis_jsh; jcont++) {
	    dcomplex cumsum(0.0);
	    for(int iprim = 0; iprim < n_prim_i; iprim++) {
	      for(int jprim = 0; jprim < n_prim_j; jprim++) {
		dcomplex s;
		s =  gto_overlap(nx_ish_iprim[ish][iprim],
				 ny_ish_iprim[ish][iprim],
				 nz_ish_iprim[ish][iprim],
				 x_ish[ish], y_ish[ish], z_ish[ish],
				 zeta_ish[ish],
				 nx_ish_iprim[jsh][jprim],
				 ny_ish_iprim[jsh][jprim],
				 nz_ish_iprim[jsh][jprim],
				 x_ish[jsh], y_ish[jsh], z_ish[jsh],
				 zeta_ish[jsh]);
		s *= coef_ish_icont_iprim[ish][icont][iprim];
		s *= coef_ish_icont_iprim[jsh][jcont][jprim];
		cumsum += s;		
	      }
	    }
	    S[idx] = cumsum;
	    idx++;	    
	  }
      }
    }    
    return S;    
  }
  dcomplex GTOs::overlap(int ish, int iprim, int jsh, int jprim) const {

    return gto_overlap(nx_ish_iprim[ish][iprim],
		       ny_ish_iprim[ish][iprim],
		       nz_ish_iprim[ish][iprim],
		       x_ish[ish], y_ish[ish], z_ish[ish],
		       zeta_ish[ish],
		       nx_ish_iprim[jsh][jprim],
		       ny_ish_iprim[jsh][jprim],
		       nz_ish_iprim[jsh][jprim],
		       x_ish[jsh], y_ish[jsh], z_ish[jsh],
		       zeta_ish[jsh]);

  }
  dcomplex GTOs::kinetic(int ish, int iprim, int jsh, int jprim) const {

    return gto_kinetic(nx_ish_iprim[ish][iprim],
		       ny_ish_iprim[ish][iprim],
		       nz_ish_iprim[ish][iprim],
		       x_ish[ish], y_ish[ish], z_ish[ish],
		       zeta_ish[ish],
		       nx_ish_iprim[jsh][jprim],
		       ny_ish_iprim[jsh][jprim],
		       nz_ish_iprim[jsh][jprim],
		       x_ish[jsh], y_ish[jsh], z_ish[jsh],
		       zeta_ish[jsh]);

  }
  dcomplex* GTOs::SMat2() const {
    
    int nb = this->size_basis();    
    dcomplex* S = new dcomplex[nb*nb];
    int idx(0);
    for(int ish = 0; ish < this->size_sh(); ish++) {
      for(int jsh = 0; jsh < this->size_sh(); jsh++) {
	int npi = this->size_prim_ish(ish);
	int npj = this->size_prim_ish(jsh);
	int idx_prim(0);	
	dcomplex* smat_prim = new dcomplex[npi*npj];	
	for(int iprim = 0; iprim < npi; iprim++) {
	  for(int jprim = 0; jprim < npj; jprim++) {
	    smat_prim[idx_prim] = this->overlap(ish, iprim, jsh, jprim);	    
	    idx_prim++;	    
	  }	      
	}
	
	for(int ibasis = 0; ibasis < this->size_basis_ish(ish); ibasis++) {
	  for(int jbasis = 0; jbasis < this->size_basis_ish(jsh); jbasis++) {
	    dcomplex s(0.0);
	    int idx_prim0(0);	    
	    for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++) {
	      for(int jprim = 0; jprim < this->size_prim_ish(jsh); jprim++) {
		s += (coef_ish_icont_iprim[ish][ibasis][iprim] * 
		      coef_ish_icont_iprim[jsh][jbasis][jprim] * 
		      smat_prim[idx_prim0]);
		idx_prim0++;		
	      }
	    }
	    S[idx] = s;
	    idx++;	    
	  }
	}
	delete smat_prim;
      }
    }
    return S;    
  }
  dcomplex* GTOs::TMat() const {
    int nb = this->size_basis();    
    dcomplex* T = new dcomplex[nb*nb];
    int idx(0);
    for(int ish = 0; ish < this->size_sh(); ish++) {
      for(int jsh = 0; jsh < this->size_sh(); jsh++) {
	int npi = this->size_prim_ish(ish);
	int npj = this->size_prim_ish(jsh);
	int idx_prim(0);	
	dcomplex* smat_prim = new dcomplex[npi*npj];	
	for(int iprim = 0; iprim < npi; iprim++) {
	  for(int jprim = 0; jprim < npj; jprim++) {
	    smat_prim[idx_prim] = this->kinetic(ish, iprim, jsh, jprim);	    
	    idx_prim++;	    
	  }	      
	}
	
	for(int ibasis = 0; ibasis < this->size_basis_ish(ish); ibasis++) {
	  for(int jbasis = 0; jbasis < this->size_basis_ish(jsh); jbasis++) {
	    dcomplex s(0.0);
	    int idx_prim0(0);	    
	    for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++) {
	      for(int jprim = 0; jprim < this->size_prim_ish(jsh); jprim++) {
		s += (coef_ish_icont_iprim[ish][ibasis][iprim] * 
		      coef_ish_icont_iprim[jsh][jbasis][jprim] * 
		      smat_prim[idx_prim0]);
		idx_prim0++;		
	      }
	    }
	    T[idx] = s;
	    idx++;	    
	  }
	}
	delete smat_prim;
      }
    }
    return T;    
  }
  void GTOs::CalcMat(dcomplex** s, dcomplex** t, dcomplex** dz, dcomplex** v) {
    
    int nb = this->size_basis();    
    dcomplex* S = new dcomplex[nb*nb];
    dcomplex* T = new dcomplex[nb*nb];
    dcomplex* Dz = new dcomplex[nb*nb];
    dcomplex* V = new dcomplex[nb*nb];    
    for(int idx=0; idx<nb*nb; idx++) {
      S[idx] = 7.7; T[idx] = 7.7; Dz[idx] = 7.7; V[idx] = 7.7;
    }
    for(int ish = 0; ish < this->size_sh(); ish++) {
      for(int jsh = 0; jsh < this->size_sh(); jsh++) {
	dcomplex zetai = zeta_ish[ish]; dcomplex zetaj = zeta_ish[jsh];	
	dcomplex zetaP = zetai + zetaj;
	dcomplex xi(x_ish[ish]); dcomplex yi(y_ish[ish]); dcomplex zi(z_ish[ish]);
	dcomplex xj(x_ish[jsh]); dcomplex yj(y_ish[jsh]); dcomplex zj(z_ish[jsh]);
	dcomplex wPx    = (zetai*xi + zetaj*xj)/zetaP;
	dcomplex wPy    = (zetai*yi + zetaj*yj)/zetaP;
	dcomplex wPz    = (zetai*zi + zetaj*zj)/zetaP;	
	dcomplex d2 = dist2(xi-xj,yi-yj,zi-zj);	
	dcomplex eAB = exp(-zetai*zetaj/zetaP*d2);
	dcomplex ce = eAB * pow(M_PI/zetaP, 1.5);

	int npi = this->size_prim_ish(ish);
	int npj = this->size_prim_ish(jsh);
	int idx_prim(0);	
	dcomplex* smat_prim = new dcomplex[npi*npj];
	dcomplex* zmat_prim = new dcomplex[npi*npj];
	dcomplex* tmat_prim = new dcomplex[npi*npj];
	dcomplex* vmat_prim = new dcomplex[npi*npj];
	
	for(int iprim = 0; iprim < npi; iprim++) {
	  for(int jprim = 0; jprim < npj; jprim++) {
	    int nxi = nx_ish_iprim[ish][iprim]; int nxj = nx_ish_iprim[jsh][jprim];
	    int nyi = ny_ish_iprim[ish][iprim]; int nyj = ny_ish_iprim[jsh][jprim];
	    int nzi = nz_ish_iprim[ish][iprim]; int nzj = nz_ish_iprim[jsh][jprim];

	    // ---- S mat ----
	    dcomplex dx00 = coef_d(zetaP, wPx, xi, xj, nxi, nxj, 0);
	    dcomplex dy00 = coef_d(zetaP, wPy, yi, yj, nyi, nyj, 0);
	    dcomplex dz00 = coef_d(zetaP, wPz, zi, zj, nzi, nzj, 0);

	    smat_prim[idx_prim] = ce * dx00 * dy00 * dz00;

	    // ---- z mat ----
	    dcomplex dz10 = coef_d(zetaP, wPz, zi, zj, nzi+1, nzj, 0);
	    zmat_prim[idx_prim] = ce*(dx00*dy00*dz10 + zi*dx00*dy00*dz00);

	    // ---- t mat ----
	    dcomplex dx02 = coef_d(zetaP, wPx, xi, xj, nxi, nxj+2, 0);
	    dcomplex dy02 = coef_d(zetaP, wPy, yi, yj, nyi, nyj+2, 0);
	    dcomplex dz02 = coef_d(zetaP, wPz, zi, zj, nzi, nzj+2, 0);
	    dcomplex tmp(0.0);
	    tmp += -2.0*zetaj * (2*nxj+2*nyj+2*nzj+3.0) * dx00*dy00*dz00;
	    tmp += 4.0*zetaj*zetaj*(dx02*dy00*dz00+dx00*dy02*dz00+dx00*dy00*dz02);
	    if(nxj > 1) {
	      dcomplex dx = coef_d(zetaP, wPx, xi, xj, nxi, nxj-2, 0);
	      tmp += 1.0*nxj*(nxj-1) * dx * dy00 * dz00;
	    }
	    if(nyj > 1) {
	      dcomplex dy = coef_d(zetaP, wPy, yi, yj, nyi, nyj-2, 0);
	      tmp += 1.0*nyj*(nyj-1) * dx00 * dy * dz00;
	    }
	    if(nzj > 1) {
	      dcomplex dz = coef_d(zetaP, wPz, zi, zj, nzi, nzj-2, 0);
	      tmp += 1.0*nzj*(nzj-1) * dx00 * dy00 * dz;
	    }
	    tmat_prim[idx_prim] = -0.5 * ce * tmp;

	    // ---- v mat ----
	    dcomplex* dxs = new dcomplex[nxi + nxj + 1];
	    dcomplex* dys = new dcomplex[nyi + nyj + 1];
	    dcomplex* dzs = new dcomplex[nzi + nzj + 1];
	    dcomplex cx(0.0); dcomplex cy(0.0); dcomplex cz(0.0);
	    for(int nx = 0; nx <= nxi+nxj; nx++)
	      dxs[nx] = coef_d(zetaP, wPx, xi, xj, nxi, nxj, nx);
	    for(int ny = 0; ny <= nyi+nyj; ny++)
	      dys[ny] = coef_d(zetaP, wPy, yi, yj, nyi, nyj, ny);
	    for(int nz = 0; nz <= nzi+nzj; nz++)
	      dzs[nz] = coef_d(zetaP, wPz, zi, zj, nzi, nzj, nz);
	    dcomplex* Fjs = IncompleteGamma(nxi+nxj+nyi+nyj+nzi+nzj,
					    zetaP * dist2(wPx-cx,wPy-cy,wPz-cz));
	    dcomplex cumsum(0);
	    for(int nx = 0; nx <= nxi + nxj; nx++)
	      for(int ny = 0; ny <= nyi + nyj; ny++)
		for(int nz = 0; nz <= nzi + nzj; nz++) {
		  cumsum += (dxs[nx] * dys[ny] * dzs[nz] *
			     coef_R(zetaP, wPx, wPy, wPz, 0.0, 0.0, 0.0,
				    nx, ny, nz, 0, Fjs));
		}
	    vmat_prim[idx_prim] = -2.0*M_PI/zetaP * eAB * cumsum;
	    delete dxs; delete dys; delete dzs;
	    delete Fjs;

	    // ---- index ----
	    idx_prim++;	    
	  }	      
	}
	for(int ibasis = 0; ibasis < this->size_basis_ish(ish); ibasis++) {
	  for(int jbasis = 0; jbasis < this->size_basis_ish(jsh); jbasis++) {
	    dcomplex s(0.0); dcomplex dz(0.0); dcomplex t(0.0);	dcomplex v(0.0);
	    int idx_prim0(0);	    
	    for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++) {
	      for(int jprim = 0; jprim < this->size_prim_ish(jsh); jprim++) {
		s += (coef_ish_icont_iprim[ish][ibasis][iprim] * 
		      coef_ish_icont_iprim[jsh][jbasis][jprim] * 
		      smat_prim[idx_prim0]);
		t += (coef_ish_icont_iprim[ish][ibasis][iprim] * 
		      coef_ish_icont_iprim[jsh][jbasis][jprim] * 
		      tmat_prim[idx_prim0]);
		dz += (coef_ish_icont_iprim[ish][ibasis][iprim] * 
		       coef_ish_icont_iprim[jsh][jbasis][jprim] * 
		       zmat_prim[idx_prim0]);
		v += (coef_ish_icont_iprim[ish][ibasis][iprim] * 
		      coef_ish_icont_iprim[jsh][jbasis][jprim] * 
		      vmat_prim[idx_prim0]);
		idx_prim0++;
	      }
	    }
	    int i = offset_ish[ish] + ibasis;
	    int j = offset_ish[jsh] + jbasis;
	    int idx = i+j*nb;	
	    S[idx] = s; Dz[idx] = dz; T[idx] = t; V[idx] = v;	    
	  }
	}
	delete smat_prim;
	delete tmat_prim;
	delete zmat_prim;
	delete vmat_prim;
      }
    }
    *s = S; *dz = Dz; *t = T; *v = V;    
  }
  /*
  MatrixSet CalcMat(const CartGTOs& a, const MolePot& v, const CartGTOs& b) {
    
    int numA = a.size();
    int numB = b.size();

    for(int A = 0; A < numA; A++)
      for(int B = 0; B < numB; B++) {

	


      }
  }
  */
}
