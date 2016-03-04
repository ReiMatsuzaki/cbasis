#include <iostream>
#include "molint.hpp"
#include "macros.hpp"


namespace l2func {

  typedef MultArray3<dcomplex> A3dc;

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

  void IncompleteGamma_F1(int max_m, dcomplex z, dcomplex* res_list) {

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

    for(int m = 0; m <= max_m; m++) 
      res_list[m] = dcomplex(fmR[m]+bmR[m], fmI[m]+bmI[m]);

    delete anR; delete anI; delete fmR; delete fmI; delete bmR; delete bmI;
  }

  void IncompleteGamma_F2(int max_m, dcomplex z, dcomplex* res) {

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
      IncompleteGamma_F2(max_m, dcomplex(x, -y), res);
      for(int m = 0; m <= max_m; m++)
	res[m] = conj(res[m]);
    }

    //    double pi(M_PI);
    //    double z2 = x*x+y*y;

    int NR(47);

    dcomplex* Bn = new dcomplex[NR+1];
    Bn[0] = 1.0;
    Bn[1] = 1.0 + 0.5*z;
    for(int n = 2; n <= NR; n++) 
      Bn[n] = Bn[n-1] + z*z/(4*(2*n-1)*(2*n-3)*1.0)*Bn[n-2];

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
      res[m] = div(1, 2*m+1) * An[NR]/Bn[NR];
    }
  }

  void IncompleteGamma(int max_m, dcomplex z, dcomplex* res_list) {
    
    double x = real(z);
    double y = imag(z);
    double eps(pow(10.0, -10.0));

    if(x < -eps) {
      std::string msg;
      SUB_LOCATION(msg);
      msg += "negative Re[z] is not implemented yet.";
      throw std::runtime_error(msg);
    } 
    
    if(x > -eps && y > -eps) {
      if(x < 21.0 && x+y < 37.0) 
	IncompleteGamma_F2(max_m, z, res_list);
      else
	IncompleteGamma_F1(max_m, z, res_list);
    } else {
      IncompleteGamma(max_m, dcomplex(x, -y), res_list);
      for(int m = 0; m <= max_m; m++)
	res_list[m] = conj(res_list[m]);
    }
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

  A3dc calc_d_coef(int max_ni, int max_nj, int max_n,
		   dcomplex zetaP, dcomplex wPx, dcomplex xi, dcomplex xj,
		   dcomplex* buffer) {

    A3dc data(buffer,
	      0, max_ni,
	      0, max_nj,
	      -max_ni-max_nj, max_ni+max_nj+max_n);

    data.set(0, 0, 0, 1.0);
    for(int n = -max_ni -max_nj; n < 0; n++) 
      data.set(0, 0, n, 0.0);
    for(int n = 1; n <= max_ni+max_nj+max_n; n++) 
      data.set(0, 0, n, 0.0);

    for(int ni_nj = 1; ni_nj <= max_ni + max_nj; ni_nj++) {
      for(int ni = 0; ni <= std::min(max_ni, ni_nj); ni++) {
	int nj = ni_nj - ni;
	if(nj <= max_nj) {
	  for(int n = -(max_ni + max_nj - ni_nj);
	      n <= max_ni + max_nj + max_n - ni_nj; n++) {
	    if(ni > 0) {
	      dcomplex v = (1.0/(2.0*zetaP) * data.get_safe(ni-1, nj, n-1) +
			    (wPx - xi)      * data.get_safe(ni-1, nj, n) + 
			    (n + 1.0)      * data.get_safe(ni-1, nj, n+1));
	      data.set(ni, nj, n, v);
	    } else if(nj > 0) {
	      dcomplex v = (1.0/(2.0*zetaP) * data.get_safe(ni, nj-1, n-1) +
			    (wPx - xj)      * data.get_safe(ni, nj-1, n) + 
			    (n + 1.0)      * data.get_safe(ni, nj-1, n+1));
	      data.set(ni, nj, n, v);	    
	    }
	  }
	}
      }
    }

    return data;
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
    
    dcomplex* Fjs = new dcomplex[nAx+nAy+nAz+nBx+nBy+nBz+1];
    IncompleteGamma(nAx+nAy+nAz+nBx+nBy+nBz,
		    zetaP * dist2(wPx-wCx, wPy-wCy, wPz-wCz), Fjs);

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
  void GTOs::AddAtom(dcomplex q, dcomplex x, dcomplex y, dcomplex z) {
    x_iat.push_back(x);
    y_iat.push_back(y);
    z_iat.push_back(z);
    q_iat.push_back(q);
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
    int np = this->size_prim();
    dcomplex* S = new dcomplex[nb*nb];
    dcomplex* T = new dcomplex[nb*nb];
    dcomplex* Dz = new dcomplex[nb*nb];
    dcomplex* V = new dcomplex[nb*nb];    
    dcomplex* smat_prim = new dcomplex[np*np];
    dcomplex* zmat_prim = new dcomplex[np*np];    
    dcomplex* vmat_prim = new dcomplex[np*np];
    dcomplex* tmat_prim = new dcomplex[np*np];
    dcomplex* dsx_buff = new dcomplex[5*7*10];
    dcomplex* dsy_buff = new dcomplex[5*7*10];
    dcomplex* dsz_buff = new dcomplex[5*7*10];
    dcomplex** Fjs_iat = new dcomplex*[this->size_atom()];
    
    for(int idx=0; idx<nb*nb; idx++) {
      S[idx] = 7.7; T[idx] = 7.7; Dz[idx] = 7.7; V[idx] = 7.7;
    }
    
    int maxn(0);
    int* maxn_ish = new int[this->size_sh()];
    for(int ish = 0; ish < this->size_sh(); ish++) {
      int sum_ni(0);
      for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++) {
	  int ni = (nx_ish_iprim[ish][iprim] +
		    ny_ish_iprim[ish][iprim] +
		    nz_ish_iprim[ish][iprim]); 
	  sum_ni = sum_ni < ni ? ni : sum_ni;
      }
      maxn_ish[ish] = sum_ni;
      if(maxn < sum_ni)
	maxn = sum_ni;
    }
    for(int iat = 0; iat < this->size_atom(); iat++) {
      Fjs_iat[iat] = new dcomplex[2*maxn+1];
    }
    
    for(int ish = 0; ish < this->size_sh(); ish++) {
      int jsh0 = ish;
      for(int jsh = jsh0; jsh < this->size_sh(); jsh++) {
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

	int mi = maxn_ish[ish]; int mj = maxn_ish[jsh];
	A3dc dxmap = calc_d_coef(mi,mj+2,mi+mj,zetaP,wPx,xi,xj,dsx_buff);
	A3dc dymap = calc_d_coef(mi,mj+2,mi+mj,zetaP,wPy,yi,yj,dsy_buff);
	A3dc dzmap = calc_d_coef(mi,mj+2,mi+mj,zetaP,wPz,zi,zj,dsz_buff);
	
	for(int iat = 0; iat < this->size_atom(); iat++) {
	  dcomplex cx(x_iat[iat]); dcomplex cy(y_iat[iat]); dcomplex cz(z_iat[iat]);
	  IncompleteGamma(maxn_ish[ish] + maxn_ish[jsh],
			  zetaP * dist2(wPx-cx,wPy-cy,wPz-cz),
			  Fjs_iat[iat]);
	}

	int npi = this->size_prim_ish(ish);
	int npj = this->size_prim_ish(jsh);
	for(int iprim = 0; iprim < npi; iprim++) {
	  int jprim0 = ish == jsh ? iprim : 0;
	  for(int jprim = jprim0; jprim < npj; jprim++) {
	    int nxi = nx_ish_iprim[ish][iprim]; int nxj = nx_ish_iprim[jsh][jprim];
	    int nyi = ny_ish_iprim[ish][iprim]; int nyj = ny_ish_iprim[jsh][jprim];
	    int nzi = nz_ish_iprim[ish][iprim]; int nzj = nz_ish_iprim[jsh][jprim];

	    // ---- S mat ----
	    dcomplex dx00 = dxmap.get(nxi, nxj, 0);
	    dcomplex dy00 = dymap.get(nyi, nyj, 0);
	    dcomplex dz00 = dzmap.get(nzi, nzj, 0);

	    // ---- z mat ----
	    dcomplex dz01 = dzmap.get(nzi, nzj+1, 0);
	    
	    // ---- t mat ----
	    dcomplex dx02 = dxmap.get(nxi, nxj+2, 0);
	    dcomplex dy02 = dymap.get(nyi, nyj+2, 0);
	    dcomplex dz02 = dzmap.get(nzi, nzj+2, 0);
	    dcomplex t_ele(0.0);
	    t_ele += -2.0*zetaj * (2*nxj+2*nyj+2*nzj+3.0) * dx00*dy00*dz00;
	    t_ele += 4.0*zetaj*zetaj*(dx02*dy00*dz00+dx00*dy02*dz00+dx00*dy00*dz02);
	    if(nxj > 1) {
	      dcomplex dx = dxmap.get(nxi, nxj-2, 0);
	      t_ele += 1.0*nxj*(nxj-1) * dx * dy00 * dz00;
	    }
	    if(nyj > 1) {
	      dcomplex dy = dymap.get(nyi, nyj-2, 0);
	      t_ele += 1.0*nyj*(nyj-1) * dx00 * dy * dz00;
	    }
	    if(nzj > 1) {
	      dcomplex dz = dzmap.get(nzi, nzj-2, 0);
	      t_ele += 1.0*nzj*(nzj-1) * dx00 * dy00 * dz;
	    }
	    
	    dcomplex v_ele(0);
	    for(int nx = 0; nx <= nxi + nxj; nx++)
	      for(int ny = 0; ny <= nyi + nyj; ny++)
		for(int nz = 0; nz <= nzi + nzj; nz++)
		  for(int iat = 0; iat < this->size_atom(); iat++) {
		    v_ele += (dxmap.get(nxi, nxj, nx) *
			      dymap.get(nyi, nyj, ny) *
			      dzmap.get(nzi, nzj, nz) *
			      coef_R(zetaP, wPx, wPy, wPz,
				     x_iat[iat], y_iat[iat], z_iat[iat],
				     nx, ny, nz, 0, Fjs_iat[iat]));
		  
		}
	    	    
	    int idx = iprim * npj + jprim;
	    smat_prim[idx] = ce * dx00 * dy00 * dz00;
	    zmat_prim[idx] = ce*(dx00*dy00*dz01 + zj*dx00*dy00*dz00);
	    tmat_prim[idx] = -0.5 * ce * t_ele;
	    vmat_prim[idx] = -2.0*M_PI/zetaP * eAB * v_ele;
	    if(ish == jsh && iprim != jprim) {
	      int idx1 = jprim * npj + iprim;
	      smat_prim[idx1] = smat_prim[idx];
	      zmat_prim[idx1] =zmat_prim[idx];
	      tmat_prim[idx1] = tmat_prim[idx];
	      vmat_prim[idx1] = vmat_prim[idx];
	    }
	  }
	}

	for(int ibasis = 0; ibasis < this->size_basis_ish(ish); ibasis++) {
	  int jbasis0 = ish == jsh ? ibasis : 0;
	  for(int jbasis = jbasis0; jbasis < this->size_basis_ish(jsh); jbasis++) {
	    dcomplex s(0.0); dcomplex dz(0.0); dcomplex t(0.0);	dcomplex v(0.0);
	    
	    for(int iprim = 0; iprim < npi; iprim++) {
	      
	      for(int jprim = 0; jprim < npj; jprim++) {
		int idx_prim0(iprim * npj + jprim);
		dcomplex cc(coef_ish_icont_iprim[ish][ibasis][iprim]*
			    coef_ish_icont_iprim[jsh][jbasis][jprim]);
		s += cc*smat_prim[idx_prim0];
		t += cc*tmat_prim[idx_prim0];
		dz += cc*zmat_prim[idx_prim0];
		v += cc*vmat_prim[idx_prim0];
		idx_prim0++;
	      }
	    }
	    int i = offset_ish[ish] + ibasis;
	    int j = offset_ish[jsh] + jbasis;
	    int idx = i+j*nb;	
	    S[idx] = s; Dz[idx] = dz; T[idx] = t; V[idx] = v;
	    if(ish != jsh || ibasis != jbasis) {
	      int idx1 = j + i*nb;
	      S[idx1] = s; Dz[idx1] = dz; T[idx1] = t; V[idx1] = v;
	    }
	  }
	}

      }
    }
    *s = S; *dz = Dz; *t = T; *v = V;    
    delete smat_prim;
    delete tmat_prim;
    delete zmat_prim;
    delete vmat_prim;
    delete dsx_buff;
    delete dsy_buff;
    delete dsz_buff;
    for(int iat = 0; iat < this->size_atom(); iat++)
      delete Fjs_iat[iat];
    delete[] Fjs_iat;
  }
  void GTOs::Show() const {
    std::cout << "Shell" << std::endl;
    for(int ish = 0; ish < this->size_sh(); ish++)
      std::cout << ish << ": " << zeta_ish[ish]
		<< x_ish[ish] << y_ish[ish] << z_ish[ish] << std::endl;
    std::cout << "Offset" << std::endl;
    for(int ish = 0; ish < this->size_sh(); ish++)
      std::cout << ish << offset_ish[ish] << std::endl;
    std::cout << "n" << std::endl;
    for(int ish = 0; ish < this->size_sh(); ish++)
      for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++) {
	std::cout << ish << iprim << " ";
	std::cout << nx_ish_iprim[ish][iprim];
	std::cout << ny_ish_iprim[ish][iprim];
	std::cout << nz_ish_iprim[ish][iprim];
	std::cout << std::endl;
      }

    std::cout << "Coef" << std::endl;
    for(int ish = 0; ish < this->size_sh(); ish++)
      for(int ibasis = 0; ibasis < this->size_basis_ish(ish); ibasis++)
	for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++)
	  std::cout << ish << ibasis << iprim << coef_ish_icont_iprim[ish][ibasis][iprim] << std::endl;
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
