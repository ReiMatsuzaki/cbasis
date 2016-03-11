#include "spec_func.hpp"

namespace l2func {

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
      double t1 = (2.0*m+1) / (2.0*m+3);
      double t2 = (2.0*m+1) / ((2*m+3)*(2*m+5));
      double t3 = (double((2*m+1)*((2*m+1)*(2*m+1)+44)) /
		   double(60*(2*m+3)*(2*m+5)*(2*m+7)));
	
      An[1] = Bn[1] - t1*z;
      An[2] = Bn[2] - t1*z - t2*z*z;
      An[3] = Bn[3] - t1*z - t2*z*z - t3*z*z*z;
      for(int n = 4; n <= NR; n++) {
	double F1 = double(2*n-2*m-5)/double(2*(2*n-3)*(2*n+2*m+1));
	double F2 = double(1) / double(4*(2*n-1)*(2*n-3));
	double F3 = double(-F1) / double(4*(2*n-3)*(2*n-5));
	double E = -F1;
	An[n] = (1.0+F1*z)*An[n-1] + (E + F2*z)*z*An[n-2] + F3*z*z*z*An[n-3];
      }
      res[m] = 1.0 / (2*m+1) * An[NR]/Bn[NR];
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

  MultArray<dcomplex, 3> calc_d_coef(int max_ni, int max_nj, int max_n,
				     dcomplex zetaP, dcomplex wPx,
				     dcomplex xi, dcomplex xj,
				     dcomplex* buffer) {

    MultArray<dcomplex, 3> data(buffer,
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

}
