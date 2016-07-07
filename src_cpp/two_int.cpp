#include <iostream>
#include "typedef.hpp"
#include "two_int.hpp"

using namespace std;

namespace l2func {

  dcomplex ERIEle(CartGTO& i, CartGTO& j, CartGTO& k, CartGTO& l) {
    dcomplex zetaP = i.zeta + j.zeta;
    dcomplex zetaPp= k.zeta + l.zeta;
    dcomplex wPx = (i.zeta*i.x + j.zeta*j.x)/zetaP;
    dcomplex wPy = (i.zeta*i.y + j.zeta*j.y)/zetaP;
    dcomplex wPz = (i.zeta*i.z + j.zeta*j.z)/zetaP;
    dcomplex wPpx= (k.zeta*k.x + l.zeta*l.x)/zetaPp;
    dcomplex wPpy= (k.zeta*k.y + l.zeta*l.y)/zetaPp;
    dcomplex wPpz= (k.zeta*k.z + l.zeta*l.z)/zetaPp;
    dcomplex d2 = dist2(i.x-j.x, i.y-j.y, i.z-j.z);
    dcomplex d2p= dist2(k.x-l.x, k.y-l.y, k.z-l.z);

    static const int num(1000);
    MultArray<dcomplex, 3> dx(num),  dy(num),  dz(num);
    MultArray<dcomplex, 3> dxp(num), dyp(num), dzp(num);
    calc_d_coef(i.nx, j.nx, i.nx+j.nx, zetaP,  wPx,  i.x, j.x, dx);
    calc_d_coef(i.ny, j.ny, i.ny+j.ny, zetaP,  wPy,  i.y, j.y, dy);
    calc_d_coef(i.nz, j.nz, i.nz+j.nz, zetaP,  wPz,  i.z, j.z, dz);
    calc_d_coef(k.nx, l.nx, k.nx+l.nx, zetaPp, wPpx, k.x, l.x, dxp);
    calc_d_coef(k.ny, l.ny, k.ny+l.ny, zetaPp, wPpy, k.y, l.y, dyp);
    calc_d_coef(k.nz, l.nz, k.nz+l.nz, zetaPp, wPpz, k.z, l.z, dzp);
    
    double delta(0.0000000000001);
    dcomplex arg = zetaP * zetaPp / (zetaP + zetaPp) * dist2(wPx-wPpx, wPy-wPpy, wPz-wPpz);
    int mm = i.nx+j.nx+k.nx+l.nx+i.ny+j.ny+k.ny+l.ny+i.nz+j.nz+k.nz+l.nz;
    dcomplex Fj_or_Gj[num];
    dcomplex c(0);
    if(real(arg)+delta > 0.0) {
      IncompleteGamma(mm, arg, Fj_or_Gj);
      c = exp(-i.zeta*j.zeta/zetaP*d2 -k.zeta*l.zeta/zetaPp*d2p);
    } else {
      ExpIncompleteGamma(mm, -arg, Fj_or_Gj);
      c = exp(-i.zeta*j.zeta/zetaP*d2 -k.zeta*l.zeta/zetaP*d2p -arg);
    }

    MultArray<dcomplex, 3> Rrs(num);
    Rrs.SetRange(0, mm, 0, mm, 0, mm);
    for(int nx = 0; nx <= mm; nx++)
      for(int ny = 0; ny <= mm; ny++)
	for(int nz = 0; nz <= mm; nz++)
	  if(nx + ny + nz <= mm) {
	    Rrs(nx, ny, nz) = c * coef_R(zetaP*zetaPp/(zetaP+zetaPp),
					 wPx, wPy, wPz, wPpx, wPpy, wPpz,
					 nx, ny, nz, 0, Fj_or_Gj);
	  }
    dcomplex cumsum(0);
    for(int Nx  = 0; Nx  <= i.nx + j.nx; Nx++)
    for(int Nxp = 0; Nxp <= k.nx + l.nx; Nxp++)
    for(int Ny  = 0; Ny  <= i.ny + j.ny; Ny++)
    for(int Nyp = 0; Nyp <= k.ny + l.ny; Nyp++)
    for(int Nz  = 0; Nz  <= i.nz + j.nz; Nz++)
    for(int Nzp = 0; Nzp <= k.nz + l.nz; Nzp++) {
      dcomplex r0;
      r0 = Rrs(Nx+Nxp, Ny+Nyp, Nz+Nzp);
      cumsum += (dx(i.nx, j.nx, Nx) * dxp(k.nx, l.nx, Nxp) *
		 dy(i.ny, j.ny, Ny) * dyp(k.ny, l.ny, Nyp) *
		 dz(i.nz, j.nz, Nz) * dzp(k.nz, l.nz, Nzp) *
		 r0 * pow(-1.0, Nxp+Nyp+Nzp));
    }
    
    dcomplex lambda = 2.0*pow(M_PI, 2.5)/(zetaP * zetaPp * sqrt(zetaP + zetaPp));
    return lambda * cumsum;
  }
  
}
