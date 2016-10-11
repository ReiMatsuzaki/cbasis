/**
   One electron integrals
*/

#include <sstream>
#include "one_int.hpp"
#include "mol_func.hpp"
#include "symmolint.hpp"

using namespace std;
using namespace Eigen;

namespace l2func {

  // ==== Slow routines ====
  dcomplex SMatEle(CartGTO& a, CartGTO& b) {

    dcomplex zetaP, wPx, wPy, wPz;
    zetaP = a.zeta + b.zeta;
    wPx = (a.zeta * a.x + b.zeta * b.x) / zetaP;
    wPy = (a.zeta * a.y + b.zeta * b.y) / zetaP;
    wPz = (a.zeta * a.z + b.zeta * b.z) / zetaP;

    MultArray<dcomplex, 3> dxmap(1000), dymap(1000), dzmap(1000);
    calc_d_coef(a.nx, b.nx, 0, zetaP, wPx, a.x, b.x, dxmap);
    calc_d_coef(a.ny, b.ny, 0, zetaP, wPy, a.y, b.y, dymap);
    calc_d_coef(a.nz, b.nz, 0, zetaP, wPz, a.z, b.z, dzmap);

    dcomplex d2, eAB, ce;
    d2 = pow(a.x-b.x, 2) + pow(a.y-b.y, 2) + pow(a.z-b.z, 2);
    eAB = exp(-a.zeta*b.zeta/zetaP*d2);
    ce = eAB * pow(M_PI/zetaP, 1.5);

    return ce * dxmap(a.nx, b.nx, 0) * dymap(a.ny, b.ny, 0) * dzmap(a.nz, b.nz, 0);
  }
  dcomplex TMatEle(CartGTO& a, CartGTO& b) {

    dcomplex zetaP, wPx, wPy, wPz;
    zetaP = a.zeta + b.zeta;
    wPx = (a.zeta * a.x + b.zeta * b.x) / zetaP;
    wPy = (a.zeta * a.y + b.zeta * b.y) / zetaP;
    wPz = (a.zeta * a.z + b.zeta * b.z) / zetaP;

    MultArray<dcomplex, 3> dxmap(1000), dymap(1000), dzmap(1000);
    calc_d_coef(a.nx, b.nx+2, 0, zetaP, wPx, a.x, b.x, dxmap);
    calc_d_coef(a.ny, b.ny+2, 0, zetaP, wPy, a.y, b.y, dymap);
    calc_d_coef(a.nz, b.nz+2, 0, zetaP, wPz, a.z, b.z, dzmap);
    dcomplex dx00, dy00, dz00, dx02, dy02, dz02;
    dx00 = dxmap(a.nx, b.nx ,0);
    dy00 = dymap(a.ny, b.ny ,0);
    dz00 = dzmap(a.nz, b.nz ,0);
    dx02 = dxmap(a.nx, b.nx+2 ,0);
    dy02 = dymap(a.ny, b.ny+2 ,0);
    dz02 = dzmap(a.nz, b.nz+2 ,0);

    dcomplex d2, eAB, ce;
    d2 = pow(a.x-b.x, 2) + pow(a.y-b.y, 2) + pow(a.z-b.z, 2);
    eAB = exp(-a.zeta*b.zeta/zetaP*d2);
    ce = eAB * pow(M_PI/zetaP, 1.5);

    dcomplex tele(0);
    tele += -2.0*(2*b.nx+2*b.ny+2*b.nz+3) * b.zeta * dx00 * dy00 * dz00;
    tele += +4.0*b.zeta*b.zeta*(dx02*dy00*dz00+dx00*dy02*dz00+dx00*dy00*dz02);
    if(b.nx > 1) {
      dcomplex dx = dxmap(a.nx, b.nx-2, 0);
      tele += 1.0*b.nx*(b.nx-1) * dx * dy00 * dz00;
    }
    if(b.ny > 1) {
      dcomplex dy = dymap(a.ny, b.ny-2, 0);
      tele += 1.0*b.ny*(b.ny-1) * dx00 * dy * dz00;
    }
    if(b.nz > 1) {
      dcomplex dz = dzmap(a.nz, b.nz-2, 0);
      tele += 1.0*b.nz*(b.nz-1) * dx00 * dy00 * dz;
    }

    return -0.5*ce*tele;
  }
  dcomplex VMatEle(CartGTO& a, Vector3cd at, CartGTO& b) {

    dcomplex zetaP, wPx, wPy, wPz, wKx, wKy, wKz, d2p, d2;
    zetaP = a.zeta + b.zeta;
    wPx = (a.zeta * a.x + b.zeta * b.x) / zetaP;
    wPy = (a.zeta * a.y + b.zeta * b.y) / zetaP;
    wPz = (a.zeta * a.z + b.zeta * b.z) / zetaP;
    wKx = at[0]; wKy = at[1]; wKz = at[2];
    d2p = dist2(wPx-wKx, wPy-wKy, wPz-wKz);
    d2  = dist2(a.x-b.x, a.y-b.y, a.z-b.z);

    MultArray<dcomplex, 3> dxmap(1000), dymap(1000), dzmap(1000);
    calc_d_coef(a.nx, b.nx, a.nx+b.nx, zetaP, wPx, a.x, b.x, dxmap);
    calc_d_coef(a.ny, b.ny, a.ny+b.ny, zetaP, wPy, a.y, b.y, dymap);
    calc_d_coef(a.nz, b.nz, a.nz+b.nz, zetaP, wPz, a.z, b.z, dzmap);

    dcomplex arg = zetaP * d2p;
    double delta(0.0000000000001);
    dcomplex Fjs[100], Gjs[100];
    if(real(arg)+delta > 0.0) {
      IncompleteGamma(a.nx+a.ny+a.nz+b.nx+b.ny+b.nz, arg, Fjs);
    } else {
      ExpIncompleteGamma(a.nx+a.ny+a.nz+b.nx+b.ny+b.nz, -arg, Gjs);
    }

    MultArray<dcomplex, 3> rmap(1000);
    rmap.SetRange(0, a.nx+b.nx, 0, a.ny+b.ny, 0, a.nz+b.nz);
    for(int nx = 0; nx <= a.nx + b.nx; nx++) {
      for(int ny = 0; ny <= a.ny + b.ny; ny++) {
	for(int nz = 0; nz <= a.nz + b.nz; nz++) {
	  if(real(arg)+delta > 0.0) {
	    rmap(nx, ny, nz) = -2.0*M_PI/zetaP*exp(-a.zeta*b.zeta/zetaP*d2) *
	      coef_R(zetaP, wPx, wPy, wPz, wKx, wKy, wKz, nx, ny, nz, 0, Fjs);
	  } else {
	    rmap(nx, ny, nz) = -2.0*M_PI/zetaP*exp(-a.zeta*b.zeta/zetaP*d2-arg) *
	      coef_R(zetaP, wPx, wPy, wPz, wKx, wKy, wKz, nx, ny, nz, 0, Gjs);
	  }
	}
      }
    }

    dcomplex vele(0);
    for(int nx = 0; nx <= a.nx + b.nx; nx++) {
      for(int ny = 0; ny <= a.ny + b.ny; ny++) {
	for(int nz = 0; nz <= a.nz + b.nz; nz++) {
	  vele += (dxmap(a.nx, b.nx, nx) *
		   dymap(a.ny, b.ny, ny) *
		   dzmap(a.nz, b.nz, nz) *
		   rmap(nx, ny, nz));
	}
      }
    }
    return vele;

  }
  dcomplex DXMatEle(CartGTO& a, CartGTO& b) {
    
    /**
       Dx [ (x-x0)^n exp(-z(x-x0)^2) ] 
       =  { n(x-x0)^(n-1) -2z(x-x0)^(n+1)} exp(-z(x-x0)^2)
     */

    dcomplex ele(0);
    
    if(b.nx != 0) {
      
      CartGTO b1 = b;
      b1.nx--;
      ele += dcomplex(b.nx) * SMatEle(a, b1);
      
    }
    
    CartGTO b2 = b;
    b2.nx++;
    ele += -2.0 * b.zeta * SMatEle(a, b2);
    
    return ele;
    
  }
  dcomplex DYMatEle(CartGTO& a, CartGTO& b) {
    
    dcomplex ele(0);
    
    if(b.ny != 0) {
      
      CartGTO b1 = b;
      b1.ny--;
      ele += dcomplex(b.ny) * SMatEle(a, b1);
      
    }
    
    CartGTO b2 = b;
    b2.ny++;
    ele += -2.0 * b.zeta * SMatEle(a, b2);
    
    return ele;

  }
  dcomplex DZMatEle(CartGTO& a, CartGTO& b) {

    dcomplex ele(0);
    
    if(b.nz != 0) {
      
      CartGTO b1 = b;
      b1.nz--;
      ele += dcomplex(b.nz) * SMatEle(a, b1);
      
    }
    
    CartGTO b2 = b;
    b2.nz++;
    ele += -2.0 * b.zeta * SMatEle(a, b2);
    
    return ele;

  }

  // ==== SymGTOs ====
  struct PrimBasis {
    A4dc s, t, v, x, y, z, dx, dy, dz;
    PrimBasis(int num):
      s(num), t(num), v(num), x(num), y(num), z(num),
      dx(num), dy(num), dz(num) {}
    void SetRange(int niat, int nipn, int njat, int njpn) {
      s.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
      t.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
      v.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
      x.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
      y.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
      z.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
      dx.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
      dy.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
      dz.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
    }
  };
  dcomplex calc_tele(SubIt isub, SubIt jsub, dcomplex zetaj, 
		     int ipn, int jpn,
		     A3dc& dxmap, A3dc& dymap, A3dc& dzmap) {
    int nxi, nxj, nyi, nyj, nzi, nzj;
    nxi = isub->nx(ipn); nxj = jsub->nx(jpn);
    nyi = isub->ny(ipn); nyj = jsub->ny(jpn);
    nzi = isub->nz(ipn); nzj = jsub->nz(jpn);
    dcomplex dx00, dy00, dz00, dx02, dy02, dz02;
    dx00 = dxmap(nxi, nxj ,0);
    dy00 = dymap(nyi, nyj ,0);
    dz00 = dzmap(nzi, nzj ,0);
    dx02 = dxmap(nxi, nxj+2 ,0);
    dy02 = dymap(nyi, nyj+2 ,0);
    dz02 = dzmap(nzi, nzj+2 ,0);

    dcomplex t_ele(0.0);		    
    t_ele += -2.0*(2*nxj+2*nyj+2*nzj+3)*zetaj*dx00*dy00*dz00;
    t_ele += 4.0*zetaj*zetaj*(dx02*dy00*dz00+dx00*dy02*dz00+dx00*dy00*dz02);
    if(nxj > 1) {
      dcomplex dx = dxmap(nxi, nxj-2, 0);
      t_ele += 1.0*nxj*(nxj-1) * dx * dy00 * dz00;
    }
    if(nyj > 1) {
      dcomplex dy = dymap(nyi, nyj-2, 0);
      t_ele += 1.0*nyj*(nyj-1) * dx00 * dy * dz00;
    }
    if(nzj > 1) {
      dcomplex dz = dzmap(nzi, nzj-2, 0);
      t_ele += 1.0*nzj*(nzj-1) * dx00 * dy00 * dz;
    }
    return t_ele;

  }
  dcomplex calc_vele(const SymGTOs gtos, SubIt isub, SubIt jsub, int ipn, int jpn,
		     A3dc& dxmap, A3dc& dymap, A3dc& dzmap, A4dc& rmap) {

    int nxi, nxj, nyi, nyj, nzi, nzj;
    nxi = isub->nx(ipn); nxj = jsub->nx(jpn);
    nyi = isub->ny(ipn); nyj = jsub->ny(jpn);
    nzi = isub->nz(ipn); nzj = jsub->nz(jpn);

    dcomplex v_ele(0.0);
    for(int nx = 0; nx <= nxi + nxj; nx++)
      for(int ny = 0; ny <= nyi + nyj; ny++)
	for(int nz = 0; nz <= nzi + nzj; nz++)
	  for(int kat = 0; kat < gtos->size_atom(); kat++) {

	    v_ele += (gtos->q_at(kat) *
		      dxmap(nxi, nxj, nx) *
		      dymap(nyi, nyj, ny) *
		      dzmap(nzi, nzj, nz) *
		      rmap(nx, ny, nz, kat));
		  }
    return v_ele;
  }
  void CalcPrim(const SymGTOs gtos, SubIt isub, SubIt jsub, int iz, int jz,
		PrimBasis& prim, bool calc_coulomb) {
    dcomplex zetai, zetaj;
    zetai = isub->zeta_iz[iz]; zetaj = jsub->zeta_iz[jz];
    dcomplex zetaP = zetai + zetaj;
    int niat(isub->size_at()); int njat(jsub->size_at());
    int nipn(isub->size_pn()); int njpn(jsub->size_pn());

    static const int num(1000);
    A3dc dxmap(num), dymap(num), dzmap(num);
    A4dc rmap(num);
    A2dc Fjs_iat(num);

    prim.SetRange(niat, nipn, njat, njpn);
    Fjs_iat.SetRange(0, gtos->size_atom(), 0, isub->maxn+jsub->maxn);

    for(int iat = 0; iat < niat; iat++) {
      // int jat0 = (isub == jsub && iz == jz ? iat : 0);
      for(int jat = 0; jat < njat; jat++) { 
	dcomplex xi, xj, yi, yj, zi, zj, wPx, wPy, wPz;
	xi = isub->x(iat); yi = isub->y(iat); zi = isub->z(iat);
	xj = jsub->x(jat); yj = jsub->y(jat); zj = jsub->z(jat);
	
	wPx = (zetai*xi+zetaj*xj)/zetaP;
	wPy = (zetai*yi+zetaj*yj)/zetaP;
	wPz = (zetai*zi+zetaj*zj)/zetaP;
	dcomplex d2 = dist2(xi-xj, yi-yj, zi-zj);
	dcomplex eAB = exp(-zetai*zetaj/zetaP*d2);
	dcomplex ce = eAB * pow(M_PI/zetaP, 1.5);
	int mi = isub->maxn; int mj = jsub->maxn;
	if(calc_coulomb) {
	  calc_d_coef(mi,mj+2,mi+mj,zetaP,wPx,xi,xj,dxmap);
	  calc_d_coef(mi,mj+2,mi+mj,zetaP,wPy,yi,yj,dymap);
	  calc_d_coef(mi,mj+2,mi+mj,zetaP,wPz,zi,zj,dzmap);
	  for(int kat = 0; kat < gtos->size_atom(); kat++) {
	    dcomplex d2p = dist2(wPx-gtos->x_at(kat), wPy-gtos->y_at(kat),
				 wPz-gtos->z_at(kat));
	    dcomplex arg = zetaP * d2p;
	    double delta(0.0000000000001);
	    if(real(arg)+delta > 0.0) {
	      IncompleteGamma(isub->maxn+jsub->maxn, arg, &Fjs_iat(kat, 0));
	    } else {
	      ExpIncompleteGamma(isub->maxn+jsub->maxn, -arg, &Fjs_iat(kat, 0));
	    }
	  }
	} else {
	  calc_d_coef(mi,mj+2,0,zetaP,wPx,xi,xj,dxmap);
	  calc_d_coef(mi,mj+2,0,zetaP,wPy,yi,yj,dymap);
	  calc_d_coef(mi,mj+2,0,zetaP,wPz,zi,zj,dzmap);	  
	}
	
	for(int ipn = 0; ipn < nipn; ipn++) {
	  for(int jpn = 0; jpn < njpn; jpn++) {
	    int nxi, nxj, nyi, nyj, nzi, nzj;
	    nxi = isub->nx(ipn); nxj = jsub->nx(jpn);
	    nyi = isub->ny(ipn); nyj = jsub->ny(jpn);
	    nzi = isub->nz(ipn); nzj = jsub->nz(jpn);

	    dcomplex s_ele = dxmap(nxi,nxj,0) * dymap(nyi,nyj,0) * dzmap(nzi,nzj,0);
	    dcomplex x_ele = dymap(nyi,nyj,0)*dzmap(nzi,nzj,0)*
	      (dxmap(nxi,nxj+1,0)+xj*dxmap(nxi,nxj,0));
	    dcomplex y_ele = dzmap(nzi,nzj,0)*dxmap(nxi,nxj,0)*
	      (dymap(nyi,nyj+1,0)+yj*dymap(nyi,nyj,0));
	    dcomplex z_ele = dxmap(nxi,nxj,0)*dymap(nyi,nyj,0)*
	      (dzmap(nzi,nzj+1,0)+zj*dzmap(nzi,nzj,0));
	    dcomplex t_ele = calc_tele(isub, jsub, zetaj, ipn, jpn,
				       dxmap, dymap, dzmap);	    
	    dcomplex dx_ele = -2.0*zetaj*dxmap(nxi,nxj+1,0);
	    if(nxj>0)
	      dx_ele += dcomplex(nxj)*dxmap(nxi,nxj-1,0);
	    dx_ele *= dymap(nyi,nyj,0) * dzmap(nzi,nzj,0);
	    dcomplex dy_ele = -2.0*zetaj*dymap(nyi,nyj+1,0);
	    if(nyj>0)
	      dy_ele += dcomplex(nyj)*dymap(nyi,nyj-1,0);
	    dy_ele *= dxmap(nxi,nxj,0) * dzmap(nzi,nzj,0);
	    dcomplex dz_ele = -2.0*zetaj*dzmap(nzi,nzj+1,0);
	    dz_ele += dcomplex(nzj)*dzmap(nzi,nzj-1,0);
	    dz_ele *= dymap(nyi,nyj,0) * dxmap(nxi,nxj,0);
	      
	    prim.s(iat, ipn, jat, jpn) =  ce * s_ele;
	    prim.t(iat, ipn, jat, jpn) =  -0.5* ce * t_ele;
	    prim.x(iat, ipn, jat, jpn) =  ce*x_ele;
	    prim.y(iat, ipn, jat, jpn) =  ce*y_ele;
	    prim.z(iat, ipn, jat, jpn) =  ce*z_ele;
	    prim.dx(iat, ipn, jat, jpn) = ce*dx_ele;
	    prim.dy(iat, ipn, jat, jpn) = ce*dy_ele;
	    prim.dz(iat, ipn, jat, jpn) = ce*dz_ele;
	    rmap.SetRange(0, nxi+nxj, 0, nyi+nyj,
			  0, nzi+nzj, 0, gtos->size_atom());

	    if(calc_coulomb) {	      
	      for(int kat = 0; kat < gtos->size_atom(); kat++) {
		dcomplex d2p = dist2(wPx-gtos->x_at(kat), wPy-gtos->y_at(kat),
				     wPz-gtos->z_at(kat));
		dcomplex arg = zetaP * d2p;
		double delta(0.0000000000001);
		dcomplex wKx = gtos->x_at(kat);
		dcomplex wKy = gtos->y_at(kat);
		dcomplex wKz = gtos->z_at(kat);
		if(real(arg)+delta > 0.0) {
		  for(int nx = 0; nx <= nxi + nxj; nx++)
		    for(int ny = 0; ny <= nyi + nyj; ny++)
		      for(int nz = 0; nz <= nzi + nzj; nz++)
			rmap(nx, ny, nz, kat) =
			  -2.0*M_PI/zetaP*eAB * coef_R(zetaP, wPx, wPy, wPz,
						       wKx, wKy, wKz, nx, ny, nz,
						       0, &Fjs_iat(kat, 0));
		} else {
		  dcomplex arg_other = -zetai * zetaj / zetaP * d2 - arg;
		  for(int nx = 0; nx <= nxi + nxj; nx++)
		    for(int ny = 0; ny <= nyi + nyj; ny++)
		      for(int nz = 0; nz <= nzi + nzj; nz++) {
			rmap(nx, ny, nz, kat) =
			  -2.0*M_PI/zetaP*exp(arg_other) *
			  coef_R(zetaP, wPx, wPy, wPz,
				 wKx, wKy, wKz, nx, ny, nz,
				 0, &Fjs_iat(kat, 0)); 
		      }
		}
		
	      }
	      prim.v(iat, ipn, jat, jpn) =  
		    calc_vele(gtos, isub, jsub, ipn, jpn, dxmap, dymap, dzmap, rmap);
	    }
	  }}}}

  }
  void CalcTrans(SubIt isub, SubIt jsub, int iz, int jz,
		 PrimBasis& prim, BMatSet mat_map) {

    int niat(isub->size_at()); int njat(jsub->size_at());
    int nipn(isub->size_pn()); int njpn(jsub->size_pn());
    //    cout << "tras before set : " << niat << nipn << njat << njpn << endl;

    for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) {
      for(RdsIt jrds = jsub->rds.begin(); jrds != jsub->rds.end();++jrds) {
	dcomplex cumsum_s(0.0), cumsum_t(0.0), cumsum_v(0.0);
	dcomplex cumsum_x(0.0), cumsum_y(0.0), cumsum_z(0.0);
	dcomplex cumsum_dx(0.0), cumsum_dy(0.0), cumsum_dz(0.0);
	for(int iat = 0; iat < niat; iat++) {
	  for(int ipn = 0; ipn < nipn; ipn++) {
	    for(int jat = 0; jat < njat; jat++) { 
	      for(int jpn = 0; jpn < njpn; jpn++) {
		dcomplex cc = 
		  irds->coef_iat_ipn(iat, ipn) *
		  jrds->coef_iat_ipn(jat, jpn) * 
		  irds->coef_iz(iz) * 
		  jrds->coef_iz(jz);
		cumsum_s += cc*prim.s(iat, ipn, jat, jpn);
		cumsum_t += cc*prim.t(iat, ipn, jat, jpn);
		cumsum_v += cc*prim.v(iat, ipn, jat, jpn);
		cumsum_x += cc*prim.x(iat, ipn, jat, jpn);
		cumsum_y += cc*prim.y(iat, ipn, jat, jpn);
		cumsum_z += cc*prim.z(iat, ipn, jat, jpn);
		cumsum_dx += cc*prim.dx(iat, ipn, jat, jpn);
		cumsum_dy += cc*prim.dy(iat, ipn, jat, jpn);
		cumsum_dz += cc*prim.dz(iat, ipn, jat, jpn);
	      }}}}
	int i(irds->offset + iz); int j(jrds->offset + jz);
	int isym(irds->irrep); int jsym(jrds->irrep);
	
	mat_map->SelfAdd("s", isym, jsym, i, j, cumsum_s);
	mat_map->SelfAdd("t", isym, jsym, i, j, cumsum_t);
	mat_map->SelfAdd("v", isym, jsym, i, j, cumsum_v);
	mat_map->SelfAdd("x", isym, jsym, i, j, cumsum_x);
	mat_map->SelfAdd("y", isym, jsym, i, j, cumsum_y);
	mat_map->SelfAdd("z", isym, jsym, i, j, cumsum_z);
	mat_map->SelfAdd("dx", isym, jsym, i, j, cumsum_dx);
	mat_map->SelfAdd("dy", isym, jsym, i, j, cumsum_dy);
	mat_map->SelfAdd("dz", isym, jsym, i, j, cumsum_dz);	
      }
    }

  }
  BMatSet CalcMat(SymGTOs a, SymGTOs b, bool calc_coulomb) {

    BMatSet bmat(new _BMatSet);
   
    if(not a->setupq)
      a->SetUp();

    if(not b->setupq)
      b->SetUp();

    int num_sym(a->sym_group->num_class());

    for(Irrep isym = 0; isym < num_sym; isym++) {
      for(Irrep jsym = 0; jsym < num_sym; jsym++) {
	int numi = a->size_basis_isym(isym);
	int numj = b->size_basis_isym(jsym);
	MatrixXcd s = MatrixXcd::Zero(numi, numj);
	bmat->SetMatrix("s", isym, jsym, s);
	MatrixXcd t = MatrixXcd::Zero(numi, numj); 
	bmat->SetMatrix("t", isym, jsym, t);
	MatrixXcd v = MatrixXcd::Zero(numi, numj); 
	bmat->SetMatrix("v", isym, jsym, v);

	MatrixXcd x = MatrixXcd::Zero(numi, numj); 
	bmat->SetMatrix("x", isym, jsym, x);
	MatrixXcd y = MatrixXcd::Zero(numi, numj); 
	bmat->SetMatrix("y", isym, jsym, y);
	MatrixXcd z = MatrixXcd::Zero(numi, numj); 
	bmat->SetMatrix("z", isym, jsym, z);
	
	MatrixXcd dx = MatrixXcd::Zero(numi, numj);
	bmat->SetMatrix("dx", isym, jsym, dx);
	MatrixXcd dy = MatrixXcd::Zero(numi, numj);
	bmat->SetMatrix("dy", isym, jsym, dy);
	MatrixXcd dz = MatrixXcd::Zero(numi, numj);
	bmat->SetMatrix("dz", isym, jsym, dz);
      }
    }

    PrimBasis prim(100);
    
    for(SubIt isub = a->subs.begin(); isub != a->subs.end(); ++isub) {
      for(SubIt jsub = b->subs.begin(); jsub != b->subs.end(); ++jsub) {
	for(int iz = 0; iz < isub->size_zeta(); iz++) {
	  for(int jz = 0; jz < jsub->size_zeta(); jz++) {
	    CalcPrim(a, isub, jsub, iz, jz, prim, calc_coulomb);
	    CalcTrans(isub, jsub, iz, jz, prim, bmat);
	  }
	}
      }
    }

    return bmat;

  }
  BMatSet CalcMat_Complex(SymGTOs g, bool calc_coulomb) {

    return CalcMat(g, g, calc_coulomb);
    
  }
  BMatSet CalcMat_Hermite(SymGTOs g, bool calc_coulomb) {
    SymGTOs c = g->Conj();
    return CalcMat(c, g, calc_coulomb);
  }

  /*
  BMatSet CalcMat_V(SymGTOs a, SymGTOs b, Eigen::Vector3cd xyz, dcomplex q) {
    MatrixXcd xyzq(4,1);
    xyzq << xyz[0], xyz[1], xyz[2], q;

    SymGTOs aa = a->Clone();
    aa->SetAtoms(xyzq);
    aa->SetUp();

  }
  */
}