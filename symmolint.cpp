#include <iostream>
#include <numeric>
#include "symmolint.hpp"
#include "spec_func.hpp"

namespace l2func {

  using namespace std;
  using namespace Eigen;
  typedef vector<SubSymGTOs>::iterator SubIt;
  typedef vector<ReductionSets>::iterator RdsIt;
  typedef vector<SubSymGTOs>::const_iterator cSubIt;
  typedef vector<ReductionSets>::const_iterator cRdsIt;
  typedef MultArray<dcomplex, 3> A3dc;

  SubSymGTOs::SubSymGTOs(Eigen::MatrixXcd xyz, Eigen::MatrixXi ns,
			 std::vector<ReductionSets> ao, Eigen::VectorXcd zs) :
    xyz_iat(xyz), ns_ipn(ns), rds(ao), zeta_iz(zs) {

    maxn = 0;
    for(int ipn = 0; ipn < ns_ipn.cols(); ipn++) {
      if(maxn < ns_ipn(0, ipn) + ns_ipn(1, ipn) + ns_ipn(2, ipn))
	maxn = ns_ipn(0, ipn) + ns_ipn(1, ipn) + ns_ipn(2, ipn);
    }

    for(RdsIt it = rds.begin(); it != rds.end(); ++it) {
      it->set_zs_size(zs.size());
    }

  }

  SubSymGTOs Sub_s(Irrep sym, Vector3cd xyz, VectorXcd zs) {

    MatrixXcd xyz_in(3, 1);  xyz_in << xyz[0] , xyz[1] , xyz[2];
    MatrixXi ns  = MatrixXi::Zero(3, 1);
    MatrixXcd cs = MatrixXcd::Ones(1,1);
    ReductionSets rds(sym, cs);
    vector<ReductionSets> rds_list; rds_list.push_back(rds);

    SubSymGTOs sub(xyz_in, ns, rds_list, zs);
    return sub;
  }
  SubSymGTOs Sub_pz(Irrep sym, Vector3cd xyz, VectorXcd zs) {

    MatrixXcd xyz_in(3, 1);  xyz_in << xyz[0] , xyz[1] , xyz[2];
    MatrixXi ns(3, 1); ns << 0, 0, 1;
    MatrixXcd cs = MatrixXcd::Ones(1,1);
    ReductionSets rds(sym, cs);
    vector<ReductionSets> rds_list; rds_list.push_back(rds);

    SubSymGTOs sub(xyz_in, ns, rds_list, zs);
    return sub;    
  }

  // ==== SymmetryGroup ====
  SymmetryGroup::SymmetryGroup(int order) {
    if(order != 2 && order != 1) {
      string msg; SUB_LOCATION(msg);
      msg += "Now, only order == 1 and 2 is supported.";
      throw runtime_error(msg);	
    }
    order_ = order;
  }
  void SymmetryGroup::CheckIrrep(Irrep a) const {
    if(a < 0 || this->order() <= a) {
      std::string msg; SUB_LOCATION(msg);
      throw runtime_error(msg);
    }
  }
  // Irrep SymmetryGroup::GetIrrep(int n) const {
  //    this->CheckIrrep(n);
  //    return n;
  //  }
  bool SymmetryGroup::IncludeScalar_2(Irrep a, Irrep b) const {
    
    this->CheckIrrep(a);
    this->CheckIrrep(b);

    if(a == b)
      return true;
    else
      return false;
  }
  bool SymmetryGroup::IncludeZ_2(Irrep a,Irrep b) const {

    this->CheckIrrep(a);
    this->CheckIrrep(b);

    if(this->order() == 1)
      return true;

    if(a != b)
      return true;
    else
      return false;
  }
  
  SymmetryGroup SymmetryGroup_Cs() {
    SymmetryGroup cs(2);
    return cs;
  }
  SymmetryGroup SymmetryGroup_C1() {
    SymmetryGroup c1(1);
    return c1;
  }

  // ==== SymGTOs ====
  // ---- Constructors ----
  SymGTOs::SymGTOs(SymmetryGroup _sym_group):
    sym_group(_sym_group), setupq(false) {}

  // ---- Accessors ----
  int SymGTOs::size_atom() const {return xyzq_iat.cols(); }
  int SymGTOs::size_basis_isym(Irrep isym) const {

    int cumsum(0);
    for(cSubIt isub = subs.begin(); isub != subs.end(); ++isub) {

      for(cRdsIt irds = isub->rds.begin(); irds != isub->rds.end(); 
	  ++irds) {

	if(irds->sym == isym) 
	  cumsum += isub->size_zeta();
      }
    }
    return cumsum;

  }

  // ---- Add ----
  void SymGTOs::SetAtoms(MatrixXcd _xyzq_iat) {

    if(_xyzq_iat.rows() != 4) {
      string msg; SUB_LOCATION(msg);
      msg += "xyzq_iat.rows() must be 4 ";
      throw runtime_error(msg);
    }

    xyzq_iat = _xyzq_iat;

  }
  void SymGTOs::AddSub(SubSymGTOs sub) {

    setupq = false;
    for(RdsIt irds = sub.rds.begin();
	irds != sub.rds.end(); ++irds) {
      irds->offset = 0;
    }
    subs.push_back(sub);

  }
  void SymGTOs::SetUp() {
    this->SetOffset();
    this->Normalize();
    setupq = true;
  }
  void SymGTOs::SetOffset() {
    
    vector<int> num_sym(10, 0);
    
    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(RdsIt irds = isub->rds.begin(), end = isub->rds.end();
	  irds != end; ++irds) {
	irds->offset = num_sym[irds->sym];
	num_sym[irds->sym] += isub->size_zeta();
      }
    }

  }
  void SymGTOs::Normalize() {
    
    dcomplex* dsx_buff = new dcomplex[100];
    dcomplex* dsy_buff = new dcomplex[100];
    dcomplex* dsz_buff = new dcomplex[100];

    // >>> Irrep Adapted GTOs >>>
    for(SubIt isub = subs.begin(), end = subs.end(); isub != end; ++isub) {
      for(int iz = 0; iz < isub->size_zeta(); iz++) {
	for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end();
	    ++irds) {
	  dcomplex zetai = isub->zeta_iz[iz];
	  dcomplex zetaP = zetai + zetai;
	  int niat(isub->size_at()); int nipn(isub->size_pn());
	  dcomplex norm2(0.0);

	  // >>> Primitive GTOs >>>
	  for(int iat = 0; iat < niat; iat++) {
	    for(int jat = 0; jat < niat; jat++) {
	      dcomplex xi, xj, yi, yj, zi, zj, wPx, wPy, wPz;
		xi = isub->xyz_iat(0, iat); xj = isub->xyz_iat(0, jat);
		yi = isub->xyz_iat(1, iat); yj = isub->xyz_iat(1, jat);
		zi = isub->xyz_iat(2, iat); zj = isub->xyz_iat(2, jat);
		wPx = (zetai*xi+zetai*xj)/zetaP;
		wPy = (zetai*yi+zetai*yj)/zetaP;
		wPz = (zetai*zi+zetai*zj)/zetaP;
		dcomplex d2 = pow(xi-xj,2) + pow(yi-yj,2) + pow(zi-zj,2);	
		dcomplex eAB = exp(-zetai*zetai/zetaP*d2);
		dcomplex ce = eAB * pow(M_PI/zetaP, 1.5);
		int mi = isub->maxn;
		A3dc dxmap = calc_d_coef(mi,mi,0, zetaP,wPx,xi,xj,dsx_buff);
		A3dc dymap = calc_d_coef(mi,mi,0, zetaP,wPy,yi,yj,dsy_buff);
		A3dc dzmap = calc_d_coef(mi,mi,0, zetaP,wPz,zi,zj,dsz_buff);

		for(int ipn = 0; ipn < nipn; ipn++) {
		  for(int jpn = 0; jpn < nipn; jpn++) {
		    int nxi, nxj, nyi, nyj, nzi, nzj;
		    nxi = isub->ns_ipn(0, ipn); nxj = isub->ns_ipn(0, jpn);
		    nyi = isub->ns_ipn(1, ipn); nyj = isub->ns_ipn(1, jpn);
		    nzi = isub->ns_ipn(2, ipn); nzj = isub->ns_ipn(2, jpn);
		    dcomplex dx00, dy00, dz00;
		    dx00 = dxmap.get_safe(nxi, nxj ,0);
		    dy00 = dymap.get_safe(nyi, nyj ,0);
		    dz00 = dzmap.get_safe(nzi, nzj ,0);
		    dcomplex s_ele = dx00 * dy00 * dz00;
		    norm2 += ce * s_ele *
		      irds->coef_iat_ipn(iat, ipn) *
		      irds->coef_iat_ipn(jat, jpn);
		  }
		}
	    }
	  }
	  // <<< Primitive GTOs <<<
	  irds->coef_iz[iz] = 1.0/sqrt(norm2);
	}
      }
    }
    // <<< Irrep Adapted GTOs <<<

    delete[] dsx_buff;
    delete[] dsy_buff;
    delete[] dsz_buff;
  }
  
  // ---- Calculation ----
  void SymGTOs::loop() {}
  BMatMap SymGTOs::STVMat() {

    if(not setupq) {
      string msg; SUB_LOCATION(msg); 
      msg += "call SetUp before calculation";
      throw std::runtime_error(msg);
    }

    int max_isym(0);
    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) 
	if(max_isym < irds->sym)
	  max_isym = irds->sym;
    }
    int max_n(0);
    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(int ipn = 0; ipn < isub->ns_ipn.cols(); ipn++)
	for(int i = 0; i < 3; i++)
	  if(max_n < isub->ns_ipn(i, ipn))
	    max_n = isub->ns_ipn(i, ipn);
    }

    BMatMap mat_map;
    mat_map["s"] = BMat(); mat_map["t"] = BMat(); mat_map["v"] = BMat();
    for(Irrep sym = 0; sym <= max_isym; sym++) {
      int num = this->size_basis_isym(sym);
      cout << "(sym, num) =  " << sym << ", " << num << endl;
      mat_map["s"][make_pair(sym, sym)] = MatrixXcd::Zero(num, num);
      mat_map["t"][make_pair(sym, sym)] = MatrixXcd::Zero(num, num);
      mat_map["v"][make_pair(sym, sym)] = MatrixXcd::Zero(num, num);
    }
    
    dcomplex* bufs = new dcomplex[100];
    dcomplex* buft = new dcomplex[100];
    dcomplex* bufv = new dcomplex[100];
    dcomplex* dsx_buff = new dcomplex[100];
    dcomplex* dsy_buff = new dcomplex[100];
    dcomplex* dsz_buff = new dcomplex[100];
    dcomplex** Fjs_iat;    

    Fjs_iat = new dcomplex*[this->size_atom()];
    for(int iat = 0; iat < this->size_atom(); iat++)
      Fjs_iat[iat] = new dcomplex[2*max_n+1];
    
    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(SubIt jsub = subs.begin(); jsub != subs.end(); ++jsub) {
	
	// -- search symmetry for isub and jsub --
	bool is_zero(true);
	for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) 
	  for(RdsIt jrds = jsub->rds.begin(); jrds != jsub->rds.end(); ++jrds) 
	    if(sym_group.IncludeScalar_2(irds->sym, jrds->sym))
	      is_zero = false;

	if(is_zero)
	  break;

	// -- loop over each zeta --
	for(int iz = 0; iz < isub->size_zeta(); iz++) {
	  for(int jz = 0; jz < jsub->size_zeta(); jz++) {
	    dcomplex zetai, zetaj;
	    zetai = isub->zeta_iz[iz]; zetaj = jsub->zeta_iz[jz];
	    dcomplex zetaP = zetai + zetaj;
	    int niat(isub->size_at()); int njat(jsub->size_at());
	    int nipn(isub->size_pn()); int njpn(jsub->size_pn());

	    // -- primitive basis --
	    MultArray<dcomplex, 4> s_prim(bufs, 0, niat, 0, nipn, 0, njat, 0, njpn);
	    MultArray<dcomplex, 4> t_prim(buft, 0, niat, 0, nipn, 0, njat, 0, njpn);
	    MultArray<dcomplex, 4> v_prim(bufv, 0, niat, 0, nipn, 0, njat, 0, njpn);

	    for(int iat = 0; iat < niat; iat++) {
	      for(int jat = 0; jat < njat; jat++) { 
		dcomplex xi, xj, yi, yj, zi, zj, wPx, wPy, wPz;
		xi = isub->xyz_iat(0, iat); xj = jsub->xyz_iat(0, jat);
		yi = isub->xyz_iat(1, iat); yj = jsub->xyz_iat(1, jat);
		zi = isub->xyz_iat(2, iat); zj = jsub->xyz_iat(2, jat);
		wPx = (zetai*xi+zetaj*xj)/zetaP;
		wPy = (zetai*yi+zetaj*yj)/zetaP;
		wPz = (zetai*zi+zetaj*zj)/zetaP;
		dcomplex d2 = pow(xi-xj,2) + pow(yi-yj,2) + pow(zi-zj,2);	
		dcomplex eAB = exp(-zetai*zetaj/zetaP*d2);
		dcomplex ce = eAB * pow(M_PI/zetaP, 1.5);
		int mi = isub->maxn; int mj = jsub->maxn;
		A3dc dxmap = calc_d_coef(mi,mj+2,mi+mj,zetaP,wPx,xi,xj,dsx_buff);
		A3dc dymap = calc_d_coef(mi,mj+2,mi+mj,zetaP,wPy,yi,yj,dsy_buff);
		A3dc dzmap = calc_d_coef(mi,mj+2,mi+mj,zetaP,wPz,zi,zj,dsz_buff);
		for(int kat = 0; kat < this->size_atom(); kat++) {
		  dcomplex dx = wPx-xyzq_iat(0, kat);
		  dcomplex dy = wPy-xyzq_iat(1, kat);
		  dcomplex dz = wPz-xyzq_iat(2, kat);
		  dcomplex d2p = dx*dx + dy*dy + dz*dz;
		  IncompleteGamma(isub->maxn+jsub->maxn, zetaP * d2p, Fjs_iat[iat]);
		}

		for(int ipn = 0; ipn < nipn; ipn++) {
		  for(int jpn = 0; jpn < njpn; jpn++) {
		    int nxi, nxj, nyi, nyj, nzi, nzj;
		    nxi = isub->ns_ipn(0, ipn); nxj = jsub->ns_ipn(0, jpn);
		    nyi = isub->ns_ipn(1, ipn); nyj = jsub->ns_ipn(1, jpn);
		    nzi = isub->ns_ipn(2, ipn); nzj = jsub->ns_ipn(2, jpn);
		    dcomplex dx00, dy00, dz00, dx02, dy02, dz02;
		    dx00 = dxmap.get_safe(nxi, nxj ,0);
		    dy00 = dymap.get_safe(nyi, nyj ,0);
		    dz00 = dzmap.get_safe(nzi, nzj ,0);
		    dx02 = dxmap.get_safe(nxi, nxj+2 ,0);
		    dy02 = dymap.get_safe(nyi, nyj+2 ,0);
		    dz02 = dzmap.get_safe(nzi, nzj+2 ,0);
		    dcomplex s_ele = dx00 * dy00 * dz00;
		    dcomplex t_ele(0.0);
		    t_ele += -2.0*(2*nxj+2*nyj+2*nzj+3)*zetaj*dx00*dy00*dz00;
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
		    dcomplex v_ele(0.0);
		    for(int nx = 0; nx <= nxi + nxj; nx++)
		      for(int ny = 0; ny <= nyi + nyj; ny++)
			for(int nz = 0; nz <= nzi + nzj; nz++)
			  for(int kat = 0; kat < size_atom(); kat++) {
			    v_ele += (xyzq_iat(3, kat) *
				      dxmap.get_safe(nxi, nxj, nx) *
				      dymap.get_safe(nyi, nyj, ny) *
				      dzmap.get_safe(nzi, nzj, nz) *
				      coef_R(zetaP, wPx, wPy, wPz,
					     xyzq_iat(0, kat),
					     xyzq_iat(1, kat),
					     xyzq_iat(2, kat),
					     nx, ny, nz, 0, Fjs_iat[kat]));
		  
			}
		    s_prim.set_safe(iat, ipn, jat, jpn, ce * s_ele);
		    t_prim.set_safe(iat, ipn, jat, jpn, -0.5* ce * t_ele);
		    v_prim.set_safe(iat, ipn, jat, jpn, -2.0*M_PI/zetaP*eAB * v_ele);
		  }}}}

	    // -- contractions --
	    for(RdsIt irds = isub->rds.begin();
		irds != isub->rds.end();
		++irds) {
	      for(RdsIt jrds = jsub->rds.begin(); jrds != jsub->rds.end();++jrds) {
		if(sym_group.IncludeScalar_2(irds->sym, jrds->sym)) {
		  dcomplex cumsum_s(0.0), cumsum_t(0.0), cumsum_v(0.0);
		  for(int iat = 0; iat < isub->size_at(); iat++) {
		    for(int ipn = 0; ipn < isub->size_pn(); ipn++) {
		      for(int jat = 0; jat < jsub->size_at(); jat++) { 
			for(int jpn = 0; jpn < jsub->size_pn(); jpn++) {
			  dcomplex cc = 
			    irds->coef_iat_ipn(iat, ipn) *
			    jrds->coef_iat_ipn(jat, jpn) * 
			    irds->coef_iz(iz) * 
			    jrds->coef_iz(jz);
			  cumsum_s += cc*s_prim.get_safe(iat, ipn, jat, jpn);
			  cumsum_t += cc*t_prim.get_safe(iat, ipn, jat, jpn);
			  cumsum_v += cc*v_prim.get_safe(iat, ipn, jat, jpn);
			}}}}
		  int i(irds->offset + iz); int j(jrds->offset + jz);
		  int sym(irds->sym);
		  //cout << sym << ", " << i << ", "<< j << endl;
		  mat_map["s"][make_pair(sym, sym)](i, j) = cumsum_s;
		  mat_map["t"][make_pair(sym, sym)](i, j) = cumsum_t;
		  mat_map["v"][make_pair(sym, sym)](i, j) = cumsum_v;
		}
	      }
	    }
	  }
	}
      }
    }
    return mat_map;
  }
}

