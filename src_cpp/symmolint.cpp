#include <iostream>
#include <numeric>
#include "angmoment.hpp"
#include "mol_func.hpp"

#include "symmolint.hpp"

namespace l2func {

  using namespace std;
  using namespace Eigen;
  typedef vector<SubSymGTOs>::iterator SubIt;
  typedef vector<Reduction>::iterator RdsIt;
  typedef vector<SubSymGTOs>::const_iterator cSubIt;
  typedef vector<Reduction>::const_iterator cRdsIt;
  typedef MultArray<dcomplex, 2> A2dc;
  typedef MultArray<dcomplex, 3> A3dc;
  typedef MultArray<dcomplex, 4> A4dc;

  // ==== coef R ====
  // -- more efficient code can be used!
  void calc_R_coef(dcomplex zetaP,
		   dcomplex wPx, dcomplex wPy, dcomplex wPz,
		   Eigen::MatrixXcd xyzq_kat, A2dc& Fjs_kat,
		   int mx, int my, int mz, int mat, A4dc& res) {
    
    res.SetRange(0, mx, 0, my, 0, mz, 0, mat);
    for(int nx = 0; nx <= mx; nx++)
      for(int ny = 0; ny <= my; ny++)
	for(int nz = 0; nz <= mz; nz++)
	  for(int kat = 0; kat < mat; kat++) {
	    dcomplex v = coef_R(zetaP, wPx, wPy, wPz,
				xyzq_kat(0, kat),
				xyzq_kat(1, kat),
				xyzq_kat(2, kat),
				nx, ny, nz, 0, &Fjs_kat(kat, 0));
	    res(nx, ny, nz, kat) = v;
	  }

  }

  // ==== SymmetryGroup ====
  SymmetryGroup::SymmetryGroup(int order, string name) {
    if(order != 2 && order != 1) {
      string msg; SUB_LOCATION(msg);
      msg += "Now, only order == 1 and 2 is supported.";
      throw runtime_error(msg);	
    }
    order_ = order;
    name_ = name;
  }
  void SymmetryGroup::CheckIrrep(Irrep a) const {
    if(a < 0 || this->order() <= a) {
      std::string msg; SUB_LOCATION(msg);
      throw runtime_error(msg);
    }
  }
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
  string SymmetryGroup::str() const {
    ostringstream oss; 
    oss << "==== SymmetryGroup ====" << endl;
    oss << "name : " << name_ << endl;
    oss << "order: " << order_ << endl;
    return oss.str();
  }
  void SymmetryGroup::Display() const {
    cout << this->str(); 
  }
  
  SymmetryGroup SymmetryGroup_Cs() {
    SymmetryGroup cs(2, "Cs");
    return cs;
  }
  SymmetryGroup SymmetryGroup_C1() {
    SymmetryGroup c1(1, "C1");
    return c1;
  }
  Irrep Cs_Ap() { return 0;}
  Irrep Cs_App() { return 1; }

  // ==== Reduction Sets ====
  string Reduction::str() const {
    ostringstream oss;
    oss << "==== RedcutionSets ====" << endl;
    oss << "irrep : " << irrep << endl;
    oss << "coef_iat_ipn: " << endl << coef_iat_ipn << endl;
    oss << "coef_iz: " << endl <<  coef_iz << endl;
    oss << "offset:  " << offset << endl;
    return oss.str();
  }
  void Reduction::Display() const {
    cout << this->str() ;
  }

  
  // ==== Sub ====
  SubSymGTOs::SubSymGTOs() {
    xyz_iat = MatrixXcd::Zero(3, 0);
    ns_ipn = MatrixXi::Zero(3, 0);
    zeta_iz = VectorXcd::Zero(0);
    setupq = false;
  }
  void SubSymGTOs::SetUp() {

    // ---- check values ----
    for(RdsIt it = rds.begin(); it != rds.end(); ++it) {
      if(it->size_at() != this->size_at() |
	 it->size_pn() != this->size_pn()) {
	string msg; SUB_LOCATION(msg);
	ostringstream oss; oss << msg;
	oss << ": size mismatch.\n"
	    << "size_at (ReductionSets, SubSymGTOs) = "
	    << it->size_at() << this->size_at()
	    << "size_pn (ReductionSets, SubSymGTOs) = "
	    << it->size_pn() << this->size_pn() << endl;
	throw runtime_error(oss.str());
      }
    }

    if(xyz_iat.rows() != 3 || ns_ipn.rows() != 3) {
      string msg; SUB_LOCATION(msg);
      throw runtime_error(msg);
    }

    // ---- compute internal values ----
    maxn = 0;
    for(int ipn = 0; ipn < ns_ipn.cols(); ipn++) {
      if(maxn < ns_ipn(0, ipn) + ns_ipn(1, ipn) + ns_ipn(2, ipn))
	maxn = ns_ipn(0, ipn) + ns_ipn(1, ipn) + ns_ipn(2, ipn);
    }
    for(RdsIt it = rds.begin(); it != rds.end(); ++it) {
      it->set_zs_size(zeta_iz.size());
    }



    setupq = true;
  }
  void SubSymGTOs::AddXyz(Vector3cd xyz) {

    setupq = false;
    int num_atom = xyz_iat.cols();
    MatrixXcd res(3, num_atom + 1);

    for(int i = 0; i < num_atom; i++) {
      for(int j = 0; j < 3; j++)
	res(j, i) = xyz_iat(j, i);
    }

    res(0, num_atom) = xyz(0);
    res(1, num_atom) = xyz(1);
    res(2, num_atom) = xyz(2);

    xyz_iat.swap(res);    
  }
  void SubSymGTOs::AddNs(Vector3i ns) {

    setupq = false;

    int num = ns_ipn.cols();
    MatrixXi res(3, num+ 1);

    for(int i = 0; i < num; i++) {
      for(int j = 0; j < 3; j++)
	res(j, i) = ns_ipn(j, i);
    }

    res(0, num) =ns(0);
    res(1, num) =ns(1);
    res(2, num) =ns(2);

    ns_ipn.swap(res);    
  }
  void SubSymGTOs::AddZeta(const VectorXcd& zs) {
    setupq = false;
    VectorXcd res(zeta_iz.size() + zs.size());
    for(int i = 0; i < zeta_iz.size(); i++)
      res(i) = zeta_iz(i);
    for(int i = 0; i < zs.size(); i++)
      res(i+zeta_iz.size()) = zs(i);
    zeta_iz.swap(res);
  }
  void SubSymGTOs::AddRds(const Reduction& _rds) {
    setupq = false;
    rds.push_back(_rds);
  }
  SubSymGTOs::SubSymGTOs(MatrixXcd xyz, MatrixXi ns,
			 vector<Reduction> ao, VectorXcd zs) :
    xyz_iat(xyz), ns_ipn(ns), zeta_iz(zs), rds(ao) {

    setupq = false;

    this->SetUp();

  }
  string SubSymGTOs::str() const {
    ostringstream oss;
    oss << "==== SubSymGTOs ====" << endl;
    oss << "xyz : " << endl <<  xyz_iat << endl;
    oss << "ns  : " << endl <<  ns_ipn << endl;
    oss << "zeta: " << endl <<  zeta_iz << endl;
    for(cRdsIt it = rds.begin(); it != rds.end(); ++it)
      oss << it->str();
    return oss.str();
  }
  void SubSymGTOs::Display() const {
    cout << this->str();
  }
  SubSymGTOs Sub_s(Irrep sym, Vector3cd xyz, VectorXcd zs) {

    MatrixXcd xyz_in(3, 1);  xyz_in << xyz[0] , xyz[1] , xyz[2];
    MatrixXi ns  = MatrixXi::Zero(3, 1);
    MatrixXcd cs = MatrixXcd::Ones(1,1);
    Reduction rds(sym, cs);
    vector<Reduction> rds_list; rds_list.push_back(rds);

    SubSymGTOs sub(xyz_in, ns, rds_list, zs);
    return sub;
  }
  SubSymGTOs Sub_pz(Irrep sym, Vector3cd xyz, VectorXcd zs) {

    MatrixXcd xyz_in(3, 1);  xyz_in << xyz[0] , xyz[1] , xyz[2];
    MatrixXi ns(3, 1); ns << 0, 0, 1;
    MatrixXcd cs = MatrixXcd::Ones(1,1);
    Reduction rds(sym, cs);
    vector<Reduction> rds_list; rds_list.push_back(rds);

    SubSymGTOs sub(xyz_in, ns, rds_list, zs);
    return sub;    
  }
  SubSymGTOs Sub_TwoSGTO(SymmetryGroup sym, Irrep irrep,
			 Vector3cd xyz, VectorXcd zs) {

    sym.CheckIrrep(irrep);

    MatrixXi ns_in   = MatrixXi::Zero(3, 1);
    if(sym.name() == "Cs") {
      MatrixXcd xyz_in(3, 2) ;
      xyz_in << xyz(0), xyz(0), xyz(1), xyz(1), xyz(2), -xyz(2);
      MatrixXcd cs(2, 1);
      if(irrep == Cs_Ap()) {
	cs << 1.0, 1.0;
      } else if (irrep == Cs_App()){
	cs << 1.0, -1.0;
      }
      vector<Reduction> rds;
      rds.push_back(Reduction(irrep, cs));
      SubSymGTOs sub(xyz_in, ns_in, rds, zs);
      return sub;

    } else {
      string msg; SUB_LOCATION(msg);
      msg += "only Cs symmetry is implemented now";
      throw runtime_error(msg);
    }
  }
  
 
  // ==== SymGTOs ====
  // ---- Constructors ----
  SymGTOs::SymGTOs(SymmetryGroup _sym_group):
    sym_group(_sym_group), setupq(false)  {
    xyzq_iat = MatrixXcd::Zero(4, 0);
  }

  // ---- Accessors ----
  int SymGTOs::size_atom() const {return xyzq_iat.cols(); }
  int SymGTOs::size_basis_isym(Irrep isym) const {

    int cumsum(0);
    for(cSubIt isub = subs.begin(); isub != subs.end(); ++isub) {

      for(cRdsIt irds = isub->rds.begin(); irds != isub->rds.end(); 
	  ++irds) {

	if(irds->irrep == isym) 
	  cumsum += isub->size_zeta();
      }
    }
    return cumsum;

  }
  string SymGTOs::str() const {
    ostringstream oss;
    oss << "==== SymGTOs ====" << endl;
    oss << "Set Up?" << (setupq ? "Yes" : "No") << endl;
    oss << sym_group.str();
    for(cSubIt it = subs.begin(); it != subs.end(); ++it) 
      oss << it->str();
    oss << "xyzq:" << endl;
    oss << xyzq_iat << endl;
    return oss.str();
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
  void SymGTOs::AddAtom(Eigen::Vector3cd _xyz, dcomplex q) {

    int num_atom = xyzq_iat.cols();
    MatrixXcd res(4, num_atom + 1);

    for(int i = 0; i < num_atom; i++) {
      for(int j = 0; j < 4; j++)
	res(j, i) = xyzq_iat(j, i);
    }

    res(0, num_atom) = _xyz(0);
    res(1, num_atom) = _xyz(1);
    res(2, num_atom) = _xyz(2);
    res(3, num_atom) = q;

    xyzq_iat.swap(res);

  }
  void SymGTOs::AddSub(SubSymGTOs sub) {

    subs.push_back(sub);

  }

  // ---- SetUp ----
  void SymGTOs::SetUp() {
    
    for(SubIt it = subs.begin(); it != subs.end(); ++it) {
      // -- setup sub --
      if(it->setupq == false)
	it->SetUp();

      // -- init offset --
      for(RdsIt irds = it->begin_rds(); irds != it->end_rds(); ++irds) {
	sym_group.CheckIrrep(irds->irrep);
	irds->offset = 0;
      }
    }

    this->SetOffset();
    this->Normalize();
    setupq = true;
  }
  void SymGTOs::SetOffset() {
    
    vector<int> num_irrep(10, 0);
    
    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(RdsIt irds = isub->rds.begin(), end = isub->rds.end();
	  irds != end; ++irds) {
	irds->offset = num_irrep[irds->irrep];
	num_irrep[irds->irrep] += isub->size_zeta();
      }
    }

  }
  void SymGTOs::Normalize() {

    A3dc dxmap(100), dymap(100), dzmap(100);

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
		calc_d_coef(mi,mi,0, zetaP,wPx,xi,xj,dxmap);
		calc_d_coef(mi,mi,0, zetaP,wPy,yi,yj,dymap);
		calc_d_coef(mi,mi,0, zetaP,wPz,zi,zj,dzmap);
		
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

	  if(abs(norm2) < pow(10.0, -14.0)) {
	    string msg; SUB_LOCATION(msg);
	    ostringstream oss; oss << msg << ": " << endl;
	    oss << "norm is too small" << endl;
	    oss << "iz  : " << iz << endl;
	    oss << "zeta: " << zetai << endl;
	    oss << "norm2: " << norm2 << endl;	    
	    throw runtime_error(oss.str());
	  }

	  // <<< Primitive GTOs <<<
	  irds->coef_iz(iz) = 1.0/sqrt(norm2);
	}
      }
    }
    // <<< Irrep Adapted GTOs <<<

  }

  // ---- utils ----
  // -- not used now --
  int SymGTOs::max_n() const {
    int max_n(0);
    for(cSubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(int ipn = 0; ipn < isub->ns_ipn.cols(); ipn++)
	for(int i = 0; i < 3; i++)
	  if(max_n < isub->ns_ipn(i, ipn))
	    max_n = isub->ns_ipn(i, ipn);
    }    
    return max_n;
  }
  
  // ---- CalcMat ----
  struct PrimBasis {
    A4dc s, t, v, z;
    PrimBasis(int num):
      s(num), t(num), v(num), z(num) {}
    void SetRange(int niat, int nipn, int njat, int njpn) {
      s.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
      t.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
      v.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
      z.SetRange(0, niat, 0, nipn, 0, njat, 0, njpn);
    }
  };
  dcomplex dist2(dcomplex x, dcomplex y, dcomplex z) {
    return x*x + y*y + z*z;
  } 
  dcomplex calc_tele(SubIt isub, SubIt jsub, dcomplex zetaj, 
		     int ipn, int jpn,
		     A3dc& dxmap, A3dc& dymap, A3dc& dzmap) {
    int nxi, nxj, nyi, nyj, nzi, nzj;
    nxi = isub->nx(ipn); nxj = jsub->nx(jpn);
    nyi = isub->ny(ipn); nyj = jsub->ny(jpn);
    nzi = isub->nz(ipn); nzj = jsub->nz(jpn);
    dcomplex dx00, dy00, dz00, dx02, dy02, dz02, dz01;
    dx00 = dxmap.get_safe(nxi, nxj ,0);
    dy00 = dymap.get_safe(nyi, nyj ,0);
    dz00 = dzmap.get_safe(nzi, nzj ,0);
    dz01 = dzmap.get_safe(nzi, nzj+1 ,0);
    dx02 = dxmap.get_safe(nxi, nxj+2 ,0);
    dy02 = dymap.get_safe(nyi, nyj+2 ,0);
    dz02 = dzmap.get_safe(nzi, nzj+2 ,0);

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
    return t_ele;

  }
  dcomplex calc_vele(const SymGTOs& gtos, SubIt isub, SubIt jsub, int ipn, int jpn,
		     A3dc& dxmap, A3dc& dymap, A3dc& dzmap, A4dc& rmap) {

    int nxi, nxj, nyi, nyj, nzi, nzj;
    nxi = isub->nx(ipn); nxj = jsub->nx(jpn);
    nyi = isub->ny(ipn); nyj = jsub->ny(jpn);
    nzi = isub->nz(ipn); nzj = jsub->nz(jpn);

    dcomplex v_ele(0.0);
    for(int nx = 0; nx <= nxi + nxj; nx++)
      for(int ny = 0; ny <= nyi + nyj; ny++)
	for(int nz = 0; nz <= nzi + nzj; nz++)
	  for(int kat = 0; kat < gtos.size_atom(); kat++) {
	    v_ele += (gtos.q_at(kat) *
		      dxmap(nxi, nxj, nx) *
		      dymap(nyi, nyj, ny) *
		      dzmap(nzi, nzj, nz) *
		      rmap(nx, ny, nz, kat));
		  }
    return v_ele;
  }

  void CalcPrim(const SymGTOs gtos, SubIt isub, SubIt jsub, int iz, int jz,
		A3dc& dxmap, A3dc& dymap, A3dc& dzmap, A4dc& rmap, A2dc& Fjs_iat,
		PrimBasis& prim, bool calc_coulomb) {
    dcomplex zetai, zetaj;
    zetai = isub->zeta_iz[iz]; zetaj = jsub->zeta_iz[jz];
    dcomplex zetaP = zetai + zetaj;
    int niat(isub->size_at()); int njat(jsub->size_at());
    int nipn(isub->size_pn()); int njpn(jsub->size_pn());

    prim.SetRange(niat, nipn, njat, njpn);
    Fjs_iat.SetRange(0, gtos.size_atom(), 0, isub->maxn+jsub->maxn);

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
	  for(int kat = 0; kat < gtos.size_atom(); kat++) {
	    dcomplex d2p = dist2(wPx-gtos.x_at(kat), wPy-gtos.y_at(kat),
				 wPz-gtos.z_at(kat));
	    IncompleteGamma(isub->maxn+jsub->maxn, zetaP * d2p, &Fjs_iat(kat, 0));
	  }
	} else {
	  calc_d_coef(mi,mj+2,0,zetaP,wPx,xi,xj,dxmap);
	  calc_d_coef(mi,mj+2,0,zetaP,wPy,yi,yj,dymap);
	  calc_d_coef(mi,mj+2,0,zetaP,wPz,zi,zj,dzmap);	  
	}
	  
	for(int ipn = 0; ipn < nipn; ipn++) {
	  for(int jpn = 0; jpn < njpn; jpn++) {
	    int nxi, nxj, nyi, nyj, nzi, nzj;
	    nxi = isub->ns_ipn(0, ipn); nxj = jsub->ns_ipn(0, jpn);
	    nyi = isub->ns_ipn(1, ipn); nyj = jsub->ns_ipn(1, jpn);
	    nzi = isub->ns_ipn(2, ipn); nzj = jsub->ns_ipn(2, jpn);

	    dcomplex s_ele = dxmap(nxi,nxj,0) * dymap(nyi,nyj,0) * dzmap(nzi,nzj,0);
	    dcomplex z_ele = dxmap(nxi,nxj,0)*dymap(nyi,nyj,0)*
	      (dzmap(nzi,nzj+1,0)+zj*dzmap(nzi,nzj,0));
	    dcomplex t_ele = calc_tele(isub, jsub, zetaj, ipn, jpn,
				       dxmap, dymap, dzmap);	    

	    prim.s(iat, ipn, jat, jpn) =  ce * s_ele;
	    prim.t(iat, ipn, jat, jpn) =  -0.5* ce * t_ele;
	    prim.z(iat, ipn, jat, jpn) =  ce*z_ele;

	    if(calc_coulomb) {
	      calc_R_coef(zetaP, wPx, wPy, wPz, gtos.xyzq_iat, Fjs_iat,
			  nxi+nxj, nyi+nyj, nzi+nzj, gtos.size_atom(), rmap);
	      dcomplex v_ele = calc_vele(gtos, isub, jsub, ipn, jpn,
					 dxmap, dymap, dzmap, rmap);
	      prim.v(iat, ipn, jat, jpn) =  -2.0*M_PI/zetaP*eAB * v_ele;
	    }
	  }}}}

  }
  void CalcTrans(SubIt isub, SubIt jsub, int iz, int jz,
		 PrimBasis& prim, BMatSet& mat_map) {

    for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) {
      for(RdsIt jrds = jsub->rds.begin(); jrds != jsub->rds.end();++jrds) {
	dcomplex cumsum_s(0.0), cumsum_t(0.0), cumsum_v(0.0);
	dcomplex cumsum_z(0.0);
	for(int iat = 0; iat < isub->size_at(); iat++) {
	  for(int ipn = 0; ipn < isub->size_pn(); ipn++) {
	    for(int jat = 0; jat < jsub->size_at(); jat++) { 
	      for(int jpn = 0; jpn < jsub->size_pn(); jpn++) {
		dcomplex cc = 
		  irds->coef_iat_ipn(iat, ipn) *
		  jrds->coef_iat_ipn(jat, jpn) * 
		  irds->coef_iz(iz) * 
		  jrds->coef_iz(jz);
		cumsum_s += cc*prim.s.get_safe(iat, ipn, jat, jpn);
		cumsum_t += cc*prim.t.get_safe(iat, ipn, jat, jpn);
		cumsum_v += cc*prim.v.get_safe(iat, ipn, jat, jpn);
		cumsum_z += cc*prim.z.get_safe(iat, jpn, jat, jpn);
	      }}}}
	int i(irds->offset + iz); int j(jrds->offset + jz);
	int isym(irds->irrep); int jsym(jrds->irrep);
		
	mat_map.SelfAdd("s", isym, jsym, i, j, cumsum_s);
	mat_map.SelfAdd("t", isym, jsym, i, j, cumsum_t);
	mat_map.SelfAdd("v", isym, jsym, i, j, cumsum_v);
	mat_map.SelfAdd("z", isym, jsym, i, j, cumsum_z);
      }
    }
  }

  void SymGTOs::loop() {

    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(SubIt jsub = subs.begin(); jsub != subs.end(); ++jsub) {

	// -- loop over each zeta --
	for(int iz = 0; iz < isub->size_zeta(); iz++) {
	  for(int jz = 0; jz < jsub->size_zeta(); jz++) {
	    dcomplex zetai, zetaj;
	    zetai = isub->zeta_iz[iz]; zetaj = jsub->zeta_iz[jz];
	    int niat(isub->size_at()); int njat(jsub->size_at());
	    int nipn(isub->size_pn()); int njpn(jsub->size_pn());

	    // -- primitive basis --
	    for(int iat = 0; iat < niat; iat++) {
	      for(int jat = 0; jat < njat; jat++) { 
		for(int ipn = 0; ipn < nipn; ipn++) {
		  for(int jpn = 0; jpn < njpn; jpn++) {

		  }}}}

	    // -- contractions --
	    for(RdsIt irds = isub->rds.begin();irds != isub->rds.end();++irds) {
	      for(RdsIt jrds = jsub->rds.begin(); jrds != jsub->rds.end();++jrds) {
	      }}

	  }}
      }}
    }
  void SymGTOs::CalcMatOther(SymGTOs& o, bool calc_coulomb, BMatSet* res) {

    if(not setupq)
      this->SetUp();

    if(not o.setupq)
      o.SetUp();

    BMatSet mat_map(sym_group.order());
    for(Irrep isym = 0; isym < sym_group.order(); isym++) {
      for(Irrep jsym = 0; jsym < o.sym_group.order(); jsym++) {
	int numi = this->size_basis_isym(isym);
	int numj = o.size_basis_isym(jsym);
	MatrixXcd s = MatrixXcd::Zero(numi, numj);
	mat_map.SetMatrix("s", isym, jsym, s);
	MatrixXcd t = MatrixXcd::Zero(numi, numj); 
	mat_map.SetMatrix("t", isym, jsym, t);
	MatrixXcd v = MatrixXcd::Zero(numi, numj); 
	mat_map.SetMatrix("v", isym, jsym, v);
	MatrixXcd z = MatrixXcd::Zero(numi, numj); 
	mat_map.SetMatrix("z", isym, jsym, z);
      }
    }

    PrimBasis prim(100);
    A4dc rmap(100);
    A3dc dxmap(100),  dymap(100),  dzmap(100);
    A2dc Fjs_iat(100);
    
    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(SubIt jsub = o.subs.begin(); jsub != o.subs.end(); ++jsub) {
	for(int iz = 0; iz < isub->size_zeta(); iz++) {
	  for(int jz = 0; jz < jsub->size_zeta(); jz++) {
	    CalcPrim(*this, isub, jsub, iz, jz, dxmap, dymap, dzmap, rmap,
		     Fjs_iat, prim, calc_coulomb);
	    CalcTrans(isub, jsub, iz, jz, prim, mat_map);
	  }
	}
      }
    }
    *res = mat_map;    

  }
  void SymGTOs::CalcMat(BMatSet* res) {
    this->CalcMatOther(*this, true, res);
/*
    if(not setupq)
      this->SetUp();

    BMatSet mat_map(sym_group.order());
    for(Irrep isym = 0; isym < sym_group.order(); isym++) {
      for(Irrep jsym = 0; jsym < sym_group.order(); jsym++) {
	int numi = this->size_basis_isym(isym);
	int numj = this->size_basis_isym(jsym);
	MatrixXcd s = MatrixXcd::Zero(numi, numj);
	mat_map.SetMatrix("s", isym, jsym, s);
	MatrixXcd t = MatrixXcd::Zero(numi, numj); 
	mat_map.SetMatrix("t", isym, jsym, t);
	MatrixXcd v = MatrixXcd::Zero(numi, numj); 
	mat_map.SetMatrix("v", isym, jsym, v);
	MatrixXcd z = MatrixXcd::Zero(numi, numj); 
	mat_map.SetMatrix("z", isym, jsym, z);
      }
    }

    PrimBasis prim(100);
    A4dc rmap(100);
    A3dc dxmap(100),  dymap(100),  dzmap(100);
    A2dc Fjs_iat(100);
    
    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(SubIt jsub = subs.begin(); jsub != subs.end(); ++jsub) {
	for(int iz = 0; iz < isub->size_zeta(); iz++) {
	  for(int jz = 0; jz < jsub->size_zeta(); jz++) {
	    CalcPrim(*this, isub, jsub, iz, jz, dxmap, dymap, dzmap, rmap,
		     Fjs_iat, prim);
	    CalcTrans(isub, jsub, iz, jz, prim, mat_map);
	  }
	}
      }
    }
    *res = mat_map;
    */
  }

  // ---- AtR ----
  bool IsCenter(SubIt isub, int iat, double eps) {
    dcomplex x  = isub->xyz_iat(0, iat);
    dcomplex y  = isub->xyz_iat(1, iat);
    dcomplex z  = isub->xyz_iat(2, iat);
    dcomplex a2 = x*x+y*y+z*z;
    dcomplex a = sqrt(a2);
    return (abs(a) < eps);
  }
  void AtR_Ylm_cen(SubIt isub, int iz,
		   int ipn, dcomplex r, int L, int M, dcomplex* v, dcomplex* dv) {

    int nx = isub->ns_ipn(0, ipn);
    int ny = isub->ns_ipn(1, ipn);
    int nz = isub->ns_ipn(2, ipn);
    int nn = nx + ny + nz;
    dcomplex zeta = isub->zeta_iz(iz);
    if(L == 0) {
      // -- s-GTO --
      if(nn == 0) {
	dcomplex c(sqrt(4.0*M_PI));
	*v = c*r*exp(-zeta*r*r);
	*dv= c*(1.0 -2.0*zeta*r*r) * exp(-zeta*r*r);
      } else {
	*v = 0.0;
	*dv= 0.0;
      }
    } else if(L == 1) {
      // -- p-GTO --
      if(nx == 0 && ny == 0 && nz == 1 && M == 0) {
	dcomplex c(sqrt(4.0*M_PI/3.0));
	*v = c*r*r*exp(-zeta*r*r);
	*dv= c*(2.0*r - 2.0*zeta*r*r*r )*exp(-zeta*r*r);
      }	  
      else if(nx == 0 && ny == 1 && nz == 0 && M == 1)
	throw runtime_error("0101 is not implemented");
      else if(nx == 1 && ny == 0 && nz == 0 && M ==-1)
	throw runtime_error("0101 is not implemented");	
      else {
	*v = 0.0;
	*dv= 0.0;
      }
    } else {
      string msg; SUB_LOCATION(msg);
      msg += "L>1 is not implemented";
      throw runtime_error(msg);
    }
    
  }
  void AtR_Ylm_noncen(SubIt isub, int iz, int iat, int ipn,
		      dcomplex r, int L, int M, dcomplex* v, dcomplex *dv) {

    dcomplex* il_ipl  = new dcomplex[2*L+2];
    dcomplex* il = &il_ipl[0];
    dcomplex* ipl= &il_ipl[L+1];
    dcomplex* ylm = new dcomplex[num_lm_pair(L)];

    dcomplex res;

    int nn = (isub->ns_ipn(0, ipn) +
	      isub->ns_ipn(1, ipn) +
	      isub->ns_ipn(2, ipn));
    dcomplex zeta = isub->zeta_iz(iz);
    dcomplex x  = isub->xyz_iat(0, iat);
    dcomplex y  = isub->xyz_iat(1, iat);
    dcomplex z  = isub->xyz_iat(2, iat);
    dcomplex xxyy = x*x+y*y;
    dcomplex a2 = xxyy+z*z;
    dcomplex a = sqrt(a2);

    if(abs(a) < 0.000001) {
      string msg; SUB_LOCATION(msg);
      msg += "(x,y,z) is not centered on origin."; 
      throw runtime_error(msg);
    }

    dcomplex expz = exp(-zeta*(r*r+a2));
    dcomplex theta = acos(z / a);
    dcomplex phi   = (abs(xxyy) < 0.00001) ? 0.0 : acos(x / sqrt(xxyy));
    ModSphericalBessel(2.0*zeta*a*r, L, il_ipl);
    RealSphericalHarmonics(theta, phi, L, ylm);
    if(nn == 0) {
      
      dcomplex c((4.0*M_PI) * pow(-1.0, M) * ylm[lm_index(L, -M)]);
      *v = c * r * il[L] * expz;
      *dv = (c*   il[L]*expz +
	     c*r*ipl[L]*expz*(+2.0*zeta*a) +
	     c*r* il[L]*expz*(-2.0*zeta*r));
    } else {
      string msg; SUB_LOCATION(msg);
      msg += ": not implemented yet for p or higher orbital";
      throw runtime_error(msg);
    }		  

    delete[] il_ipl;
    delete[] ylm;
    
  }      
  void SymGTOs::AtR_Ylm(int L, int M, int irrep, const VectorXcd& cs_ibasis,
			const VectorXcd& rs,
			VectorXcd* res_vs, VectorXcd* res_dvs ) {

    if(not setupq) 
      this->SetUp();

    if(!is_lm_pair(L, M)) {
      string msg; SUB_LOCATION(msg);
      msg += ": invalid L,M pair";
      throw runtime_error(msg);
    }

    try{
      sym_group.CheckIrrep(irrep);
    } catch(const runtime_error& e) {
      string msg; SUB_LOCATION(msg);
      msg += ": Invalid irreq.\n";
      msg += e.what(); throw runtime_error(msg);
    }

    if(cs_ibasis.size() != this->size_basis_isym(irrep)) {
      string msg; SUB_LOCATION(msg);
      ostringstream oss;
      oss << msg << endl << ": size of cs must be equal to basis size" << endl;
      oss << "(L, M, irrep)     = " << L << M << irrep << endl;
      oss << "csibasis.size()   = " << cs_ibasis.size() << endl;
      oss << "size_basis_isym() = " << this->size_basis_isym(irrep) << endl;
      oss << endl;
      oss << "print SymGTOs:" << endl;
      oss << this->str();
      
      throw runtime_error(oss.str());
    }

    VectorXcd vs  = VectorXcd::Zero(rs.size());   // copy
    VectorXcd dvs = VectorXcd::Zero(rs.size());   // copy
    dcomplex* ylm = new dcomplex[num_lm_pair(L)];
    dcomplex* il  = new dcomplex[2*L+2];
    double eps(0.0000001);

    // Y00   = 1/sqrt(4pi)
    // r Y10 = sqrt(3/4pi) z
    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) {
	if(irds->irrep != irrep) 
	  continue;
	
	for(int iz = 0; iz < isub->size_zeta(); iz++) {
	  for(int iat = 0; iat < isub->size_at(); iat++) {
	    for(int ipn = 0; ipn < isub->size_pn(); ipn++) {
	      for(int ir = 0; ir < rs.size(); ir++) {
		int ibasis = irds->offset + iz;
		dcomplex c = (cs_ibasis(ibasis) *
			      irds->coef_iat_ipn(iat, ipn) *
			      irds->coef_iz(iz));
		dcomplex v, dv;
		if(IsCenter(isub, iat, eps)) {
		  AtR_Ylm_cen(isub, iz, ipn, rs[ir], L, M, &v, &dv);
		} else {
		  AtR_Ylm_noncen(isub, iz, iat, ipn, rs[ir], L, M, &v, &dv);
		}
		vs[ir] += c * v; dvs[ir]+= c * dv;
	      }
	    }
	  }
	}
      }
    }
    delete[] ylm;
    delete[] il;
    res_vs->swap(vs);
    res_dvs->swap(dvs);
  }

  // ---- Correct Sign ----
  void SymGTOs::CorrectSign(int L, int M, int irrep, Eigen::VectorXcd& cs) {
			    

    if(not setupq) {
      string msg; SUB_LOCATION(msg);
      msg += ": call SetUp() before calculation.";
      throw runtime_error(msg);
    }

    if(!is_lm_pair(L, M)) {
      string msg; SUB_LOCATION(msg);
      msg += ": invalid L,M pair";
      throw runtime_error(msg);
    }

    try{
      sym_group.CheckIrrep(irrep);
    } catch(const runtime_error& e) {
      string msg; SUB_LOCATION(msg);
      msg += ": Invalid irreq.\n";
      msg += e.what(); throw runtime_error(msg);
    }

    if(cs.size() != this->size_basis_isym(irrep)) {
      string msg; SUB_LOCATION(msg);
      msg += ": size of cs must be equal to basis size";
      throw runtime_error(msg);
    }

    VectorXcd rs(1); rs << 0.0;
    VectorXcd vs, ds;
    this->AtR_Ylm(L, M, irrep, cs, rs, &vs, &ds);

    if(ds(0).real() < 0) {
      cs = -cs;
    }
  }
  
}

