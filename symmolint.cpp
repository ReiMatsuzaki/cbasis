#include <iostream>
#include <numeric>
#include "angmoment.hpp"
#include "mol_func.hpp"

#include "symmolint.hpp"

namespace l2func {

  using namespace std;
  using namespace Eigen;
  typedef vector<SubSymGTOs>::iterator SubIt;
  typedef vector<ReductionSets>::iterator RdsIt;
  typedef vector<SubSymGTOs>::const_iterator cSubIt;
  typedef vector<ReductionSets>::const_iterator cRdsIt;
  typedef MultArray<dcomplex, 3> A3dc;
  typedef MultArray<dcomplex, 4> A4dc;

  // ==== Data structure ====
  MultArray<dcomplex, 4> calc_R_coef(dcomplex zetaP,
				     dcomplex wPx, dcomplex wPy, dcomplex wPz,
				     Eigen::MatrixXcd xyzq_kat, dcomplex** Fjs_kat,
				     int mx, int my, int mz, int mat, dcomplex* buf) {

    MultArray<dcomplex, 4> res(buf, 0, mx, 0, my, 0, mz, 0, mat);
    for(int nx = 0; nx <= mx; nx++)
      for(int ny = 0; ny <= my; ny++)
	for(int nz = 0; nz <= mz; nz++)
	  for(int kat = 0; kat < mat; kat++) {
	    dcomplex v = coef_R(zetaP, wPx, wPy, wPz,
				xyzq_kat(0, kat),
				xyzq_kat(1, kat),
				xyzq_kat(2, kat),
				nx, ny, nz, 0, Fjs_kat[kat]);
	    res.set(nx, ny, nz, kat, v);
	  }
    return res;
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
  string ReductionSets::str() const {
    ostringstream oss;
    oss << "==== RedcutionSets ====" << endl;
    oss << "irrep : " << irrep << endl;
    oss << "coef_iat_ipn: " << endl << coef_iat_ipn << endl;
    oss << "coef_iz: " << endl <<  coef_iz << endl;
    oss << "offset:  " << offset << endl;
    return oss.str();
  }
  void ReductionSets::Display() const {
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
  void SubSymGTOs::AddRds(const ReductionSets& _rds) {
    setupq = false;
    rds.push_back(_rds);
  }
  SubSymGTOs::SubSymGTOs(MatrixXcd xyz, MatrixXi ns,
			 vector<ReductionSets> ao, VectorXcd zs) :
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
      vector<ReductionSets> rds;
      rds.push_back(ReductionSets(irrep, cs));
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

	  if(abs(norm2) < pow(10.0, -7.0)) {
	    string msg; SUB_LOCATION(msg);
	    msg += ": Norm is too small";
	    throw runtime_error(msg);
	  }

	  // <<< Primitive GTOs <<<
	  irds->coef_iz(iz) = 1.0/sqrt(norm2);
	}
      }
    }
    // <<< Irrep Adapted GTOs <<<

    delete[] dsx_buff;
    delete[] dsy_buff;
    delete[] dsz_buff;
  }
  
  // ---- Calculation ----
  // -- to be removed
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
  void SymGTOs::CalcMat(BMatSet* res) {

    if(not setupq)
      this->SetUp();

    int max_n(0);
    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(int ipn = 0; ipn < isub->ns_ipn.cols(); ipn++)
	for(int i = 0; i < 3; i++)
	  if(max_n < isub->ns_ipn(i, ipn))
	    max_n = isub->ns_ipn(i, ipn);
    }

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
    
    dcomplex* bufs = new dcomplex[100];    
    dcomplex* buft = new dcomplex[100];
    dcomplex* bufv = new dcomplex[100];
    dcomplex* bufz = new dcomplex[100];
    dcomplex* dsx_buff = new dcomplex[100];
    dcomplex* dsy_buff = new dcomplex[100];
    dcomplex* dsz_buff = new dcomplex[100];
    dcomplex** Fjs_iat;    

    Fjs_iat = new dcomplex*[this->size_atom()];
    for(int iat = 0; iat < this->size_atom(); iat++)
      Fjs_iat[iat] = new dcomplex[2*max_n+1];
    
    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(SubIt jsub = subs.begin(); jsub != subs.end(); ++jsub) {

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
	    MultArray<dcomplex, 4> z_prim(bufz, 0, niat, 0, nipn, 0, njat, 0, njpn);

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
		  IncompleteGamma(isub->maxn+jsub->maxn, zetaP * d2p, Fjs_iat[kat]);
		}

		for(int ipn = 0; ipn < nipn; ipn++) {
		  for(int jpn = 0; jpn < njpn; jpn++) {
		    int nxi, nxj, nyi, nyj, nzi, nzj;
		    nxi = isub->ns_ipn(0, ipn); nxj = jsub->ns_ipn(0, jpn);
		    nyi = isub->ns_ipn(1, ipn); nyj = jsub->ns_ipn(1, jpn);
		    nzi = isub->ns_ipn(2, ipn); nzj = jsub->ns_ipn(2, jpn);
		    dcomplex dx00, dy00, dz00, dx02, dy02, dz02, dz01;
		    dx00 = dxmap.get_safe(nxi, nxj ,0);
		    dy00 = dymap.get_safe(nyi, nyj ,0);
		    dz00 = dzmap.get_safe(nzi, nzj ,0);
		    dz01 = dzmap.get_safe(nzi, nzj+1 ,0);
		    dx02 = dxmap.get_safe(nxi, nxj+2 ,0);
		    dy02 = dymap.get_safe(nyi, nyj+2 ,0);
		    dz02 = dzmap.get_safe(nzi, nzj+2 ,0);
		    dcomplex s_ele = dx00 * dy00 * dz00;
		    dcomplex z_ele = dx00 * dy00 * (dz01 + zj *dz00);
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
		    z_prim.set_safe(iat, ipn, jat, jpn, ce*z_ele);
		  }}}}

	    // -- contractions --
	    for(RdsIt irds = isub->rds.begin();
		irds != isub->rds.end();
		++irds) {
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
			cumsum_s += cc*s_prim.get_safe(iat, ipn, jat, jpn);
			cumsum_t += cc*t_prim.get_safe(iat, ipn, jat, jpn);
			cumsum_v += cc*v_prim.get_safe(iat, ipn, jat, jpn);
			cumsum_z += cc*z_prim.get_safe(iat, jpn, jat, jpn);
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
	}
      }
    }
    delete[] bufs;
    delete[] buft;
    delete[] bufv;
    delete[] bufz;
    delete[] dsx_buff;
    delete[] dsy_buff;
    delete[] dsz_buff;
    for(int iat = 0; iat < this->size_atom(); iat++)
      delete[] Fjs_iat[iat];
    delete[] Fjs_iat;
    *res = mat_map;
  }
  bool IsCenter(SubIt isub, int iat, double eps) {
    dcomplex x  = isub->xyz_iat(0, iat);
    dcomplex y  = isub->xyz_iat(1, iat);
    dcomplex z  = isub->xyz_iat(2, iat);
    dcomplex a2 = x*x+y*y+z*z;
    dcomplex a = sqrt(a2);
    return (abs(a) < eps);
  }
  dcomplex AtR_Ylm_cen(SubIt isub, int iz,
		       int ipn, dcomplex r, int L, int M) {

    int nx = isub->ns_ipn(0, ipn);
    int ny = isub->ns_ipn(1, ipn);
    int nz = isub->ns_ipn(2, ipn);
    int nn = nx + ny + nz;
    dcomplex zeta = isub->zeta_iz(iz);
    if(L == 0) {
      // -- s-GTO --
      if(nn == 0) {
	return sqrt(4.0*M_PI)*r*exp(-zeta*r*r);
      } else {
	return 0.0;
      }
    }
    if(L == 1) {
      // -- p-GTO --
      if(nx == 0 && ny == 0 && nz == 1 && M == 0)
	return sqrt(4.0*M_PI/3.0)*r*r*exp(-zeta*r*r);	  
      if(nx == 0 && ny == 1 && nz == 0 && M == 1)
	throw runtime_error("0101 is not implemented");
      if(nx == 1 && ny == 0 && nz == 0 && M ==-1)
	throw runtime_error("0101 is not implemented");	
      return 0.0;
    }

    string msg; SUB_LOCATION(msg);
    msg += "L>1 is not implemented";
    throw runtime_error(msg);
    
  }
  dcomplex AtR_Ylm_noncen(SubIt isub, int iz, int iat, int ipn,
			  dcomplex r, int L, int M) {

    dcomplex* il  = new dcomplex[L+1];
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
    ModSphericalBessel(2.0*zeta*a*r, L, il);
    RealSphericalHarmonics(theta, phi, L, ylm);
    if(nn == 0) {
      res = r * (4.0*M_PI) * il[L] * expz * pow(-1.0, M) * ylm[lm_index(L, -M)];
    } else {
      string msg; SUB_LOCATION(msg);
      msg += ": not implemented yet for p or higher orbital";
      throw runtime_error(msg);
    }		  

    delete[] il;
    delete[] ylm;
    return res;
    
  }
			      
  void SymGTOs::AtR_Ylm(int L, int M, int irrep, const VectorXcd& cs_ibasis,
			const VectorXcd& rs, VectorXcd* res ) {

    if(not setupq) 
      this->SetUp();

    try{
      sym_group.CheckIrrep(irrep);
    } catch(const runtime_error& e) {
      string msg; SUB_LOCATION(msg);
      msg += ": Invalid irreq.\n";
      msg += e.what(); throw runtime_error(msg);
    }

    if(cs_ibasis.size() != this->size_basis_isym(irrep)) {
      string msg; SUB_LOCATION(msg);
      msg += ": size of cs must be equal to basis size";
      throw runtime_error(msg);
    }

    VectorXcd vs = VectorXcd::Zero(rs.size()); // copy
    dcomplex* ylm = new dcomplex[num_lm_pair(L)];
    dcomplex* il  = new dcomplex[L+1];
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
		if(IsCenter(isub, iat, eps)) {
		  vs[ir] += c * AtR_Ylm_cen(isub, iz, ipn, rs[ir], L, M);
		} else {
		  vs[ir] += c * AtR_Ylm_noncen(isub, iz, iat, ipn, rs[ir], L, M);
		}
	      }
	    }
	  }
	}
      }
    }
    delete[] ylm;
    delete[] il;
    res->swap(vs);
  }
  
}

