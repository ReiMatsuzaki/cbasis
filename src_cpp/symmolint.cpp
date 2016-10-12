#include <iostream>
#include <stdexcept>
#include <numeric>
#include "cfunc.hpp"
#include "angmoment.hpp"
#include "mol_func.hpp"
#include "one_int.hpp"
#include "two_int.hpp"
#include "symmolint.hpp"

namespace l2func {

  using namespace std;
  using namespace Eigen;
  typedef vector<SubSymGTOs>::iterator SubIt;
  typedef vector<Reduction>::iterator RdsIt;
  typedef vector<SubSymGTOs>::const_iterator cSubIt;
  typedef vector<Reduction>::const_iterator cRdsIt;
  typedef MultArray<dcomplex, 1> A1dc;
  typedef MultArray<dcomplex, 2> A2dc;
  typedef MultArray<dcomplex, 3> A3dc;
  typedef MultArray<dcomplex, 4> A4dc;


  // ==== ERI method ====
  ERIMethod::ERIMethod(): symmetry(0), coef_R_memo(0), perm(0) {}
  void ERIMethod::set_symmetry(int s) {symmetry = s; }
  void ERIMethod::set_coef_R_memo(int s) {coef_R_memo = s; }
  void ERIMethod::set_perm(int s) {perm = s; }

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
    zeta_iz = VectorXcd::Zero(0);
    setupq = false;
  }
  void SubSymGTOs::SetUp() {

    // ---- check values ----
    if(sym_group.get() == NULL) {
      string msg; SUB_LOCATION(msg);
      msg += ": sym_group is not set.";
      throw runtime_error(msg);
    }

    if(this->size_at() * this->size_pn() == 0) {
      string msg; SUB_LOCATION(msg);
      msg += ": primitive GTO is not set.";
      throw runtime_error(msg);
    }
    if(this->zeta_iz.size() == 0) {
      string msg; SUB_LOCATION(msg);
      msg += ": zeta is not set.";
      throw runtime_error(msg);
    }
    if(this->rds.size() == 0) {
      string msg; SUB_LOCATION(msg);
      msg += ": ReductionSet is not set.";
      throw runtime_error(msg);
    }
    for(RdsIt it = rds.begin(); it != rds.end(); ++it) {
      if(it->size_at() != this->size_at() ||
	 it->size_pn() != this->size_pn()) {
	string msg; SUB_LOCATION(msg);
	ostringstream oss; oss << msg;
	oss << ": size mismatch."
	    << endl
	    << "size_at (Reduction, xyz) = "
	    << "(" << it->size_at() << ", " << this->size_at() << ")"
	    << endl
	    << "size_pn (Reduction, ns) = "
	    << "(" << it->size_pn() << "," << this->size_pn() << ")"
	    << endl;
	throw runtime_error(oss.str());
      }
    }

    // ---- compute internal values ----
    ip_iat_ipn = MatrixXi::Zero(this->size_at(), this->size_pn());
    int ip(0);
    for(int iat = 0; iat < this->size_at(); ++iat)
      for(int ipn = 0; ipn <  this->size_pn(); ++ipn) {
	ip_iat_ipn(iat, ipn) = ip;
	ip++;
      }
    
    maxn = 0;
    for(int ipn = 0; ipn < size_pn(); ipn++) {
      if(maxn < nx(ipn) + ny(ipn) + nz(ipn))
	maxn = nx(ipn) + ny(ipn) + nz(ipn);
    }
    for(RdsIt it = rds.begin(); it != rds.end(); ++it) {
      it->set_zs_size(zeta_iz.size());
    }

    // -- build PrimGTO set.
    int nat(this->size_at());
    int npn(this->size_pn());
    vector<PrimGTO> gtos(nat * npn);
    for(int iat = 0; iat < nat; iat++) 
      for(int ipn = 0; ipn < npn; ipn++)  {
	int ip = this->ip_iat_ipn(iat, ipn);
	gtos[ip] = PrimGTO(this->nx(ipn),
			   this->ny(ipn),
			   this->nz(ipn),
			   this->x(iat),
			   this->y(iat),
			   this->z(iat));
      }
    // -- compute symmetry operation and primitive GTO relations --
    sym_group->CalcSymMatrix(gtos, this->ip_jg_kp, this->sign_ip_jg_kp);

    // -- set flag --
    setupq = true;
  }
  void SubSymGTOs::AddXyz(Vector3cd xyz) {

    setupq = false;
    this->x_iat.push_back(xyz[0]);
    this->y_iat.push_back(xyz[1]);
    this->z_iat.push_back(xyz[2]);

  }
  void SubSymGTOs::AddNs(Vector3i ns) {
    setupq = false;
    this->nx_ipn.push_back(ns[0]);
    this->ny_ipn.push_back(ns[1]);
    this->nz_ipn.push_back(ns[2]);
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
  string SubSymGTOs::str() const {
    ostringstream oss;
    oss << "==== SubSymGTOs ====" << endl;
    //oss << "sym : " << sym_group->str() << endl;
    oss << "xyz : " << endl;
    for(int iat = 0; iat < size_at(); iat++)
      oss <<  x(iat) << y(iat) << z(iat) << endl;
    oss << "ns  : " << endl;
    for(int ipn = 0; ipn < size_pn(); ipn++) 
      oss <<  nx(ipn) << ny(ipn) << nz(ipn) << endl;
    oss << "zeta: " << endl <<  zeta_iz << endl;
    oss << "maxn: " << maxn << endl;
    for(cRdsIt it = rds.begin(); it != rds.end(); ++it)
      oss << it->str();
    oss << "ip_iat_ipn: " << endl << ip_iat_ipn << endl;
    oss << "ip_jg_kp: " << endl << ip_jg_kp << endl;
    oss << "sign_ip_jg_kp: " << endl << sign_ip_jg_kp << endl;
    return oss.str();
  }
  void SubSymGTOs::Display() const {
    cout << this->str();
  }
  SubSymGTOs Sub_s(Irrep irrep, Vector3cd xyz, VectorXcd zs) {

    MatrixXcd cs = MatrixXcd::Ones(1,1);
    Reduction rds(irrep, cs);
    MatrixXi symmat = MatrixXi::Ones(1, 1);
    MatrixXi signmat= MatrixXi::Ones(1, 1);

    SubSymGTOs sub;
    sub.AddXyz(xyz);
    sub.AddNs(Vector3i(0, 0, 0));
    sub.AddZeta(zs);
    sub.AddRds(rds);

    return sub;
  }
  SubSymGTOs Sub_pz(Irrep irrep, Vector3cd xyz, VectorXcd zs) {

    MatrixXcd cs = MatrixXcd::Ones(1,1);
    Reduction rds(irrep, cs);

    SubSymGTOs sub;
    sub.AddXyz(xyz);
    sub.AddNs(Vector3i(0, 0, 1));
    sub.AddRds(rds);
    sub.AddZeta(zs);

    return sub;    
  }
  SubSymGTOs Sub_TwoSGTO(pSymmetryGroup sym, Irrep irrep,
			 Vector3cd xyz, VectorXcd zs) {

    sym->CheckIrrep(irrep);
    SubSymGTOs sub;
    
    sub.SetSym(sym);

    if(sym->name() == "Cs") {
      dcomplex x = xyz[0];
      dcomplex y = xyz[1];
      dcomplex z = xyz[2];
      sub.AddXyz(xyz);
      sub.AddXyz(Vector3cd(x, y, -z));
      sub.AddNs(Vector3i(0, 0, 0));
      sub.AddZeta(zs);

      MatrixXcd cs(2, 1);
      if(irrep == sym->GetIrrep("A'")) {
	cs << 1.0, 1.0;
      } else if (irrep == sym->GetIrrep("A''")){
	cs << 1.0, -1.0;
      }
      Reduction rds(irrep, cs);
      sub.AddRds(rds);
      
      sub.SetUp();
      return sub;

    } else {
      string msg; SUB_LOCATION(msg);
      msg += "only Cs symmetry is implemented now";
      throw runtime_error(msg);
    }
  }
  SubSymGTOs Sub_mono(Irrep irrep,
		      Vector3cd xyz, Vector3i ns, VectorXcd zs) {

    SubSymGTOs sub;

    sub.AddXyz(xyz);
    sub.AddNs(ns);
    sub.AddZeta(zs);
    sub.AddRds(Reduction(irrep, MatrixXcd::Ones(1, 1)));

    return sub;
  }
  
  // ==== SymGTOs ====
  // ---- Constructors ----
  _SymGTOs::_SymGTOs():
    setupq(false)  {
    xyzq_iat = MatrixXcd::Zero(4, 0);
  }

  // ---- Accessors ----
  int _SymGTOs::size_atom() const {return xyzq_iat.cols(); }
  int _SymGTOs::size_basis() const {
    int cumsum(0);
    for(cSubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      cumsum += isub->size_zeta() * isub->rds.size();
    }
    return cumsum;
  }
  int _SymGTOs::size_basis_isym(Irrep isym) const {

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
  string _SymGTOs::str() const {
    ostringstream oss;
    oss << "==== SymGTOs ====" << endl;
    oss << "Set Up?" << (setupq ? "Yes" : "No") << endl;
    oss << sym_group->str();
    for(cSubIt it = subs.begin(); it != subs.end(); ++it) 
      oss << it->str();
    oss << "xyzq:" << endl;
    oss << xyzq_iat << endl;
    return oss.str();
  }

  // ---- Add ----
  void _SymGTOs::SetSym(pSymmetryGroup sym) {
    sym_group = sym;
  }
  void _SymGTOs::SetAtoms(MatrixXcd _xyzq_iat) {

    if(_xyzq_iat.rows() != 4) {
      string msg; SUB_LOCATION(msg);
      msg += "xyzq_iat.rows() must be 4 ";
      throw runtime_error(msg);
    }

    if(_xyzq_iat.cols() < 1) {
      string msg; SUB_LOCATION(msg);
      msg += "xyzq_iat.cols() must be gerater or equal to 1 ";
      throw runtime_error(msg);
    }    

    xyzq_iat = _xyzq_iat;

  }
  void _SymGTOs::AddAtom(Eigen::Vector3cd _xyz, dcomplex q) {

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
  void _SymGTOs::AddSub(SubSymGTOs sub) {

    subs.push_back(sub);

  }
  
  // ---- Other ----
  SymGTOs _SymGTOs::Clone() const {

    SymGTOs gtos = SymGTOs(new _SymGTOs());
    gtos->SetSym(this->sym_group);
    gtos->xyzq_iat = this->xyzq_iat;
    gtos->subs = this->subs;
    gtos->SetUp();
    return gtos;

  }
  SymGTOs _SymGTOs::Conj() const {

    SymGTOs gtos = this->Clone();
    for(SubIt isub = gtos->subs.begin(); isub != gtos->subs.end(); ++isub) {
      isub->zeta_iz = isub->zeta_iz.conjugate();
      for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) {
	irds->coef_iat_ipn = irds->coef_iat_ipn.conjugate();
	irds->coef_iz = irds->coef_iz.conjugate();
      }
    }
    gtos->SetUp();
    return gtos;
  }

  // ---- SetUp ----
  void _SymGTOs::SetUp() {

    if(sym_group.get() == NULL) {
      string msg; SUB_LOCATION(msg);
      msg += ": symmetry_group is not set";
      throw runtime_error(msg);
    }
    
    for(SubIt it = subs.begin(); it != subs.end(); ++it) {

      // -- setup sub --
      it->sym_group = this->sym_group;
      if(it->setupq == false)
	it->SetUp();

      // -- check each sub --
      if(it->ip_jg_kp.rows() != sym_group->order()) {
	string msg; SUB_LOCATION(msg); 
	msg += ": size of symmetry transformation matrix and num_class of symmetry group is not matched.";
	cout << "ip_jg_kp:" << endl;
	cout << it->ip_jg_kp << endl;
	throw runtime_error(msg);
      }
      
      // -- init offset --
      for(RdsIt irds = it->begin_rds(); irds != it->end_rds(); ++irds) {
	sym_group->CheckIrrep(irds->irrep);
	irds->offset = 0;
      }
    }

    this->SetOffset();
    this->Normalize();
    setupq = true;
  }
  void _SymGTOs::SetOffset() {
    
    vector<int> num_irrep(10, 0);
    
    for(SubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(RdsIt irds = isub->rds.begin(), end = isub->rds.end();
	  irds != end; ++irds) {
	irds->offset = num_irrep[irds->irrep];
	num_irrep[irds->irrep] += isub->size_zeta();
      }
    }

  }
  void _SymGTOs::Normalize() {

    /*
  A3dc dxmap(100), dymap(100), dzmap(100);
    for(SubIt isub = subs.begin(), end = subs.end(); isub != end; ++isub) {
      for(int iz = 0; iz < isub->size_zeta(); iz++) {

	dcomplex zetai = isub->zeta_iz[iz];
	dcomplex zetaP = zetai + zetai;

	int mi = isub->maxn;
	calc_d_coef(mi,mi,0, zetaP, 0.0, 0.0, 0.0, dxmap);
	calc_d_coef(mi,mi,0, zetaP, 0.0, 0.0, 0.0, dymap);
	calc_d_coef(mi,mi,0, zetaP, 0.0, 0.0, 0.0, dzmap);

	int nipn(isub->size_pn());
	for(int ipn = 0; ipn < nipn; ipn++) {
	  int nxi, nyi, nzi;
	  nxi = isub->nx(ipn);
	  nyi = isub->ny(ipn);
	  nzi = isub->nz(ipn);
	  dcomplex dx00, dy00, dz00;
	  dx00 = dxmap(nxi, nxi ,0);
	  dy00 = dymap(nyi, nyi ,0);
	  dz00 = dzmap(nzi, nzi ,0);
	  dcomplex s_ele = dx00 * dy00 * dz00;
	  dcomplex ce = pow(M_PI/zetaP, 1.5);	  
	  for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) {
	    irds->coef_iz(iz) = 1.0/sqrt(ce*s_ele);
	  }
	}
      }
    }
    */

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
	      xi = isub->x(iat); xj = isub->x(jat);
	      yi = isub->y(iat); yj = isub->y(jat);
	      zi = isub->z(iat); zj = isub->z(jat);
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
		  nxi = isub->nx(ipn); nxj = isub->nx(jpn);
		  nyi = isub->ny(ipn); nyj = isub->ny(jpn);
		  nzi = isub->nz(ipn); nzj = isub->nz(jpn);
		  dcomplex dx00, dy00, dz00;
		  dx00 = dxmap(nxi, nxj ,0);
		  dy00 = dymap(nyi, nyj ,0);
		  dz00 = dzmap(nzi, nzj ,0);
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
	    oss << "isub:"  << distance(subs.begin(), isub) << endl;
	    oss << "iz  : " << iz << endl;
	    oss << "zeta: " << zetai << endl;
	    oss << "irds:" << distance(isub->rds.begin(), irds) << endl;
	    oss << "nipn: " << nipn << endl;
	    oss << "niat: " << niat << endl;
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
  int _SymGTOs::max_n() const {
    int max_n(0);
    for(cSubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      for(int ipn = 0; ipn < isub->size_pn(); ipn++) {
	int nx(isub->nx(ipn));
	int ny(isub->ny(ipn));
	int nz(isub->nz(ipn));
	if(max_n < nx)
	  max_n = nx;
	if(max_n < ny)
	  max_n = ny;
	if(max_n < nz)
	  max_n = nz;
      }
    }    
    return max_n;
  }
  
  // ---- CalcMat ---- 
  void _SymGTOs::loop() {

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
  int _SymGTOs::max_num_prim() const {
    
    int max_num_prim(0);
    for(cSubIt isub = subs.begin(); isub != subs.end(); ++isub) {
      max_num_prim = std::max(max_num_prim, isub->size_prim());
    }
    return max_num_prim;

  }

  // ---- AtR ----
  bool IsCenter(SubIt isub, int iat, double eps) {
    dcomplex x  = isub->x(iat);
    dcomplex y  = isub->y(iat);
    dcomplex z  = isub->z(iat);
    dcomplex a2 = x*x+y*y+z*z;
    dcomplex a = sqrt(a2);
    return (abs(a) < eps);
  }
  void AtR_Ylm_cen(SubIt isub, int iz,
		   int ipn, dcomplex r, int L, int M, dcomplex* v, dcomplex* dv) {

    int nx = isub->nx(ipn);
    int ny = isub->ny(ipn);
    int nz = isub->nz(ipn);
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
      if((nx == 0 && ny == 0 && nz == 1 && M == 0) ||
	 (nx == 0 && ny == 1 && nz == 0 && M ==-1) ||
	 (nx == 1 && ny == 0 && nz == 0 && M ==+1)) {
	dcomplex c(sqrt(4.0*M_PI/3.0));
	*v = c*r*r*exp(-zeta*r*r);
	*dv= c*(2.0*r - 2.0*zeta*r*r*r )*exp(-zeta*r*r);
      }	else {
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

    int nn = (isub->nx(ipn) +
	      isub->ny(ipn) +
	      isub->nz(ipn));
    dcomplex zeta = isub->zeta_iz(iz);
    dcomplex x  = isub->x(iat);
    dcomplex y  = isub->y(iat);
    dcomplex z  = isub->z(iat);
    dcomplex xxyy = x*x+y*y;
    dcomplex a2 = xxyy+z*z;
    dcomplex a = sqrt(a2);

    if(abs(a) < 0.000001) {
      string msg; SUB_LOCATION(msg);
      msg += "(x,y,z) is not centered on origin."; 
      throw runtime_error(msg);
    }

    dcomplex expz = exp(-zeta*(r*r+a2));
    dcomplex theta = cacos(z / a);
    dcomplex phi   = (abs(xxyy) < 0.00001) ? 0.0 : cacos(x / sqrt(xxyy));
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
  void _SymGTOs::AtR_Ylm(int L, int M, int irrep, const VectorXcd& cs_ibasis,
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
      sym_group->CheckIrrep(irrep);
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
  void _SymGTOs::CorrectSign(int L, int M, int irrep, Eigen::VectorXcd& cs) {

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
      sym_group->CheckIrrep(irrep);
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

  SymGTOs CreateSymGTOs() {
    SymGTOs ptr(new _SymGTOs);
    return ptr;
  }  
  
}

