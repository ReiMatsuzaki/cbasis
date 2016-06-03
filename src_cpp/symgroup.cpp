#include <ostream>
#include "symgroup.hpp"
#include "macros.hpp"

using namespace std;

namespace l2func {

  // ==== Utilities ====
  bool IsNear(const dcomplex& a, const dcomplex& b) {
    double eps(0.00000001);
    return abs(a-b) < eps;
  }
  
  // ==== primitive GTO ====
  PrimGTO::PrimGTO():
    nx(0), ny(0), nz(0), x(0), y(0), z(0) {}
  PrimGTO::PrimGTO(int _nx, int _ny, int _nz,
		   dcomplex _ax, dcomplex _ay, dcomplex _az):
    nx(_nx), ny(_ny), nz(_nz), x(_ax), y(_ay), z(_az) {}
  string PrimGTO::str() const {
    ostringstream oss;
    oss << *this;
    return oss.str();
  }
  ostream& operator<< (std::ostream& oss, const PrimGTO& o) {
    oss << "GTO(" << o.nx << o.ny << o.nz << ": " <<
      o.x << ", " << o.y << ", " << o.z << ")";
    return oss;
  }
  bool IsNear(const PrimGTO& a, const PrimGTO& b) {
    return (a.nx == b.nx &&
	    a.ny == b.ny &&
	    a.nz == b.nz &&
	    IsNear(a.x, b.x) &&
	    IsNear(a.y, b.y) &&
	    IsNear(a.z, b.z));
  }
  
  // ==== Symmetry operation ====
  // ---- Interface ----
  int ISymOp::Op(const PrimGTO& a, const PrimGTO& b) const {
    PrimGTO op_a(a); // copy a
    bool is_prim;
    int sig;
    this->getOp(a, &op_a, &sig, &is_prim);

    if(!is_prim)
      return 0;

    if(op_a.nx == b.nx &&
       op_a.ny == b.ny &&
       op_a.nz == b.nz &&
       IsNear(op_a.x, b.x) &&
       IsNear(op_a.y, b.y) &&
       IsNear(op_a.z, b.z))
      return sig;
    else
      return 0;
  }

  // ---- Mult ----
  Mult::Mult(ISymOp *_sym_op, int _mult) {
    sym_op = _sym_op->Clone();
    mult = _mult;
    if(mult < 2) {
      string msg; SUB_LOCATION(msg);
      msg += ": mult must be bigger than 1";
    }
  }
  Mult::~Mult() {
    delete sym_op;
  }
  ISymOp* Mult::Clone() const {
    ISymOp *ptr = new Mult(*this);
    return ptr;
  }
  void Mult::getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *prim) const {
    PrimGTO tmp0(a);
    PrimGTO tmp1(a);
    *sig = 1;
    int tmp_sig = 1;
    bool tmp_prim;
    for(int i = 0; i < this->mult; i++) {
      sym_op->getOp(tmp0, &tmp1, &tmp_sig, &tmp_prim);
      if(!tmp_prim) {
	*prim = false;
	return;
      }
      *sig *= tmp_sig;
      tmp0 = tmp1;
    }
    *b = tmp1;
    *prim = true;
  }
  string Mult::str() const {
    ostringstream oss;
    oss << this->mult << sym_op->str();
    return oss.str();
  }

  // ---- Product ----
  Prod::Prod(ISymOp* _a, ISymOp* _b) {
    a = _a->Clone();
    b = _b->Clone();
  }
  Prod::~Prod() {
    delete a;
    delete b;
  }
  ISymOp* Prod::Clone() const {
    ISymOp* ptr = new Prod(this->a, this->b);
    return ptr;
  }
  void Prod::getOp(const PrimGTO& x, PrimGTO *y, int *sig, bool *prim) const {
    PrimGTO bx(x);
    int sig_bx;
    bool prim_bx;
    b->getOp(x, &bx, &sig_bx, &prim_bx);
    if(!prim_bx) {
      *prim = false;
      return;
    }

    int sig_a;
    bool prim_a;
    a->getOp(bx, y, &sig_a, &prim_a);

    if(!prim_a) {
      *prim = false;
      return;
    }

    *sig = sig_a * sig_bx;
    *prim = true;
  }
  std::string Prod::str() const {
    ostringstream oss;
    oss << "prod(" << a->str() << ", " << b->str() << ")";
    return oss.str();
  }

  // ---- Id ----  
  ISymOp* Id::Clone() const {
    ISymOp *ptr = new Id(*this);
    return ptr;
  }
  void Id::getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *prim) const {
    *b = a;
    *sig = 1;    
    *prim = true;
  }
  string Id::str() const {
    return "E";
  }

  // ---- Cyclic ---- 
  Cyclic::Cyclic(Coord _coord, int _n): coord(_coord), n(_n) {
    if(n < 2) {
      string msg; SUB_LOCATION(msg);
      msg += ": n must be positive integer greater than 1";
      throw runtime_error(msg);
    }
  }
  ISymOp* Cyclic::Clone() const {
    ISymOp *ptr = new Cyclic(*this);
    return ptr;
  }
  void Cyclic::getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *prim) const {

    if(coord == CoordX && n == 2) {
      *b = a;
      b->y = -a.y;
      b->z = -a.z;
      *sig = pow(-1, a.ny + a.nz);
      *prim = true;
    } else if(coord == CoordY && n == 2) {
      *b = a;
      b->x = -a.x;
      b->z = -a.z;
      *sig = pow(-1, a.nx + a.nz);
      *prim = true;
    } else if(coord == CoordZ && n == 2) {
      *b = a;
      b->x = -a.x;
      b->y = -a.y;
      *sig = pow(-1, a.nx + a.ny);
      *prim = true;
    } else if(coord == CoordX && n == 4) {
      *b = a;
      b->y = -a.z;
      b->z = +a.y;
      b->ny = a.nz;
      b->nz = a.ny;
      *sig = pow(-1, a.nz);
      *prim = true;      
    } else if(coord == CoordY && n == 4) {
      *b = a;
      b->z = -a.x;
      b->x = +a.z;
      b->nz = a.nx;
      b->nx = a.nz;
      *sig = pow(-1, a.nx);
      *prim = true;
    } else if(coord == CoordZ && n == 4) {
      *b = a;
      b->x = -a.y;
      b->y = +a.x;
      b->nx = a.ny;
      b->ny = a.nx;
      *sig = pow(-1, a.ny);
      *prim = true;
    } else {
      string msg; SUB_LOCATION(msg);
      msg += ": not supported yet";
      throw runtime_error(msg);
    }
  }
  string Cyclic::str() const {
    ostringstream oss;
    oss << "C";
    if(this->coord == CoordX)
      oss << "x";
    if(this->coord == CoordY)
      oss << "y";
    if(this->coord == CoordZ)
      oss << "z";
    oss << "(" << this->n << ")";
    return oss.str();    
  }

  // ---- Reflection ----
  Reflect::Reflect(Coord _coord): coord(_coord) {}
  ISymOp* Reflect::Clone() const {
    ISymOp* ptr = new Reflect(*this);
    return ptr;
  }
  void Reflect::getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *prim) const {

    *b = a;
    *prim = true;

    if(coord == CoordX) {
      b->x = -a.x;
      *sig = pow(-1, a.nx);
    } 
    else if(coord == CoordY) {
      b->y = -a.y;
      *sig = pow(-1, a.ny);
    } else if(coord == CoordZ) {
      b->z = -a.z;
      *sig = pow(-1, a.nz);
    }
    
  }
  string Reflect::str() const {
    ostringstream oss;
    oss << "SIG" ;
    if(this->coord == CoordX)
      oss << "x";
    if(this->coord == CoordY)
      oss << "y";
    if(this->coord == CoordZ)
      oss << "z";
    return oss.str();
  }

  // ---- Inversion Center ----
  ISymOp* InvCent::Clone() const {
    ISymOp *ptr = new InvCent(*this);
    return ptr;
  }
  void InvCent::getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *prim) const {
    *b = a;
    b->x = -a.x;
    b->y = -a.y;
    b->z = -a.z;
    *sig = pow(-1, a.nx + a.ny + a.nz);
    *prim = true;    
  }
  string InvCent::str() const {
    return "i";
  }
}
