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
  PrimGTO::PrimGTO(int _nx, int _ny, int _nz,
		   dcomplex _ax, dcomplex _ay, dcomplex _az):
    nx(_nx), ny(_ny), nz(_nz), x(_ax), y(_ay), z(_az) {}
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
  SymOpMult::SymOpMult(ISymOp *_sym_op, int _mult) {
    sym_op = _sym_op;
    mult = _mult;
    if(mult < 2) {
      string msg; SUB_LOCATION(msg);
      msg += ": mult must be bigger than 1";
    }
  }
  void SymOpMult::getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *prim) const {
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
  string SymOpMult::str() const {
    ostringstream oss;
    oss << this->mult << sym_op->str();
    return oss.str();
  }

  // ---- Id ----  
  void SymOpId::getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *prim) const {
    *b = a;
    *sig = 1;    
    *prim = true;
  }
  string SymOpId::str() const {
    return "E";
  }

  // ---- Cyclic ---- 
  SymOpCyclic::SymOpCyclic(Axis _axis, int _n): axis(_axis), n(_n) {
    if(n < 2) {
      string msg; SUB_LOCATION(msg);
      msg += ": n must be positive integer greater than 1";
      throw runtime_error(msg);
    }
  }
  void SymOpCyclic::getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *prim) const {

    if(axis == AxisX && n == 2) {
      *b = a;
      b->y = -a.y;
      b->z = -a.z;
      *sig = pow(-1, a.ny + a.nz);
      *prim = true;
    } else if(axis == AxisY && n == 2) {
      *b = a;
      b->x = -a.x;
      b->z = -a.z;
      *sig = pow(-1, a.nx + a.nz);
      *prim = true;
    } else if(axis == AxisZ && n == 2) {
      *b = a;
      b->x = -a.x;
      b->y = -a.y;
      *sig = pow(-1, a.nx + a.ny);
      *prim = true;
    } else if(axis == AxisX && n == 4) {
      *b = a;
      b->y = -a.z;
      b->z = +a.y;
      b->ny = a.nz;
      b->nz = a.ny;
      *sig = pow(-1, a.nz);
      *prim = true;      
    } else if(axis == AxisY && n == 4) {
      *b = a;
      b->z = -a.x;
      b->x = +a.z;
      b->nz = a.nx;
      b->nx = a.nz;
      *sig = pow(-1, a.nx);
      *prim = true;
    } else if(axis == AxisZ && n == 4) {
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
  string SymOpCyclic::str() const {
    ostringstream oss;
    oss << "C";
    if(this->axis == AxisX)
      oss << "x";
    if(this->axis == AxisY)
      oss << "y";
    if(this->axis == AxisZ)
      oss << "z";
    oss << "(" << this->n << ")";
    return oss.str();    
  }

}
