#ifndef SYMGROUP_H
#define SYMGROUP_H

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "typedef.hpp"
#include "mult_array.hpp"

namespace l2func {
  // ==== primitive GTO ====
  class PrimGTO {
  public:
    int nx, ny, nz;
    dcomplex x, y, z;
    PrimGTO();
    PrimGTO(int _nx, int _ny, int _nz, dcomplex _ax, dcomplex _ay, dcomplex _az);
    std::string str() const;
  };
  std::ostream& operator<< (std::ostream& oss, const PrimGTO& o);
  bool IsNear(const PrimGTO& a, const PrimGTO& b);

  enum Coord {
    CoordX,
    CoordY,
    CoordZ
  };

  // ==== Symmetry operation ====
  // ---- Interface ----
  class ISymOp {
  public:
    /**
       return relation of two PrimGTO a and b.
       1    (Op(a) == +b)
       -1   (Op(a) == -b)
       0    (otherwise)
    */
    virtual ~ISymOp() {}
    int Op(const PrimGTO& a, const PrimGTO& b) const;
    /**
       a :  operated GTO
       b : resultant primitive GTO if exist
       sig : 1 or -1
       is_prim : true => resultant GTO is primitive
     */
    virtual void getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *is_prim) const = 0;
    virtual std::string str() const = 0;
  };
  typedef boost::shared_ptr<ISymOp> SymOp;

  // ---- multiple ----
  class Mult : public ISymOp {
    int mult;
    SymOp sym_op;
  public:
    Mult(SymOp sym_op, int mult);
    ~Mult();
    void getOp(const PrimGTO& a, PrimGTO *b, int *sig, bool *is_prim) const ;
    std::string str() const;
  };
  SymOp mult(SymOp a, int n);
  
  // ---- product ----
  class Prod : public ISymOp {
    SymOp a;
    SymOp b;
  public:
    Prod(SymOp _a, SymOp _b);
    ~Prod();
    void getOp(const PrimGTO& x, PrimGTO *y, int *sig, bool *prim) const;
    std::string str() const;
  };
  SymOp prod(SymOp a, SymOp b);
  
  // ---- Identity ----
  class Id : public ISymOp {
  public:
    ~Id();
    void getOp(const PrimGTO& a, PrimGTO *b, int *sig, bool *is_prim) const ;
    std::string str() const;
  };
  SymOp id();

  // ---- Cyclic ----  
  class Cyclic : public ISymOp {
  public:
    Coord coord;
    int n;
    Cyclic(Coord _axis, int _n);
    void getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *is_prim) const;
    std::string str() const;
  };
  SymOp cyclic(Coord coord, int n);

  // ---- Reflection ----
  class Reflect : public ISymOp {
  public:
    Coord coord;
    Reflect(Coord _axis);
    void getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *is_prim) const;
    std::string str() const;
  };
  SymOp reflect(Coord coord);

  // ---- Inversion Center ----
  class InvCent : public ISymOp {
  public:
    void getOp(const PrimGTO& a, PrimGTO* b, int *sig, bool *is_prim) const;
    std::string str() const;
  };
  SymOp inv();
  
  // ==== Symmetry Group ====
  typedef int Irrep;

  class SymmetryGroup {
  public:
    int order_;
    std::string name_;
    std::vector<SymOp> sym_op_;    
    std::vector<std::string> irrep_name_;
    MultArray<int, 2> character_table_;
    MultArray<bool, 3> prod_table_;
    Irrep irrep_s;
    Irrep irrep_x;
    Irrep irrep_y;
    Irrep irrep_z;    
  public:
    // ---- typedef ----
    typedef std::vector<SymOp>::const_iterator ItSymOp;

    // ---- Constructors ----
    SymmetryGroup(int order, std::string name);

    // ---- Accessor ----
    int order()  const { return order_; }
    std::string name() const { return name_; }
    std::string str() const;
    void Display() const;
    
    // ---- Calculation ----
    void CheckIrrep(Irrep a);
    bool Non0_Scalar(Irrep a, Irrep b);
    bool Non0_Z(Irrep a, Irrep b);
    bool Non0_4(Irrep a, Irrep b, Irrep c, Irrep d);

    // ---- Specific symmetry group ----
    static SymmetryGroup C1();
    static SymmetryGroup Cs();
    static SymmetryGroup C2h();
    static SymmetryGroup D2h();
  };
  SymmetryGroup SymmetryGroup_Cs();
  Irrep Cs_Ap();
  Irrep Cs_App();
  SymmetryGroup SymmetryGroup_C1();
  SymmetryGroup SymmetryGroup_SO3(int maxl);
  Irrep SO3_LM(int L, int M);


}
#endif
