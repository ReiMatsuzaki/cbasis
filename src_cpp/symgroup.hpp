#ifndef SYMGROUP_H
#define SYMGROUP_H

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "typedef.hpp"
#include "mult_array.hpp"
#include <Eigen/Core>

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
  
  // ==== Symmetry operation class ====
  typedef std::vector<SymOp> SymOpClass;
  SymOpClass ClassMono(SymOp o);

  // ==== Symmetry Group ====
  typedef int Irrep;
  //typedef Eigen::Matrix<std::complex<int>, Eigen::Dynamic, Eigen::Dynamic> MatrixXci;
  //typedef Eigen::Matrix<std::complex<int>, 1, Eigen::Dynamic> VectorXci;

  class SymmetryGroup;
  typedef boost::shared_ptr<SymmetryGroup> pSymmetryGroup;
  class SymmetryGroup {
  public:
    std::string name_;
    int id_num_; 
    std::vector<SymOp> sym_op_;    
    std::vector<SymOpClass> sym_op_class_;    
    std::vector<std::string> irrep_name_;
    Eigen::MatrixXi character_table_;
    MultArray<bool, 3> prod_table_;
    Irrep irrep_s;
    Irrep irrep_x;
    Irrep irrep_y;
    Irrep irrep_z;    
  public:
    // ---- typedef ----
    typedef std::vector<SymOp>::const_iterator ItSymOp;

    // ---- Constructors ----
    SymmetryGroup(int num_class, std::string name, int id_num);

    // ---- Accessor ----
    int order()  const { return sym_op_.size(); }
    int num_class()  const { return sym_op_class_.size(); }
    int id_num() const { return id_num_;}
    std::string name() const { return name_; }
    std::string str() const;
    void Display() const;
    Irrep GetIrrep(std::string name) const;
    std::string GetIrrepName(Irrep irrep) const;
    
    // ---- Calculation ----
    void setProdTable();
    void setSymOp();
    void CheckIrrep(Irrep a);
    bool IsSame(pSymmetryGroup o);
    bool Non0_Scalar(Irrep a, Irrep b);
    bool Non0_Z(Irrep a, Irrep b);
    bool Non0_4(Irrep a, Irrep b, Irrep c, Irrep d);
    /**
       Build matrix representing symmetry operation for gtos.
       a(I, i) = j  : Ith symmetry operation for ith GTO  is equal to jth GTO
     */
    void CalcSymMatrix(const std::vector<PrimGTO>& gtos,
		       Eigen::MatrixXi& a, Eigen::MatrixXi& sig);

    // ---- Specific symmetry group ----
    
    static pSymmetryGroup C1();
    static pSymmetryGroup Cs();
    static pSymmetryGroup C2h();
    static pSymmetryGroup C2v();
    static pSymmetryGroup D2h();
    static pSymmetryGroup C4();
    //    static SymmetryGroup C4v();
  };

}
#endif
