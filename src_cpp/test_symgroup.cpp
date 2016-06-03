#include <iostream>
#include <gtest/gtest.h>
#include "gtest_plus.hpp"
#include "symgroup.hpp"
#include "macros.hpp"

using namespace std;
using namespace Eigen;
using namespace l2func;

dcomplex random_complex() {

  double r0 = (double)rand() / ((double)RAND_MAX+1);
  double i0 = (double)rand() / ((double)RAND_MAX+1);
  return dcomplex(r0, i0);

}
void RandomPrim(PrimGTO* g) {
  int maxn(4);
  g->nx = rand() % maxn;
  g->ny = rand() % maxn;
  g->nz = rand() % maxn;
  g->x = random_complex();
  g->y = random_complex();
  g->z = random_complex();
}
void ExpectOpEq(ISymOp* a, ISymOp* b, int num = 5) {
  
  PrimGTO x0, ay, by;
  for(int i = 0; i < num; i++) {
    RandomPrim(&x0);

    bool prim_a, prim_b;
    int  sig_a,  sig_b;
    a->getOp(x0, &ay, &sig_a, &prim_a);
    b->getOp(x0, &by, &sig_b, &prim_b);
    if(!prim_a || !prim_b) {
      string msg; SUB_LOCATION(msg);
      msg += "result is not primitive";
      throw runtime_error(msg);
    }
    EXPECT_TRUE(IsNear(ay, by))
      << a->str() << endl
      << b->str() << endl
      << ": " << x0.str() << endl << ay.str() << endl << by.str();
    EXPECT_EQ(sig_a, sig_b)
      << a->str() << endl
      << b->str() << endl
      << ": " << x0.str() << endl << ay.str() << endl << by.str();
  }
}

TEST(first, first) {
  EXPECT_EQ(2, 1+1);
}
TEST(SymOp, Id) {

  PrimGTO a(2, 1, 3, 0.1, 0.2, 0.3);
  PrimGTO b(a);
  PrimGTO c(2, 1, 3, 0.0, 0.1, 0.3);
  ISymOp *id = new Id();
  //  bool is_prim;
  //  int sig;
  
  EXPECT_EQ(1, id->Op(a, b));
  EXPECT_EQ(0, id->Op(a, c));

  delete id;
  
}
TEST(SymOp, C2) {

  ISymOp *cx2 = new Cyclic(CoordX, 2);
  ISymOp *cy2 = new Cyclic(CoordY, 2);
  ISymOp *cz2 = new Cyclic(CoordZ, 2);

  PrimGTO s(0, 0, 0, 0.1, 0.2, 0.3);
  PrimGTO px(1, 0, 0, 0.1, 0.2, 0.3);
  PrimGTO px_x(1, 0, 0, 0.1, -0.2, -0.3);
  PrimGTO px_y(1, 0, 0, -0.1, +0.2, -0.3);
  PrimGTO px_z(1, 0, 0, -0.1, -0.2, 0.3);

  EXPECT_EQ(+1, cx2->Op(px, px_x));
  EXPECT_EQ(-1, cy2->Op(px, px_y));
  EXPECT_EQ(-1, cz2->Op(px, px_z));

  PrimGTO py(0, 1, 0, 0.1, 0.2, 0.3);
  PrimGTO pz(0, 0, 1, 0.1, 0.2, 0.3);
  
  delete cx2;
  delete cy2;
  delete cz2;

}
TEST(SymOp, C2_C4) {

  Coord axis_list[3] = {CoordX, CoordY, CoordZ};

  for(int i = 0; i < 3; i++) {
    Coord axis = axis_list[i];
    ISymOp *id   = new Id();
    ISymOp *C2   = new Cyclic(axis, 2);
    ISymOp *C2_2 = new Mult(C2, 2);
    ISymOp *C2_C2= new Prod(C2, C2);

    ISymOp *C4   = new Cyclic(axis, 4);
    ISymOp *C4_C4 = new Prod(C4, C4);
    ISymOp *C4_2 = new Mult(C4, 2);
    ISymOp *C4_4 = new Mult(C4, 4);

    ExpectOpEq(id, C2_2);
    ExpectOpEq(id, C2_C2);
    ExpectOpEq(id, C4_4);
    ExpectOpEq(C2, C4_2);
    ExpectOpEq(C2, C4_C4);

    delete id;
    delete C2;
    delete C2_2;
    delete C2_C2;
    delete C4;
    delete C4_C4;
    delete C4_2;    
    delete C4_4;
  }
}  
TEST(SymOp, Sig) {
  Coord axis_list[3] = {CoordX, CoordY, CoordZ};
  for(int i = 0; i < 3; i++) {
    Coord coord = axis_list[i];
    ISymOp *id   = new Id();
    ISymOp *SIG  = new Reflect(coord);
    ISymOp *SIG2  = new Prod(SIG, SIG);
    ExpectOpEq(id, SIG2, 1);
  }
}
TEST(SymOp, SigC2) {

  Coord axis_list[3] = {CoordX, CoordY, CoordZ};
  for(int i = 0; i < 3; i++) {
    Coord axis = axis_list[i];
    ISymOp *inv   = new InvCent();
    ISymOp *C2   = new Cyclic(axis, 2);
    ISymOp *SIG  = new Reflect(axis);
    ISymOp *C2_SIG = new Prod(C2, SIG);
    ISymOp *SIG_C2 = new Prod(SIG, C2);
    
    ExpectOpEq(inv, C2_SIG);
    ExpectOpEq(inv, SIG_C2);

    delete inv;
    delete C2;
    delete SIG;
    delete C2_SIG;
    delete SIG_C2;
  }  
}
TEST(SymOp, D2h) {
  ISymOp *C2z = new Cyclic(CoordZ, 2);
  ISymOp *C2x = new Cyclic(CoordX, 2);
  ISymOp *C2y = new Cyclic(CoordY, 2);
  ISymOp *SIGxy = new Reflect(CoordZ);
  ISymOp *SIGzx = new Reflect(CoordY);
  ISymOp *SIGyz = new Reflect(CoordX);
  ISymOp *Inv = new InvCent();

  ISymOp *t = new Prod(C2z, C2x);
  ExpectOpEq(C2y, t);
  delete t;
  
  t = new Prod(C2z, C2y);
  ExpectOpEq(C2x, t);
  delete t;

  t = new Prod(C2z, SIGxy);
  ExpectOpEq(Inv, t);
  delete t;

  t = new Prod(C2x, C2y);
  ExpectOpEq(C2z, t);
  delete t;

  t = new Prod(C2x, SIGxy);
  ExpectOpEq(SIGzx, t);
  delete t;

  t = new Prod(C2y, SIGxy);
  ExpectOpEq(SIGyz, t);
  delete t;

  delete C2z; 
  delete C2x;
  delete C2y;
  delete SIGxy;
  delete SIGzx;
  delete SIGyz; 
  delete Inv;
}

int main (int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
