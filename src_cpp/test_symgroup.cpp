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
void ExpectOpEq(SymOp a, SymOp b, int num = 5) {
  
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
    SymOp C2 = cyclic(axis, 2);
    SymOp C4 = cyclic(axis, 4);

    ExpectOpEq(id(), mult(C2, 2));
    ExpectOpEq(id(), prod(C2, C2));
    ExpectOpEq(id(), mult(C4, 4));
    ExpectOpEq(C2, mult(C4, 2));
    ExpectOpEq(C2, prod(C4, C4));

  }
}
TEST(SymOp, Sig) {
  Coord axis_list[3] = {CoordX, CoordY, CoordZ};
  for(int i = 0; i < 3; i++) {
    Coord coord = axis_list[i];
    SymOp SIG  = reflect(coord);
    ExpectOpEq(id(), prod(SIG, SIG));
  }
}
TEST(SymOp, SigC2) {

  Coord axis_list[3] = {CoordX, CoordY, CoordZ};
  for(int i = 0; i < 3; i++) {
    Coord axis = axis_list[i];
    SymOp C2   = cyclic(axis, 2);
    SymOp SIG  = reflect(axis);
    ExpectOpEq(inv(), prod(C2, SIG));
    ExpectOpEq(inv(), prod(SIG, C2));
  }  
}
TEST(SymOp, D2h) {
  SymOp C2z = cyclic(CoordZ, 2);
  SymOp C2x = cyclic(CoordX, 2);
  SymOp C2y = cyclic(CoordY, 2);
  SymOp SIGxy = reflect(CoordZ);
  SymOp SIGzx = reflect(CoordY);
  SymOp SIGyz = reflect(CoordX);

  ExpectOpEq(C2y, prod(C2z, C2x));
  ExpectOpEq(C2x, prod(C2z, C2y));
  ExpectOpEq(inv(), prod(C2z, SIGxy));
  ExpectOpEq(C2z, prod(C2x, C2y));
  ExpectOpEq(SIGzx, prod(C2x, SIGxy));
  ExpectOpEq(SIGyz, prod(C2y, SIGxy));

}

int main (int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
