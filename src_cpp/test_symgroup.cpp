#include <iostream>
#include <gtest/gtest.h>
#include "gtest_plus.hpp"
#include "symgroup.hpp"

using namespace std;
using namespace Eigen;
using namespace l2func;

TEST(first, first) {
  EXPECT_EQ(2, 1+1);
}
TEST(SymOp, Id) {

  PrimGTO a(2, 1, 3, 0.1, 0.2, 0.3);
  PrimGTO b(a);
  PrimGTO c(2, 1, 3, 0.0, 0.1, 0.3);
  ISymOp *id = new SymOpId();
  //  bool is_prim;
  //  int sig;
  
  EXPECT_EQ(1, id->Op(a, b));
  EXPECT_EQ(0, id->Op(a, c));

  delete id;
  
}
TEST(SymOp, C2) {

  ISymOp *cx2 = new SymOpCyclic(AxisX, 2);
  ISymOp *cy2 = new SymOpCyclic(AxisY, 2);
  ISymOp *cz2 = new SymOpCyclic(AxisZ, 2);

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

  Axis axis_list[3] = {AxisX, AxisY, AxisZ};

  for(int i = 0; i < 3; i++) {
    Axis axis = axis_list[i];
    ISymOp *C2   = new SymOpCyclic(axis, 2);
    ISymOp *C4   = new SymOpCyclic(axis, 4);
    ISymOp *C4_2 = new SymOpMult(C4, 2);
    ISymOp *C4_4 = new SymOpMult(C4, 4);

    PrimGTO a(1, 2, 3, 0.1, 0.2, 0.3);
    PrimGTO C2_a(a);
    PrimGTO C4_2_a(a);
    int sig1, sig2;
    bool prim1, prim2;

    C2->getOp(  a, &C2_a,   &sig1, &prim1);
    C4_2->getOp(a, &C4_2_a, &sig2, &prim2);
    
    EXPECT_EQ(sig1, sig2);
    EXPECT_TRUE(prim1);
    EXPECT_TRUE(prim2);
    EXPECT_TRUE(IsNear(C2_a, C4_2_a));

    PrimGTO C4_4_a(a);
    C4_4->getOp(a, &C4_4_a, &sig1, &prim1);
    EXPECT_EQ(1, sig1);
    EXPECT_TRUE(prim1);
    EXPECT_TRUE(IsNear(a, C4_4_a));

    delete C2;
    delete C4;
    delete C4_2;    
    delete C4_4;
  }
}  


int main (int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
