#include <iostream>
#include <gtest/gtest.h>

#include "r1_lc.hpp"
#include "r1basis.hpp"

using namespace std;
using namespace cbasis;

TEST(EigenPlus, Canonical) {
  EXPECT_EQ(2, 1+1);
}
TEST(TestR1LC, TestSTO) {

  LC_STOs s = Create_LC_STOs();
  s->Add(1.1, 2, 1.1);
  s->Add(1.1, 3, 1.3);
  
  cout << s->str() << endl;

}
TEST(TestR1Basis, TestSTO) {

  STOs s = Create_STOs();
  s->AddPrim(2, 1.1);
  s->AddPrim(3, 1.3);
  s->SetUp();
  
  cout << s->str() << endl;

}

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}

