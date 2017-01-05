#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include <gtest/gtest.h>
#include "../src_cpp/bmatset.hpp"
#include "read_aoint_wrapper.hpp"

using namespace std;
using boost::format;
using namespace cbasis;
using namespace Eigen;


TEST(TestFirst, First) {
  /*
  int aoint_ifile;
  char aoint_filename[20] = "out/AOINTS";
  //long aoint_filename_length = 20;
  long aoint_filename_length = strlen(aoint_filename)+1;
  bool succ;

  aoint_filename[aoint_filename_length-1] = ' ';
  open_file_binary_read_(&aoint_ifile, &succ, aoint_filename, aoint_filename_length);
  
  char blabel[80];
  double repfunc;
  int    nst;
  short int nd[8];
  
  read_header_(&aoint_ifile, blabel, &repfunc, &nst, nd);
  cout << "blabel: " << blabel << endl;
  cout << "repfunc: " << repfunc << endl;
  cout << "nst: " << nst << endl;
  cout << "nd: " << nd[0] << nd[1] << endl;
  close_file_(&aoint_ifile);
  EXPECT_TRUE(succ);
  EXPECT_EQ(2, 1+1);
  */
}
TEST(TestFirst, third) {

  AoIntsHeader ao;
  BMat s, t, v;
  char filename[20] = "out/AOINTS";
  ReadAOINTS(filename, &ao, &s, &t, &v);
  cout << s(0, 0)(0, 0) << endl;
  cout << s(0, 0)(0, 1) << endl;
  cout << s(0, 0)(1, 0) << endl;
}
TEST(TestFirst, Second) {
  int aoint_ifile;
  char aoint_filename[20] = "out/AOINTS";
  //long aoint_filename_length = 20;
  long aoint_filename_length = strlen(aoint_filename)+1;
  bool succ;

  aoint_filename[aoint_filename_length-1] = ' ';
  open_file_binary_read_(&aoint_ifile, &succ, aoint_filename, aoint_filename_length);

  AoIntsHeader ao;
  aoints_read_header_(&aoint_ifile, &ao);

  cout << ao.blabel << endl;
  cout << ao.repfunc << endl;
  cout << ao.nst << endl;
  cout << ao.nd[0] << endl;
  cout << ao.ityp[0] << endl;
  cout << "ns: " << ao.ns << endl;
  cout << "zscale: " << ao.zscale.re << "+" << ao.zscale.im << endl;

  // -- BMat --
  BMat smat, tmat, vmat;
  int num_isym[10];
  aoints_read_mat_structure_(&aoint_ifile, num_isym);
  cout << format("num_isym = ");
  for(int isym = 0; isym < 8; isym++) {
    int ni = num_isym[isym];
    cout << format(" %d") % ni;
    smat(isym, isym) = MatrixXcd::Zero(ni, ni);
    tmat(isym, isym) = MatrixXcd::Zero(ni, ni);
    vmat(isym, isym) = MatrixXcd::Zero(ni, ni);
  }
  cout << endl;
  close_file_(&aoint_ifile);

  // re read
  open_file_binary_read_(&aoint_ifile, &succ, aoint_filename, aoint_filename_length);
  aoints_read_header_(&aoint_ifile, &ao);
  int is[1080];
  int js[1080];
  int isyms[1080];
  for_complex vs[1080];
  int num;
  bool is_end_block = false;
  for(bool is_end_block = false; not is_end_block; ) {
    aoints_read_mat_value_block_(&aoint_ifile, &num, &is_end_block, vs, is, js, isyms);
    cout << "num : " << num << endl;
    cout << "is end : " << (is_end_block ? "yes" : "no") << endl;
    for(int i = 0; i < num; i++) {
      smat(isyms[i]-1, isyms[i]-1)(is[i]-1, js[i]-1) = dcomplex(vs[i].re, vs[i].im);
    }
  }
  for(bool is_end_block = false; not is_end_block; ) {
    aoints_read_mat_value_block_(&aoint_ifile, &num, &is_end_block, vs, is, js, isyms);
    for(int i = 0; i < num; i++) {
      tmat(isyms[i]-1, isyms[i]-1)(is[i]-1, js[i]-1) = dcomplex(vs[i].re, vs[i].im);
    }
  }
  for(bool is_end_block = false; not is_end_block; ) {
    aoints_read_mat_value_block_(&aoint_ifile, &num, &is_end_block, vs, is, js, isyms);
    for(int i = 0; i < num; i++) {
      vmat(isyms[i]-1, isyms[i]-1)(is[i]-1, js[i]-1) = dcomplex(vs[i].re, vs[i].im);
    }
  }    
  close_file_(&aoint_ifile);
}

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
