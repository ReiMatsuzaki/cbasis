#include <fstream>
#include <iostream>
#include <stdexcept>
#include <boost/format.hpp>
#include "../external/picojson/picojson.h"
#include "../utils/eigen_plus.hpp"
#include "../utils/timestamp.hpp"
#include "../src_cpp/symmolint.hpp"
#include "../src_cpp/one_int.hpp"
#include "../src_cpp/two_int.hpp"
#include "../src_cpp/mo.hpp"
#include "../src_cpp/read_json.hpp"

using namespace std;
using namespace cbasis;

// -- Input --
string comment, in_json, out_eigvecs, out_eigvals;
SymmetryGroup sym;
Molecule mole;
SymGTOs gtos;
int num_ele;
int max_iter;
double tol;
ERIMethod eri_method;

// -- Results --
MO mo;
bool conv;

void Parse() {
  PrintTimeStamp("Parse", NULL);
  
  ifstream f;
  picojson::value json;

  f.open(in_json.c_str(), ios::in);
  if(f.fail()) {
    cerr << "opening input json file failed." << endl;
    cerr << "filename : " << in_json << endl;
    exit(1);
  }

  try {
    f >> json;
    if(not json.is<picojson::object>()) {
      throw runtime_error("json is not object");
    }
    picojson::object& obj = json.get<picojson::object>();
    comment = ReadJson<string>(obj, "comment");
    sym = ReadJson<SymmetryGroup>(obj, "sym");
    mole = NewMolecule(sym);
    ReadJson_Molecule(obj, "molecule", mole);
    gtos = NewSymGTOs(mole);
    ReadJson_SymGTOs_Subs(obj, "basis", gtos);
    num_ele   = ReadJson<int>(obj, "num_ele");
    max_iter  = ReadJson<int>(obj, "max_iter");
    tol       = ReadJson<double>(obj, "tol");
    eri_method = ReadJson<ERIMethod>(obj, "eri_method");
    //    in_eigvecs = ReadJson<string>(obj, "in_eigvecs");
    out_eigvecs = ReadJson<string>(obj, "out_eigvecs");
    out_eigvals = ReadJson<string>(obj, "out_eigvals");
    
  } catch(exception& e) {
    cerr << "error on parse json" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  f.close();

  try {
    gtos->SetUp();
  } catch(exception& e) {
    cerr << "error on setup of SymGTOs" << endl;
    cerr << e.what() << endl;
    cerr << "gtos:" << endl;
    cerr << gtos->str() << endl;
    exit(1);
  }

}
void PrintIn() {
  PrintTimeStamp("PrintIn", NULL);
  cout << "comment: " << comment << endl;
  cout << "in_json: " << in_json << endl;
  cout << "out_eigvecs: " << out_eigvecs << endl;
  cout << "out_eigvals: " << out_eigvals << endl;
  cout << "ERIMethod_use_symmetry: " << eri_method.symmetry << endl;
  cout << "ERIMethod_use_memo: " << eri_method.coef_R_memo << endl;
  cout << "ERIMethod_use_perm: " << eri_method.perm << endl;
  cout << "symmetry: " << sym->name() << endl;
  cout << "molecule: " << endl << mole->show() << endl;
  cout << "num_ele: " << num_ele << endl;
  cout << "gtos: " << endl << gtos->show() << endl;
}
void CalcMat() {
  
  //PrintTimeStamp("Mat", NULL);
    /*
    try {
      E.Read(in_eigvals);
      C.Read(in_eigvecs);
    } catch(exception& e) {
      cerr << "error on reading initial eigvals or eigvecs" << endl;
      cerr << e.what() << endl;
      exit(1);
    }
    vector<int> occ_num = CalcOccNum(E, sym->order(), num_ele/2);    
    CalcSTVMat(gtos, gtos, &S, &T, &V);

    ERIMethod eri_method;
    eri = CalcERI_Complex(gtos, eri_method);
    */
    
}
void CalcMain() {
  
  PrintTimeStamp("Calc", NULL);
  BMatSet mat_set;

  try {
    mat_set = CalcMat_Complex(gtos, true);
  } catch(exception& e) {
    cerr << "error on calculating mat" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  B2EInt  eri;
  try {
    eri = CalcERI_Complex(gtos, eri_method);
  } catch(exception& e) {
    cerr << "error on calculating eri" << endl;
    cerr << e.what() << endl;
    exit(1);
  } 

  try {
    mo = CalcRHF(sym, mat_set, eri, num_ele, max_iter, tol, &conv, 1);
  } catch(exception& e) {
    cerr << "error on RHF" << endl;
    cerr << e.what() << endl;
  }
  //  for(Irrep irrep = 0; irrep < sym->order(); irrep++) {
  //  }
}
void PrintOut() {

  PrintTimeStamp("PrintOut", NULL);

  try {
    mo->C.Write(out_eigvecs);
  } catch(exception& e) {
    cerr << "error on writing coef matrix" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  
  try {
    mo->eigs.Write(out_eigvals);
  } catch(exception& e) {
    cerr << "error on writing eig vector" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  
  cout << "convergence: " << (conv ? "Yes" : "No") << endl;
  cout << "orbital_energy: " << mo->eigs(0)(0) << endl;
  cout << "ele_energy: " << mo->energy << endl;
  cout << "mole_energy: " << mo->energy + mole->NucEnergy() << endl;

  /*
  picojson::object out;
  out["comment"] = picojson::value(comment);
  out["in_json"] = picojson::value(in_json);
  picojson::array ary;
  ary.push_back(picojson::value(mo->energy.real()));
  ary.push_back(picojson::value(mo->energy.imag()));
  out["E0"] = picojson::value(ary);
  ofstream ofs(out_json.c_str())
  */
}
int main (int argc, char *argv[]) {
  cout << ">>>> rhf >>>>" << endl;
  if(argc != 2) {
    cerr << "need input json file" << endl;
    exit(1);
  }
  in_json = argv[1];
  Parse();
  PrintIn();
  CalcMat();
  CalcMain();
  PrintOut();
  cout << "<<<< rhf <<<<" << endl;
}
