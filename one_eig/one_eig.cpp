#include <fstream>
#include "../utils/eigen_plus.hpp"
#include "../utils/timestamp.hpp"
#include "../src_cpp/symmolint.hpp"
#include "../src_cpp/one_int.hpp"
#include "../src_cpp/read_json.hpp"
#include "../external/picojson/picojson.h"

using namespace std;
using namespace Eigen;
using namespace cbasis;
using namespace picojson;

void PrintHelp() {
  cout << "one_eig" << endl;
}
void CheckInput(picojson::value& json) {

  if(not json.is<picojson::object>()) {
    cerr << "json is not object" << endl;
    exit(1);
  }

  // ==== Key check ====
  picojson::object& obj = json.get<picojson::object>();

  if(obj.find("sym") == obj.end()) {
    cerr << "sym not found" << endl; exit(1);
  }
  if(obj.find("basis") == obj.end()) {
    cerr << "basis not found" << endl; exit(1);
  }
  if(obj.find("molecule") == obj.end()) {
    cerr << "molecule not found" << endl; exit(1);
  }
  if(obj.find("out") == obj.end()) {
    cerr << "out not found" << endl; exit(1);
  }
  if(not obj["out"].is<picojson::object>()) {
    cerr << "out must be object" << endl; exit(1);
  }
  picojson::object& out_obj = obj["out"].get<picojson::object>();
  if(out_obj.find("json") == out_obj.end()) {
    cerr << "out must contain json" << endl; exit(1);
  }
  if(out_obj.find("prefix") == out_obj.end()) {
    cerr << "out must contain prefix" << endl; exit(1);
  } 
}
int main (int argc, char *argv[]) {

  if(argc == 1) {
    PrintHelp();
    exit(1);
  }
  
  ifstream f;
  f.open(argv[1], ios::in);
  picojson::value json;  
  
  // ==== parse json ====
  PrintTimeStamp("Parse", NULL);
  pSymmetryGroup sym;
  SymGTOs gtos(new _SymGTOs());
  Molecule mole;
  string json_out_name, eigvecs_out_name, eigvals_out_name;

  try {
    f >> json;
    CheckInput(json);
    picojson::object& obj = json.get<picojson::object>();

    // -- basis --
    sym = ReadJson<pSymmetryGroup>(obj["sym"]);
    mole = ReadJson<Molecule>(obj["molecule"]);
    ReadJson_SymGTOs_Subs(obj["basis"], gtos);
    gtos->SetSym(sym);
    gtos->SetMolecule(mole);
    gtos->SetUp();

    // -- out --
    object& out_obj = obj["out"].get<object>();
    json_out_name = out_obj["json"].get<string>();
    string prefix = out_obj["prefix"].get<string>();
    eigvecs_out_name = prefix + "eigvecs.bin";
    eigvals_out_name = prefix + "eigvals.bin";    
    
  } catch(exception& e) {
    string msg = "error on parsing json\n";
    msg += e.what();
    cerr << msg << endl;
    exit(1);
  }
  cout << gtos->str() << endl;

  // ==== calculation ====
  PrintTimeStamp("Calc", NULL);
  BMat S, T, V, C;
  BVec E;
  gtos->InitBMat(0, &S); gtos->InitBMat(0, &T); gtos->InitBMat(0, &V);
  gtos->InitBMat(0, &C); gtos->InitBVec(&E);
  
  CalcSTVMat(gtos, gtos, &S, &T, &V);  

  for(int irrep = 0; irrep < sym->order(); irrep++) {
    MatrixXcd h = -0.5 * T(irrep, irrep) + V(irrep, irrep);
    MatrixXcd s = S(irrep, irrep);
    SymGenComplexEigenSolver solver(h, s);
    E(irrep) = solver.eigenvalues();
    C(irrep, irrep) = solver.eigenvectors();
  }
  
  // ==== output ====
  PrintTimeStamp("Out", NULL);
  cout << "E0 = " << E(0)(0) << endl;
  E.Write(eigvals_out_name);
  C.Write(eigvecs_out_name);
  picojson::object out;
  out["in_json"] = picojson::value(argv[1]);
  out["eigvals"] = picojson::value(eigvals_out_name);
  out["eigvecs"] = picojson::value(eigvecs_out_name);
  picojson::value out_val(out);
  ofstream of(json_out_name.c_str(), ios::out);
  of << out_val.serialize();
  
}
