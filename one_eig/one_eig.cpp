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
int main (int argc, char *argv[]) {

  if(argc == 1) {
    PrintHelp();
    exit(1);
  }
  
  ifstream f;
  f.open(argv[1], ios::in);
  picojson::value json;  
  cout << ">>>> one_eig >>>>" << endl;
  // ==== parse json ====
  PrintTimeStamp("Parse", NULL);
  SymmetryGroup sym;
  SymGTOs gtos;
  Molecule mole;
  string comment, out_json, out_eigvecs, out_eigvals;

  try {
    f >> json;
    picojson::object& obj = json.get<picojson::object>();

    // -- Comment --
    comment = ReadJson<string>(obj, "comment");

    // -- Basis --
    sym = ReadJson<SymmetryGroup>(obj, "sym");
    
    mole = NewMolecule(sym);
    ReadJson_Molecule(obj, "molecule", mole);
    
    gtos = NewSymGTOs(mole);
    ReadJson_SymGTOs_Subs(obj, "basis", gtos);    

    // -- out --
    out_json = ReadJson<string>(obj, "out_json");
    out_eigvecs =ReadJson<string>(obj, "out_eigvecs");
    out_eigvals =ReadJson<string>(obj, "out_eigvals");    
    
  } catch(exception& e) {
    string msg = "error on parsing json\n";
    msg += e.what();
    cerr << msg << endl;
    exit(1);
  }
  try {
    gtos->SetUp();
  } catch(exception& e) {
    cerr << "error on setup" << endl;
    cerr << "error message:" << endl;
    cerr << e.what() << endl;
    cerr << "gtos::" << endl;
    cerr << gtos->str() << endl;
    exit(1);
  }
  
  cout << "comment: " << comment << endl;
  cout << "in_json: " << argv[1] << endl;
  cout << "out_json: " << out_json << endl;
  cout << "out_eigvecs: " << out_eigvecs << endl;
  cout << "out_eigvals: " << out_eigvals << endl;
  cout << "symmetry: " << sym->name() << endl;
  cout << "molecule: " << endl << mole->show() << endl;
  cout << "gtos: " << endl << gtos->show() << endl;

  // ==== calculation ====
  PrintTimeStamp("Calc", NULL);
  BMat S, T, V, C;
  BVec E;
  
  CalcSTVMat(gtos, gtos, &S, &T, &V);

  for(int irrep = 0; irrep < sym->order(); irrep++) {
    if(gtos->size_basis_isym(irrep) != 0) {
      MatrixXcd h = T(irrep, irrep) + V(irrep, irrep);
      MatrixXcd s = S(irrep, irrep);
      SymGenComplexEigenSolver solver(h, s);
      E(irrep) = solver.eigenvalues();
      C(irrep, irrep) = solver.eigenvectors();
      /*
      MatrixXcd Ctmp;
      VectorXcd Etmp;
      generalizedComplexEigenSolve(h, s, &Ctmp, &Etmp);
      cout << Etmp << endl;
      */

      ComplexEigenSolver<MatrixXcd> s_solver(s);
      cout << "eig_of_s: " << s_solver.eigenvalues()[0] << endl;
    }
  }
  
  // ==== output ====
  PrintTimeStamp("Out", NULL);
  cout << "E0 = " << E(0)(0) << endl;

  cout << "writing E" << endl;
  E.Write(out_eigvals);
  cout << "writing C" << endl;
  C.Write(out_eigvecs);
  cout << "out.json" << endl;
  
  picojson::object out;
  out["comment"] = picojson::value(comment);
  out["in_json"] = picojson::value(argv[1]);
  out["out_eigvals"] = picojson::value(out_eigvals);
  out["out_eigvecs"] = picojson::value(out_eigvecs);
  out["E0"] = ToJson(E(0)(0));
  int irrep0(0);
  out["irrep0"] = ToJson(irrep0);
  picojson::value out_val(out);
  ofstream of(out_json.c_str(), ios::out);
  of << out_val.serialize();
  cout << "<<<< one_eig <<<<" << endl;
}
