#include <stdexcept>
#include <fstream>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/QR>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>
#include "../utils/timestamp.hpp"
#include "../utils/eigen_plus.hpp"
#include "../utils/read_json.hpp"
#include "../r1basis/r1basis.hpp"

using namespace std;
using namespace boost;
using namespace Eigen;
using namespace picojson;
using namespace cbasis;

void DrivenTerm(string channel, string dipole, LC_STOs stos) {

  bool is_length = dipole == "length";
  bool is_velocity = dipole == "velocity";
  bool is_1skp = channel == "1s->kp";
  if(is_1skp && is_length) {
    stos->Add(2.0, 2, 1.0);
  } else if(is_1skp && is_velocity) {
    stos->Add(-2.0, 1, 1.0);
  } else {
    throw runtime_error("not supported combination of channel and dipole.");
  }
}

template<int MB, int MD>
class R1Driv {
public:
  typedef typename _EXPs<MB>::EXPs BasisEXPs;
  typedef typename _EXPs<MD>::LC_EXPs DrivLC_EXPs;
  object& obj;
  string in_file;
  string out_file;
  string comment;

  LinearSolver linear_solver;
  
  dcomplex Z;
  int L;
  dcomplex E;
  
  BasisEXPs basis;
  DrivLC_EXPs driv;

  string write_psi_filename;
  VectorXd write_psi_rs;

  R1Driv(object& _obj): obj(_obj), basis(new _EXPs<MB>()), driv(new _LC_EXPs<MD>()) {}
  void Parse() {
    PrintTimeStamp("Parse", NULL);

    // -- parse simple --            
    try {
      this->comment = ReadJson<string>(obj, "comment");
      this->linear_solver = ReadJsonWithDefault
	<LinearSolver>(obj, "linear_solver", LinearSolver());
      cout << format("comment: %s\n") % this->comment;
    } catch(std::exception& e) {
      cerr << "error on parsing json\n";
      cerr << e.what() << endl;
      exit(1);
    }
    
    // -- wave func --
    try {
      CheckObject<object>(obj, "write_psi");
      object& wavefunc_obj = obj["write_psi"].get<object>();
      this->write_psi_filename = ReadJson<string>(wavefunc_obj, "file");
      this->write_psi_rs = ReadJson<VectorXd>(wavefunc_obj, "rs");
    } catch(std::exception& e) {
      cerr << "error on parsing wavefunc\n";
      cerr << e.what() << endl;
      exit(1);
    }
    
    // -- parse target --
    try {      
      CheckObject<object>(obj, "target");
      object& target_obj = obj["target"].get<object>();
      string target_type = ReadJson<string>(target_obj, "type");
      if(target_type == "h_pi") {
	cout << "target_type: h_pi\n";
	string channel = ReadJson<string>(target_obj, "channel");
	cout << format("channel: %s\n") % channel;
	string dipole  = ReadJson<string>(target_obj, "dipole");
	cout << format("dipole: %s\n") % dipole;
	DrivenTerm(channel, dipole, this->driv);	
      } else if(target_type == "custom") {
	cout << format("target_type: %s\n") % target_type;

	// - driv -
	CheckObject<object>(target_obj, "driv");
	object& driv_obj = target_obj["driv"].get<object>();
	//string driv_type = ReadJso<string>(driv_obj, "type");
	CheckObject<picojson::array>(driv_obj, "value");
	picojson::array& vals = driv_obj["value"].get<picojson::array>();	  
	BOOST_FOREACH(value& val, vals) {
	  CheckValue<object>(val);
	  object& obj = val.get<object>();
	  dcomplex c = ReadJson<dcomplex>(obj, "c");
	  int      n = ReadJson<int>(obj, "n");
	  dcomplex z = ReadJson<dcomplex>(obj, "z");
	  this->driv->Add(c, n, z);
	}

	// - other -
	Z = ReadJson<dcomplex>(target_obj, "Z");
	cout << "Z: " << Z << endl;
	L = ReadJson<int>(target_obj, "L");
	cout << "L: " << L << endl;
	E = ReadJson<dcomplex>(target_obj, "E");
	cout << "E: " << E << endl;
	
      } else {
	throw runtime_error("not supported type in target");
      }
    } catch(std::exception& e) {
      cerr << "error on parsing target\n";
      cerr << e.what() << endl;
      exit(1);
    }

    // -- parse basis --
    try {
      CheckObject<object>(obj, "basis");
      object& basis_obj = obj["basis"].get<object>();
      CheckObject<picojson::array>(basis_obj, "value");
      picojson::array& basis_vals = basis_obj["value"].get<picojson::array>();
      BOOST_FOREACH(value& val, basis_vals) {
	CheckValue<object>(val);
	object& obj0 = val.get<object>();
	int n = ReadJson<int>(obj0, "n");
	cout << "n: " << n << endl;
	if(obj0.find("z") != obj0.end()) {
	  dcomplex z = ReadJson<dcomplex>(obj0, "z");
	  this->basis->AddPrim(n, z);
	  cout << "z: " << z << endl;
	} else if(obj0.find("czs") != obj0.end()) {
	  MatrixXcd czs = ReadJson<MatrixXcd>(obj0, "czs");
	  typename _EXPs<MB>::LC_EXPs lc(new _LC_EXPs<MB>());
	  for(int i = 0; i < czs.rows(); i++) {
	    lc->Add(czs(i, 0), n, czs(i, 1));
	  }
	  this->basis->AddLC(lc);
	  cout << "czs: " << endl << czs << endl;
	} else {
	  throw runtime_error("key \"z\" or \"czs\" are necessary in value in basis");
	}
      }
      
    } catch(std::exception& e) {
      cerr << "error on parsing basis\n";
      cerr << e.what() << endl;
      exit(1);
    }
    
  }
  void PrintIn() {
    PrintTimeStamp("PrintIn", NULL);
    cout << "linear_solver: " << linear_solver.show() << endl;
    cout << "write_psi_filename: " << this->write_psi_filename << endl;
    VectorXd& rs = this->write_psi_rs;
    cout << format("write_psi_rs: [%5.3f, %5.3f, ..., %5.3f]\n")% rs[0] % rs[1] % rs[rs.size()-1];    
    cout << "Z: " << Z << endl;
    cout << "L: " << L << endl;
    cout << "E: " << E << endl;
    cout << format("driven_term: %s\n") % driv->str();
    cout << format("basis:\n");
    cout << basis->str() << endl;
  }
  void Calc() {
    PrintTimeStamp("Calc", NULL);
    
    // -- set up basis --
    this->basis->SetUp();

    // -- build matrix --
    // D = E - T = E + 1/2 D2 - L(L+1)/2r2 +Z/r 
    MatrixXcd D, R2, R1, S;
    this->basis->InitMat(D);
    this->basis->InitMat(R2);
    this->basis->InitMat(R1);
    this->basis->InitMat(S);
    
    this->basis->CalcD2Mat(D);
    D *= 0.5;
    this->basis->CalcRmMat(-2, R2);
    D += -0.5*(L*(L+1))  * R2;
    this->basis->CalcRmMat(-1, R1);
    D += Z * R1;
    this->basis->CalcRmMat(0, S);
    D += E *S;

    // -- build vector --
    VectorXcd m; this->basis->InitVec(m);
    this->basis->CalcVec(driv, m);

    // -- solve linear problem --
    VectorXcd c;
    this->linear_solver.Solve(D, m, &c);

    // -- alpha --
    dcomplex alpha = TDot(m, c);
    cout << format("alpha: %20.15f, %20.15f\n") % alpha.real() % alpha.imag();    

    // -- wave function --
    VectorXd& rs = this->write_psi_rs;
    VectorXcd ys = this->basis->AtR(rs.cast<dcomplex>(), c);
    int num = ys.size();
    ofstream f(this->write_psi_filename.c_str(), ios::out);
    f << "r,re_y,im_y\n";
    for(int i = 0; i < num; i++) {
      f << format("%f,%f,%f\n") % rs[i] % ys[i].real() % ys[i].imag();
    }
  }  
  void Run() {
    this->Parse();
    this->PrintIn();
    this->Calc();
  }
};
int main(int argc, char *argv[]) {
  
  cout << ">>>> r1_driv >>>>\n";

  // -- argument number --
  if(argc != 2) {
    cout << "one argument is necessary\n";
    exit(1);
  }

  // -- open and read file --  
  ifstream f(argv[1]);
  if(f.fail()) {
    throw runtime_error("open file failed");
  }
  value json; f >> json;
  if(not json.is<object>()) {
    throw runtime_error("input json file is not object");
  }
  object& obj = json.get<object>();
  cout << format("in_file: %s\n") % argv[1];
  
  // -- check function type --
  try {
    CheckObject<object>(obj, "target");
    object& target_obj = obj["target"].get<object>();
    CheckObject<object>(obj, "basis");
    object& basis_obj = obj["basis"].get<object>();
    int MS = 1;
    string target_type = ReadJson<string>(target_obj, "type");
    if(target_type == "h_pi") {
      MS = 1;
    } else if(target_type == "custom") {
      CheckObject<object>(target_obj, "driv");
      object& driv_obj = target_obj["driv"].get<object>();
      string driv_type = ReadJson<string>(driv_obj, "type");
      if(driv_type == "STO")
	MS = 1;
      else if(driv_type == "GTO")
	MS = 2;
      else
	throw runtime_error("unsupported type in driv in target");
    }
    string basis_type = ReadJson<string>(basis_obj, "type");
    int MB = 1;
    if(basis_type == "STO") {
      MB = 1;
    } else if(basis_type == "GTO") {
      MB = 2;
    } else {
      cerr << "unsupported type in basis\n";
      exit(1);
    }
    
    // -- start main --
    if(MB == 1 && MS == 1) {
      R1Driv<1,1> r1driv(obj);
      r1driv.Run();
    }
  } catch(std::exception& e) {
    cerr << "error on calculation\n";
    cerr << e.what() << endl;
    exit(1);
  }
  
  //  r1driv.in_file = argv[1];
  //  
  //  r1driv.Parse(obj);
  //  cout << "<<<< r1_driv <<<<\n";
}

