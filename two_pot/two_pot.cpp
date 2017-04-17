#include <fstream>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <Eigen/Core>
#include <Eigen/QR>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include "../utils/eigen_plus.hpp"
#include "../utils/timestamp.hpp"
#include "../utils/read_json.hpp"
#include "../src_cpp/symmolint.hpp"
#include "../src_cpp/angmoment.hpp"
#include "../src_cpp/mo.hpp"
#include "../src_cpp/one_int.hpp"
#include "../src_cpp/two_int.hpp"
#include "../src_cpp/symmol_read_json.hpp"

using namespace std;
using namespace Eigen;
using namespace cbasis;
using namespace picojson;
using namespace boost::assign;
using boost::format;

// -- const --
double au2ev = 27.2114;
double c_light = 137.035999139;
double au2mb = 5.291772 * 5.291772;

// -- Input --
string comment;
string in_json, out_json, cs_csv;

// -- Output --
ofstream f_para;

// -- target --
SymmetryGroup sym;
Molecule mole;
dcomplex Z;
Molecule mole0;
vector<double> w_list;
int ne;
string calc_type;
//bool use_stex;
enum ECalcTerm {ECalcTerm_One, ECalcTerm_Full };
ECalcTerm calc_term;
ERIMethod eri_method;

// -- Initial state --
SymGTOs basis0;
Irrep irrep0;
int i0;
dcomplex E0;
VectorXcd c0;

// -- psi1 --
SymGTOs basis1;

// -- psi0/chi0 --
vector<int> Ls;
map<int, SymGTOs> basis_psi0_L;
map<int, SymGTOs> basis_c_psi0_L;
map<int, SymGTOs> basis_chi0_L;

// -- solver --
string solve_driv_type;
LinearSolver linear_solver;

// -- write_psi --
string write_psi_filename;
VectorXd write_psi_rs;
int   write_psi_lmax;

// -- intermediate --
// ---- for psi1 ----
BMat S1, T1, V1, L1;
BVec sX1, sY1, sZ1, sDX1, sDY1, sDZ1;
BVec cX1, cY1, cZ1, cDX1, cDY1, cDZ1;

// ---- for psi0_L --
map<int, BMat> S0L, T0L, V0L, L0L;
map<int, BMat> V0L1, HV0L1;
map<int, BVec> sX0L, sY0L, sZ0L, sDX0L, sDY0L, sDZ0L;
map<int, BVec> c0L, Hc0L;     // coef for psi0_p
map<int, BVec> s0L_chi; // driven term for psi_p

MultArray<dcomplex, 2> impsi0_muphi(100),  impsi0_muphi_v(100);
MultArray<dcomplex, 2> impsi0_v_psi1(100), impsi0_v_psi1_v(100);
MultArray<dcomplex, 2> impsi0_chi(100);

// -- results --
MatrixXd result;
int num_header = 14;
int idx_cs(0), idx_cs_sigu(1), idx_cs_piu(2), idx_beta(3);
int idx_cs_alpha(4), idx_cs_sigu_alpha(5), idx_cs_piu_alpha(6);
int idx_cs_v(7), idx_cs_sigu_v(8), idx_cs_piu_v(9), idx_beta_v(10);
int idx_cs_alpha_v(11), idx_cs_sigu_alpha_v(12), idx_cs_piu_alpha_v(13);

void s2y(MultArray<dcomplex, 1>& sM,
	 MultArray<dcomplex, 1> *yM, int L) {
  /**
     compute spherical harmonics values from solid spherical harmonics values.
     
     case M > 0
     <0 | Y_LM> = (-)^M / sqrt(2) * (<0|S_LM + j S_L-M>)
     .          = (-)^M / sqrt(2) * (<0|S_LM >+ j <0|S_L-M>)
     <0 | Y_L-M> =   1.0 / sqrt(2) * (<0|S_LM > - j <0|S_L-M>)
   */

  dcomplex ii(0, 1);
  (*yM)(0) = sM(0);
  for(int M = 1; M <= L; M++) {
    (*yM)(M) = pow(-1, M) * (sM(M) + ii*sM(-M)) / sqrt(2.0);
    (*yM)(-M)=              (sM(M) - ii*sM(-M)) / sqrt(2.0);
  }
}
void s2y_bra(MultArray<dcomplex, 1>& sM,
	 MultArray<dcomplex, 1> *yM, int L) {
  /**
     <Y_LM | 0> = <0^* |Y_LM^*> = 
   */

  dcomplex ii(0, 1);
  (*yM)(0) = sM(0);
  for(int M = 1; M <= L; M++) {
    (*yM)(M) = pow(-1, M) * (sM(M) - ii*sM(-M)) / sqrt(2.0);
    (*yM)(M) =              (sM(M) + ii*sM(-M)) / sqrt(2.0);
  }

}
void ss0_2_yy0(MultArray<dcomplex, 1>& sM,
	       MultArray<dcomplex, 1> *yM, int L) {
  /**
     <S_L0 | S_L'0 | 0> = <Y_L0 | Y_L'0 | 0>
     <S_L+M | S_L'0 | 0> = 0
     <S_L+M | S_L'M | 0> 
     .= 2^(-0.5) { (-)^M <Y_L+M|S_L'M|0> -j(-)^M<Y_L-M|S_L'M|0> }
     .= 1/2 {    <Y_L+M|Y_L'+M|0> + j <Y_L+M|Y_L'-M|0>
     .        -j (-)^M <Y_L-M|Y_L'M|0> -1(-)^M}
     <S_L-M | S_L'0 | 0> 
     .= 2^(-0.5)      { <Y_L+M|S_L'0|0> +j <Y_L-M|S_L'0|0> }
     <S_L+M | S_L' | 0> 
     .= 2^(-0.5)(-)^M { <Y_L+M|S_L'0|0> -j <Y_L-M|S_L'0|0> }
     <S_L+M | S_L'0 | 0> = 2^(-0.5)(-)^M {  <Y_L+M | S_L'0 | 0> 
     .                                    -j<Y_L-M | S_L'0 | 0> }


     <Y_+M | Y_M | 0> 
     .  = 2^(-1) { <S_+M | S_M | 0> + <S_-M | S_-M | 0> }
     <Y_-M | Y_-M | 0> 
     .  = 2^(-1) { <S_+M | S_M | 0> + <S_-M | S_-M | 0> }

   */
}
double CoulombShift(dcomplex eta, int L) {
  // -- eta = Z1*Z2/k
  dcomplex ii(0,1);
  dcomplex argment = 1.0 + L + ii * eta;
  double re = argment.real();
  double im = argment.imag();
  gsl_sf_result lnr, arg;
  gsl_sf_lngamma_complex_e(re, im, &lnr, &arg);
  return arg.val;
}
dcomplex CoefAzeta0(dcomplex w, MultArray<dcomplex, 2>& Alm,
		    int zeta,
		    vector<int>& L_list, vector<int>& Ms) {
  /**
     see
     J. C. Tully, R. S. Berry, and B. J. Dalton, 
     Physical Review 176, 95 (1968)
  */
  
  dcomplex c0 = 4.0 * M_PI * M_PI / (c_light*w);
  dcomplex cumsum(0);
  BOOST_FOREACH(int L1, L_list) {
    BOOST_FOREACH(int L2, L_list) {
      dcomplex c1 = sqrt(1.0 * (2*L1+1)*(2*L2+1)) / (4*M_PI*(2*zeta+1));
	
      BOOST_FOREACH(int M1, Ms) {
	BOOST_FOREACH(int M2, Ms) {
	  int mzeta = -M1+M2;
	  dcomplex c3 = (cg_coef(1,  1,-M1,M2,  zeta,mzeta) *
			 cg_coef(1,  1,0,0,     zeta,0) *
			 cg_coef(L2,L1,+M2,-M1, zeta,mzeta) *
			 cg_coef(L2,L1,0,0,     zeta,0));
	  dcomplex c4 = Alm(L1,M1) * conj(Alm(L2,M2));
	  cumsum += c0*c1*c3*c4;
	}
      }
    }
  }
  
  return cumsum;
}
double Coef_FixedMole(dcomplex w, MultArray<dcomplex, 2>& Alm,
		      vector<int>& L_list, int L, int M ) {
  /**
     Gives coefficient of differential cross section for 
     fixed molecule axis with paralled electronic field.
     see note @ 2017/3/15

     Inputs
     ------
     Alm(L, M) = w i sqrt(2/3) (i^-L) e^{i eta_l} <phi^-, mu phi_i>

     d sigma/d Omega = sum_L C_L P_L

   */

  dcomplex cumsum(0);
  BOOST_FOREACH(int L1, L_list) {    
    BOOST_FOREACH(int L2, L_list) {
      dcomplex c1 = sqrt(1.0*(2*L1+1)*(2*L2+1));
      dcomplex c2 = cg_coef(L1,L2,0,0, L,0) * cg_coef(L1,L2,M,-M, L,0);
      dcomplex c3 = Alm(L1, M) * conj(Alm(L2, M));
      cumsum += c1 * c2 * c3;
    }
  }
  return (cumsum).real();
}
void Parse() {
  PrintTimeStamp("Parse", NULL);  
  try {
    ifstream f(in_json.c_str());  
    value json; f >> json;
    if(not json.is<object>()) 
      throw runtime_error("invalid json file");
    object& obj = json.get<object>();
    comment = ReadJson<string>(obj, "comment");
    calc_type = ReadJson<string>(obj, "calc_type");

    // -- write coeff for molecule fixed with parallel electronic field --
    string para_csv = ReadJsonWithDefault<string>(obj, "para_csv", "");
    cout << "para_csv: " << para_csv << endl;
    if(para_csv != "") {
      f_para.open(para_csv.c_str());
      if(f_para.fail()) {
	throw runtime_error("failed to open file.");
      }
    }
    
    // -- write wave function --
    CheckObject<object>(obj, "write_psi");
    object& write_psi_obj = obj["write_psi"].get<object>();
    write_psi_filename = ReadJson<string>(write_psi_obj, "file");
    write_psi_rs = ReadJson<VectorXd>(write_psi_obj, "rs");
    write_psi_lmax = ReadJson<int>(write_psi_obj, "lmax");

    // -- number of calculation term --
    string str_calc_term = ReadJsonWithDefault<string>(obj, "calc_term", "full");
    if(str_calc_term == "one") {
      calc_term = ECalcTerm_One;
    } else if (str_calc_term == "full") {
      calc_term = ECalcTerm_Full;
    } else {
      throw runtime_error("calc_term must be \"one\" or \"full\"");
    }

    // solver for driven equation
    solve_driv_type = ReadJsonWithDefault
      <string>(obj, "solve_driv", "linear_solve");

    if(solve_driv_type == "linear_solve") {
      linear_solver = ReadJsonWithDefault
	<LinearSolver>(obj, "linear_solver", LinearSolver());
    } else if(solve_driv_type == "eigen_value") {

    } else {
      THROW_ERROR("unsupported solve_driv");
    }
    ne = ReadJson<int>(obj, "num_ele");

    if(calc_type == "STEX" or calc_type == "RPA") {
      eri_method = ReadJson<ERIMethod>(obj, "eri_method");
    }
    sym     = ReadJson<SymmetryGroup>(obj, "sym");
    out_json = ReadJson<string>(obj, "out_json");
    cs_csv   = ReadJson<string>(obj, "cs_csv");
    mole = NewMolecule(sym); ReadJson_Molecule(obj, "molecule", mole);
    Z = ReadJson<dcomplex>(obj, "Z");
    int lmax = ReadJsonWithDefault<int>(obj, "lmax", 3);
    for(int L = 1; L <= lmax; L += 2) {
      Ls.push_back(L);
    }
    cout << "Ls: ";
    BOOST_FOREACH(int L, Ls) {
      cout << L << " " ;
    }
    cout << endl;
    cout << "Ls0 = " << Ls[0] << endl;
    cout << "Ls1 = " << Ls[Ls.size()-1] << endl;
    impsi0_muphi.SetRange( Ls[0], Ls[Ls.size()-1], -1,1);
    impsi0_muphi_v.SetRange( Ls[0], Ls[Ls.size()-1], -1,1);

    impsi0_v_psi1.SetRange(Ls[0], Ls[Ls.size()-1], -1,1);
    impsi0_v_psi1_v.SetRange(Ls[0], Ls[Ls.size()-1], -1,1);
    impsi0_chi.SetRange(   Ls[0], Ls[Ls.size()-1], -1,1);

    cout << "mole" << endl;
    mole0 = NewMolecule(sym);
    Atom cen0 = NewAtom("CEN0", Z); cen0->Add(0, 0, 0);
    mole0->Add(cen0);
    
    basis0 = NewSymGTOs(mole); ReadJson_SymGTOs_Subs(obj, "basis0", basis0);
    basis1 = NewSymGTOs(mole); ReadJson_SymGTOs_Subs(obj, "basis1", basis1);
    ReadJson_Orbital(obj, "orbital0", sym, &E0, &c0, &irrep0, &i0);

    VectorXi Ms(3); Ms << -1, 0, 1;
    VectorXcd zeta_chi(1); zeta_chi << ReadJson<dcomplex>(obj, "zeta_chi");
    BOOST_FOREACH(int L, Ls) {
      string name;
      if(L == 1) 
	name = "p";
      else if(L == 3)
	name = "f";
      else if(L == 5)
	name = "h";
      else
	throw runtime_error("unsupported L");
      cout << "start raeding psi0 for " << name << endl;
      VectorXcd zeta0 = ReadJson<VectorXcd>(obj, "zeta0_" + name);
      basis_psi0_L[L] = NewSymGTOs(mole0);
      basis_psi0_L[L]->NewSub("CEN0").SolidSH_Ms(L, Ms, zeta0);
      basis_chi0_L[L] = NewSymGTOs(mole0);
      basis_chi0_L[L]->NewSub("CEN0").SolidSH_Ms(L, Ms, zeta_chi);
      cout << "end raeding psi0 for " << name << endl;
    }

    VectorXd _ws = ReadJson<VectorXd>(obj, "ws");
    for(int i = 0; i < _ws.size(); i++) {
      w_list.push_back(_ws(i));
    }
    result = MatrixXd::Zero(w_list.size(), num_header);

    
  } catch(exception& e) {
    cerr << "error on parsing json" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  cout << "basis0 and basis1 setup" << endl;
  try {
    basis0->SetUp();    
    basis1->SetUp();
  } catch(exception& e) {
    cerr << "error on SetUp basis: " << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  cout << "basis_psi0_L setup" << endl;
  try {
    BOOST_FOREACH(int L, Ls) {
      basis_psi0_L[L]->SetUp();
    }
  } catch(exception& e) {
    cerr << "error on SetUp basis: " << endl;
    cerr << e.what() << endl;
    cerr << basis_psi0_L[1]->str() << endl;
    cerr << basis_psi0_L[3]->str() << endl;
    exit(1);
  }
  cout << "basis_chi setup" << endl;
  try {
    BOOST_FOREACH(int L, Ls) {
      basis_chi0_L[L]->SetUp();
    }
  } catch(exception& e) {
    cerr << "error on SetUp basis_chi0_L" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  cout << "basis_c_psi0 conj" << endl;
  BOOST_FOREACH(int L, Ls) {
    basis_c_psi0_L[L] = basis_psi0_L[L]->Conj();
  }
  
}
void PrintIn() {
  PrintTimeStamp("PrintIn", NULL);
  cout << "comment: " << comment << endl;
  cout << "in_json: " << in_json << endl;
  cout << "write_psi_filename: " << write_psi_filename << endl;
  VectorXd& rs = write_psi_rs;
  cout << format("write_psi_rs: [%5.3f, %5.3f, ..., %5.3f]\n") % rs[0] % rs[1] % rs[rs.size()-1];
  cout << "write_psi_lmax: " << write_psi_lmax << endl;
  cout << "cs_csv: "  << cs_csv;
  cout << "out_json: " << out_json << endl;
  cout << "calc_type: " << calc_type << endl;
  cout << "calc_term: " << (calc_term == ECalcTerm_One ? "one" : "full") << endl;
  cout << "linear_solver: " << linear_solver.str() << endl;
  cout << "ERIMethod_use_symmetry: " << eri_method.symmetry << endl;
  cout << "ERIMethod_use_memo: " << eri_method.coef_R_memo << endl;
  cout << "ERIMethod_use_perm: " << eri_method.perm << endl;  
  cout << "Ne: " << ne << endl;
  cout << "E0: " << E0 << endl;
  cout << "Z: " << Z << endl;
  cout << "symmetry: " << sym->name() << endl;
  cout << "molecule:" << endl << mole->show() << endl;
  cout << "basis0:" << endl << basis0->show() << endl;  
  cout << "basis1:" << endl << basis1->show() << endl;
  BOOST_FOREACH(int L, Ls) {
    cout << "basis0_:" << L << endl << basis_psi0_L[L]->show() << endl;
  }
}
void CalcMat() {

  S1.set_name("S1"); T1.set_name("T1"); V1.set_name("V1");
  
  Irrep x = sym->irrep_x(); Irrep y = sym->irrep_y(); Irrep z = sym->irrep_z();

  PrintTimeStamp("psi1", NULL);
  CalcSTVMat(basis1, basis1, &S1, &T1, &V1);
  
  PrintTimeStamp("psi1/init", NULL);
  BMat X1i, DX1i, Y1i, DY1i, Z1i, DZ1i;
  CalcDipMat(basis1, basis0, &X1i, &Y1i, &Z1i, &DX1i, &DY1i, &DZ1i);
  sX1(x) = X1i(x, 0) * c0; sDX1(x) = DX1i(x, 0) * c0; 
  sY1(y) = Y1i(y, 0) * c0; sDY1(y) = DY1i(y, 0) * c0; 
  sZ1(z) = Z1i(z, 0) * c0; sDZ1(z) = DZ1i(z, 0) * c0;  

  PrintTimeStamp("psi0,chi0", NULL);
  BOOST_FOREACH(int L, Ls) {
    cout << "L = " << L << endl;
    SymGTOs psi0   = basis_psi0_L[L];
    SymGTOs c_psi0 = basis_c_psi0_L[L];
    SymGTOs chi0 = basis_chi0_L[L];
    CalcSTVMat(psi0, psi0, &S0L[L], &T0L[L], &V0L[L]);

    BMat S0L_chi; CalcSMat(psi0, chi0, &S0L_chi);
    s0L_chi[L](x) = S0L_chi(x, x).col(0);
    s0L_chi[L](y) = S0L_chi(y, y).col(0);
    s0L_chi[L](z) = S0L_chi(z, z).col(0);

    BMat X0i, DX0i, Y0i, DY0i, Z0i, DZ0i;
    CalcDipMat(psi0, basis0, &X0i, &Y0i, &Z0i, &DX0i, &DY0i, &DZ0i);
    sX0L[L](x) = X0i(x,0) * c0; sDX0L[L](x) = DX0i(x,0) * c0;
    sY0L[L](y) = Y0i(y,0) * c0; sDY0L[L](y) = DY0i(y,0) * c0;
    sZ0L[L](z) = Z0i(z,0) * c0; sDZ0L[L](z) = DZ0i(z,0) * c0;

    BMat V01_full, HV01_full, V01_0th, HV01_0th;
    CalcVMat(psi0,   mole,  basis1, &V01_full);
    CalcVMat(c_psi0, mole,  basis1, &HV01_full);
    CalcVMat(psi0,   mole0, basis1, &V01_0th);
    CalcVMat(c_psi0, mole0, basis1, &HV01_0th);

    V0L1[L](x,x) = V01_full(x,x) - V01_0th(x,x);
    V0L1[L](y,y) = V01_full(y,y) - V01_0th(y,y);
    V0L1[L](z,z) = V01_full(z,z) - V01_0th(z,z);
    
    HV0L1[L](x,x) = HV01_full(x,x) - HV01_0th(x,x);
    HV0L1[L](y,y) = HV01_full(y,y) - HV01_0th(y,y);
    HV0L1[L](z,z) = HV01_full(z,z) - HV01_0th(z,z);
  }
}
void CalcMatSTEX() {
  PrintTimeStamp("MatSTEX_1", NULL);

  B2EInt eri_J_11 = CalcERI(basis1, basis1, basis0, basis0, eri_method);
  B2EInt eri_K_11 = CalcERI(basis1, basis0, basis0, basis1, eri_method);
  AddJ(eri_J_11, c0, irrep0, 1.0, V1); AddK(eri_K_11, c0, irrep0, 1.0, V1);

  PrintTimeStamp("MatSTEX_01", NULL);
  BOOST_FOREACH(int L, Ls) {
    cout << "L = " << L << endl;
    SymGTOs psi0   = basis_psi0_L[L];
    SymGTOs c_psi0 = basis_c_psi0_L[L];
    B2EInt eri_JC = CalcERI(psi0,   basis1, basis0, basis0, eri_method);
    B2EInt eri_JH = CalcERI(c_psi0, basis1, basis0, basis0, eri_method);
    B2EInt eri_KC = CalcERI(psi0,   basis0, basis0, basis1, eri_method);
    B2EInt eri_KH = CalcERI(c_psi0, basis0, basis0, basis1, eri_method);
    AddJ(eri_JC, c0, irrep0, 1.0, V0L1[L]);  AddK(eri_KC, c0, irrep0, 1.0, V0L1[L]);
    AddJ(eri_JH, c0, irrep0, 1.0, HV0L1[L]); AddK(eri_KH, c0, irrep0, 1.0, HV0L1[L]);
  }
}
void CalcMatRPA() {
  PrintTimeStamp("MatRPA_1", NULL);

  // ==== psi1 ====
  B2EInt eri_J_11 = CalcERI(basis1, basis1, basis0, basis0, eri_method);
  B2EInt eri_K_11 = CalcERI(basis1, basis0, basis0, basis1, eri_method);
  
  // -- compute A matrix 
  // -- A = T + V_ne + J + K - epsilon0  
  CalcSTVMat(basis1, basis1, &S1, &T1, &V1);  
  BMat A("A"); Copy(T1, A);    
  try {
    A.Add(1.0, V1);
  } catch(exception& e) {
    cout << "error on A+=V1\n";
    cout << e.what() << endl;
    THROW_ERROR("ERROR");
  } 
  AddJ(eri_J_11, c0, irrep0, 1.0, A);
  AddK(eri_K_11, c0, irrep0, 1.0, A);
  try {
    A.Add(-E0, S1);
  } catch(exception& e) {
    cout << e.what() << endl;
    THROW_ERROR("error on A+=-E0 S1")
      }

  // -- compute B matrix --
  // -- B = K
  BMat B("B"); Copy(T1, B); B.SetZero();  
  AddK(eri_K_11, c0, irrep0, 1.0, B);

  // -- compute V0 --
  BMat V1_0th("V1_0th");
  CalcVMat(basis1, mole0, basis1, &V1_0th);
  
  // -- compute interaction for psi1 --
  // -- (H-T) = sqrt((A+B)(A-B)) - T + epsilon0
  BMat ApB("ApB"); Copy(A, ApB); ApB.Add(+1.0, B);  
  BMat AmB("AmB"); Copy(A, AmB); AmB.Add(-1.0, B);
  BMat H2("H2");  Multi(ApB, AmB, H2);
  BMatSqrt(H2, V1);
  V1.Add(+E0, S1);

  try {
    V1.Add(-1.0, T1);
  } catch(exception& e) {
    cout << e.what() << endl;
    THROW_ERROR("V1+=-T1");
  }
  
  // ==== psi0/psi1 ====  
  PrintTimeStamp("MatRPA_2", NULL);
  BOOST_FOREACH(int L, Ls) {
    cout << "L = " << L << endl;
    SymGTOs psi0   = basis_psi0_L[L];
    SymGTOs c_psi0 = basis_c_psi0_L[L];

    BMat T01, HT01, S01, HS01, V01, HV01;
    BMat V01_0th, HV01_0th;
    CalcSTMat(psi0,   basis1, &S01,  &T01);
    CalcSTMat(c_psi0, basis1, &HS01, &HT01);
    CalcVMat(psi0,    mole,  basis1, &V01);
    CalcVMat(c_psi0,  mole,  basis1, &HV01);
    CalcVMat(psi0,    mole0, basis1, &V01_0th);
    CalcVMat(c_psi0,  mole0, basis1, &HV01_0th);

    // -- ERI --
    B2EInt eri_JC = CalcERI(psi0,   basis1, basis0, basis0, eri_method);
    B2EInt eri_JH = CalcERI(c_psi0, basis1, basis0, basis0, eri_method);
    B2EInt eri_KC = CalcERI(psi0,   basis0, basis0, basis1, eri_method);
    B2EInt eri_KH = CalcERI(c_psi0, basis0, basis0, basis1, eri_method);

    // -- A matrix --
    // -- A = T + V_ne + J + K - epsilon0  
    BMat CA; Copy(T01,  CA); CA.Add(1.0, V01);
    BMat HA; Copy(HT01, HA); HA.Add(1.0, HV01);    
    AddJ(eri_JC, c0, irrep0, 1.0, CA); AddK(eri_KC, c0, irrep0, 1.0, CA);
    AddJ(eri_JH, c0, irrep0, 1.0, HA); AddK(eri_KH, c0, irrep0, 1.0, HA);
    CA.Add(-E0, S01);
    HA.Add(-E0, HS01);

    // -- B matrix --
    BMat CB; Copy(T1, CB); CB.SetZero();
    BMat HB; Copy(T1, HB); HB.SetZero();
    AddK(eri_KC, c0, irrep0, 1.0, CB);
    AddK(eri_KH, c0, irrep0, 1.0, HB);

    // -- interaction --
    // V = H-H_0 = sqrt((A+B)(A-B)) + epsilon - T - V0
    BMat C_ApB; Copy(CA, C_ApB); C_ApB.Add(+1.0, CB);
    BMat C_AmB; Copy(CA, C_AmB); C_AmB.Add(-1.0, CB);
    BMat C_H2; Multi(C_ApB, C_AmB, C_H2);
    BMatSqrt(C_H2, V0L1[L]);
    V0L1[L].Add(E0, S01);
    V0L1[L].Add(-1.0, T01);
    V0L1[L].Add(-1.0, V01_0th);    
    BMat H_ApB; Copy(HA, H_ApB); H_ApB.Add(+1.0, HB);
    BMat H_AmB; Copy(HA, H_AmB); H_AmB.Add(-1.0, HB);
    BMat H_H2; Multi(H_ApB, H_AmB, H_H2);
    BMatSqrt(H_H2, HV0L1[L]);
    HV0L1[L].Add(E0,   HS01);
    HV0L1[L].Add(-1.0, HT01);
    HV0L1[L].Add(-1.0, HV01_0th);
    

  }
  
}
void CalcDriv(int iw) {
  /**
     inputs:
     .  S1, T1, V1,  S0L, T0L, V0L

     outputs:
     . c[D][XYZ]1, c0L
   */
  
  PrintTimeStamp("CalcDriv", NULL);
  if(solve_driv_type == "linear_solve") {
    CalcDriv_linear_solve(iw);
  } else if(solve_driv_type == "eigen_value") {
    CalcDriv_eigen_value(iw);
  }
}
void CalcDriv_linear_solve(int iw) {
  PrintTimeStamp("CalcDriv_linear_solve", NULL);
  double w = w_list[iw];
  Irrep x = sym->irrep_x();
  Irrep y = sym->irrep_y();
  Irrep z = sym->irrep_z(); 

  dcomplex ene = E0 + w;
  
  // -- Compute psi1 --
  L1(x,x) = S1(x,x) * ene - T1(x,x) - V1(x,x);
  ColPivHouseholderQR<MatrixXcd> pivx = L1(x,x).colPivHouseholderQr();
  cX1(x)  = pivx.solve(sX1(x));
  cDX1(x) = pivx.solve(sDX1(x));
  
  L1(y,y) = S1(y,y) * ene - T1(y,y) - V1(y,y);
  ColPivHouseholderQR<MatrixXcd> pivy = L1(y,y).colPivHouseholderQr();
  cY1(y) = pivy.solve(sY1(y));
  cDY1(y) = pivy.solve(sDY1(y));    
  
  L1(z,z) = S1(z,z) * ene - T1(z,z) - V1(z,z);
  ColPivHouseholderQR<MatrixXcd> pivz = L1(z,z).colPivHouseholderQr();
  cZ1(z)  = pivz.solve(sZ1(z));
  cDZ1(z) = pivz.solve(sDZ1(z));

  // -- Compute psi0_p --
  BOOST_FOREACH(int L, Ls) {
    L0L[L](x,x) = S0L[L](x,x) * ene - T0L[L](x,x) - V0L[L](x,x);
    c0L[L](x) = L0L[L](x,x).colPivHouseholderQr().solve(s0L_chi[L](x));
    Hc0L[L](x) = c0L[L](x).conjugate();

    L0L[L](y,y) = S0L[L](y,y) * ene - T0L[L](y,y) - V0L[L](y,y);
    c0L[L](y)   = L0L[L](y,y).colPivHouseholderQr().solve(s0L_chi[L](y));
    Hc0L[L](y)  = c0L[L](y).conjugate();

    L0L[L](z,z) = S0L[L](z,z) * ene - T0L[L](z,z) - V0L[L](z,z);
    c0L[L](z)   = L0L[L](z,z).colPivHouseholderQr().solve(s0L_chi[L](z));
    Hc0L[L](z)  = c0L[L](z).conjugate();
  }
}
void CalcDriv_eigen_value(int iw) {
  PrintTimeStamp("CalcDriv_linear_solve", NULL);
  double w = w_list[iw];
  Irrep x = sym->irrep_x();
  Irrep y = sym->irrep_y();
  Irrep z = sym->irrep_z(); 

  dcomplex ene = E0 + w;
  
}
void CalcBraket() {

  // <ImPsi|0> = <(y-Hy)/2j|0>
  //           = -1/2j*(<y|0>-<Hy|0>)
  //           = -1/2j*((Hy|0)-(y|0))
  //           = -1/2j*((Hy|0)-(y|0))  
  
  //  PrintTimeStamp("calc_braket", NULL);
  Irrep x = sym->irrep_x();
  Irrep y = sym->irrep_y();
  Irrep z = sym->irrep_z();
  dcomplex m2(0, 2);
  
  BOOST_FOREACH(int L, Ls) {
    BVec& c0 = c0L[L];
    BVec& Hc0 = Hc0L[L];

    // -- <ImPsi0 | mu | PhiInit> --  
    impsi0_muphi(L, +1) = TDot(c0(x), sX0L[L](x)).imag();
    impsi0_muphi(L, -1) = TDot(c0(y), sY0L[L](y)).imag();
    impsi0_muphi(L,  0) = TDot(c0(z), sZ0L[L](z)).imag();
    
    impsi0_muphi_v(L, +1) = TDot(c0(x), sDX0L[L](x)).imag();
    impsi0_muphi_v(L, -1) = TDot(c0(y), sDY0L[L](y)).imag();
    impsi0_muphi_v(L,  0) = TDot(c0(z), sDZ0L[L](z)).imag();

    BMat& V  = V0L1[L];
    BMat& HV = HV0L1[L];
    //cout << "L = " << L << endl;
    //cout << "V mat: \n" << V.full_print() << endl;
    //cout << "VH mat: \n" << HV.full_print() << endl;
    
    impsi0_v_psi1(L, +1) = (TDot(c0(x),   V( x,x)*cX1(x))
			    -TDot(Hc0(x), HV(x,x)*cX1(x)))/m2;
    impsi0_v_psi1(L, -1) = (TDot(c0(y),   V( y,y)*cY1(y))
			    -TDot(Hc0(y), HV(y,y)*cY1(y)))/m2;
    impsi0_v_psi1(L,  0) = (TDot(c0(z),   V( z,z)*cZ1(z))
			    -TDot(Hc0(z), HV(z,z)*cZ1(z)))/m2;

    impsi0_v_psi1_v(L, +1) = (TDot(c0(x),   V( x,x)*cDX1(x))
			      -TDot(Hc0(x), HV(x,x)*cDX1(x)))/m2;
    impsi0_v_psi1_v(L, -1) = (TDot(c0(y),   V( y,y)*cDY1(y))
			      -TDot(Hc0(y), HV(y,y)*cDY1(y)))/m2;
    impsi0_v_psi1_v(L,  0) = (TDot(c0(z),   V( z,z)*cDZ1(z))
			      -TDot(Hc0(z), HV(z,z)*cDZ1(z)))/m2;

    
    impsi0_chi(L, +1) = TDot(c0(x), s0L_chi[L](x)).imag();
    impsi0_chi(L, -1) = TDot(c0(y), s0L_chi[L](y)).imag();
    impsi0_chi(L,  0) = TDot(c0(z), s0L_chi[L](z)).imag();    
    
  }

}
void CalcMain_alpha(int iw) {
  
  //  PrintTimeStamp("calc_alpha", NULL);  
  double w = w_list[iw];
  Irrep x = sym->irrep_x();
  Irrep y = sym->irrep_y();
  Irrep z = sym->irrep_z();
  
  dcomplex alpha_x = TDot(cX1(x), sX1(x))/3.0*(1.0*ne);
  dcomplex alpha_y = TDot(cY1(y), sY1(y))/3.0*(1.0*ne);
  dcomplex alpha_z = TDot(cZ1(z), sZ1(z))/3.0*(1.0*ne);
  dcomplex alpha = (alpha_x+alpha_y+alpha_z);
  double cs = 4.0*M_PI*w/ c_light * alpha.imag()* au2mb;
  double cs_sigu = 4.0*M_PI*w/ c_light * alpha_z.imag()* au2mb;
  double cs_piu  = 4.0*M_PI*w/ c_light * (alpha_x+alpha_y).imag()* au2mb;  
  result(iw, idx_cs_alpha) = cs;
  result(iw, idx_cs_sigu_alpha) = cs_sigu;
  result(iw, idx_cs_piu_alpha) = cs_piu;  

  dcomplex alpha_dx = TDot(cDX1(x), sDX1(x))/3.0*(1.0*ne);
  dcomplex alpha_dy = TDot(cDY1(y), sDY1(y))/3.0*(1.0*ne);
  dcomplex alpha_dz = TDot(cDZ1(z), sDZ1(z))/3.0*(1.0*ne);  
  dcomplex alpha_d = (alpha_dx+alpha_dy+alpha_dz);
  double cs_v      = 4.0*M_PI/(w*c_light) * alpha_d.imag() * au2mb;
  double cs_sigu_v = 4.0*M_PI/(w*c_light) * alpha_dz.imag() * au2mb;
  double cs_piu_v  = 4.0*M_PI/(w*c_light) * (alpha_dx+alpha_dy).imag() * au2mb;
  result(iw, idx_cs_alpha_v)      = cs_v;
  result(iw, idx_cs_sigu_alpha_v) = cs_sigu_v;
  result(iw, idx_cs_piu_alpha_v)  = cs_piu_v;

  cout << format("Cs(total,alpha): %20.10f, %20.10f\n") % cs % cs_v;
  cout << format("Cs(sigu,alpha):  %20.10f, %20.10f\n") % cs_sigu % cs_sigu_v;
  cout << format("Cs(piu,alpha):   %20.10f, %20.10f\n") % cs_piu % cs_piu_v;
  
}
void CalcMain(int iw) {
  // -- element of dl_ss is solid spherical, but
  // -- in the special case, it is equivalent to spherical haromonics
  
  // -- total --
  vector<int> Ms_sigu; Ms_sigu += 0;
  vector<int> Ms_piu; Ms_piu += 1,-1;
  vector<int> Ms; Ms += -1,0,+1;
  
  double w = w_list[iw];
  dcomplex k = sqrt(2.0*(w + E0));
  MultArray<dcomplex, 2> Alm(100);   Alm.SetRange(1,3,  -1,1);
  MultArray<dcomplex, 2> Alm_v(100); Alm_v.SetRange(1,3,  -1,1);
  dcomplex c = (sqrt(k/2.0) * dcomplex(0,2) / sqrt(k) * sqrt(3.0/(4.0*M_PI))
		* sqrt(1.0*ne));
  BOOST_FOREACH(int L, Ls) {
    dcomplex ii(0, 1);    
    BOOST_FOREACH(int M, Ms) {
      dcomplex sign = impsi0_chi(L,M)/abs(impsi0_chi(L,M));
      dcomplex etal = exp(ii*CoulombShift(-Z/k, L));
      dcomplex coef_L = sign * w * ii * sqrt(2.0/3.0) * pow(ii, -L) * etal * c; 
      dcomplex psi0_other, psi0_other_v;
      if(calc_term == ECalcTerm_One) {
	psi0_other   = impsi0_muphi(L,M);
	psi0_other_v = impsi0_muphi_v(L,M);
      } else if(calc_term == ECalcTerm_Full) {
	psi0_other   = impsi0_muphi(L,M)   + impsi0_v_psi1(L,M);
	psi0_other_v = impsi0_muphi_v(L,M) + impsi0_v_psi1_v(L,M);
      }
      Alm(L, M)   = coef_L * psi0_other   / sqrt(sign * impsi0_chi(L,M));
      Alm_v(L, M) = coef_L * psi0_other_v / sqrt(sign * impsi0_chi(L,M));
    }
  }

  // Length form, molecular axis averaged
  dcomplex A00 = CoefAzeta0(w, Alm, 0, Ls, Ms);
  dcomplex A20 = CoefAzeta0(w, Alm, 2, Ls, Ms);
  dcomplex A00_sigu = CoefAzeta0(w, Alm, 0, Ls, Ms_sigu);
  dcomplex A00_piu = CoefAzeta0(w, Alm, 0, Ls, Ms_piu);
  double cs   =    (4.0 * M_PI * A00      * au2mb).real(); 
  double cs_sigu = (4.0 * M_PI * A00_sigu * au2mb).real();
  double cs_piu =  (4.0 * M_PI * A00_piu  * au2mb).real();
  double beta = (A20/A00).real();    
  result(iw, idx_cs) = cs;
  result(iw, idx_beta) = beta;
  result(iw, idx_cs_sigu) = cs_sigu;
  result(iw, idx_cs_piu) = cs_piu;

  cout << "coefficient A_lm\n";
  BOOST_FOREACH(int L, Ls) {
    BOOST_FOREACH(int M, Ms) {
      if(L==0 && M!=0) 
	break;
      cout << format("%d, %d, %20.10f, %20.10f\n") % L % M
	% Alm(L, M).real() % Alm(L, M).imag();
    }
  }

  // Velocity form, molecular axis averaged
  dcomplex A00_v = CoefAzeta0(w, Alm_v, 0, Ls, Ms);
  dcomplex A20_v = CoefAzeta0(w, Alm_v, 2, Ls, Ms);
  dcomplex A00_piu_v  = CoefAzeta0(w, Alm_v, 0, Ls, Ms_piu);
  dcomplex A00_sigu_v = CoefAzeta0(w, Alm_v, 0, Ls, Ms_sigu);
  double cs_v      = (4 * M_PI * A00_v     /(w*w) * au2mb).real();
  double cs_piu_v  = (4 * M_PI * A00_piu_v /(w*w) * au2mb).real();
  double cs_sigu_v = (4 * M_PI * A00_sigu_v/(w*w) * au2mb).real();
  double beta_v    = (A20_v/A00_v).real();
  result(iw, idx_cs_v) = cs_v;
  result(iw, idx_beta_v) = beta_v; 
  result(iw, idx_cs_sigu_v) = cs_sigu_v;
  result(iw, idx_cs_piu_v)  = cs_piu_v;
  //  cout << "cross sections (length form, velocity form)" << endl;
  cout << format("Cs(total): %20.10f, %20.10f\n") % cs % cs_v;
  cout << format("Cs(sig_u): %20.10f, %20.10f\n") % cs_sigu % cs_sigu_v;
  cout << format("Cs(pi_u) : %20.10f, %20.10f\n") % cs_piu % cs_piu_v;
  cout << format("beta     : %20.10f, %20.10f\n") % beta % beta_v;

  // fixed molecule
  double coeff =  4.0*M_PI*M_PI / (c_light*w) * 1.0/(4.0*M_PI) * au2mb;
  cout << format("coef: %20.10f\n") % coeff;
  int lmax = 2*Ls[Ls.size()-1];
  if(f_para.is_open()) {
    f_para << "w";
    for(int L = 0; L <= lmax; L++) {
      f_para << format(",%d") % L;
    }
    f_para << ",cs_sigu\n";
  }

  for(int M = 0; M <= 1; M++) {
    cout << format("Coef(M=%d) = \n") % M;
    for(int L = 0; L <= lmax; L++) {      
      double c_l = Coef_FixedMole(w, Alm,   Ls, L, M)      *coeff;
      double c_v = Coef_FixedMole(w, Alm_v, Ls, L, M)/(w*w)*coeff;
      cout << format("%d: %20.10f, %20.10f\n") % L % c_l % c_v; 
      if(f_para.is_open()) {
	f_para << format(",%20.10f") % c;
      }
    }
    
    double c0_l = Coef_FixedMole(w, Alm,   Ls, 0, M)       * coeff * 4.0*M_PI/3.0;
    double c0_v = Coef_FixedMole(w, Alm_v, Ls, 0, M)/(w*w) * coeff * 4.0*M_PI/3.0;
    if(M==1) {
      c0_l *= 2.0; c0_v *= 2.0;
    }
    cout << format("Cs(M=%d,fixed): %20.10f, %20.10f\n") % M % c0_l % c0_v;
    if(f_para.is_open()) {
      f_para << format(",%20.10f\n") % c0;
    }
  }
}
void PrintOut() {

  int num(w_list.size());
  
  ofstream f(cs_csv.c_str(), ios::out);
  f << "w,ene,cs,cs_sigu,cs_piu,beta,cs_alpha,cs_sigu_alpha,cs_piu_alpha,cs_v,cs_sigu_v,cs_piu_v,beta_v,cs_alpha_v,cs_sigu_alpha_v,cs_piu_alpha_v" << endl;
  for(int i = 0; i < num; i++) {
    f << w_list[i] << "," << w_list[i] + E0.real();
    for(int idx = 0; idx < num_header; idx++) {
      f << "," << result(i, idx);
    }
    f << endl;
  }
  f.close();

  
  if(write_psi_filename != "") {
    VectorXd& rs = write_psi_rs;
    int irr = sym->irrep_z();
    ofstream f(write_psi_filename.c_str(), ios::out);
    f << "r";
    for(int L = 0; L <= write_psi_lmax; L++) {
      f << format(",re_%d,im_%d") % L % L;
    }
    f << endl;
    
    vector<VectorXcd> ys(write_psi_lmax+1), dys(write_psi_lmax+1);
    for(int L = 0; L <= write_psi_lmax; L++) {
      basis1->AtR_Ylm(L, 0, irr, cZ1(irr), rs.cast<dcomplex>(), &ys[L], &dys[L]);
    }
    
    for(int i = 0; i < rs.size(); i++) {
      f << format("%f") % rs[i];
      for(int L = 0; L <= write_psi_lmax; L++) {
	dcomplex y = ys[L][i];
	f << format(",%f,%f") % y.real() % y.imag();
      }
      f << endl;
    }
  }
}
int main(int argc, char *argv[]) {
  cout << ">>>> two_pot >>>>" << endl;
  if(argc == 1) {
    cerr << "need one argument" << endl;
    exit(1);
  }
  in_json = argv[1];
  Parse();
  PrintIn();
  CalcMat();
  if(calc_type == "STEX") {
    CalcMatSTEX();
  } else if(calc_type == "RPA") {
    try {
      CalcMatRPA();
    } catch(exception& e) {      
      cout << e.what() << endl;
      cout << "error on CalcMatRPA\n";
      return 1;
    }
  }
  PrintTimeStamp("Calc", NULL);
  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    double w = w_list[iw];
    cout << "w_eV: " << w * au2ev << endl;
    cout << "w_au: " << w << endl;
    cout << "E_eV: " << (w + E0) * au2ev << endl;
    cout << "E_au: " << w + E0 << endl;
    cout << "k_au: " << sqrt(2.0*(w + E0)) << endl;
    CalcDriv(iw);
    CalcBraket();
    CalcMain_alpha(iw);
    CalcMain(iw);
  }
  PrintOut();
  cout << "<<<< two_pot <<<<" << endl;
  return 0;
}


