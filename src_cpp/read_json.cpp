#include <fstream>
#include "../utils/macros.hpp"
#include "read_json.hpp"

using namespace std;
using namespace picojson;
using namespace Eigen;

namespace cbasis {

  template<> void CheckValue<object>(picojson::value& val, int n, int m) {
    
    if(not val.is<object>()) {
      throw(runtime_error("value is not object. "));
    }
    if(n > 0) {
      if((int)val.get<object>().size() != n)
	throw(runtime_error("invalid size object."));
    }
  }
  template<> void CheckValue<string>(picojson::value& val, int n, int m) {
    if(not val.is<string>()) {
      throw(runtime_error("value is not string. "));
    }    
  }
  template<> void CheckValue<double>(picojson::value& val, int n, int m) {
    
    if(not val.is<double>()) {
      throw(runtime_error("value is not double."));
    }    
  }
  template<> void CheckValue<dcomplex>(picojson::value& val, int n, int m) {
    string msg = "value is not compelx (double or two element array of double).";
    if(val.is<array>()) {
      array& ary = val.get<array>();
      if(ary.size() != 2) {
	throw(runtime_error(msg));
      }
      if(not ary[0].is<double>())
	throw(runtime_error(msg));
      if(not ary[1].is<double>())
	throw(runtime_error(msg));
    } else if(not val.is<double>()) {
      throw runtime_error(msg);
    }
  }
  template<> void CheckValue<int>(picojson::value& val, int n, int m) {
    string msg = "value is not int. ";
    if(not val.is<double>()) {
      throw runtime_error(msg);
    }
    double x(val.get<double>());
    double y((double)(int)x);
    double eps(0.0000000001);
    if(abs(x-y) > eps) {
      throw runtime_error(msg);
    }
  }
  template<> void CheckValue<array>(picojson::value& val, int n, int m) {
    if(not val.is<array>()) {
      throw(runtime_error("value  is not array. "));
    }

    if(n > 0) {
      if((int)val.get<array>().size() != n) {
	throw(runtime_error("invalid size of array. "));
      }
    }
  }
  template<> void CheckValue<VectorXcd>(picojson::value& val, int n, int m) {

    try {
      CheckValue<array>(val, n);
    } catch(exception& e) {
      throw(runtime_error("value is not array for VectorXcd."));
    }
    array& ary = val.get<array>();
    for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
      try {
	CheckValue<dcomplex>(*it);
      } catch(...) {
	throw(runtime_error("element of value is not dcomplex for VectorXcd."));
      }
    }
    
  }
  template<> void CheckValue<VectorXi>(picojson::value& val, int n, int m) {

    try {
      CheckValue<array>(val, n);
    } catch(...) {
      throw(runtime_error("value is not array for VectorXi."));
    }
    array& ary = val.get<array>();
    for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
      try {
	CheckValue<int>(*it);
      } catch(...) {
	throw(runtime_error("element is not int for VectorXi."));
      }
    }
    
  }
  template<> void CheckValue<MatrixXi>(picojson::value& val, int n, int m) {
    
    try {
      CheckValue<array>(val, -1);
    } catch(exception& e) {
      throw(runtime_error("value is not array for MatrixXi"));
    }
    array& ary = val.get<array>();
    
    if(n > 0) 
      if(m != (int)ary.size()) 
	throw(runtime_error("invalid rows size for MatrixXi"));
	
    for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
      try {
	CheckValue<array>(*it, n);
      } catch(...) {
	throw(runtime_error("element of value is not array for MatrixXi."));
      }
      array& aryary = it->get<array>();
      if(m > 0) {
	if(m != (int)aryary.size()) {
	  throw(runtime_error("invalid cols size of element for MatrixXi."));
	}
      }
      for(array::iterator jt = aryary.begin(); jt != aryary.end(); ++jt) {
	try {
	  CheckValue<int>(*jt);
	} catch(...) {
	  throw(runtime_error("element of element of value is not int for MatrixXi"));
	}
      }
    }
        
  }
  template<> void CheckValue<MatrixXcd>(picojson::value& val, int n, int m) {
    
    try {
      CheckValue<array>(val, -1);
    } catch(exception& e) {
      throw(runtime_error("value is not array for MatrixXcd"));
    }
    array& ary = val.get<array>();
    for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
      try {
	CheckValue<array>(*it, n);
      } catch(...) {
	throw(runtime_error("element of value is not array for MatrixXcd."));
      }
      array& aryary = it->get<array>();
      if(n > 0) {
	if(n != (int)aryary.size()) {
	  throw(runtime_error("invalid size of element for MatrixXcd."));
	}
      }
      for(array::iterator jt = aryary.begin(); jt != aryary.end(); ++jt) {
	try {
	  CheckValue<dcomplex>(*jt);
	} catch(...) {
	  throw(runtime_error("element of element of value is not dcomplex for MatrixXcd"));
	}
      }
    }  
  }
  template<class T>
  void CheckObject(picojson::object& obj, std::string k, int n, int m) {
    
    string key = "key \"" + k + "\"";
    if(obj.find(k) == obj.end()) {
      throw(runtime_error(key + " not found."));
    }
    
    try {
      CheckValue<T>(obj[k], n, m);
    } catch(exception& e) {
      string msg = "error on parsing value of " + key + "\n";
      msg += e.what();
      throw(runtime_error(msg));
    }
  }

  template<> string ReadJson<string>(value& json, int n, int m) {
    CheckValue<string>(json);
    return json.get<string>();
  }
  template<> double ReadJson<double>(value& json, int n, int m) {
    CheckValue<double>(json);
    return json.get<double>();
  }  
  template<> int ReadJson<int>(value& json, int n, int m) {
    CheckValue<int>(json, n, m);
    return (int)json.get<double>();;
  }  
  template<> dcomplex ReadJson<dcomplex>(value& json, int n, int m) {
    CheckValue<dcomplex>(json);
    dcomplex x(0);
    if(json.is<double>()) {
      x = json.get<double>();
    } else if(json.is<array>()) {
      array& ary = json.get<array>();
      x = dcomplex(ary[0].get<double>(), ary[1].get<double>());
    }   
    return x;
  }  
  template<> VectorXcd ReadJson<VectorXcd>(value& json, int n, int m) {

    CheckValue<VectorXcd>(json, n);
    
    array& ary = json.get<array>();
    VectorXcd vec = VectorXcd::Zero(ary.size());

    for(int i = 0; i <(int) ary.size(); i++) {
      vec.coeffRef(i) = ReadJson<dcomplex>(ary[i]);
    }
    return vec;
  }
  template<> MatrixXcd ReadJson<MatrixXcd>(value& json, int _n, int _m) {
    CheckValue<MatrixXcd>(json, _n, _m);
    array& ary = json.get<array>();
    int n(0), m(0);
    n = ary.size();
    m = ary[0].get<array>().size();

    MatrixXcd mat = MatrixXcd::Zero(n, m);

    for(int i = 0; i < n; i++) {
      array& aryary = ary[i].get<array>();
      for(int j = 0; j < m; j++) {
	mat.coeffRef(i, j) = ReadJson<dcomplex>(aryary[j]);
      }
    }
    
    return mat;
  }
  template<> VectorXi ReadJson<VectorXi>(value& json, int n, int m) {
    CheckValue<VectorXi>(json, n);
    array& ary = json.get<array>();
    VectorXi vec = VectorXi::Zero(ary.size());

    for(int i = 0; i < (int)ary.size(); i++) {
      int x(ary[i].get<double>());
      vec.coeffRef(i) = x;
    }
    
    return vec;
  }  
  template<> MatrixXi ReadJson<MatrixXi>(value& json, int _n, int _m) {
    CheckValue<MatrixXi>(json, _n, _m);
    array& ary = json.get<array>();
    int n(0), m(0);
    n = ary.size();
    m = ary[0].get<array>().size();
    MatrixXi mat = MatrixXi::Zero(n, m);

    for(int i = 0; i < n; i++) {
      array& aryary = ary[i].get<array>();
      for(int j = 0; j < m; j++) {
	mat.coeffRef(i, j) = (int)aryary[j].get<double>();
      }
    }
    
    return mat;
  }
  
  template<class T>
  T ReadJson(object& json, string key, int n, int m) {
    if(json.find(key) == json.end()) {
      throw runtime_error("key \"" + key + "\" not found");
    }
    return ReadJson<T>(json[key], n, m);
  }

  template<>
  value ToJson<dcomplex>(dcomplex& x) {
    double eps(0.000000000000001);
    if(abs(x.imag()) < eps ) {
      double re_x(x.real());
      return value(re_x);
    } else {
      array xs(2);
      xs[0] = value(x.real());
      xs[1] = value(x.imag());
      return value(xs);
    }
  }
  template<>
  value ToJson<int>(int& x) {
    double y(x);
    return value(y);
  }  
  
  template<> pSymmetryGroup ReadJson<pSymmetryGroup>( value& v, int n, int m) {

    string sym;
    try {
      sym = ReadJson<string>(v);
    } catch(exception& e) {
      string msg = "error on parsing SymmetryGroup.\n";
      msg += e.what();
      throw(runtime_error(msg));
    }
    
    if(sym == "C1")
      return SymmetryGroup::C1();
    else if(sym == "Cs")
      return SymmetryGroup::Cs();
    else if(sym == "C2h")
      return SymmetryGroup::Cs();
    else if(sym == "D2h")
      return SymmetryGroup::D2h();
    else
      throw(runtime_error("error on parsing SymmetryGroup.\nsym must be C1,Cs,C2h or D2h"));
  }  
  template<> Molecule ReadJson<Molecule>(value& v, int n, int m) {
    
    Molecule mole(new _Molecule());
    try {
      CheckValue<array>(v);
      array& atoms = v.get<array>();
      for(array::iterator it = atoms.begin(); it != atoms.end(); ++it) {

	// -- build --
	CheckValue<object>(*it);
	object& atom = it->get<object>();
	string name = ReadJson<string>(atom, "name");
	VectorXcd xyz = ReadJson<VectorXcd>(atom, "xyz");
	dcomplex q = ReadJson<dcomplex>(atom, "q");
	mole->Add(name, xyz, q);
      }
      
    } catch(exception& e) {
      string msg = "error on parsing Molecule\n";
      msg += e.what();
      throw(runtime_error(msg));
    }
    
    return mole;
  }
  template<> Reduction ReadJson<Reduction>(picojson::value& json, int n, int m) {

    try {
      CheckValue<object>(json);
      object& obj = json.get<object>();
      int irrep = ReadJson<int>(obj, "irrep");
      MatrixXcd coef = ReadJson<MatrixXcd>(obj, "coef");
      Reduction rds(irrep, coef);
      return rds;

    } catch(exception& e) {
      string msg = "error on parsing Reduction.\n";
      msg += e.what();
      throw(runtime_error(msg));      
    }
    
  }
  void ReadJson_SymGTOs_Subs_cart( object basis, SymGTOs gtos) {

    try {
      Molecule mole = gtos->GetMolecule();
      
      // -- pn --
      Vector3i ns(ReadJson<VectorXi>(basis, "ns", -1, 3));

      // -- xyz --
      string atom_name = ReadJson<string>(basis["atom"]);
      
      Vector3cd xyz(ReadJson<VectorXcd>(basis["xyz"], -1, 3));

      // -- zeta --
      VectorXcd zeta =  ReadJson<VectorXcd>(basis["zeta"]);

      // -- build --    
      SubSymGTOs sub(Sub_mono(0, xyz, ns, zeta));
      gtos->AddSub(sub);

    } catch(exception& e) {
      string msg = "error on parsing SymGTOs_Subs_cart.\n";
      msg += e.what();
      throw(runtime_error(msg));      
    }    

  }
  void ReadJson_SymGTOs_Subs_full( object basis, SymGTOs gtos) {
    
    try {
      SubSymGTOs sub;
      
      CheckObject<array>(basis, "rds");
      array& rds_list = basis["rds"].get<array>();
      for(array::iterator it = rds_list.begin(); it != rds_list.end(); ++it) {
	Reduction rds = ReadJson<Reduction>(*it);
	sub.AddRds(rds);
      }

      MatrixXi ns = ReadJson<MatrixXi>(basis["ns"]);
      for(int i = 0; i < ns.rows(); i++) {
	sub.AddNs(ns(i, 0), ns(i, 1), ns(i, 2));
      }

      MatrixXcd xyz = ReadJson<MatrixXcd>(basis["xyz"]);
      for(int i = 0; i < xyz.rows(); i++) {
	sub.AddXyz(xyz(i, 0), xyz(i, 1), xyz(i, 2));
      }

      VectorXcd zeta = ReadJson<VectorXcd>(basis["zeta"]);
      sub.AddZeta(zeta);

      gtos->AddSub(sub);
      
    } catch(exception& e) {
      string msg = "error on parsing SymGTOs_Subs_cart.\n";
      msg += e.what();
      throw(runtime_error(msg));      
    }    
  }
  void ReadJson_SymGTOs_Subs(picojson::value& json, SymGTOs gtos) {

    try {
      CheckValue<array>(json);
      array& ary = json.get<array>();
      
      for(array::iterator it = ary.begin(); it != ary.end(); ++it) {

	CheckValue<object>(*it);
	object& basis = it->get<object>();

	CheckObject<string>(basis, "type");
	string type = ReadJson<string>(basis, "type");

	if(type == "cart")
	  ReadJson_SymGTOs_Subs_cart(basis, gtos);
	else if(type == "full")
	  ReadJson_SymGTOs_Subs_full(basis, gtos);
	else 
	  throw(runtime_error("unsupported type: " + type));
      }
    } catch(exception& e) {
      string msg = "error on parsing SymGTOs_Subs\n";
      msg += e.what();
      throw runtime_error(msg);
    }
  }


}
