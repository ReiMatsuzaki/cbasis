#include <fstream>
#include "../utils/macros.hpp"
#include "read_json.hpp"

using namespace std;
using namespace picojson;
using namespace Eigen;

namespace cbasis {

  template<> void CheckValue<object>(picojson::value& val, int n, string l) {
    
    if(not val.is<object>()) {
      throw(runtime_error("value is not string. " + l));
    }
    if(n > 0) {
      if((int)val.get<object>().size() != n)
	throw(runtime_error("invalid size object. " + l));
    }
  }
  template<> void CheckValue<string>(picojson::value& val, int n, string l) {
    if(not val.is<string>()) {
      throw(runtime_error("value is not string. " + l));
    }    
  }
  template<> void CheckValue<double>(picojson::value& val, int n, string l) {
    if(not val.is<double>()) {
      throw(runtime_error("value is not string. " + l));
    }    
  }
  template<> void CheckValue<dcomplex>(picojson::value& val, int n, string l) {
    string msg = "value  must be complex (double or two element array of double)." +l;
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
  template<> void CheckValue<int>(picojson::value& val, int n, string l) {
    string msg = "value  must be int. " +l;
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
  template<> void CheckValue<array>(picojson::value& val, int n, string l) {
    if(not val.is<array>()) {
      throw(runtime_error("value  is not array. " + l));
    }

    if(n > 0) {
      if((int)val.get<array>().size() != n) {
	throw(runtime_error("invalid size of array. "+l));
      }
    }
  }
  template<> void CheckValue<VectorXcd>(picojson::value& val, int n, string l) {

    try {
      CheckValue<array>(val, n, l);
    } catch(...) {
      throw(runtime_error("val must be array for VectorXcd. "+l));
    }
    array& ary = val.get<array>();
    for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
      try {
	CheckValue<dcomplex>(*it);
      } catch(...) {
	throw(runtime_error("element must be dcomplex for VectorXcd. "+l));
      }
    }
    
  }
  template<> void CheckValue<VectorXi>(picojson::value& val, int n, string l) {

    try {
      CheckValue<array>(val, n, l);
    } catch(...) {
      throw(runtime_error("val must be array for VectorXcd. "+l));
    }
    array& ary = val.get<array>();
    for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
      try {
	CheckValue<int>(*it);
      } catch(...) {
	throw(runtime_error("element must be dcomplex for VectorXcd. "+l));
      }
    }
    
  }
  template<> void CheckValue<MatrixXi>(picojson::value& val, int n, string l) {
    try {
      CheckValue<array>(val, -1);
    } catch(exception& e) {
      string msg = e.what();
      msg += "\nval must be array for MatrixXi. "+l;
      throw(runtime_error(msg));
    }
    array& ary = val.get<array>();
    for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
      try {
	CheckValue<array>(*it, n);
      } catch(...) {
	throw(runtime_error("element of val must be array for MatrixXi. "+l));
      }
      array& aryary = it->get<array>();
      if(n > 0) {
	if(n != (int)aryary.size()) {
	  throw(runtime_error("invalid size of element for MatrixXi. " + l));
	}
      }
      for(array::iterator jt = aryary.begin(); jt != aryary.end(); ++jt) {
	try {
	  CheckValue<int>(*jt);
	} catch(...) {
	  throw(runtime_error("element of element of val must be int for MatrixXi" + l));
	}
      }
    }
        
  }
  template<> void CheckValue<MatrixXcd>(picojson::value& val, int n, string l) {
    try {
      CheckValue<array>(val, -1);
    } catch(exception& e) {
      string msg = e.what();
      msg += "\nval must be array for MatrixXcd. "+l;
      throw(runtime_error(msg));
    }
    array& ary = val.get<array>();
    for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
      try {
	CheckValue<array>(*it, n);
      } catch(...) {
	throw(runtime_error("element of val must be array for MatrixXi. "+l));
      }
      array& aryary = it->get<array>();
      if(n > 0) {
	if(n != (int)aryary.size()) {
	  throw(runtime_error("invalid size of element for MatrixXi. " + l));
	}
      }
      for(array::iterator jt = aryary.begin(); jt != aryary.end(); ++jt) {
	try {
	  CheckValue<dcomplex>(*jt);
	} catch(...) {
	  throw(runtime_error("element of element of val must be int for MatrixXi" + l));
	}
      }
    }  
  }  
  template<> dcomplex ReadJson<dcomplex>(value& json) {
    dcomplex x(0);
    if(json.is<double>()) {
      x = json.get<double>();
    } else if(json.is<array>()) {
      array& ary = json.get<array>();
      x = dcomplex(ary[0].get<double>(), ary[1].get<double>());
    }   
    return x;
  }  
  template<> VectorXcd ReadJson<VectorXcd>(value& json) {
    
    array& ary = json.get<array>();
    VectorXcd vec = VectorXcd::Zero(ary.size());

    for(int i = 0; i <(int) ary.size(); i++) {
      vec.coeffRef(i) = ReadJson<dcomplex>(ary[i]);
    }
    return vec;
  }
  template<> MatrixXcd ReadJson<MatrixXcd>(value& json) {

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
  template<> VectorXi ReadJson<VectorXi>(value& json) {
    
    array& ary = json.get<array>();
    VectorXi vec = VectorXi::Zero(ary.size());

    for(int i = 0; i < (int)ary.size(); i++) {
      int x(ary[i].get<double>());
      vec.coeffRef(i) = x;
    }
    
    return vec;
  }  
  template<> MatrixXi ReadJson<MatrixXi>(value& json) {

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
  void CheckObject(picojson::object& obj, std::string k, int n, string l) {
    if(obj.find(k) == obj.end()) {
      throw(runtime_error("key \"" + k + "\" not found"));
    }
    CheckValue<T>(obj[k], n, l);
  }
    
  template<> pSymmetryGroup ReadJson<pSymmetryGroup>( value& v) {
    if(not v.is<string>())
      throw(runtime_error("sym is string"));
    string sym = v.get<string>();
    
    if(sym == "C1")
      return SymmetryGroup::C1();
    else if(sym == "Cs")
      return SymmetryGroup::Cs();
    else if(sym == "C2h")
      return SymmetryGroup::Cs();
    else if(sym == "D2h")
      return SymmetryGroup::D2h();
    else
      throw(runtime_error("sym must be C1,Cs,C2h or D2h"));
  }  
  template<> Molecule ReadJson<Molecule>(value& v) {

    Molecule mole(new _Molecule());
    if(v.is<array>()) {

      const array& atoms = v.get<array>();
      for(array::const_iterator it = atoms.begin();
	  it != atoms.end(); ++it) {

	// -- each object --
	if(not it->is<object>())
	  throw(runtime_error("element of molecule object must be object"));
	object atom = it->get<object>();

	// -- xyz --
	CheckObject<string>(atom, "name", 1, "molecule");
	string name = atom["name"].get<string>();

	CheckObject<VectorXcd>(atom, "xyz", 3, "molecule");
	VectorXcd xyz = ReadJson<VectorXcd>(atom["xyz"]);

	CheckObject<dcomplex>(atom, "q", 1, "molecule");
	dcomplex q = ReadJson<dcomplex>(atom["q"]);

	// -- build --
	mole->Add(name, xyz, q);
      }
    }
    return mole;
  }
  template<> Reduction ReadJson<Reduction>(picojson::value& json) {
    
    if(not json.is<object>()) 
      throw(runtime_error("Reduction must be object"));

    object& obj = json.get<object>();
    
    if(obj.find("irrep") == obj.end())
      throw(runtime_error("irrep not found in Reduction"));

    if(obj.find("coef") == obj.end())
      throw(runtime_error("coef not found in Reduction"));

    int irrep = (int)obj["irrep"].get<double>();
    MatrixXcd coef = ReadJson<MatrixXcd>(obj["coef"]);
    
    Reduction rds(irrep, coef);
    return rds;
    
  }
  void ReadJson_SymGTOs_Subs_cart( object basis, SymGTOs gtos) {

    string l("SymGTOs_cart");
    
    // -- pn --
    CheckObject<VectorXi>(basis, "ns", 3, l+" ns");
    VectorXi _ns =  ReadJson<VectorXi>(basis["ns"]);
    Vector3i ns(_ns);
    
    // -- xyz --
    CheckObject<VectorXcd>(basis, "xyz", 3, l + " xyz");
    VectorXcd _xyz = ReadJson<VectorXcd>(basis["xyz"]);
    Vector3cd xyz(_xyz);

    // -- zeta --
    CheckObject<VectorXcd>(basis, "zeta", -1, l + " zeta");
    VectorXcd zeta =  ReadJson<VectorXcd>(basis["zeta"]);

    // -- build --    
    SubSymGTOs sub(Sub_mono(0, xyz, ns, zeta));
    gtos->AddSub(sub);

  }
  void ReadJson_SymGTOs_Subs_full( object basis, SymGTOs gtos) {
    
    string l("SymGTOs_Subs_full");

    SubSymGTOs sub;
    
    CheckObject<array>(basis, "rds", -1, l);
    array& rds_list = basis["rds"].get<array>();
    for(array::iterator it = rds_list.begin(); it != rds_list.end(); ++it) {
      Reduction rds = ReadJson<Reduction>(*it);
      sub.AddRds(rds);
    }

    CheckObject<MatrixXi>(basis, "ns", 3, "SymGTOs_Subs_full");
    MatrixXi ns = ReadJson<MatrixXi>(basis["ns"]);
    for(int i = 0; i < ns.rows(); i++) {
      sub.AddNs(ns(i, 0), ns(i, 1), ns(i, 2));
    }

    CheckObject<MatrixXcd>(basis, "xyz", 3, "SymGTOs_Subs_full");
    MatrixXcd xyz = ReadJson<MatrixXcd>(basis["xyz"]);
    for(int i = 0; i < ns.rows(); i++) {
      sub.AddXyz(xyz(i, 0), xyz(i, 1), xyz(i, 2));
    }

    CheckObject<VectorXcd>(basis, "zeta", -1, "SymGTOs_Subs_full");
    VectorXcd zeta = ReadJson<VectorXcd>(basis["zeta"]);
    sub.AddZeta(zeta);
    
    gtos->AddSub(sub);
  }
  void ReadJson_SymGTOs_Subs(picojson::value& json, SymGTOs gtos) {

    if(not json.is<array>())
      throw runtime_error("json is not array");
    const array& ary = json.get<array>();
    
    for(array::const_iterator it = ary.begin(); it != ary.end(); ++it) {

      // -- check object --
      if(not it->is<object>())
	throw(runtime_error("element of basis must be object"));
      object basis = it->get<object>();

      // -- type --
      CheckObject<string>(basis, "type");
      string type = basis["type"].get<string>();

      if(type == "cart")
	ReadJson_SymGTOs_Subs_cart(basis, gtos);
      else if(type == "full")
	ReadJson_SymGTOs_Subs_full(basis, gtos);
      else 
	throw(runtime_error("now only type==cart is supported"));
    }
  }


}
