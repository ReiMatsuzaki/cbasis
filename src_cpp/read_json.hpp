#include "symmolint.hpp"
#include "../external/picojson/picojson.h"

namespace cbasis {

  template<class T>
  void CheckValue(picojson::value& val, int n=-1, int m=-1);
  
  template<class T>
  void CheckObject(picojson::object& obj, std::string k, int n=-1, int m=-1);
  
  void CheckJson_string(picojson::object& obj, std::string k);
  void CheckJson_double(picojson::object& obj, std::string k);
  void CheckJson_complex(picojson::object& obj, std::string k);
  void CheckJson_array(picojson::object& obj, std::string k);
  void CheckJson_VectorXcd(picojson::value& json);


  template<class T>
  T ReadJson(picojson::value& json, int n=-1, int m=-1);
  template<class T>
  T ReadJson(picojson::object& obj, std::string, int n=-1, int m=-1);
  template<class T>
  T ReadJsonWithDefault(picojson::object& obj, std::string,
			T t, int n=-1, int m=-1);  

  template<class T>
  picojson::value ToJson(T& t);

  Reduction ReadJson_Reduction(picojson::value& json, SymmetryGroup sym);
  void ReadJson_Molecule(picojson::value& json, Molecule mole);
  void ReadJson_Molecule(picojson::object& obj, std::string k, Molecule mole);
  void ReadJson_SymGTOs_Subs(picojson::value& json, SymGTOs gtos);
  void ReadJson_SymGTOs_Subs(picojson::object& obj, std::string k, SymGTOs gtos);

  void ReadJson_Orbital_file(picojson::object& initstate,
			     SymmetryGroup _sym, dcomplex *_E0,
			     Eigen::VectorXcd *_c0, Irrep *_irrep0, int *_i0);
  void ReadJson_Orbital(picojson::object& obj, string k, SymmetryGroup _sym,
			dcomplex *_E0, Eigen::VectorXcd *_c0,
			Irrep *_irrep0, int *_i0);
}

