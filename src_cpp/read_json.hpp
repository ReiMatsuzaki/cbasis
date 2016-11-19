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
  void ReadJson_SymGTOs_Subs(picojson::value& json, SymGTOs gtos);

  template<class T>
  T ReadJson(picojson::value& json);

}

