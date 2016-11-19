#include "symmolint.hpp"
#include "../external/picojson/picojson.h"

namespace cbasis {

  template<class T>
  void CheckValue(picojson::value& val, int n=-1, std::string label="");
  
  template<class T>
  void CheckObject(picojson::object& obj, std::string k,
		   int n=-1, std::string l="");
  
  void CheckJson_string(picojson::object& obj, std::string k);
  void CheckJson_double(picojson::object& obj, std::string k);
  void CheckJson_complex(picojson::object& obj, std::string k);
  void CheckJson_array(picojson::object& obj, std::string k);
  //void CheckJson_array(picojson::object& json, int num);
  void CheckJson_VectorXcd(picojson::value& json);
  //  void CheckJson_VectorXcd(picojson::value& json, int num);
  //  void CheckJson_MatrixXcd(picojson::value& json);
  
  //  void ReadJson_SymmetryGroup(picojson::value& json, pSymmetryGroup sym);
  //  void ReadJson_Molecule(picojson::value& json, Molecule mole);
  //  void ReadJson_Reduction(picojson::value& json, Reduction *rds);
  void ReadJson_SymGTOs_Subs(picojson::value& json, SymGTOs gtos);

  template<class T>
  T ReadJson(picojson::value& json);
  //  void ReadJson_Complex(picojson::value& json, dcomplex *x);
  //  void ReadJson_VectorXcd(picojson::value& json, Eigen::VectorXcd *vec);
  //  void ReadJson_MatrixXcd(picojson::value& json, Eigen::MatrixXcd *mat);
  //  void ReadJson_VectorXi(picojson::value& json, Eigen::VectorXi *vec);

  
  //  void SymGTOsReadJson(SymGTOs gtos, std::string json);
  //  void SymGTOsReadFile(SymGTOs gtos, std::string filename);
}

