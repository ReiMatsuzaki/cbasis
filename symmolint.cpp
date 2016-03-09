#include <iostream>
#include "symmolint.hpp"

namespace l2func {

  using namespace std;
  using namespace Eigen;

  void SymGTOs::AddSymmetry(Symmetry sym) {
    symmetries.push_back(sym);
  }
  /*
  void SymGTOs::AddShell(const Shell& shell) {

    shells_.push_back(shell);
    if(shells_.size() == 1)
      shells_[0].offset = 0;
    else {
      Shell& prev = shells_[shells_.size()-2];
      shells_[shells_.size()-1].offset = prev.offset + prev.size_contractions();
    }
  }
  */
  void SymGTOs::loop() {

    typedef vector<SubSymGTOs>::iterator SubIt;
    typedef vector<Contraction>::iterator ContIt;
    
    MatrixXcd ss(5, 5);

    for(SubIt isub = sub_sym_gtos.begin(); isub != sub_sym_gtos.end(); ++isub) {
      for(SubIt jsub = sub_sym_gtos.begin(); jsub != sub_sym_gtos.end(); ++jsub) {
	
	// -- search symmetry for isub and jsub --


	// loop over each zeta
	for(int iz = 0; iz < isub->size_zeta(); iz++) {
	  for(int jz = 0; jz < jsub->size_zeta(); jz++) {

	    // -- primitive basis --
	    MatrixXcd s(isub->size_at() * isub->size_pn(),
			jsub->size_at() * jsub->size_pn());
	    int iprim(0);
	    for(int iat = 0; iat < isub->size_at(); iat++) {
	      for(int ipn = 0; ipn < isub->size_pn(); ipn++) {
		int jprim(0);
		for(int jat = 0; jat < jsub->size_at(); jat++) { 
		  for(int jpn = 0; jpn < jsub->size_pn(); jpn++) {
		    //s(iprim, jprim) = 1.0;
		    s(0, 0) = 1.0;
		    jprim++;
		  }
		}
		iprim++;
	      }
	    }

	    // -- contractions --
	    for(int icont = 0; icont < isub->size_cont(); icont++) {
	      for(int jcont = 0; jcont < jsub->size_cont(); jcont++) {
		dcomplex cumsum(0.0);
		int iprim2(0);
		for(int iat = 0; iat < isub->size_at(); iat++) {
		  for(int ipn = 0; ipn < isub->size_pn(); ipn++) {
		    int jprim2(0);
		    for(int jat = 0; jat < jsub->size_at(); jat++) { 
		      for(int jpn = 0; jpn < jsub->size_pn(); jpn++) {
			cout << iat << ipn << jat << jpn << endl;
			dcomplex cc = 
			  isub->contraction_icont[icont].coef_iat_ipn(iat, ipn) *
			  jsub->contraction_icont[jcont].coef_iat_ipn(jat, jpn);
			cumsum += cc*s(iprim2, jprim2);
			jprim2++;
		      }
		    }
		  }
		}
		iprim2++;
		ss(0, 0) = cumsum;
	      }
	    }
	  }
	}
      }
    }
  }
}

