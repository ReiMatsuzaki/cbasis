#include <Eigen/Core>
#include <Eigen/Dense>
#include "mo.hpp"
#include "eigen_plus.hpp"

using namespace Eigen;

namespace l2func {
  _MO::_MO() {}

  typedef pair<Irrep, dcomplex> IrrepComplex;
  struct Compare_IrrepEig {
    bool operator()(const IrrepComplex& a, const IrrepComplex& b) {
      return real(a.second) < real(b.second);
    }
  };
  vector<int> CalcOccNum(const BVec& eigs, int num_sym, int total_occ_orb) {
  /**
     Give occupation number from orbital energies for each symmetry.
     
     Inputs 
     ------
     const BVec& eigs : orbital energies for each symmetries.
     int num_sym      : number of symmetry
     int num_orb      : total number of occupied orbitals
     
     Outputs
     -------
     vector<int> res : number of occupied orbital for each symmetry
   */

    // -- build vector<IrrepIndexEig> --
    vector<IrrepComplex> irrep_eig_list;
    for(BVec::const_iterator it = eigs.begin(); it != eigs.end(); ++it) {

      Irrep irrep = it->first;
      VectorXcd eig = it->second;

      for(int index = 0; index < eig.size(); index++) {
	IrrepComplex ic(irrep, eig[index]);
	irrep_eig_list.push_back(ic);
      }
    }

    // -- sort --
    sort(irrep_eig_list.begin(), irrep_eig_list.end(), Compare_IrrepEig());

    // -- extract results --
    vector<int> res(num_sym, 0);
    for(int i = 0; i < total_occ_orb; i++) {
      res[irrep_eig_list[i].first]++;
    }
    return res;
  }
  vector<Irrep> CalcIrrepList(BMatSet& mat_set) {
    vector<Irrep> irrep_list;
    int num_block(mat_set.block_num());
    for(Irrep irrep = 0; irrep < num_block; irrep++) {
      if(mat_set.Exist("t", irrep, irrep)) {
	irrep_list.push_back(irrep);
      }
    }    
    return irrep_list;
  }
  MO CalcOneEle(pSymmetryGroup sym, BMatSet& mat_set, int) {

    MO mo(new _MO);

    // ---- get non0 symmetry ----
    mo->irrep_list = CalcIrrepList(mat_set);
    typedef vector<Irrep>::iterator It;
    for(It it = mo->irrep_list.begin(), end = mo->irrep_list.end();
	it != end; ++it) {
      Irrep irrep(*it);
      pair<Irrep, Irrep> ii(irrep, irrep);
      mo->H[ii] = (mat_set.GetMatrix("t", irrep, irrep) +
		   mat_set.GetMatrix("v", irrep, irrep));
      mo->S[ii] = mat_set.GetMatrix("s", irrep, irrep);
      mo->C[ii] = MatrixXcd::Zero(1, 1);
      mo->eigs[irrep] = VectorXcd::Zero(1);
      generalizedComplexEigenSolve(mo->H[ii], mo->S[ii], &mo->C[ii],
				   &mo->eigs[irrep]);
    }

    // ---- other value ----
    mo->sym = sym;
    mo->num_occ_irrep = CalcOccNum(mo->eigs, sym->num_class(), 1);
    for(It it = mo->irrep_list.begin(), end = mo->irrep_list.end();
	it != end; ++it) {
      if(mo->num_occ_irrep[*it] == 1)
	mo->energy = mo->eigs[*it][0];
    }
    return mo;
  }
  MO CalcRHF(SymGTOs& gtos, int nele, int max_iter, double eps,
	     bool *is_conv, int debug_lvl) {
    int num_basis(gtos.size_basis());
    BMatSet mat_set; 
    IB2EInt *eri = new B2EIntMem(pow(num_basis, 4));

    gtos.CalcMat(&mat_set);
    gtos.CalcERI(eri, 2);

    return CalcRHF(gtos.sym_group, mat_set, eri, nele, max_iter,
		   eps, is_conv, debug_lvl);

  }
  MO CalcRHF(pSymmetryGroup sym, BMatSet& mat_set, IB2EInt* eri,
	     int nele, int max_iter, double eps, bool *is_conv, int debug_lvl) {
    
    if(nele == 1) {
      *is_conv = true;
      return CalcOneEle(sym, mat_set, debug_lvl);
    }

    if(nele % 2 == 1) {
      string msg; SUB_LOCATION(msg); msg += "nele must be even integer.";
      throw runtime_error(msg);
    }
    
    int nocc = nele/2;
    *is_conv = false;
    MO mo(new _MO);
    mo->sym = sym;

    // ---- get non0 symmetry ----
    mo->irrep_list = CalcIrrepList(mat_set);

    // ---- initilize ----
    BMat FOld;    
    vector<int> num_irrep(sym->num_class(), 0);
    typedef vector<Irrep>::iterator It;
    for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it ) {
      Irrep irrep = *it;      
      pair<Irrep, Irrep> ii(irrep, irrep);
      mo->H[ii] = (mat_set.GetMatrix("t", irrep, irrep) +
		   mat_set.GetMatrix("v", irrep, irrep));
      mo->S[ii] = mat_set.GetMatrix("s", irrep, irrep);
      mo->F[ii] = mo->H[ii];
      int n(mo->H[ii].rows());
      num_irrep[irrep] = n;
      mo->C[ii] = MatrixXcd::Zero(n, n);
      mo->P[ii] = MatrixXcd::Zero(n, n);
      FOld[ii] = MatrixXcd::Zero( n, n);
      mo->eigs[*it] = VectorXcd::Zero(n);
    }

    // ---- SCF calculation ----
    for(int iter = 0; iter < max_iter; iter++) {
      
      // -- solve --
      for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	pair<Irrep, Irrep> ii(make_pair(*it, *it));
	generalizedComplexEigenSolve(mo->F[ii], mo->S[ii],
				     &mo->C[ii], &mo->eigs[*it]);
      }

      // -- number of occupied orbitals --
      mo->num_occ_irrep = CalcOccNum(mo->eigs, sym->num_class(), nocc);      

      // -- density matrix --
      for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	pair<Irrep, Irrep> ii(make_pair(*it, *it));
	int n_irrep(num_irrep[*it]);
	MatrixXcd& P_ii = mo->P[ii];
	MatrixXcd& C_ii = mo->C[ii];
	for(int i = 0; i < mo->num_occ_irrep[*it]; i++) {
	  for(int k = 0; k < n_irrep; k++)
	    for(int l = 0; l < n_irrep; l++)
	      P_ii(k, l) = 2.0 * C_ii(k, i) * C_ii(l, i);
	}
      }
      
      // -- update Fock matrix --
      for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	Irrep irrep(*it);
	pair<Irrep, Irrep> ii(make_pair(irrep, irrep));
	FOld[ii].swap(mo->F[ii]);
	mo->F[ii] = mo->H[ii];
      }
      int ib,jb,kb,lb,i,j,k,l,t;
      dcomplex v;
      eri->Reset();
      while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
	if(ib == jb && kb == lb) {
	  pair<Irrep, Irrep> ii(make_pair(ib, jb));
	  pair<Irrep, Irrep> kk(make_pair(kb, lb));
	  mo->F[ii](i, j) += mo->P[kk](k, l) * v;
	}
	if(ib == lb && kb == jb) {
	  pair<Irrep, Irrep> il(make_pair(ib, lb));
	  pair<Irrep, Irrep> jk(make_pair(jb, kb));
	  mo->F[il](i, l) -= 0.5 * mo->P[jk](k, j) * v;
	}
      }

      // -- Convergence check --
      bool conv(true);
      for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	pair<Irrep, Irrep> ii(make_pair(*it, *it));
	if(abs((FOld[ii] - mo->F[ii]).mean()) > eps) {
	  conv = false;
	}
      }
      FOld = mo->F; // copy

      // -- calculate total energy --
      mo->energy = 0.0;
      for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	pair<Irrep, Irrep> ii(make_pair(*it, *it));
	MatrixXcd& P(mo->P[ii]);
	MatrixXcd& H(mo->H[ii]);
	MatrixXcd& F(mo->F[ii]);
	int n(num_irrep[*it]);
	for(int i = 0; i < n; i++)
	  for(int j = 0; j < n; j++) 
	    mo->energy += 0.5 * P(i, j) * (H(i, j) + F(i, j));
      }
      
      // -- print current status --
      if(debug_lvl > 0) {
	for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	  Irrep irrep = *it;
	  cout << irrep << ": " << mo->num_occ_irrep[irrep];
	  for(int i = 0; i < mo->num_occ_irrep[irrep]; i++) {
	    cout << mo->eigs[irrep][i] << ", ";
	  }
	  cout << endl;	
	}
	cout << endl;       
      }

      // -- if convergence, break loop. --
      if(conv) {
	*is_conv = true;
	break;
      }
    }
    return mo;
  }
  void CalcSEHamiltonian(MO mo, IB2EInt* eri, Irrep I0, int i0, BMat* h_stex) {

    typedef vector<Irrep>::iterator It;
    BMat res;
    
    // set res as Fock matrix
    for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
      Irrep irrep(*it);
      pair<Irrep, Irrep> ii(make_pair(irrep, irrep));
      res[ii] = mo->H[ii];
    }

    // loop ERI and add to J and K
    int ib,jb,kb,lb,i,j,k,l,t;
    dcomplex v;
    eri->Reset();
    while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
      if(ib == jb && kb == I0 && lb == I0) {
	// Add J
	pair<Irrep, Irrep> ij(ib, jb), kl(I0, I0);
	MatrixXcd& C = mo->C[kl];
	res[ij](i, j) += C(k, i0) * C(l, i0) * v;
      }
      
      if(ib == lb && jb == I0 && kb == I0) {
	// Add K
	pair<Irrep, Irrep> il(ib, lb), jk(I0, I0);
	MatrixXcd& C = mo->C[jk];
	res[il](i, l) += C(j, i0) * C(k, i0) * v;
      }
    }  
    h_stex->swap(res);
  }
  dcomplex CalcAlpha(MO mo, BMatSet& mat_set, Irrep I0, int i0, BMat& h_stex,
		     double w) {
    dcomplex a(0);
    typedef vector<Irrep>::const_iterator It;
    dcomplex ene = w + mo->eigs[I0](i0);
    for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {      
      Irrep irrep = *it;
      if(mat_set.Exist("z", irrep, I0)) {
	pair<Irrep, Irrep> ii(irrep, irrep);
	const MatrixXcd& S = mat_set.GetMatrix("s", irrep, irrep);
	const MatrixXcd& H = h_stex[ii];
	const MatrixXcd& Z = mat_set.GetMatrix("z", irrep, I0);
	MatrixXcd L = S * ene - H;
	VectorXcd m = Z * mo->C[make_pair(I0, I0)].col(i0);
	VectorXcd c = L.colPivHouseholderQr().solve(m);
	int num(m.size());
	for(int i = 0; i < num; i++) {
	  a += m[i] * c[i];
	}
      }
    }
    return a;
  }
  double PITotalCrossSection(dcomplex alpha, double w, int num_occ_ele) {
    double au2mb(pow(5.291772, 2));
    double c(137.0);
    double c0 = -4.0 * M_PI * w / c * imag(alpha) * au2mb;
    if(num_occ_ele == 1)
      return c0;
    else
      return 2.0 * c0;
  }
}
