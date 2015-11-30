# Overview

This library support numerical calculations using l2 basis
functions. Now only one dimendional Slater Type Orbital(STO) 
and Gauss Type Orbitals are supported. Complex version
of these functions are also supported.

# Source files
## Utilities

	- macros.hpp : simple macro is defined


## Math

	- fact.hpp : calculation of factorial
	- erfc.hpp : evaluations of real and complex error functions
	- lgamma.hpp : Lower Gamma function


## Linear Space

	-func.hpp
	
	- exp_func.hpp : STO and GTO
	- cut_exp.hpp : finite ranged STO/GTO
	- lin_func.hpp : linear combination of same type functions

	- op.hpp: linear operator object	

	- cip.hpp  : complex symmetric inner product


## hydrogen

	- hatom.hpp : hydrogen atom eigen solutions
	- hatom_h.hpp: header version of hatom.hpp


# todo 

 - refactoring for linear operators
 - python bindings only LinearComb?
 - consider for create new operator from existing one (e.g. Hamiltonian of Hydrogen atom)




