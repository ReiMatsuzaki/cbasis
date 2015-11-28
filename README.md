# Overview

This library support numerical calculations using l2 basis
functions. Now only one dimendional Slater Type Orbital(STO) 
and Gauss Type Orbitals are supported. Complex version
of these functions are also supported.

# Source files

- macros.hpp : simple macro is defined
- fact.hpp : calculation of factorial
- erfc.hpp : evaluations of real and complex error functions
- exp_func.hpp : STO and GTO
- cut_exp.hpp : finite ranged STO/GTO
- cip.hpp  : complex symmetric inner product
- lin_func.hpp : linear combination of same type functions
- prim.hpp : to be delete
- lcomb.hpp: to be delete
- op.hpp: linear operator object
- hatom.hpp: hydrogen atom eigen solutions

# todo 

 - refactoring for linear operators
 - python bindings only LinearComb?
 - consider for create new operator from existing one (e.g. Hamiltonian of Hydrogen atom)




