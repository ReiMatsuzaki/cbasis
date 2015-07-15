# Overview

This library support numerical calculations using l2 basis
functions. Now only one dimendional Slater Type Orbital(STO) 
and Gauss Type Orbitals are supported. Complex version
of these functions are also supported.

# Source files

- macros.hpp : simple macro is defined
- fact.hpp : calculation of factorial
- erfc.hpp : evaluations of real and complex error functions
- prim.hpp : STO and GTO
- lcomb.hpp: linear combinations of primitive basis
- op.hpp: linear operator object
- hatom.hpp: hydrogen atom eigensolutions


# todo 

 - refactoring for linear operators
 - python bindings only LinearComb?
 - consider for create new operator from existing one (e.g. Hamiltonian of Hydrogen atom)




