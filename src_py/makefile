CXXFLAGS=${INC_PATH} -Wall -O3  
#CXXFLAGS=${INC_PATH} -Wall
#CXXFLAGS=${INC_PATH} -Wall -g
#CXXFLAGS=${INC_PATH} -Wall -pg -g -fno-inline
PYBIND=

#OBJS= ../cip_exp.o ../cut_exp.o ../exp_func.o  ../delta.o ../fact.o ../erfc.o ../math_utils.o
#OBJS3D = ../angmoment.o ../cints.o ../gto3dset.o
#OBJMOLINT = ../molint.o ../angmoment.o ../math_utils.o

# ==== L2(R1) ====
l2func_bind.so: l2func_bind.cpp ../fact.o ../erfc.o ../lgamma.o ../math_utils.o ../exp_func.o ../cut_exp.o ../cip_exp.o ../delta.o
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ $^ ${OBJS} ${OBJS3D} ${CXXFLAGS} -lboost_python -lboost_numpy -lgsl -lblas -framework Python

.PHONY: check
check: l2func_bind.so
	python test.py


# ==== SymMolInt ====
symmolint_bind.so: symmolint_bind.cpp ../symmolint.o ../fact.o ../mol_func.o ../bmatset.o ../angmoment.o ../eigen_plus.o 
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ $^  ${CXXFLAGS} -lboost_python -lboost_numpy -lgsl -lblas -framework Python
check_symmolint_bind: symmolint_bind.so
	python test_symmolint.py

# ==== R1GTOs ====
r1gtoint_bind.so: r1gtoint_bind.cpp ../fact.o ../erfc.o ../lgamma.o ../math_utils.o ../exp_func.o ../cut_exp.o ../cip_exp.o ../delta.o ../r1gtoint.o
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ $^  ${CXXFLAGS} -lboost_python -lboost_numpy -lgsl -lblas -framework Python
check_r1gtoint: r1gtoint_bind.so
	python test_r1gtoint.py
install_r1gtoint: r1gtoint_bind.so
	cp r1gtoint_bind.so ~/local/Python/
	cp r1gtoint.py ~/local/Python/

# ==== operation ====
.PHONY: clean
clean:
	rm -f *.o
	rm -f *.so



# ==== old ====
check_gto3d: l2func_bind.so
	python test_gto3d.py

l2func_molint_bind.o: l2func_molint_bind.cpp
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ -c $< ${CXXFLAGS} -framework Python

l2func_molint_bind.so: l2func_molint_bind.cpp ${OBJMOLINT} ../eigen_plus.o
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ $^  ${CXXFLAGS} -lboost_python -lboost_numpy -lgsl -lblas -framework Python

check_molint_bind: l2func_molint_bind.so
	python test_molint.py
