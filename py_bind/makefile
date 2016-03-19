OBJS= ../cip_exp.o ../cut_exp.o ../exp_func.o  ../delta.o ../fact.o ../erfc.o ../math_utils.o
OBJS3D = ../angmoment.o ../cints.o ../gto3dset.o
OBJMOLINT = ../molint.o ../angmoment.o ../math_utils.o

l2func_bind.so: l2func_bind.cpp ${OBJS} ${OBJS3D}
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ l2func_bind.cpp ${OBJS} ${OBJS3D} ${CXXFLAGS} -lboost_python -lboost_numpy -lgsl -lblas -framework Python

check_gto3d: l2func_bind.so
	python test_gto3d.py

l2func_molint_bind.o: l2func_molint_bind.cpp
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ -c $< ${CXXFLAGS} -framework Python

l2func_molint_bind.so: l2func_molint_bind.cpp ${OBJMOLINT} ../eigen_plus.o
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ $^  ${CXXFLAGS} -lboost_python -lboost_numpy -lgsl -lblas -framework Python

check_molint_bind: l2func_molint_bind.so
	python test_molint.py


# ==== SymMolInt ====
symmolint_bind.so: symmolint_bind.cpp ../symmolint.o ../fact.o ../mol_func.o ../bmatset.o ../angmoment.o ../eigen_plus.o 
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ $^  ${CXXFLAGS} -lboost_python -lboost_numpy -lgsl -lblas -framework Python

check_symmolint_bind: symmolint_bind.so
	python test_symmolint.py



check: l2func_bind.so
	python test.py

clean:
	rm -f *.o
	rm -f *.so
