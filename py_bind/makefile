OBJS= ../cip_exp.o ../cut_exp.o ../exp_func.o  ../delta.o ../fact.o ../erfc.o ../math_utils.o
OBJS3D = ../angmoment.o ../cints.o ../gto3dset.o

l2func_bind.so: l2func_bind.cpp ${OBJS} ${OBJS3D}
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ l2func_bind.cpp ${OBJS} ${OBJS3D} ${CXXFLAGS} -lboost_python -lboost_numpy -lgsl -lblas -framework Python

check_gto3d: l2func_bind.so
	python test_gto3d.py

check: l2func_bind.so
	python test.py
