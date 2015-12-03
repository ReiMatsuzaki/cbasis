OBJS=../cip.o ../fact.o ../erfc.o ../exp_func.o ../math_utils.o
l2func_bind.so: l2func_bind.cpp ${OBJS}
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ l2func_bind.cpp ${OBJS} ${CXXFLAGS} -lboost_python  -framework Python

test: l2func_bind.so
	python test.py
