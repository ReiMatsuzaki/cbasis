PY_FLAGS= -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -lboost_python -DPIC -bundle -fPIC -framework Python 
#PY_FLAGS= -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -framework Python 

basic.so: basic.cpp
	${CXX} -o $@ $^ ${CXXFLAGS} ${PY_FLAGS}

check:basic.so
	python -c "from basic import *; print add(1,1)"
	rm basic.so
