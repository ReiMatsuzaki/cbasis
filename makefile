include local.mk
#CXXFLAGS=${INC_PATH} -Wall -O3  
CXXFLAGS=${INC_PATH} -Wall -g

MATH_OBJS=erfc.o math_utils.o
FUNC_OBJS=cut_exp.o exp_func.o delta.o
CIP_OBJS=cip_exp.o

OBJS=${CIP_OBJS} ${FUNC_OBJS} ${MATH_OBJS} 

# this command run correctly on mac
l2.a: ${OBJS}
	ar r $@ ${OBJS}

#fact.o:  fact.cpp  fact.hpp 
erfc.o:  erfc.cpp  erfc.hpp 
math_utils.o: math_utils.cpp math_utils.hpp
cut_exp.o: cut_exp.hpp exp_func.hpp linspace.hpp
exp_func.o: exp_func.cpp exp_func.hpp linspace.hpp
cip_exp.o: cip_exp.cpp
angmoment.o: angmoment.cpp math_utils.hpp  angmoment.hpp 
molint.o: molint.cpp molint.hpp math_utils.hpp macros.hpp

gto3dset.o: gto3dset.cpp gto3dset.hpp

test.o: test.cpp
test: test.o l2.a
	${CXX} -o $@ ${CXXFLAGS} ${LIBGTEST} test.o l2.a

test_gto3d.o: test_gto3d.cpp gto3d.hpp 
test_gto3d: test_gto3d.o cints.o angmoment.o gto3dset.o math_utils.o molint.o
	${CXX} -o $@ ${CXXFLAGS} ${LIBGTEST} $^ -lgsl
.PHONY: check_gto3d
check_gto3d: test_gto3d
	./test_gto3d

l2func_bind.so: wrapper.cpp ${OBJS}
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ wrapper.cpp ${OBJS} ${CXXFLAGS} -lboost_python  -framework Python

utest_py: l2func_bind.so utest.py
	python utest.py


.PHONY: check
check: test
	./test

clean:
	rm -f *.o
	rm -f *.a
	rm -f utest
