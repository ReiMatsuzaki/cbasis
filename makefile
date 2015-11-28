include local.mk
CXXFLAGS=${INC_PATH} -Wall -O3 

#OBJS=hatom.o op.o lcomb.o cut_prim.o delta.o prim.o erfc.o fact.o
MATH_OBJS=fact.o erfc.o
FUNC_OBJS=exp_func.o 
CIP_OBJS=cip.o

OBJS=${CIP_OBJS} ${FUNC_OBJS} ${MATH_OBJS} 

# this command run correctly on rcclsc
#l2.a: ${OBJS}
#	ar -cr -o $@ ${OBJS}

# this command run correctly on mac
l2.a: ${OBJS}
	ar r $@ ${OBJS}


exp_func.o: exp_func.hpp func.hpp
erfc.o:  erfc.cpp  erfc.hpp
fact.o:  fact.cpp  fact.hpp
test.o: test.cpp

test: test.o l2.a
	${CXX} -o $@ ${CXXFLAGS} ${LIBGTEST} test.o l2.a

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
