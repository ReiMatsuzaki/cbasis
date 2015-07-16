include local.mk
CXXFLAGS=${INC_PATH} -Wall -O3 

OBJS=hatom.o op.o lcomb.o prim.o erfc.o fact.o

# this command run correctly on rcclsc
#l2.a: ${OBJS}
#	ar -cr -o $@ ${OBJS}

# this command run correctly on mac
l2.a: ${OBJS}
	ar r $@ ${OBJS}

hatom.o: hatom.cpp hatom.hpp
op.o: op.cpp op.hpp
lcomb.o: lcomb.cpp lcomb.hpp
prim.o:  prim.cpp  prim.hpp
erfc.o:  erfc.cpp  erfc.hpp
fact.o:  fact.cpp  fact.hpp
utest.o: utest.cpp

utest: utest.o l2.a
	${CXX} -o $@ ${CXXFLAGS} ${LIBGTEST} utest.o l2.a
	./utest

l2func_bind.so: wrapper.cpp lcomb.o prim.o erfc.o fact.o 
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ wrapper.cpp lcomb.o prim.o erfc.o fact.o ${CXXFLAGS} -lboost_python  -framework Python

utest_py: l2func_bind.so utest.py
	python utest.py

clean:
	rm -f *.o
	rm -f *.a
	rm -f utest
