include local.mk

CXXFLAGS+=-Wall -PIC
#ifeq (\$(OS),Linux)
PY_FLAGS=-I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -shared -fPIC -lboost_python 

BINORIG=bin
ARCH=fast
BINDIR=${BINORIG}/${ARCH}

$(shell mkdir -p ${BINDIR})

ifeq (${ARCH}, fast)
  CXXFLAGS=${INC_PATH} -Wall -O3 -DARG_NO_CHECK -DEIGEN_NO_DEBUG  -fPIC
endif


SRCS:=$(wildcard *.cpp)
DEPS:=$(SRCS:%.cpp=${BINDIR}/%.d)
-include ${DEPS}

## ==== basic relation ====
${BINDIR}/%.o: %.cpp
	${CXX} -c -o $@ -MMD ${CXXFLAGS} $<


## ==== Main build ====
R1BASIS_OBJS=r1_lc.o r1basis.o
OBJS=$(foreach o, ${R1BASIS_OBJS}, ${BINDIR}/$o)
${BINDIR}/r1basis_bind.so: r1basis_bind.cpp ${OBJS}
	${CXX} -o $@ $^ ${CXXFLAGS} ${PY_FLAGS}
	cp $@ .


## ==== launch unit test ====
check: ${BINDIR}/r1basis_bind.so
	python test.py

clean:
	rm -f *.so
	rm -fr ${BINDIR}
