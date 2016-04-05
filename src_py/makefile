# ==== Common settings ====
CXXFLAGS=-Wall
PY_FLAGS=-I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -lboost_python -lboost_numpy  -DPIC -bundle -fPIC -framework Python
LDFLAGS=-lgsl -lblas

BIN_ORIG=bin
ARCH=debug
BINDIR=bin/${ARCH}
$(shell mkdir -p ${BINDIR})

ifeq (${ARCH},fast)
  # calculation speed
  CXXFLAGS=${INC_PATH} -Wall -O3 -DARG_NO_CHECK -DEIGEN_NO_DEBUG
endif
ifeq (${ARCH},debug)
  # debug mode
  CXXFLAGS=${INC_PATH} -Wall -g 
endif
ifeq (${ARCH},mem)
  # for memory check. (very slow)
  CXXFLAGS=${INC_PATH} -Wall -pg -g -fno-inline -DARG_NO_CHECK -DEIGEN_NO_DEBUG
endif


SRCS:=$(wildcard *.cpp)
DEPS:=$(SRCS:%.cpp=${BINDIR}/%.d)
-include ${DEPS}

# -- basic relation --
${BINDIR}/%.o: ../src_cpp/%.cpp
	${CXX} -c -o $@ -MMD ${CXXFLAGS} $<

# ==== L2(R1) ====
L2_OBJS=fact.o erfc.o lgamma.o math_utils.o exp_func.o cut_exp.o cip_exp.o delta.o
${BINDIR}/l2func_bind.so: l2func_bind.cpp $(foreach o, ${L2_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${PY_FLAGS} ${LDFLAGS}
.PHONY: check_l2
check_l2: ${BINDIR}/l2func_bind.so
	ln -s ${BINDIR}/l2func_bind.so .
	python test.py
	rm l2func_bind.so
install_l2: ${BINDIR}/l2func_bind.so
	cp $< ~/local/Python

# ==== SymMolInt ====
SYMMOLINT_OBJS=symmolint.o fact.o mol_func.o bmatset.o angmoment.o eigen_plus.o 
${BINDIR}/symmolint_bind.so: symmolint_bind.cpp $(foreach o, ${SYMMOLINT_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${PY_FLAGS} ${LDFLAGS}
check_symmolint_bind: ${BINDIR}/symmolint_bind.so
	ln -s ${BINDIR}/symmolint_bind.so .
	python test_symmolint.py
	rm symmolint_bind.so
install_symmolint: ${BINDIR}/symmolint_bind.so
	cp $< ~/local/Python/
	cp symmolint_plus.py ~/local/Python/
	cp symmolint.py ~/local/Python/

# ==== R1GTOs ====
R1GTOINT_OBJS=fact.o erfc.o lgamma.o math_utils.o exp_func.o \
	cut_exp.o cip_exp.o delta.o eigen_plus.o opt_alpha.o r1gtoint.o
${BINDIR}/r1gtoint_bind.so: r1gtoint_bind.cpp $(foreach o, ${R1GTOINT_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${PY_FLAGS} ${LDFLAGS}
check_r1gtoint: ${BINDIR}/r1gtoint_bind.so
	rm -f ./r1gtoint_bind.so
	ln -s ${BINDIR}/r1gtoint_bind.so .
	python test_r1gtoint.py
	rm r1gtoint_bind.so
install_r1gtoint: ${BINDIR}/r1gtoint_bind.so
	cp ${BINDIR}/r1gtoint_bind.so ~/local/Python/
	cp r1gtoint.py ~/local/Python/

# ==== operation ====
.PHONY: clean
clean:
	rm -r -f ${BINDIR}

# ==== old ====
check_gto3d: l2func_bind.so
	python test_gto3d.py

l2func_molint_bind.o: l2func_molint_bind.cpp
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ -c $< ${CXXFLAGS} -framework Python

l2func_molint_bind.so: l2func_molint_bind.cpp ${OBJMOLINT} ../eigen_plus.o
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ $^  ${CXXFLAGS} -lboost_python -lboost_numpy -lgsl -lblas -framework Python

check_molint_bind: l2func_molint_bind.so
	python test_molint.py
