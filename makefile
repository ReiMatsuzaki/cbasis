include local.mk
ARCH=debug
ifeq (${ARCH},debug)
  CXXFLAGS=${INC_PATH} -Wall -pg -g -fno-inline
else
  CXXFLAGS=${INC_PATH} -Wall -O3
endif
BINDIR=${ARCH}

VALGRIND=valgrind --error-limit=no --tool=memcheck --leak-check=full --show-reachable=no

SRCS:=$(wildcard *.cpp)
DEPS:=$(SRCS:%.cpp=${BINDIR}/%.d)

-include ${DEPS}

# -- basic relation --
${BINDIR}/%.o: %.cpp
	@[ -d ${BINDIR} ] || mkdir -p ${BINDIR}
	${CXX} -c -o $@ -MMD ${CXXFLAGS} $< 

# -- mathematics --
fact.o: fact.cpp fact.hpp
erfc.o:  erfc.cpp  erfc.hpp 
lgamma.o: lgamma.cpp lgamma.hpp typedef.hpp
mol_func.o: mol_func.cpp mol_func.hpp typedef.hpp mult_array.hpp
angmoment.o: angmoment.cpp angmoment.hpp typedef.hpp macros.hpp fact.hpp 

# -- Eigen --
eigen_plus.o: eigen_plus.cpp eigen_plus.hpp typedef.hpp macros.hpp
bmatset.o: bmatset.cpp bmatset.hpp

# -- L2(R1) --
cut_exp.o: cut_exp.hpp exp_func.hpp linspace.hpp
exp_func.o: exp_func.cpp exp_func.hpp linspace.hpp
cip_exp.o: cip_exp.cpp

# -- SymGTOs --
symmolint.o: symmolint.cpp symmolint.hpp bmatset.hpp angmoment.hpp mol_func.hpp

# -- test math --
test_math.o: test_math.cpp math_utils.hpp erfc.hpp lgamma.hpp mol_func.hpp angmoment.hpp eigen_plus.hpp bmatset.hpp
MATH_OBJS = \
	test_math.o fact.o erfc.o lgamma.o mol_func.o \
	angmoment.o eigen_plus.o bmatset.o
${BINDIR}/test_math: $(foreach o, ${MATH_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^  ${CXXFLAGS} ${LIBGTEST} -lgsl
.PHONY: check_math
check_math: ${BINDIR}/test_math
	${VALGRIND} ./$<

# -- test symmolint --
test_symmolint.o: test_symmolint.cpp symmolint.hpp
test_symmolint: test_symmolint.o symmolint.o bmatset.o angmoment.o eigen_plus.o fact.o mol_func.o
	${CXX} -o $@ $^ ${CXXFLAGS} -lgtest -lgsl
.PHONY: check_symmolint
check_symmolint: test_symmolint
	valgrind --error-limit=no --tool=memcheck --leak-check=full --show-reachable=no ./$<

# -- test r1gto --
r1gtoint.o: r1gtoint.cpp r1gtoint.hpp 
test_r1gtoint.o: test_r1gtoint.cpp r1gtoint.hpp cip_exp.hpp
test_r1gtoint: test_r1gtoint.o r1gtoint.o l2.a eigen_plus.o
	${CXX} -o $@ $^  ${CXXFLAGS} ${LIBGTEST} -lgsl
.PHONY: check_r1gtoint
check_r1gtoint: test_r1gtoint
	echo ${ARCH}
	valgrind --error-limit=no --tool=memcheck --leak-check=full --show-reachable=no ./$<


# -- test L2(R1) --
l2.a: fact.o erfc.o lgamma.o math_utils.o exp_func.o cut_exp.o cip_exp.o delta.o
	ar r $@ $^

test.o: test.cpp
test: test.o l2.a
	${CXX} -o $@ ${CXXFLAGS} $^ -lgtest

.PHONY: check
check: test
	./test

clean:
	rm -f *.o
	rm -f *.a
	rm -f utest

# ==== old ====

# -- not used now --
gto3dset.o: gto3dset.cpp gto3dset.hpp

molint.o: molint.cpp molint.hpp math_utils.hpp macros.hpp

# -- not used now --

test_gto3d.o: test_gto3d.cpp gto3d.hpp 
#test_gto3d: test_gto3d.o cints.o angmoment.o gto3dset.o math_utils.o molint.o cip_exp.o exp_func.o eigen_plus.o
test_gto3d: test_gto3d.o cints.o angmoment.o gto3dset.o molint.o spec_func.o eigen_plus.o ${OBJS} 
	${CXX} -o $@ $^  ${CXXFLAGS} ${LIBGTEST} -lgsl
.PHONY: check_gto3d
check_gto3d: test_gto3d
	./test_gto3d

test_gto3d_time.o: test_gto3d_time.cpp
test_gto3d_time:  test_gto3d_time.o cints.o angmoment.o gto3dset.o molint.o spec_func.o eigen_plus.o ${OBJS} timer.o
	${CXX} -o $@ ${CXXFLAGS} ${LIBGTEST} $^ -lgsl
.PHONY: check_gto3d_time
check_gto3d_time: test_gto3d_time
	./test_gto3d_time
profile_gto3d_time: test_gto3d_time
	iprofiler -timeprofiler ./test_gto3d_time

check_gto3d_time_10: test_gto3d_time
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 
	./test_gto3d_time; ./test_gto3d_time; ./test_gto3d_time; 

inc_gamma_grid.o: inc_gamma_grid.cpp
inc_gamma_grid: inc_gamma_grid.o molint.o
	${CXX} -o $@ ${CXXFLAGS} $^ -lgsl
.PHONY: run_inc_gamma_grid
run_inc_gamma_grid: inc_gamma_grid
	./inc_gamma_grid 

l2func_bind.so: wrapper.cpp ${OBJS}
	${CXX} -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o $@ wrapper.cpp ${OBJS} ${CXXFLAGS} -lboost_python  -framework Python

utest_py: l2func_bind.so utest.py
	python utest.py


