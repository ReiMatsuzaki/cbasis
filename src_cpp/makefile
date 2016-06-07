include local.mk
ARCH=debug
$(shell mkdir -p ${BINDIR})
VALGRIND=valgrind --error-limit=no --tool=memcheck --leak-check=full --show-reachable=no

ifeq (${ARCH},fast)
  # calculation speed
  CXXFLAGS=${INC_PATH} -Wall -O3 -DARG_NO_CHECK -DEIGEN_NO_DEBUG
  RUN=
endif
ifeq (${ARCH},debug)
  # debug mode
  CXXFLAGS=${INC_PATH} -Wall -g
  RUN=VALGRIND
endif
ifeq (${ARCH},mem)
  # for memory check. (very slow)
  CXXFLAGS=${INC_PATH} -Wall -pg -g -fno-inline -DARG_NO_CHECK -DEIGEN_NO_DEBUG
  RUN=VALGRIND
endif
ifeq (${ARCH},time)
  # for time profile
  CXXFLAGS=${INC_PATH} -Wall -pg -g -DARG_NO_CHECK -DEIGEN_NO_DEBUG
  RUN=iprofiler -timeprofiler -T 30
endif
BINDIR=bin/${ARCH}


SRCS:=$(wildcard *.cpp)
DEPS:=$(SRCS:%.cpp=${BINDIR}/%.d)
-include ${DEPS}

# -- basic relation --
${BINDIR}/%.o: %.cpp
	${CXX} -c -o $@ -MMD ${CXXFLAGS} $<

# -- test math --
MATH_OBJS = \
	test_math.o fact.o erfc.o lgamma.o mol_func.o \
	angmoment.o eigen_plus.o bmatset.o 
${BINDIR}/test_math: $(foreach o, ${MATH_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^  ${CXXFLAGS} ${LIBGTEST} -lgsl
.PHONY: check_math
check_math: ${BINDIR}/test_math
	${RUN} ./$<

# -- test symgroup --
SYMGROUP_OBJS=test_symgroup.o symgroup.o
${BINDIR}/test_symgroup: $(foreach o, ${SYMGROUP_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} -lgtest
.PHONY: check_symgroup
check_symgroup: ${BINDIR}/test_symgroup
	${RUN} ./$<

# -- test symmolint --
SYMMOLINT_OBJS = \
	test_symmolint.o symmolint.o symgroup.o bmatset.o angmoment.o \
	eigen_plus.o fact.o mol_func.o
${BINDIR}/test_symmolint: $(foreach o, ${SYMMOLINT_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} -lgtest -lgsl
.PHONY: check_symmolint
check_symmolint: ${BINDIR}/test_symmolint
	${RUN} ./$<

# -- test 2e int --
TWOINT_OBJS = test_2eint.o b2eint.o symmolint.o symgroup.o bmatset.o angmoment.o \
	eigen_plus.o fact.o mol_func.o timer.o
${BINDIR}/test_2eint:  $(foreach o, ${TWOINT_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} -lgtest -lgsl
.PHONY: check_2eint
check_2eint: ${BINDIR}/test_2eint
	${RUN} ./$<
time_2eint: ${BINDIR}/test_2eint
	iprofiler -timeprofiler ./$<

# -- test hf --
HF_OBJS = test_hf.o b2eint.o symmolint.o symgroup.o bmatset.o angmoment.o eigen_plus.o \
	fact.o mol_func.o
${BINDIR}/test_hf: $(foreach o, ${HF_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} -lgtest -lgsl
.PHONY: check_hf
check_hf: ${BINDIR}/test_hf
	${RUNE} ./$<

# -- test L2(R1) --
L2OBJS=fact.o erfc.o lgamma.o math_utils.o exp_func.o cut_exp.o exp_int.o delta.o
${BINDIR}/l2.a: $(foreach o, ${L2OBJS}, ${BINDIR}/$o)
	ar r $@ $^
${BINDIR}/test_l2: ${BINDIR}/test_l2.o ${BINDIR}/l2.a
	${CXX} -o $@ ${CXXFLAGS} $^ -lgtest
.PHONY: /check_l2
check_l2: ${BINDIR}/test_l2
	${VALGRIND} ./$<


# -- test r1gto --
${BINDIR}/r1gtoint.a: $(foreach o, r1gtoint.o op_driv.o eigen_plus.o, ${BINDIR}/$o)
	ar r $@ $^
R1GTOINT_OBJS=test_r1gtoint.o  l2.a r1gtoint.a
${BINDIR}/test_r1gtoint: $(foreach o, ${R1GTOINT_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^  ${CXXFLAGS} ${LIBGTEST} -lgsl
.PHONY: check_r1gtoint
check_r1gtoint: ${BINDIR}/test_r1gtoint
	${RUN} ./$<

# -- test r1gto derivative --
R1GTOINT_DERIV_OBJS=test_r1gtoint_deriv.o r1gtoint.o l2.a eigen_plus.o opt_alpha.o l_algebra.o op_driv.o
${BINDIR}/test_r1gtoint_deriv: $(foreach o, ${R1GTOINT_DERIV_OBJS}, ${BINDIR}/$o)
	${CXX} -o $@ $^ ${CXXFLAGS} ${LIBGTEST}
check_r1gtoint_deriv: ${BINDIR}/test_r1gtoint_deriv
	${RUN} ./$<

clean:
	rm -f ${BINDIR}/*.o
	rm -f ${BINDIR}*.a
	rm -f ${BINDIR}*.d

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


