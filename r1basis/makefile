include local.mk

VPATH= ../utils ../math
CPPFLAGS+=-Wall -fPIC
#ifeq (\$(OS),Linux)
#PY_FLAGS=-I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -shared -fPIC -lboost_python 

BINORIG=bin
ARCH=fast
BINDIR=${BINORIG}/${ARCH}

$(shell mkdir -p ${BINDIR})

ifeq (${ARCH}, fast)
  CXXFLAGS+= -O3 -DARG_NO_CHECK -DEIGEN_NO_DEBUG
endif

SRCS:=$(wildcard *.cpp)
DEPS:=$(SRCS:%.cpp=${BINDIR}/%.d)
-include ${DEPS}

## ==== basic relation ====
${BINDIR}/%.o: %.cpp
	${CXX} -c -o $@ -MMD ${CPPFLAGS} ${CXXFLAGS} $<

# ==== google test ====
# read README in googletest
CPPFLAGS += -isystem ${GTEST_DIR}/include
CXXFLAGS += -pthread
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
${BINDIR}/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc -o $@
${BINDIR}/gtest.a : ${BINDIR}/gtest-all.o
	$(AR) $(ARFLAGS) $@ $^


## ==== Main build ====
R1BASIS_OBJS=r1_lc.o r1basis.o erfc.o eigen_plus.o int_exp.o fact.o
OBJS=$(foreach o, ${R1BASIS_OBJS}, ${BINDIR}/$o)
${BINDIR}/r1basis_bind.so: r1basis_bind.cpp ${OBJS}
	${CXX} -o $@ $^ ${CXXFLAGS} ${PY_FLAGS}
	cp $@ .

${BINDIR}/r1basis_test: ${BINDIR}/test.o ${OBJS} ${BINDIR}/gtest.a
	${CXX} -o $@ $^ ${CPPFLAGS} ${CXXFLAGS} ${GTESTFLAGS}

## ==== launch unit test ====
check_erfc: ${BINDIR}/r1basis_bind.so
	python test_erfc.py

check: ${BINDIR}/r1basis_test
	./$^



check_py: ${BINDIR}/r1basis_bind.so 
	python test.py

clean:
	rm -f *.so
	rm -fr ${BINDIR}
