include ../local.mk

VPATH=../src_cpp ../utils
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

## ==== Google test ====
### read README in googletest
CPPFLAGS += -isystem ${GTEST_DIR}/include
CXXFLAGS += -pthread
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
#                $(GTEST_DIR)/include/gtest/internal/*.h

GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
${BINDIR}/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc -o $@
${BINDIR}/gtest.a : ${BINDIR}/gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

## ==== basic relation ====
${BINDIR}/%.o: %.cpp
	${CXX} -c -o $@ -MMD ${CPPFLAGS} ${CXXFLAGS} $<

## ==== Main build ====
OBJS=test_math.o erfc.o 
${BINDIR}/test_math: $(foreach o, ${OBJS}, ${BINDIR}/$o) ${BINDIR}/gtest.a
	${CXX} -o $@ $^ ${CPPFLAGS} ${CXXFLAGS} ${GTESTFLAGS}

## ==== launch unit test ====
check: ${BINDIR}/test_math
	./$^

clean:
	rm -fr ${BINDIR}
