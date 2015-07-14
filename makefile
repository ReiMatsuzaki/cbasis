include local.mk
CXXFLAGS=${INC_PATH} -Wall -O3

OBJS=hatom.o lcomb.o prim.o erfc.o fact.o

# this command run correctly on rcclsc
#l2.a: ${OBJS}
#	ar -cr -o $@ ${OBJS}

# this command run correctly on mac
l2.a: ${OBJS}
	ar r $@ ${OBJS}

hatom.o: hatom.cpp hatom.hpp
lcomb.o: lcomb.cpp lcomb.hpp
prim.o:  prim.cpp  prim.hpp
erfc.o:  erfc.cpp  erfc.hpp
fact.o:  fact.cpp  fact.hpp
utest.o: utest.cpp

utest: utest.o l2.a
	${CXX} -o $@ ${CXXFLAGS} ${LIBGTEST} utest.o l2.a
	./utest

clean:
	rm -f *.o
	rm -f *.a
	rm -f utest
