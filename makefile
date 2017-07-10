
#Make file for Phase Error for fixed lambda
CXX=g++
RM=rm -f
CPPFLAGS=-g -fopenmp -std=c++11

SRCS=CCOST_phasedrift_parallel.cpp
OBJS=$(subst .cpp,.o, $(SRCS))

all: phase_drift.exe

phase_drift.exe: $(OBJS)
	$(CXX) $(CPPFLAGS) -o phase_drift.exe $(OBJS)

ccost.o: CCOST_phasedrift_parallel.cpp
	

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) phase_drift.exe
