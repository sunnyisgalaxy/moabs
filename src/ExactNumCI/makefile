##to check: it looks like the boost libs are not used and hence not needed.
SRC=numCI.cpp types.cpp stats.cpp 
PRG=numCI

##The corresponding files are $(MOABS_ROOT)/src/bsmap/boost/boost/*hpp
ifndef BOOSTINC
	BOOSTINC=$(MOABS_ROOT)/src/bsmap/boost
endif

##The corresponding files are $(MOABS_ROOT)/src/bsmap/boost/lib/libboost*
ifndef BOOSTLIB
	BOOSTLIB=$(MOABS_ROOT)/src/bsmap/boost/lib
endif


## compiler etc settings used in default make rules
CXX := 			g++

CPPFLAGS := 		-Wall -w -I/usr/local/include  -I${BOOSTINC}

CXXFLAGS := 		-g -O2

LDLIBS := 		-L${BOOSTLIB} -lboost_program_options -Wl,-rpath=${BOOSTLIB} -L${BOOSTLIB} -lboost_system -L${BOOSTLIB} -lboost_filesystem -L${BOOSTLIB} -lboost_thread

all:	$(PRG)

OBJ=numCI.o types.o stats.o  

numCI: $(OBJ)
	$(CXX) -o $(PRG) $(OBJ) $(CPPFLAGS) $(CXXFLAGS) $(LDLIBS) 
	rm -f *.o

numCI.o: numCI.cpp
	$(CXX) -c numCI.cpp $(CPPFLAGS) $(CXXFLAGS)

types.o: types.cpp
	$(CXX) -c types.cpp $(CPPFLAGS) $(CXXFLAGS)

stats.o: stats.cpp
	$(CXX) -c stats.cpp $(CPPFLAGS) $(CXXFLAGS)


clean:
	rm -f *.o numCI

install: $(PRG)
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PRG) $(SRC_ROOT)/bin
