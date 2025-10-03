CXX = g++
CXXFLAGS = -std=c++11 -O3 -DNDEBUG -I$(HOME)/include
LDFLAGS = -L$(HOME)/lib -lsdsl -ldivsufsort -ldivsufsort64

ALL_TARGETS = simplex

all: $(ALL_TARGETS)

simplex: simplex.cpp linearalgebra.cpp
	$(CXX) $(CXXFLAGS) simplex.cpp linearalgebra.cpp -o $@ $(LDFLAGS)

clean:
	rm -f $(ALL_TARGETS)
