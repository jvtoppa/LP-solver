CXX = g++
CXXFLAGS = -std=c++11 -O3 -DNDEBUG -I$(HOME)/include

ALL_TARGETS = simplex

all: $(ALL_TARGETS)

simplex: simplex.cpp linearalgebra.cpp
	$(CXX) $(CXXFLAGS) simplex.cpp linearalgebra.cpp -o $@ 

clean:
	rm -f $(ALL_TARGETS)
