EXE=asr_spatial.exe
CPP=simulation.cpp main.cpp individual.cpp patch.cpp
HPP=simulation.hpp parameters.hpp individual.hpp patch.hpp
CXX=g++
CXXFLAGS=-Wall -ggdb -std=c++14 -O3


$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)


clean :
	rm -rf $(EXE)
