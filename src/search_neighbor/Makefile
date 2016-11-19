PS_PATH = -I../../../../fdps/src/

CXX = time CC
CXXFLAGS = -O3 
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
CXXFLAGS += -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -std=c++11
CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION

all:main_cont.out 
#all:main.out

main.out:main.cc Makefile
	$(CXX) $(PS_PATH) $(CXXFLAGS) -o $@ $<

main_cont.out:main.cc Makefile
	$(CXX) -DREAD_FILE $(PS_PATH) $(CXXFLAGS) -o $@ $<

main_org2_cont.out:main.org2.cc Makefile
	$(CXX) -DREAD_FILE $(PS_PATH) $(CXXFLAGS) -o $@ $<
