PS_PATH = -I../../../../fdps/src/

use_xc30_naoj = yes

ifeq ($(use_xc30_naoj),yes)
CXX = time CC
CXXFLAGS = -O3 
CXXFLAGS += -Wall
CXXFLAGS += -march=core-avx2
CXXFLAGS += -ffast-math -funroll-loops
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -std=c++11
CXXFLAGS += -DINTRINSIC_X86
endif

#all:make_ic.out
all:remove_particles.out

make_ic.out:make_ic.cc Makefile
	$(CXX) $(PS_PATH) $(CXXFLAGS) -o $@ $<

remove_particles.out:remove_particles.cc Makefile
	$(CXX) $(PS_PATH) $(CXXFLAGS) -o $@ $<

clean:
	rm *.out