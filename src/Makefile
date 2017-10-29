PS_PATH = -I../../FDPS/src/

#use_k_computer = yes
use_xc30_naoj = yes

ifeq ($(use_k_computer),yes)
CXX = time mpiFCCpx
CXXFLAGS = -Kfast
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -Kopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -x32
CXXFLAGS += -Xg
CXXFLAGS += -DFAST_ALL_TO_ALL_FOR_K
CXXFLAGS += -DFAST_WALK_K
CXXFLAGS += -std=c++11
#CXXFLAGS += -NRtrap
CXXFLAGS += -Nfjcex
CXXFLAGS += -Krestp=all
CXXFLAGS += -DINTRINSIC_K
endif

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
#CXXFLAGS += -DUSE_GNU_PARALLEL_SORT
endif

SRC = main_kepler.cc hard.hpp class.hpp force.hpp io.hpp kepler.hpp phantomquad_for_p3t_x86.hpp domain.hpp profile.hpp

all:main_nobarrier.out
all:main_nobarrier_cont.out
all:main_nobarrier_gas.out
all:main_nobarrier_gas_cont.out

main_nobarrier_gas_cont.out:$(SRC) Makefile
	$(CXX) -DREAD_FILE -DAERO_DRAG -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

main_nobarrier_gas.out:$(SRC) Makefile
	$(CXX) -DAERO_DRAG -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

main_nobarrier_cont.out:$(SRC) Makefile
	$(CXX) -DREAD_FILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

main_nobarrier.out:$(SRC) Makefile
	$(CXX) -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<


clean:
	rm *.out *.o
cleanall:
	rm *.out *.hpp~ *.cc~ *.h~
