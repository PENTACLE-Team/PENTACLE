#PS_PATH = -I../../fdps/src/
#PS_PATH = -I../fdps/
PS_PATH = -I../../../fdps/src/

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

#CXXFLAGS += -DONEBODY 
#CXXFLAGS += -DDEBUG_PRINT_PLANET
#CXXFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY
#CXXFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION
#CXXFLAGS += -DUSE_STD_SORT
#CXXFLAGS += -DFORDEBUG
#CXXFLAGS += -DFREE_MALLOC_TREE 
#CXXFLAGS += -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE
#CXXFLAGS += -DORIGINAL
#CXXFLAGS += -DP3T_64BIT
#CXXFLAGS += -DENERGY_CHECK 
#CXXFLAGS += -DFORCE_CHECK
#CXXFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
#CXXFLAGS += -DREAD_FILE

SRC = main_kepler.cc hard.hpp class.hpp force.hpp io.hpp kepler.hpp phantomquad_for_p3t_x86.hpp domain.hpp profile.hpp

#all:p3t_kepler_quad_merge_phantom_fix_64bit_cont.out
#all:profile_64bit_cont.out
#all:profile_nobarrier.out 
#all:profile_barrier.out 
#all:profile_barrier_64bit.out 
#all:main_nobarrier.out 
#all:main_nobarrier_gas.out 
#all:main_nobarrier_gas_cont.out 
#all:profile_nobarrier_cont.out 
#all:profile_barrier_cont.out 
#all:remove_body.out
#all:main_no_merge.out
all:eng_chk_64bit_cont.out

eng_chk_64bit_cont.out:$(SRC) Makefile
	$(CXX) -DUSE_QUAD -DDIV_FIX -DENERGY_CHECK -DP3T_64BIT -DREAD_FILE $(PS_PATH) $(CXXFLAGS) -o $@ $<

main_no_merge.out:$(SRC) Makefile
	$(CXX) -DUSE_QUAD -DDIV_FIX -DUSE_QUAD $(PS_PATH) $(CXXFLAGS) -o $@ $<

remove_body.out:$(SRC) Makefile
	$(CXX) -DREMOVE_BODY -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

main_nobarrier_gas_cont.out:$(SRC) Makefile
	$(CXX) -DREAD_FILE -DAERO_DRAG -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

main_nobarrier_gas.out:$(SRC) Makefile
	$(CXX) -DAERO_DRAG -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

main_nobarrier.out:$(SRC) Makefile
	$(CXX) -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_nobarrier_cont.out:$(SRC) Makefile
	$(CXX) -DREAD_FILE -DPROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_barrier_cont.out:$(SRC) Makefile
	$(CXX) -DREAD_FILE -DPROFILE -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

debug.out:$(SRC) Makefile
	$(CXX) -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_barrier_64bit.out:$(SRC) Makefile
	$(CXX) -DP3T_64BIT -DPROFILE -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_64bit_cont.out:$(SRC) Makefile
	$(CXX) -DREAD_FILE -DP3T_64BIT -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_64bit.out:$(SRC) Makefile
	$(CXX) -DP3T_64BIT -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile.out:$(SRC) Makefile
	$(CXX) -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_nobarrier_gas.out:$(SRC) Makefile
	$(CXX) -DAERO_DRAG -DPROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_nobarrier.out:profile_nobarrier.o Makefile
	$(CXX) -DPROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_nobarrier.o:$(SRC) Makefile
	$(CXX) -c -DPROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_barrier.out:profile_barrier.o Makefile
	$(CXX) -DPROFILE -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_barrier.o:$(SRC) Makefile
	$(CXX) -c -DPROFILE -DPARTICLE_SIMULATOR_BARRIER_FOR_PROFILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_nobarrier_gnu.out:profile_nobarrier_gnu.o Makefile
	$(CXX) -DUSE_GNU_PARALLEL_SORT -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

profile_nobarrier_gnu.o:$(SRC) Makefile
	$(CXX) -c -DUSE_GNU_PARALLEL_SORT -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

force_chk_quad.out:$(SRC) Makefile
	$(CXX) -DP3T_64BIT   -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

p3t_kepler_quad_merge_phantom_fix.out:$(SRC) Makefile
	$(CXX)  -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<

p3t_kepler_quad_merge_phantom_fix_64bit_cont.out:$(SRC) Makefile
	$(CXX) -DP3T_64BIT -DREAD_FILE -DUSE_QUAD -DDIV_FIX -DMERGE $(PS_PATH) $(CXXFLAGS) -o $@ $<


clean:
	rm *.out *.o
cleanall:
	rm *.out *.hpp~ *.cc~ *.h~
