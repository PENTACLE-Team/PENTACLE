#include<iostream>
#include<particle_simulator.hpp>
#include"phantomquad_for_p3t_x86.hpp"

int main(int argc, char *argv[]){
    PS::Initialize(argc, argv);
    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank() );

    const int nloop = 100;
    const int ni = 25;
    const int nj = 400;

    float xi[nloop][ni][3];
    float xj[nloop][nj][3];
    float mj = 1e-10;

    for(int l=0; l<nloop; l++){
	for(int i=0; i<ni; i++){
	    for(int k=0; k<3; k++){
		xi[l][i][k] = mt.genrand_res53();
	    }
	}
	for(int j=0; j<ni; j++){
	    for(int k=0; k<3; k++){
		xj[l][j][k] = mt.genrand_res53();
	    }
	}
    }

    const double eps2 = 0.0;
    const double r_crit2 = 0.1;
    PhantomGrapeQuad pg;

    double tcal = PS::GetWtime();
    for(int l=0; l<nloop; l++){
        pg.set_eps2(eps2);
        pg.set_r_crit2(r_crit2);
	for(int i=0; i<ni; i++){
	    pg.set_xi_one(i, xi[l][i][0], xi[l][i][1], xi[l][i][2]);
	}
	for(int j=0; j<nj; j++){
	    pg.set_epj_one(j, xj[l][j][0], xj[l][j][1], xj[l][j][2], mj);
	}
	pg.run_epj_for_p3t(ni, nj);
	//pg.run_epj(ni, nj);
    }
    tcal = PS::GetWtime() - tcal;

    std::cout<<"nloop="<<nloop
	     <<" ni="<<ni
	     <<" nj="<<nj
	     <<" tcal="<<tcal
	     <<" (double)(ni*nj*nloop*24)/tcal="<<(double)(ni*nj*nloop*24)/tcal
	     <<" (double)(ni*nj*nloop*30)/tcal="<<(double)(ni*nj*nloop*30)/tcal<<std::endl;

    PS::Finalize();

    return 0;
}
