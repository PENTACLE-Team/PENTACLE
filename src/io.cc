#include<iostream>
#include<fstream>
#include<iomanip>
#include<unistd.h>
#include<unordered_map>
#include<particle_simulator.hpp>
#include"kepler.hpp"
#include"io.hpp"


int main(int argc, char *argv[]){
    PS::Initialize(argc, argv);

    double a_in = 0.1;
    double a_out = 1000.0;
    double a_ice = 10.0;
    double f_ice = 4.0;
    const int n = 1000000;
    double * ax = new double[n];
    for(int i=0; i<n; i++){
	ax[i] = HayashiDistributionWithIceLine(a_in, a_out, a_ice, f_ice);
    }
    std::sort(ax, ax+n);
    /*
    for(int i=0; i<n; i++){
	std::cout<<i<<"   "<<ax[i]<<std::endl;
    }
    */
    const int n_i = 1000;
    for(int i=n_i; i<n; i+=n_i){
	double da = ax[i] - ax[i-n_i];
	std::cout<<ax[i-n_i/2]<<"    "<<(double)n_i/(2.0*M_PI*ax[i-n_i/2]*da)<<std::endl;
    }
    PS::Finalize();
    return 0;
}
