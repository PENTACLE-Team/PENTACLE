#include<iostream>
#include<cmath>
#include<particle_simulator.hpp>
#include<unordered_map>
#include"../kepler.hpp"
#include"../io.hpp"

class FP{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot_tot;
    PS::S32 n_ngb;
    PS::S64 id_ngb;
    PS::S32 nptcl_cluster; // the number of particle in reffered cluster
    PS::S32 id_cluster; // id of the cluster
    PS::S32 rank_org;
    PS::S32 adr; // point to first address of pair-interaction array
    FP * next;
    static PS::F64 r_search;
    FP(){
	id = -1;
	mass = -1.0;
	pos = vel = acc = 99999.9;
	pot_tot = 999999.9;
	n_ngb = -1;
	id_ngb = -1;
	nptcl_cluster = -1;
	id_cluster = -1;
	next = nullptr;
	rank_org = -1;
	adr = -1;
    }
    PS::F64vec getPos() const {
	return pos;
    }
    PS::F64 getCharge() const {
	return mass;
    }
    PS::F64 getRSearch() const {
	return r_search;
    }
    void setPos(const PS::F64vec & _pos) {
	pos = _pos;
    }
    void copyFromForce(const FP & p){
	n_ngb = p.n_ngb;
	id_ngb = p.id_ngb;
    }
    void copyFromFP(const FP & p){
	id = p.id;
	mass = p.mass;
	pos = p.pos;
	vel = p.vel;
	n_ngb = p.n_ngb;
	id_ngb = p.id_ngb;
	rank_org = p.rank_org;
    }
    void clear(){
	n_ngb = 0;
	id_ngb = -1;
    }
    void readAscii(FILE* fp) {
        fscanf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d \n", 
                &this->id, &this->mass, 
		&this->pos.x, &this->pos.y, &this->pos.z,  // 3-5
		&this->vel.x, &this->vel.y, &this->vel.z,  // 6-8
		&this->acc.x, &this->acc.y, &this->acc.z,  // 9-11
		&this->pot_tot, &this->n_ngb);
    }

    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %d \n", 
                this->id, this->mass, 
		this->pos.x, this->pos.y, this->pos.z,  // 3-5
		this->vel.x, this->vel.y, this->vel.z,  // 6-8
		this->acc.x, this->acc.y, this->acc.z,  // 9-11
		this->pot_tot, this->n_ngb);
    }

    void dump()const{
	std::cout<<"id= "<<id<<std::endl;
	std::cout<<"mass= "<<mass<<std::endl;
	std::cout<<"pos= "<<pos<<std::endl;
	std::cout<<"n_ngb= "<<n_ngb<<std::endl;
    }
};

template<class Tptcl0, class Tptcl1>
PS::F64 CalcEng(const Tptcl0 & p0,
		const Tptcl1 & p1){
    PS::F64vec dr = p0.pos - p1.pos;
    PS::F64vec dv = p0.vel - p1.vel;
    PS::F64 m_tot = p0.mass + p1.mass;
    PS::F64 r_inv = 1.0 / sqrt(dr*dr);
    return 0.5*dv*dv - m_tot*r_inv;
}

void CountNngb(const FP * fp_i,
	       const PS::S32 n_ip,
	       const FP * fp_j,
	       const PS::S32 n_jp,
	       FP * force){
    const PS::F64 r_search_sq = FP::r_search * FP::r_search;
    for(PS::S32 i=0; i<n_ip; i++){
	PS::F64vec pos_i = fp_i[i].pos;
	for(PS::S32 j=0; j<n_jp; j++){
	    if(fp_i[i].id == fp_j[j].id) continue;
	    PS::F64vec pos_j = fp_j[j].pos;
	    PS::F64vec rij = pos_i - pos_j;
	    if(rij*rij <= r_search_sq){
		force[i].n_ngb++;
		force[i].id_ngb = fp_j[j].id;
	    }
	}
    }
}

int main(int argc, char *argv[]){

    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);

    PS::S32 n_glb = 125000;
    PS::F64 ax_in = 0.95; //[AU]
    PS::F64 ax_out = 1.05; //[AU]
    PS::F64 ecc_sigma_hill = 2.0;
    PS::F64 inc_sigma_hill = 1.0;
    PS::F64 dens = 10.0; //[g/cm^2]
    PS::S32 seed = 0;
    PS::F64 f_ice = 1.0;
    PS::F64 ax_ice = 0.0;
    PS::F64 power = -1.5;
    while((c=getopt(argc,argv,"a:A:i:e:o:d:D:E:t:r:T:n:N:s:S:l:R:r:x:f:I:p:H:P:X:h")) != -1){
        switch(c){
        case 'a':
            ax_in = atof(optarg);
	    std::cerr<<"ax_in="<<ax_in<<std::endl;
            break;
        case 'A':
            ax_out = atof(optarg);
            std::cerr<<"ax_out="<<ax_out<<std::endl;
            break;
        case 'i':
            inc_sigma_hill = atof(optarg);
            std::cerr<<"inc_sigma_hill="<<inc_sigma_hill<<std::endl;
            break;
        case 'e':
            ecc_sigma_hill = atof(optarg);
            std::cerr<<"ecc_sigma_hill="<<ecc_sigma_hill<<std::endl;
            break;
        case 'o':
            sprintf(dir_name,optarg);
            std::cerr<<"dir_name="<<dir_name<<std::endl;
            break;
        case 'd':
            dt_soft = 1.0 / atof(optarg);
            std::cerr<<"dt_soft="<<dt_soft<<std::endl;
            break;
        case 'D':
            dt_snp = atof(optarg);
            std::cerr<<"dt_snp="<<dt_snp<<std::endl;
            break;
        case 'E':
            eta = atof(optarg);
            std::cerr<<"eta="<<eta<<std::endl;
            break;
        case 't':
            theta = atof(optarg);
            std::cerr<<"theta="<<theta<<std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            std::cerr<<"time_end="<<time_end<<std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
            break;
        case 'N':
            n_glb = atol(optarg);
            std::cerr<<"n_glb="<<n_glb<<std::endl;
            break;
        case 's':
            n_smp_ave = atoi(optarg);
            std::cerr<<"n_smp_ave="<<n_smp_ave<<std::endl;
            break;
        case 'S':
            search_factor = atoi(optarg);
            std::cerr<<"search_factor="<<search_factor<<std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr<<"n_leaf_limit="<<n_leaf_limit<<std::endl;
            break;
        case 'r':
            r_merge_factor = atof(optarg);
            std::cerr<<"r_merge_factor="<<r_merge_factor<<std::endl;
            break;
        case 'R':
            r_cut_factor = atof(optarg);
            std::cerr<<"r_cut_factor="<<r_cut_factor<<std::endl;
            break;
        case 'x':
            seed = atoi(optarg);
            std::cerr<<"seed="<<seed<<std::endl;
            break;
        case 'f':
            f_ice = atof(optarg);
            std::cerr<<"f_ice="<<f_ice<<std::endl;
            break;
        case 'I':
            ax_ice = atof(optarg);
            std::cerr<<"ax_ice="<<ax_ice<<std::endl;
            break;
        case 'p':
            power = atof(optarg);
            std::cerr<<"power="<<power<<std::endl;
	    break;
        case 'H':
            dens = atof(optarg);
            std::cerr<<"dens="<<dens<<std::endl;
	    break;
        case 'P':
            dens_planet = atof(optarg);
            std::cerr<<"dens_planet="<<dens_planet<<std::endl;
            break;
        case 'X':
	    dt_limit_hard_factor = atof(optarg);
            std::cerr<<"dt_limit_hard_factor="<<dt_limit_hard_factor<<std::endl;
            break;
        case 'h':
            std::cerr<<"a: ax_in [AU] (dafult 0.98AU)"<<std::endl;
            std::cerr<<"A: ax_out [AU] (dafult 1.02AU)"<<std::endl;
            std::cerr<<"i: inc_sigma_hill (dafult 1.0)"<<std::endl;
            std::cerr<<"e: ecc_sigma_hill (dafult 2.0)"<<std::endl;
            std::cerr<<"o: dir name of output"<<std::endl;
            std::cerr<<"d: inv_dt (dafult 16 ~ 0.01yr  )"<<std::endl;
            std::cerr<<"D: dt_snp (dafult 100 ~ 16yr  )"<<std::endl;
            std::cerr<<"E: eta (dafult 0.1)"<<std::endl;
            std::cerr<<"t: theta (dafult: 0.5)"<<std::endl;
            std::cerr<<"T: time_end (dafult: 8000 ~ 1yr)"<<std::endl;
            std::cerr<<"n: n_group_limit (dafult: 64.0)"<<std::endl;
            std::cerr<<"N: n_glb (dafult: 16384)"<<std::endl;
            std::cerr<<"s: n_smp_ave (dafult: 100)"<<std::endl;
            std::cerr<<"S: search_factor (dafult: 3.0)"<<std::endl;
            std::cerr<<"l: n_leaf_limit (dafult: 8)"<<std::endl;
            std::cerr<<"r: r_merge_factor (dafult: 1.0)"<<std::endl;
            std::cerr<<"R: r_cut_factor (dafult: 3.0)"<<std::endl;
            std::cerr<<"x: seed (dafult: 0)"<<std::endl;
	    std::cerr<<"f: f_ice (dafult: 1.0)"<<std::endl;
	    std::cerr<<"I: ax_ice (dafult: 0.0[AU])"<<std::endl;
	    std::cerr<<"p: power (dafult: -1.5)"<<std::endl;
	    std::cerr<<"H: dens [g/cm^2]@1[AU] (dafult: 10.0)"<<std::endl;
            std::cerr<<"P: dens_planet [g/cm^3](dafult: 2.0)"<<std::endl;
            std::cerr<<"X: dt_limit_hard_factor(dafult: 4.0 -> dt_limit_hard = dt_soft/4.0)"<<std::endl;
            PS::Finalize();
            return 0;
        }
    }
    PS::ParticleSystem<FPSoft> system_soft;
    system_soft.initialize();
    PS::S32 n_loc;
    SetParticleKeplerDisk(system_soft, n_glb, n_loc, time_sys,
			  ax_in, ax_out, ecc_sigma_hill, inc_sigma_hill, dens, mass_sun, ax_ice, f_ice, power, seed);

    PS::F64 r_out, r_in, vel_disp, r_hill_max_glb, mass_planet_max_glb, mass_planet_tot_glb;
    PS::F64 ecc_rms, inc_rms;
    PS::F64 Tkep_min;
    SetRoutRinRhillVdisp(system_soft,           pos_sun,
                         vel_sun,               mass_sun,
                         ratio_r_cut,           r_cut_factor,
                         r_out,                 r_in, 
                         r_hill_max_glb,        vel_disp,
			 ecc_rms,               inc_rms,
                         mass_planet_max_glb,   mass_planet_tot_glb,
			 Tkep_min);

    PS::DomainInfo dinfo;
    dinfo.initialize();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    PS::S32 nx = sqrt((PS::F64)n_proc-0.000001)+1;
    while( n_proc % nx != 0) nx++;
    PS::S32 ny = n_proc / nx;
    PS::S32 nz = 1;
    if(PS::Comm::getRank() == 0){
        fout_log<<"nx,ny,nz= "<<nx<<" "<<ny<<" "<<nz<<std::endl;
    }
    dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);
    dinfo.decomposeDomainAll(system_soft);

    system_soft.exchangeParticle(dinfo); 

    const PS::F64 theta = 0.5;
    const PS::F64 n_leaf_limit = 8;
    const PS::F64 n_group_limit = 512;
    tree_soft.initialize(n_glb, theta, n_leaf_limit, n_group_limit);
    tree_soft.calcForceAllAndWriteBack(CountNngb, system_soft, dinfo);

    n_loc = system_soft.getNumberOfParticleLocal();

    std::unordered_map<PS::S32, PS::S32> idx_to_adr_loc; // point to adr of pbase
    for(PS::S32 i=0; i<n_loc; i++){
	idx_to_adr_loc.insert( std::pair<PS::S32, PS::S32>(system[i].id, i));
    }
    PS::S32 n_loc_new = n_loc;
    for(PS::S32 i=0; i<n_loc; i++){
	FP * nbl = NULL;
	PS::S32 n_ngb = tree_soft.getNeighborListOneParticle(system[i], nbl); // self is included
	if(n_ngb <= 1) break;
	n_ngb_tot_loc += (n_ngb-1);
	PS::S32 id_new = nbl->id;
	for(PS::S32 j=1; j<n_ngb_tot_loc; j++){
	    PS::F64 eng = CalcEng(system[i], nbl[j]);
	    if(eng < 0.0){
		n_loc_new--;
	    }
	}

    }

    return 0;
}
