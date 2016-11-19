#include<iostream>
#include<cmath>
#include<particle_simulator.hpp>
#include<unordered_map>
#include"../kepler.hpp"
#include"../io.hpp"

void GatherV(){
}

struct Pbase{
    PS::S32 id;
    PS::S32 adr;
    //PS::S32 cluster_id;
    PS::S32 n_ngb;
    bool flag_search;
    Pbase * next;
    Pbase(){
	id = adr = -1;
	n_ngb = 0;
	flag_search = false;
	next = nullptr;
    }
};

template<class T>
PS::S32 GetNUnique(T * val, const PS::S32 n){
    PS::S32 n_cnt = 0;
    T ref = -9999;
    for(PS::S32 i=0; i<n; i++){
	if(ref != val[i].n_ptcl){
	    ref = val[i].n_ptcl;
	    n_cnt++;
	}
    }
    return n_cnt;
}

void MakeUniformBox(const PS::F64 mass_glb,
		    const PS::S64 n_glb,
		    const PS::S64 n_loc,
		    PS::F64 * mass,
		    PS::F64vec * pos,
		    PS::F64vec * vel,
		    const PS::S32 seed = 0){

    PS::MTTS mt;
    mt.init_genrand( PS::Comm::getRank() + PS::Comm::getNumberOfProc()*seed );
    for(PS::S32 i=0; i<n_loc; i++){
        mass[i] = mass_glb / n_glb;
	pos[i].x = mt.genrand_res53();
	pos[i].y = mt.genrand_res53();
	pos[i].z = mt.genrand_res53();
	vel[i].x = vel[i].y = vel[i].z = 0.0;
    }
}

template<class Tpsys>
void SetParticleUniformBox(Tpsys & psys,
			   const PS::S64 n_glb,
			   PS::S32 & n_loc,
			   PS::F64 & t_sys,
			   const PS::S32 seed = 0){
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
#if 1
    n_loc = n_glb / n_proc; 
    if( n_glb % n_proc > my_rank) n_loc++;
    PS::S64 i_h = n_glb/n_proc*my_rank;
    if( n_glb % n_proc  > my_rank) i_h += my_rank;
    else i_h += n_glb % n_proc;
#else
    PS::S64 i_h = 0;
    if(PS::Comm::getRank()==0){
	n_loc = n_glb;
    }
    else{
	n_loc = 0;
    }
#endif
    psys.setNumberOfParticleLocal(n_loc);

    PS::F64 * mass = new PS::F64[n_loc];
    PS::F64vec * pos = new PS::F64vec[n_loc];
    PS::F64vec * vel = new PS::F64vec[n_loc];
    t_sys = 0.0;
    const PS::F64 m_tot = 1.0;
    MakeUniformBox(m_tot, n_glb, n_loc, mass, pos, vel, seed);
    for(PS::S32 i=0; i<n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos = pos[i];
        psys[i].vel = vel[i];
        psys[i].id = i_h + i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

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
    void dump()const{
	std::cout<<"id= "<<id<<std::endl;
	std::cout<<"mass= "<<mass<<std::endl;
	std::cout<<"pos= "<<pos<<std::endl;
	std::cout<<"n_ngb= "<<n_ngb<<std::endl;
    }
};

PS::F64 FP::r_search;

struct Cluster{
    PS::S32 id;
    PS::S32 n_ptcl;
    Cluster(): id(0), n_ptcl(0){}
};

void SearchClusterList(FP * fp_target,
                       FP * cluster_top,
                       FP *& rest_first,
                       Cluster & cluster,
                       const PS::F64 r_search_sq){
    fp_target->next = nullptr;
    fp_target->id_cluster = cluster.id;
    FP * curr = rest_first;
    FP * prev = nullptr;
    while(curr != nullptr){
        FP * next = curr->next;
        if(fp_target != curr){
	    PS::F64 r_sq = (fp_target->pos - curr->pos)*(fp_target->pos - curr->pos);
            if(r_sq < r_search_sq){
                cluster.n_ptcl++;
                if(prev==nullptr) rest_first = next;
                else prev->next = next;
                cluster_top->next = curr;
                cluster_top = curr;
                SearchClusterList(curr, cluster_top, rest_first, cluster, r_search_sq);
                curr = rest_first;
                prev = nullptr;
            }
            else{
                prev = curr;
                curr = next;
            }
        }
        else{
            prev = curr;
            curr = next;
        }
    }
}

void SearchClusterSlow(FP * fp_target, 
		       FP fp[],
		       Cluster & cluster,
		       const int n, 
		       const PS::F64 r_search_sq){
    fp_target->id_cluster = cluster.id;
    for(PS::S32 i=0; i<n; i++){
	if(fp[i].id_cluster != -1) continue;
	PS::F64 r_sq = (fp_target->pos - fp[i].pos)*(fp_target->pos - fp[i].pos);
	if(r_sq <= r_search_sq){
	    cluster.n_ptcl++;
	    SearchClusterSlow(fp+i, fp, cluster, n, r_search_sq);
	}
    }
}

void SearchClusterNgb(Pbase * target, 
		      Pbase * p_first, 
		      Cluster & cluster,
		      std::unordered_map<PS::S32, PS::S32> & idx_to_adr,
		      std::pair<PS::S32, PS::S32> * idx){
    target->flag_search = true;
    cluster.n_ptcl++;
    PS::S32 target_adr = target->adr;
    PS::S32 n_ngb = target->n_ngb;
    for(PS::S32 i=0; i<n_ngb; i++){
	PS::S32 id_tmp = idx[target_adr+i].second;
	PS::S32 adr = idx_to_adr[ id_tmp ];
	//std::cout<<"adr="<<adr<<" id_tmp="<<id_tmp<<std::endl;
	if( (p_first+adr)->flag_search == false){
	    SearchClusterNgb(p_first+adr, p_first, cluster, idx_to_adr, idx);
	}
    }
}


void SearchClusterNgbWithCheck(Pbase * target, 
			       Pbase * p_first,
			       Pbase *& p_top,
			       Cluster & cluster,
			       std::unordered_map<PS::S32, PS::S32> & idx_to_adr,
			       std::vector< std::pair<PS::S32, PS::S32> > & idx,
			       bool & fg_isolated, 
			       PS::S32 & list_len,
			       std::ofstream & fout){
    target->flag_search = true;
    target->next = nullptr;
    cluster.n_ptcl++;
    PS::S32 target_adr = target->adr;
    PS::S32 n_ngb = target->n_ngb;
    for(PS::S32 i=0; i<n_ngb; i++){
	PS::S32 id_tmp = idx[target_adr+i].second;
	auto itr = idx_to_adr.find(id_tmp);
	if(itr != idx_to_adr.end()){
	    PS::S32 adr = idx_to_adr[ id_tmp ];
	    if( (p_first+adr)->flag_search == false){
		//fout<<"id_tmp="<<id_tmp<<" p_top="<<p_top<<std::endl;
		p_top->next = p_first+adr;
		p_top = p_first+adr;
		list_len++;
		SearchClusterNgbWithCheck(p_first+adr, p_first, p_top, cluster, idx_to_adr, idx, fg_isolated, list_len, fout);
	    }
	}
	else{
	    // neighbor is allocated in another node
	    fg_isolated = false;
	}
    }
}


void ClearFlag(Pbase * target, 
	       Pbase * p_first,
	       Pbase * p_top,
	       Cluster & cluster,
	       std::unordered_map<PS::S32, PS::S32> & idx_to_adr,
	       std::vector< std::pair<PS::S32, PS::S32> > & idx,
	       bool & fg_isolated){
    target->flag_search = false;
    target->next = nullptr;
    cluster.n_ptcl--;
    PS::S32 target_adr = target->adr;
    PS::S32 n_ngb = target->n_ngb;
    for(PS::S32 i=0; i<n_ngb; i++){
	PS::S32 id_tmp = idx[target_adr+i].second;
	auto itr = idx_to_adr.find(id_tmp);
	if(itr != idx_to_adr.end()){
	    PS::S32 adr = idx_to_adr[ id_tmp ];
	    if( (p_first+adr)->flag_search == true){
		p_top->next = p_first+adr;
		p_top = p_first+adr;
		ClearFlag(p_first+adr, p_first, p_top, cluster, idx_to_adr, idx, fg_isolated);
	    }
	}
	else{
	    // neighbor is allocated in another node
	    fg_isolated = false;
	}
    }
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

void CalcNothing(const FP * fp_i,
		 const PS::S32 n_ip,
		 const FP * fp_j,
		 const PS::S32 n_jp,
		 FP * force){
}

class Energy{
public:
    PS::F64 kin;
    PS::F64 pot;
    PS::F64 pot_planet;
    PS::F64 tot;
    PS::F64 disp_merge;
    void clear(){
	kin = pot = tot = 0.0;
    }
};

class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    Energy eng_init;
    Energy eng_now;
    FileHeader(){
        n_body = 0;
        time = 0.0;
        eng_init.clear();
        eng_now.clear();
    }
    FileHeader(const PS::S64 n, const PS::F64 t, const Energy & e_i, const Energy & e_n){
        n_body = n;
        time = t;
	eng_init = e_i;
	eng_now = e_n;
    }
    PS::S32 readAscii(FILE * fp){
        fscanf(fp, "%lld%lf %lf%lf%lf%lf%lf %lf%lf%lf%lf%lf\n", 
	       &n_body, &time, 
	       &eng_init.kin, &eng_init.pot, &eng_init.pot_planet, &eng_init.tot, &eng_init.disp_merge,
	       &eng_now.kin,  &eng_now.pot,  &eng_now.pot_planet,  &eng_now.tot,  &eng_now.disp_merge);
	std::cout<<"n_body="<<n_body<<" time="<<time<<std::endl;
        return n_body;
    }
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%lf\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", 
		n_body, time, 
		eng_init.kin, eng_init.pot, eng_init.pot_planet, eng_init.tot, eng_init.disp_merge,
		eng_now.kin,  eng_now.pot,  eng_now.pot_planet,  eng_now.tot,  eng_now.disp_merge);
    }

};

int main(int argc, char *argv[]){
    PS::Initialize(argc, argv);
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    char sout[1024];
    sprintf(sout, "n_partilce.dat");
    std::ofstream fout;
    if(PS::Comm::getRank() == 0){
	fout.open(sout);
    }

    const PS::F64 PI = 4.0*atan(1.0);
    PS::F64 r_search_factor = 0.01875;
    PS::F64 dt = 1.0/1024.0;
    //const PS::F64 r_search_factor = 0.3;
    //const PS::F64 dt = 1.0/64.0;
    //const PS::F64 r_search_factor = 0.6;
    //const PS::F64 dt = 1.0/32.0;

    PS::ParticleSystem<FP> system;
    system.initialize();
    PS::S32 n_loc = 0;
    PS::F64 t_sys = 0.0;
#if 0
    PS::S64 n_glb_1d = 500;
    PS::S64 n_glb = n_glb_1d*n_glb_1d*n_glb_1d;
    PS::F64 vol = 1.0*1.0*1.0;
    FP::r_search = 1.0 / n_glb_1d * 0.3;
    SetParticleUniformBox(system, n_glb, n_loc, t_sys);
#else
    const PS::F64 mass_sun = 1.0;
    const PS::F64vec pos_sun = 0.0;
    const PS::F64vec vel_sun = 0.0;
#ifdef READ_FILE
    assert(argc == 2);
    char sinput[2048];
    sprintf(sinput, argv[1]);
    if(PS::Comm::getRank() == 0){
	std::cerr<<"sinput="<<sinput<<std::endl;
    }
    FileHeader file_header;
    system.readParticleAscii(sinput, file_header);
    PS::Comm::broadcast(&file_header, 1, 0);
    t_sys = file_header.time;
    //PS::Comm::broadcast(&t_sys, 1, 0);
    PS::S64 n_glb = system.getNumberOfParticleGlobal();
    n_loc = system.getNumberOfParticleLocal();
    std::cout<<"n_glb="<<n_glb<<" n_loc="<<n_loc<<std::endl;
    for(PS::S32 i=0; i<n_loc; i++){
	assert(system[i].mass > 0.0);
	assert(system[i].id >= 0);
	//system[i].n_ngb = 0; // this is critical . but why
    }
    PS::F64 scale_hight_loc = 0.0;
    PS::F64 ax_out_loc = 0.0;
    PS::F64 ax_in_loc = 999999.9; 
    PS::F64 ecc_rms_loc = 0.0;
    PS::F64 inc_rms_loc = 0.0;
    PS::F64 mass_max_loc = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
	PS::F64 ax, ecc, inc, OMG, omg, tperi;
	PosVel2OrbParam(ax, ecc, inc, OMG, omg, tperi,
			pos_sun, system[i].pos,
			vel_sun, system[i].vel,
			mass_sun, system[i].mass);
	ecc_rms_loc += ecc * ecc;
	inc_rms_loc += inc * inc;
	scale_hight_loc += system[i].pos.z*system[i].pos.z;
	if( ax_out_loc*ax_out_loc  < ax*ax) ax_out_loc = ax;
	if( ax_in_loc*ax_in_loc  > ax*ax)   ax_in_loc = ax;
	if(mass_max_loc < system[i].mass) mass_max_loc = system[i].mass;
    }
    const PS::F64 mass_max_glb = PS::Comm::getMaxValue(mass_max_loc);
    const PS::F64 ecc_rms = sqrt(PS::Comm::getSum(ecc_rms_loc)/n_glb);
    const PS::F64 inc_rms = sqrt(PS::Comm::getSum(inc_rms_loc)/n_glb);
    const PS::F64 scale_hight = sqrt( PS::Comm::getSum(scale_hight_loc) / n_glb);
    const PS::F64 ax_out = PS::Comm::getMaxValue(ax_out_loc);
    const PS::F64 ax_in = PS::Comm::getMinValue(ax_in_loc);
    const PS::F64 vel_kep = sqrt(ax_in*ax_in*ax_in/mass_sun);
    const PS::F64 vel_disp = vel_kep*sqrt(ecc_rms*ecc_rms+inc_rms*inc_rms);
    //const PS::F64 ratio_hill = cbrt( (2.0*system[0].mass) / (3.0*mass_sun) );
    const PS::F64 ratio_hill = cbrt( (2.0*mass_max_glb) / (3.0*mass_sun) );
    FP::r_search = r_search_factor*(ax_out*ratio_hill) + vel_disp*dt*3.0;
    PS::F64 vol = PI*(ax_out*ax_out-ax_in*ax_in)*scale_hight*2.0;
    if(PS::Comm::getRank() == 0){
	std::cout<<"scale_hight= "<<scale_hight<<std::endl;
	std::cout<<"ax_out= "<<ax_out<<" ax_in= "<<ax_in<<std::endl;
	std::cout<<"ecc_rms= "<<ecc_rms<<" inc_rms= "<<inc_rms<<std::endl;
	std::cout<<"ecc_rms/ratio_hill= "<<ecc_rms/ratio_hill
		 <<" inc_rms/ratio_hill= "<<inc_rms/ratio_hill<<std::endl;
	std::cout<<"FP::r_search= "<<FP::r_search<<" ax_out*ratio_hill= "<<ax_out*ratio_hill<<" vel_disp*dt*3.0= "<<vel_disp*dt*3.0<<std::endl;
	std::cout<<"vol= "<<vol<<std::endl;
    }




#else //no READ_FILE
    const PS::S64 n_glb = 1000000;
    const PS::F64 ax = 1.0;
    const PS::F64 ax_in = 0.95;
    const PS::F64 ax_out = 1.05;
    const PS::F64 ecc_sigma_hill = 2.0;
    const PS::F64 inc_sigma_hill = ecc_sigma_hill*0.5;
    if(PS::Comm::getRank() == 0){
	std::cout<<"ax_in="<<ax_in<<" ax_out="<<ax_out<<std::endl;
	std::cout<<"ecc_sigma_hill="<<ecc_sigma_hill<<" inc_sigma_hill="<<inc_sigma_hill<<std::endl;
    }
    SetParticleKeplerDisk(system, n_glb, n_loc, t_sys, ax_in, ax_out, ecc_sigma_hill, inc_sigma_hill);
    assert(system[0].mass>=0.0);
    const PS::F64 ratio_hill = cbrt( (2.0*system[0].mass) / (3.0*mass_sun) );
    if(PS::Comm::getRank() == 0){
	std::cout<<"system[0].mass="<<system[0].mass<<std::endl;
	std::cout<<"ratio_hill= "<<ratio_hill<<std::endl;
    }
    PS::F64 scale_hight_loc = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
	scale_hight_loc += system[i].pos.z*system[i].pos.z;
    }
    PS::F64 scale_hight = sqrt( PS::Comm::getSum(scale_hight_loc) / n_glb);
    if(PS::Comm::getRank() == 0){
	std::cout<<"scale_hight= "<<scale_hight<<std::endl;
	//std::cout<<"inc_sigma_hill*ratio_hill*2.0= "<<inc_sigma_hill*ratio_hill*2.0<<std::endl;
    }
    const PS::F64 vel_kep = sqrt(ax_in*ax_in*ax_in/mass_sun);
    PS::F64 ecc_rms_loc = 0.0;
    PS::F64 inc_rms_loc = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
	PS::F64 ax, ecc, inc, OMG, omg, tperi;
	PosVel2OrbParam(ax, ecc, inc, OMG, omg, tperi,
			pos_sun, system[i].pos,
			vel_sun, system[i].vel,
			mass_sun, system[i].mass);
	ecc_rms_loc += ecc * ecc;
	inc_rms_loc += inc * inc;
    }
    const PS::F64 ecc_rms = sqrt(PS::Comm::getSum(ecc_rms_loc)/n_glb);
    const PS::F64 inc_rms = sqrt(PS::Comm::getSum(inc_rms_loc)/n_glb);
    if(PS::Comm::getRank() == 0){
	std::cout<<"ecc_rms= "<<ecc_rms<<" inc_rms= "<<inc_rms<<std::endl;
	std::cout<<"ecc_rms/ratio_hill= "<<ecc_rms/ratio_hill
		 <<" inc_rms/ratio_hill= "<<inc_rms/ratio_hill<<std::endl;
    }
    const PS::F64 vel_disp = vel_kep*sqrt(ecc_rms*ecc_rms+inc_rms*inc_rms);
    FP::r_search = r_search_factor*(ax_out*ratio_hill) + vel_disp*dt*3.0;
    //if(PS::Comm::getRank() == 0){
    std::cout<<"FP::r_search= "<<FP::r_search<<" ax_out*ratio_hill= "<<ax_out*ratio_hill<<" vel_disp*dt*3.0= "<<vel_disp*dt*3.0<<std::endl;

    PS::F64 vol = PI*(ax_out*ax_out-ax_in*ax_in)*scale_hight*2.0;
#endif //READ_FILE

#endif // 0
    PS::DomainInfo dinfo;
    dinfo.initialize();
#if 0
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain(PS::F64vec(0.0, 0.0, 0.0), PS::F64vec(1.0, 1.0, 1.0));
#endif
    dinfo.decomposeDomainAll(system);

    PS::TreeForForceShort<FP, FP, FP>::Scatter tree;
    //PS::TreeForForceLong<FP, FP, FP>::MonopoleWithScatterSearch tree;
    tree.initialize(n_glb, 0.5, 8, 512);

    // NOTE: loop should be start 0, because r_search is expnaded twice each step
    PS::S32 loop_max = 10;
    if(n_glb < 1000000){ loop_max = 11; }
    for(PS::S32 loop=0; loop<loop_max; loop++){
	std::ofstream fout_log;
	char sout_log[1024];
	sprintf(sout_log, "log_%d_%d.dat", loop, PS::Comm::getRank());
	fout_log.open(sout_log);
	fout_log<<std::endl;
	fout_log<<"-------------------------"<<std::endl;
	fout_log<<"loop="<<loop<<std::endl;
	FP::r_search = r_search_factor*(ax_out*ratio_hill) + vel_disp*dt*3.0;
	fout_log<<"ax_out="<<ax_out<<std::endl;
	fout_log<<"ratio_hill="<<ratio_hill<<std::endl;
	fout_log<<"vel_disp="<<vel_disp<<std::endl;
	fout_log<<"FP::r_search="<<FP::r_search<<std::endl;
	system.exchangeParticle(dinfo);
	n_loc = system.getNumberOfParticleLocal();
	tree.calcForceAllAndWriteBack(CountNngb, system, dinfo, true);
	//tree.calcForceAllAndWriteBack(CalcNothing, system, dinfo, false); // just only copy the results in CountNngb to epj

#if 0
	if(loop == 9){
	    PS::S32 * n_ngb_di = new PS::S32[n_loc];
	    for(PS::S32 i=0; i<n_loc; i++){
		n_ngb_di[i] = 0;
		for(PS::S32 j=0; j<n_loc; j++){
		    if(i == j) continue;
		    if( (system[i].pos-system[j].pos)*(system[i].pos-system[j].pos) <= FP::r_search*FP::r_search){
			n_ngb_di[i]++;
		    }
		}
		if(n_ngb_di[i] > 1){
		    fout_log<<"system[i].id="<<system[i].id<<" n_ngb_di[i]="<<n_ngb_di[i]<<" system[i].n_ngb="<<system[i].n_ngb<<std::endl;
		}
	    }
	}
#endif
#if 0
	std::vector< std::pair<PS::S32, PS::S32> >idx_multi_loc;
	PS::S32 n_ptcl_loc = 0;
	PS::S32 n_cluster_0_loc = 0;
	PS::S32 n_cluster_1_loc = 0;
	for(PS::S32 i=0; i<n_loc; i++){
	    FP * nbl = NULL;
	    PS::S32 n_ngb = tree.getNeighborListOneParticle(system[i], nbl); // self is included
	    PS::S32 n_ngb_new = 0;
	    for(PS::S32 ii=0; ii<n_ngb; ii++){
		if( (nbl+ii)->id == system[i].id ) continue;
		else if( ((nbl+ii)->pos - system[i].pos)*((nbl+ii)->pos - system[i].pos) < FP::r_search*FP::r_search ){
		    n_ngb_new++;
		    idx_multi_loc.push_back( std::pair<PS::S32, PS::S32>(system[i].id, (nbl+ii)->id) );
		}
	    }
	    if(n_ngb_new == 0) n_cluster_0_loc++;
	    else if(n_ngb_new > 0) n_ptcl_loc++;
	}
#else
	// new
	std::vector< std::pair<PS::S32, PS::S32> >idx_multi_loc;
	std::unordered_map<PS::S32, PS::S32> idx_to_adr_loc; // point to adr of pbase
	PS::S32 n_ptcl_loc = 0;
	PS::S32 n_cluster_0_loc = 0;
	PS::S32 n_tmp_loc = 0;
	PS::S32 n_tmp2_loc = 0;
	std::vector<PS::S32> idx_no_neighbor;
	PS::S32 n_ngb_tot_loc = 0;
	for(PS::S32 i=0; i<n_loc; i++){
	    FP * nbl = NULL;
	    PS::S32 n_ngb = tree.getNeighborListOneParticle(system[i], nbl); // self is included
	    n_ngb_tot_loc += (n_ngb-1);
	    if(system[i].n_ngb == 0){
		n_tmp2_loc++;
	    }
	    if(n_ngb == 1){
		n_tmp_loc++;
	    }
	    PS::S32 n_ngb_new = 0;
	    for(PS::S32 ii=0; ii<n_ngb; ii++){
		if( (nbl+ii)->id == system[i].id ) continue;
		else if( ((nbl+ii)->pos - system[i].pos)*((nbl+ii)->pos - system[i].pos) <= FP::r_search*FP::r_search ){
		    if(n_ngb_new == 0){
			system[i].adr = idx_multi_loc.size();
			system[i].n_ngb = 1;
			idx_no_neighbor.push_back(system[i].id);
		    }
		    //system[i].n_ngb++;
		    n_ngb_new++;
		    idx_multi_loc.push_back( std::pair<PS::S32, PS::S32>(system[i].id, (nbl+ii)->id) );
		}
	    }
	    system[i].n_ngb = n_ngb_new;
	    if(n_ngb_new == 0) n_cluster_0_loc++;
	    else if(n_ngb_new > 0){
		//idx_to_adr_loc.insert( std::pair<PS::S32, PS::S32>(system[i].id, i) );
		idx_to_adr_loc.insert( std::pair<PS::S32, PS::S32>(system[i].id, n_ptcl_loc) );
		n_ptcl_loc++;
	    }
	}
	PS::S32 n_ngb_tot_glb = PS::Comm::getSum(n_ngb_tot_loc);
	fout_log<<"n_ngb_tot_glb="<<n_ngb_tot_glb<<std::endl;
	fout_log<<"average of # of ngbs (PS::F64)n_ngb_tot_glb/n_glb="<<(PS::F64)n_ngb_tot_glb/n_glb<<std::endl;

	fout_log<<"idx_multi_loc.size()="<<idx_multi_loc.size()<<std::endl;
	fout_log<<"idx_to_adr_loc.size()="<<idx_to_adr_loc.size()<<std::endl;
	fout_log<<"n_cluster_0_loc="<<n_cluster_0_loc<<std::endl;
	PS::S32 n_cluster_0_glb = PS::Comm::getSum(n_cluster_0_loc);
	fout_log<<"n_cluster_0_glb="<<n_cluster_0_glb<<std::endl;
	PS::S32 n_tmp_glb = PS::Comm::getSum(n_tmp_loc);
	fout_log<<"n_tmp_glb="<<n_tmp_glb<<std::endl;
	PS::S32 n_tmp2_glb = PS::Comm::getSum(n_tmp2_loc);
	fout_log<<"n_tmp2_glb="<<n_tmp2_glb<<std::endl;

	Pbase * pbase_loc = new Pbase[n_ptcl_loc];

	PS::S32 n_cnt = 0;
	PS::S32 ref_id = -1000;
	for(PS::S32 i=0; i<idx_multi_loc.size(); i++){
	    if(idx_multi_loc[i].first == ref_id){
		pbase_loc[n_cnt-1].n_ngb++;
	    }
	    if(idx_multi_loc[i].first != ref_id){
		pbase_loc[n_cnt].next = nullptr;
		pbase_loc[n_cnt].id = idx_multi_loc[i].first;
		pbase_loc[n_cnt].adr = i;
		pbase_loc[n_cnt].n_ngb = 1;
		pbase_loc[n_cnt].flag_search = false;
		ref_id = idx_multi_loc[i].first;
		n_cnt++;
	    }
	}
	fout_log<<"n_cnt="<<n_cnt<<" n_ptcl_loc="<<n_ptcl_loc<<std::endl;
	assert(n_cnt == n_ptcl_loc);

	Cluster * cluster_loc = new Cluster[n_ptcl_loc];
	PS::S32 n_cluster_loc = 0;
	std::vector< std::pair<PS::S32, PS::S32> >idx_multi_loc_send;
	std::unordered_map<PS::S32, PS::S32> idx_to_adr_loc_send;
	std::vector<Pbase> pbase_send;
	n_cnt = 0;
	for(PS::S32 i=0; i<n_ptcl_loc; i++){
	    bool fg_isolated = true;
	    if(pbase_loc[i].flag_search == false){
		PS::S32 list_len = 0;
		cluster_loc[n_cluster_loc].n_ptcl = 0;
		cluster_loc[n_cluster_loc].id = n_cluster_loc;
		//fout_log<<"pbase_loc+i->id="<<(pbase_loc+i)->id<<std::endl;
		//fout_log<<"idx_to_adr_loc[(pbase_loc+i)->id]="<<idx_to_adr_loc[(pbase_loc+i)->id]<<std::endl;
		//SearchClusterNgbWithCheck(pbase_loc+i, pbase_loc, pbase_loc+i, cluster_loc[n_cluster_loc], idx_to_adr_loc, idx_multi_loc, fg_isolated, list_len, fout_log);
		Pbase * p_top = pbase_loc + i;
		SearchClusterNgbWithCheck(pbase_loc+i, pbase_loc, p_top, cluster_loc[n_cluster_loc], idx_to_adr_loc, idx_multi_loc, fg_isolated, list_len, fout_log);
		//fout_log<<"list_len="<<list_len<<std::endl;
		//fout_log<<"fg_isolated="<<fg_isolated<<std::endl;
		if(fg_isolated){
		    //fout_log<<"isolated"<<std::endl;
		    n_cnt++;
		    n_cluster_loc++;
		}
		else{
		    //fout_log<<"conected"<<std::endl;
		    //fout_log<<"cluster_loc[n_cluster_loc].n_ptcl="<<cluster_loc[n_cluster_loc].n_ptcl<<std::endl;
#if 1
		    PS::S32 n_tmp2 = 0;
		    Pbase * pb = pbase_loc+i;
		    while(pb != nullptr){
			//fout_log<<"pb->id="<<pb->id<<std::endl;
			n_cnt++;
			n_tmp2++;
			PS::S32 adr = pb->adr;
			for(PS::S32 ii=0; ii<pb->n_ngb; ii++){
			    idx_multi_loc_send.push_back( idx_multi_loc[adr+ii] );
			}
			pb = pb->next;
		    }
		    //fout_log<<"n_tmp2="<<n_tmp2<<" cluster_loc[n_cluster_loc].n_ptcl="<<cluster_loc[n_cluster_loc].n_ptcl<<std::endl;
		    cluster_loc[n_cluster_loc].n_ptcl = 0;
#else
		    ClearFlag(pbase_loc+i, pbase_loc, pbase_loc+i, cluster_loc[n_cluster_loc], idx_to_adr_loc, idx_multi_loc, fg_isolated);
		    //fout_log<<"cluster_loc[n_cluster_loc].n_ptcl="<<cluster_loc[n_cluster_loc].n_ptcl<<std::endl;
		    PS::S32 adr = pbase_loc[i].adr;
		    for(PS::S32 ii=0; ii<pbase_loc[i].n_ngb; ii++){
			idx_multi_loc_send.push_back( idx_multi_loc[adr+ii] );
		    }
#endif
		}
		//fout_log<<std::endl;
	    }
	}
	fout_log<<"n_ptcl_loc="<<n_ptcl_loc<<" n_cnt="<<n_cnt<<std::endl;
	fout_log<<"idx_multi_loc_send.size()="<<idx_multi_loc_send.size()<<std::endl;
	fout_log<<"cluster_loc[0].n_ptcl="<<cluster_loc[0].n_ptcl<<std::endl;
	std::sort(cluster_loc, cluster_loc+n_cluster_loc, 
		  [](const Cluster & l, const Cluster & r)->bool{return l.n_ptcl < r.n_ptcl;} );
	/*
	for(PS::S32 i=1; i<n_cluster_loc; i++){
	    if(cluster_loc[i].n_ptcl == cluster_loc[i-1].n_ptcl){
	    }
	}
	*/

	std::vector< std::pair<PS::S32, PS::S32> > nptcl_ncluster;
	PS::S32 ref = -99999;
	for(PS::S32 i=0; i<n_cluster_loc; i++){
	    if(ref != cluster_loc[i].n_ptcl){
		nptcl_ncluster.push_back( std::pair<PS::S32, PS::S32>(cluster_loc[i].n_ptcl, 1) );
		ref = cluster_loc[i].n_ptcl;
	    }
	    else{
		nptcl_ncluster.back().second++;
	    }
	}
	fout_log<<"nptcl_ncluster.size()="<<nptcl_ncluster.size()<<std::endl;

	delete [] pbase_loc;
	delete [] cluster_loc;
	PS::S32 n_send = nptcl_ncluster.size();
	PS::S32 * n_recv_array = new PS::S32[PS::Comm::getNumberOfProc()];
	PS::Comm::gather(&n_send, 1, n_recv_array);
	PS::S32 * n_recv_disp = new PS::S32[PS::Comm::getNumberOfProc()+1];



	n_recv_disp[0] = 0;
	for(PS::S32 i=0; i<PS::Comm::getNumberOfProc(); i++){
	    n_recv_disp[i+1] = n_recv_disp[i] + n_recv_array[i];
	}
	PS::S32 n_recv_glb = n_recv_disp[PS::Comm::getNumberOfProc()];
	std::vector< std::pair<PS::S32, PS::S32> > nptcl_ncluster_array;
	nptcl_ncluster_array.resize(n_recv_glb);
	PS::Comm::gatherV(&nptcl_ncluster[0], n_send, &nptcl_ncluster_array[0], n_recv_array, n_recv_disp);
	std::sort(nptcl_ncluster_array.begin(), nptcl_ncluster_array.end(),
		  [](const std::pair<PS::S32, PS::S32> & l, const std::pair<PS::S32, PS::S32> & r)->bool{return l.first < r.first;} );
	/*
	if(PS::Comm::getRank() == 0){
	    for(const auto & i : nptcl_ncluster_array){
		fout_log<<"n_ptcl="<<i.first<<" n_cluster="<<i.second<<std::endl;
	    }
	}
	*/
	// here, we have coneted cluster isolated in own node.




	// send pair-index array
	/*
	for(const auto & i: idx_multi_loc_send){
	    fout_log<<"i.first="<<i.first<<" i.second="<<i.second<<std::endl;
	}
	*/
	n_send = idx_multi_loc_send.size();
	PS::Comm::gather(&n_send, 1, n_recv_array);
	n_recv_disp[0] = 0;
	for(PS::S32 i=0; i<PS::Comm::getNumberOfProc(); i++){
	    n_recv_disp[i+1] = n_recv_disp[i] + n_recv_array[i];
	}
	n_recv_glb = n_recv_disp[PS::Comm::getNumberOfProc()];
	std::vector< std::pair<PS::S32, PS::S32> > idx_multi_glb;
	idx_multi_glb.resize(n_recv_glb);
	PS::Comm::gatherV(&idx_multi_loc_send[0], n_send, &idx_multi_glb[0], n_recv_array, n_recv_disp);


	if(PS::Comm::getRank() == 0){
	    std::vector<Pbase> pbase_glb;
	    std::unordered_map<PS::S32, PS::S32>  idx_to_adr_glb;
	    n_cnt = 0;
	    ref = -10000;
	    for(PS::S32 i=0; i<n_recv_glb; i++){
		/*
		fout_log<<" idx_multi_glb[i].first="<<idx_multi_glb[i].first
			<<" idx_multi_glb[i].second="<<idx_multi_glb[i].second
			<<std::endl;
		*/
		if(idx_multi_glb[i].first == ref){
		    pbase_glb[n_cnt-1].n_ngb++;
		}
		if(idx_multi_glb[i].first != ref){
		    idx_to_adr_glb.insert( std::pair<PS::S32, PS::S32>(idx_multi_glb[i].first, n_cnt) );
		    /*
		    fout_log<<"n_cnt="<<n_cnt
			    <<" idx_multi_glb[i].first="<<idx_multi_glb[i].first
			    <<" idx_multi_glb[i].second="<<idx_multi_glb[i].second
			    <<std::endl;
		    */
		    pbase_glb.push_back( Pbase() );
		    pbase_glb[n_cnt].id = idx_multi_glb[i].first;
		    pbase_glb[n_cnt].adr = i;
		    pbase_glb[n_cnt].n_ngb = 1;
		    pbase_glb[n_cnt].flag_search = false;
		    ref = idx_multi_glb[i].first;
		    n_cnt++;
		}
	    }
	    std::vector<Cluster> cluster_glb;
	    PS::S32 n_ptcl_glb = pbase_glb.size();
	    PS::S32 n_cluster_glb = 0;
	    for(PS::S32 i=0; i<n_ptcl_glb; i++){
		bool fg_dummy = true;
		if(pbase_glb[i].flag_search == false){
		    cluster_glb.push_back(Cluster());
		    cluster_glb[n_cluster_glb].n_ptcl = 0;
		    cluster_glb[n_cluster_glb].id = n_cluster_glb;
		    PS::S32 list_len = 0;
		    Pbase * p_top = &pbase_glb[0];
		    //SearchClusterNgbWithCheck(&pbase_glb[i], &pbase_glb[0], &pbase_glb[i], cluster_glb[n_cluster_glb], idx_to_adr_glb, idx_multi_glb, fg_dummy, list_len, fout_log);
		    SearchClusterNgbWithCheck(&pbase_glb[i], &pbase_glb[0], p_top, cluster_glb[n_cluster_glb], idx_to_adr_glb, idx_multi_glb, fg_dummy, list_len, fout_log);
		    /*
		    if(cluster_glb[n_cluster_glb].n_ptcl < 2){
			fout_log<<"pbase_glb[i].id="<<pbase_glb[i].id
				<<" idx_to_adr_glb[pbase_glb[i].id]="<<idx_to_adr_glb[pbase_glb[i].id]
				<<" idx_multi_glb[pbase_glb[i].adr].first="<<idx_multi_glb[pbase_glb[i].adr].first
				<<" idx_multi_glb[pbase_glb[i].adr].second="<<idx_multi_glb[pbase_glb[i].adr].second
				<<" cluster_glb[n_cluster_glb].n_ptcl="<<cluster_glb[n_cluster_glb].n_ptcl<<std::endl;
		    }
		    */
		    Pbase * pb = &pbase_glb[i];
		    while(pb != nullptr){
			pb = pb->next;
		    }
		    n_cluster_glb++;
		}
	    }
	    fout_log<<"n_cluster_glb="<<n_cluster_glb<<std::endl;
	    /*
	    for(PS::S32 i=0; i<n_cluster_glb; i++){
		fout_log<<"cluster_glb[i].id="<<cluster_glb[i].id<<" cluster_glb[i].n_ptcl="<<cluster_glb[i].n_ptcl<<std::endl;
	    }
	    */
	    std::sort(cluster_glb.begin(), cluster_glb.end(),
		      [](const Cluster & l, const Cluster & r)->bool{return l.n_ptcl < r.n_ptcl;} );

	    std::vector< std::pair<PS::S32, PS::S32> > nptcl_ncluster_glb;
	    ref = -99999;
	    for(PS::S32 i=0; i<n_cluster_glb; i++){
		if(ref != cluster_glb[i].n_ptcl){
		    nptcl_ncluster_glb.push_back( std::pair<PS::S32, PS::S32>(cluster_glb[i].n_ptcl, 1) );
		    ref = cluster_glb[i].n_ptcl;
		}
		else{
		    nptcl_ncluster_glb.back().second++;
		}
	    }
	    for(const auto & i: nptcl_ncluster_array){
		nptcl_ncluster_glb.push_back(i);
	    }

	    std::sort(nptcl_ncluster_glb.begin(), nptcl_ncluster_glb.end(), 
		      [](const std::pair<PS::S32, PS::S32> & l, const std::pair<PS::S32, PS::S32> & r)->bool{return l.first < r.first;} );

	    ref = -1000;
	    PS::S32 n_ptcl_tot = 0;
	    std::vector< std::pair<PS::S32, PS::S32> > nptcl_ncluster_new;
	    nptcl_ncluster_new.push_back( std::pair<PS::S32, PS::S32>(1, n_cluster_0_glb) );
	    n_ptcl_tot += 1*n_cluster_0_glb;
	    for(const auto & i: nptcl_ncluster_glb){
		//fout_log<<"i.first="<<i.first<<" i.second="<<i.second<<std::endl;
		n_ptcl_tot += i.first*i.second;
		if(i.first != ref){
		    nptcl_ncluster_new.push_back( std::pair<PS::S32, PS::S32>(i.first, i.second) );
		    ref = i.first;
		}
		else{
		    nptcl_ncluster_new.back().second += i.second;
		}
	    }
	    fout_log<<"n_ptcl_tot="<<n_ptcl_tot<<" n_glb="<<n_glb<<std::endl;

	    std::ofstream fout_cluster;
	    char sout_cluster[1024];
	    sprintf(sout_cluster, "cluster_%d.dat", loop);
	    fout_cluster.open(sout_cluster);

	    n_ptcl_tot = 0;
	    PS::S32 n_cluster_tot = 0;
	    for(const auto & i: nptcl_ncluster_new){
		n_ptcl_tot += i.first * i.second;
		n_cluster_tot += i.second;
	    }
	    fout_log<<"n_ptcl_tot="<<n_ptcl_tot<<std::endl;
	    assert(n_ptcl_tot == n_glb);
	    //PS::S32 n_ptcl_max = std::max(nptcl_ncluster_new.back().first, 1); // 1 means no neighbors. in this case, nptcl_ncluster_new dose not used.
	    PS::S32 n_ptcl_max = nptcl_ncluster_new.back().first;
	    fout_log<<"n_ptcl_max="<<n_ptcl_max<<std::endl;

	    n_cnt = 0;

	    PS::S32 n_ptcl_cum = 0;
	    PS::S32 n_cluster_cum = 0;

#if 0
	    // fill 0 for clusters with no particles
	    for(PS::S32 i=0; i<n_ptcl_max+1; i++){
		PS::S32 n_cluster_tmp = 0;
		PS::S32 n_ptcl_tmp = 0;
		if(nptcl_ncluster_new[n_cnt].first == i){
		    n_cluster_tmp = nptcl_ncluster_new[n_cnt].second;
		    n_ptcl_tmp = nptcl_ncluster_new[n_cnt].second*i;
		    n_cnt++;
		}
		n_ptcl_cum += n_ptcl_tmp;
		n_cluster_cum += n_cluster_tmp;
		fout_cluster<<i<<"  "<<n_ptcl_tmp<<"  "<<n_ptcl_cum<<"  "<<(PS::F64)n_ptcl_tmp/n_ptcl_tot<<"  "<<(PS::F64)n_ptcl_cum/n_ptcl_tot
			    <<"  "<<n_cluster_tmp<<"   "<<n_cluster_cum<<"  "<<(PS::F64)n_cluster_tmp/n_cluster_tot<<"   "<<(PS::F64)n_cluster_cum/n_cluster_tot<<std::endl;
	    }
#else
	    fout_cluster<<std::scientific<<std::setprecision(6);
	    for(const auto & i: nptcl_ncluster_new){
		PS::S32 n_cluster_tmp = i.second;
		PS::S32 n_ptcl_tmp = i.first*i.second;
		if(n_cluster_tmp < 1) continue;
		n_cluster_cum += n_cluster_tmp;
		n_ptcl_cum += n_ptcl_tmp;
		fout_cluster<<i.first<<"  "<<n_ptcl_tmp<<"  "<<n_ptcl_cum<<"  "<<(PS::F64)n_ptcl_tmp/n_ptcl_tot<<"  "<<(PS::F64)n_ptcl_cum/n_ptcl_tot
			    <<"  "<<n_cluster_tmp<<"   "<<n_cluster_cum<<"  "<<(PS::F64)n_cluster_tmp/n_cluster_tot<<"   "<<(PS::F64)n_cluster_cum/n_cluster_tot<<std::endl;
	    }
#endif
	    fout_cluster.close();
	    PS::S32 n_ptcl_1_cluster = 0;
	    PS::S32 n_ptcl_2_cluster = 0;
	    for(const auto & i: nptcl_ncluster_new){
		if(i.first == 1) n_ptcl_1_cluster = i.first*i.second;
		else if(i.first == 2) n_ptcl_2_cluster = i.first*i.second;
	    }
#if 0
	    PS::F64 dens = n_glb / vol;
	    PS::F64 n_ngb_ave = dens * 4.0*PI/3.0 * (FP::r_search*FP::r_search*FP::r_search);
#else
	    PS::F64 n_ngb_ave = (PS::F64)n_ngb_tot_glb/n_glb;
	    PS::F64 dens = n_ngb_ave / (4.0*PI/3.0 * (FP::r_search*FP::r_search*FP::r_search));
	    vol = n_glb / dens;
	    PS::F64 scale_hight_prev = scale_hight;
	    PS::F64 scale_hight_new = vol / (PI*(ax_out*ax_out-ax_in*ax_in)*2.0);
	    fout<<"scale_hight_prev="<<scale_hight_prev<<" scale_hight_new="<<scale_hight_new<<std::endl;
#endif
#if 0
	    // 2D 
	    PS::F64 area = PI*(ax_out*ax_out-ax_in*ax_in);
	    PS::F64 dens_2d = n_glb / area;
	    n_ngb_ave = dens_2d * PI * FP::r_search * FP::r_search;
#endif
	    //fout<<"n_ngb_ave="<<n_ngb_ave<<"  (PS::F64)n_ngb_tot_glb/n_glb="<<(PS::F64)n_ngb_tot_glb/n_glb<<std::endl;
	    PS::F64 lambda = dens*PI*FP::r_search*FP::r_search*FP::r_search;
	    //PS::F64 poisson_1 = ((PS::F64)n_ngb_tot_glb/n_glb)*exp(-((PS::F64)n_ngb_tot_glb/n_glb));
	    PS::F64 poisson_0 = exp(-n_ngb_ave);
	    PS::F64 poisson_1 = n_ngb_ave*exp(-n_ngb_ave);
	    PS::F64 p2_pred = poisson_1 * 3.0/(lambda*lambda*lambda)*(2.0-exp(-lambda)*(lambda*lambda+2.0*lambda+2.0));
	    if(loop == 0){
		fout<<"#1:r_search(hill) 2:dt 3:r_search 4:n_glb 5:p1 6:p2 7:pm 8:vol(disk) 9:vol(sphere not considered overlap) 10: faction of vol";
		fout<<" 11:n_max 12:ave n_ngb 13:p1_p 14:p2_p"<<std::endl;
	    }
	    fout<<std::scientific<<std::setprecision(6);
	    fout<<r_search_factor  //1 
		<<"  "<<dt
		<<"  "<<FP::r_search
		<<"  "<<n_glb
		<<"  "<<(PS::F64)(n_ptcl_1_cluster)/n_ptcl_tot // fraction of particles in 1-cluster // line 5
		<<"  "<<(PS::F64)(n_ptcl_2_cluster)/n_ptcl_tot // fraction of particles in 2-cluster
		<<"  "<<(PS::F64)(n_ptcl_tot - n_ptcl_1_cluster - n_ptcl_2_cluster)/n_ptcl_tot // fraction of particles in other-cluster
		<<"  "<<vol // volume of disk
		<<"  "<<4.0/3.0*PI*FP::r_search*FP::r_search*FP::r_search*n_ptcl_tot // N*vol_neighbor_sphere
		<<"  "<<(4.0/3.0*PI*FP::r_search*FP::r_search*FP::r_search*n_ptcl_tot) / vol // volume fraction
		<<"  "<<n_ptcl_max
		<<"  "<<n_ngb_ave
		<<"  "<<poisson_0 // pred n_1
		<<"  "<<poisson_1 // pred n_2
		<<"  "<<p2_pred // pred n_2
		<<std::endl;

	}
	delete [] n_recv_array;
	delete [] n_recv_disp;
#endif

	fout_log.close();
	r_search_factor *= 2.0;
	dt *= 2.0;
    }

    PS::Finalize();
    return 0;
}
