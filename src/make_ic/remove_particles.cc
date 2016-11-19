#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<unistd.h>
#include<unordered_map>
#include<particle_simulator.hpp>

bool GetFlagSnpWithEnergy(const char sinput[]){
    bool flag = true;
    if(PS::Comm::getRank() == 0){
	std::ifstream fin;
	fin.open(sinput);
	std::cout<<"sinput:"<<sinput<<std::endl;
	std::string line;
	getline(fin, line);
	std::stringstream ss(line);
	std::string str;
	PS::S32 ncnt = 0;
	while(ss>>str){
	    std::cout<<"str: "<<str<<std::endl;
	    ncnt++;
	}
	std::cout<<"line="<<line<<std::endl;
	std::cout<<"ncnt="<<ncnt<<std::endl;
	if(ncnt == 2){ flag = false; }
	else if(ncnt == 12){ flag = true; }
	else{
	    std::cerr<<"input file is wrong format"<<std::endl;
	    PS::Abort();
	}
    }
    PS::Comm::broadcast(&flag, 1, 0);
    return flag;
}

class Energy{
public:
    PS::F64 kin;
    PS::F64 pot;
    PS::F64 pot_planet;
    PS::F64 tot;
    PS::F64 disp_merge;
    PS::F64 disp_aero;
    Energy(){
        kin = pot = tot = disp_merge = pot_planet = disp_aero = 0.0;
    }
    void clear(){
        kin = pot = tot = disp_merge = pot_planet = disp_aero = 0.0;
    }
    void dump(std::ostream & fout=std::cout){
        fout<<"tot="<<tot<<" kin+pot="<<kin+pot<<" kin="<<kin<<" pot="<<pot
	    <<" disp_merge="<<disp_merge<<" pot_planet="<<pot_planet
	    <<" disp_aero"<<disp_aero<<std::endl;
    }
};

class FP{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64 pot;
    PS::S32 n_ngb;
    PS::S32 id_ngb;
    static PS::F64 r_search;
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & p) { pos = p; }
    void copyFromForce(const FP & force){
        pot = force.pot;
	n_ngb = force.n_ngb;
	id_ngb = force.id_ngb;
    }
    void copyFromFP(const FP & p){
	id = p.id;
	mass = p.mass;
	pos = p.pos;
	vel = p.vel;
	n_ngb = p.n_ngb;
	id_ngb = p.id_ngb;
    }
    void clear(){
	n_ngb = 0;
	id_ngb = -1;
    }
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %d \n", 
                this->id, this->mass, 
		this->pos.x, this->pos.y, this->pos.z,  // 3-5
		this->vel.x, this->vel.y, this->vel.z,  // 6-8
		this->acc.x, this->acc.y, this->acc.z,  // 9-11
		this->pot, this->n_ngb);
    }

    void readAscii(FILE* fp) {
        fscanf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d \n", 
                &this->id, &this->mass, 
		&this->pos.x, &this->pos.y, &this->pos.z,  // 3-5
		&this->vel.x, &this->vel.y, &this->vel.z,  // 6-8
		&this->acc.x, &this->acc.y, &this->acc.z,  // 9-11
		&this->pot, &this->n_ngb);
    }

    void dump(std::ofstream & fout){
	fout<<"id= "<<id<<std::endl;
	fout<<"mass= "<<mass<<std::endl;
	fout<<"pos= "<<pos<<std::endl;
	fout<<"vel= "<<vel<<std::endl;
	fout<<"acc= "<<acc<<std::endl;
    }
    PS::F64 getRSearch() const {
	return r_search;
    }
};

PS::F64 FP::r_search;

class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    FileHeader(){
        n_body = 0;
        time = 0.0;
    }
    FileHeader(const PS::S64 n, const PS::F64 t){
        n_body = n;
        time = t;
    }
    PS::S32 readAscii(FILE * fp){
        fscanf(fp, "%lld\t%lf\n", &n_body, &time);
	std::cout<<"n_body="<<n_body<<" time="<<time<<std::endl;
        return n_body;
    }
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%lf\n", n_body, time);
    }
};

class FileHeaderWithEnergy{
public:
    PS::S64 n_body;
    PS::F64 time;
    Energy eng_init;
    Energy eng_now;
    FileHeaderWithEnergy(){
        n_body = 0;
        time = 0.0;
        eng_init.clear();
        eng_now.clear();
    }
    FileHeaderWithEnergy(const PS::S64 n, const PS::F64 t, const Energy & e_i, const Energy & e_n){
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

template<class Tsys>
void ReadFile(const char * sinput, 
	      Tsys & system){
    const bool flag_snp_with_energy = GetFlagSnpWithEnergy(sinput);
    if(flag_snp_with_energy){
	FileHeaderWithEnergy file_header_read;
	system.readParticleAscii(sinput, file_header_read);
    }
    else{
	FileHeader file_header_read;
	system.readParticleAscii(sinput, file_header_read);
    }
}

template<class Tsys>
void WriteFile(const char * soutput,
	       Tsys & system){
    FileHeader header(system.getNumberOfParticleGlobal(), 0.0);
    system.writeParticleAscii(soutput, header);
}

template<class Tsys>
PS::F64 GetRhill(Tsys & system){
    PS::F64 mass_pla = system[0].mass;
    PS::F64 mass_sun = 1.0;
    PS::F64 r_hill = cbrt( (2.0*mass_pla)/(3.0*mass_sun) );
    return r_hill;
}

// INTERACTION FUNCTION //
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
		PS::F64vec vij = fp_i[i].vel - fp_j[j].vel;
		PS::F64 eng = 0.5*vij*vij - fp_j[j].mass/sqrt(rij*rij);
		if(eng < 0.0){
		    force[i].n_ngb++;
		    force[i].id_ngb = fp_j[j].id;
		}
	    }
	}
    }
}

template<class Tsys>
void RearrngeParticle(Tsys & system){
    PS::S32 n_loc_new = system.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc_new; i++){
	if(system[i].n_ngb > 0){
	    system[i] = system[n_loc_new-1];
	    n_loc_new--;
	    i--;
	}
    }
    system.setNumberOfParticleLocal(n_loc_new);
}

int main(int argc, char *argv[]){
    PS::Initialize(argc, argv);
    char sinput[2048];
    char soutput[2048];
    int c;
    while((c=getopt(argc,argv,"i:o:h")) != -1){
        switch(c){
        case 'i':
            sprintf(sinput, optarg);
            std::cerr<<"sinput="<<sinput<<std::endl;
            break;
        case 'o':
            sprintf(soutput, optarg);
            std::cerr<<"soutput="<<soutput<<std::endl;
            break;
	}
    }
    PS::ParticleSystem<FP> system;
    system.initialize();
    const PS::S32 n_smp_ave = 100;
    system.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);
    ReadFile(sinput, system);
    FP::r_search = GetRhill(system);
    std::cout<<"FP::r_search="<<FP::r_search<<std::endl;
    PS::DomainInfo dinfo;
    dinfo.initialize();
    dinfo.decomposeDomainAll(system);

    system.exchangeParticle(dinfo);

    PS::TreeForForceShort<FP, FP, FP>::Scatter tree;
    tree.initialize(system.getNumberOfParticleGlobal(), 1.0, 8, 512);
    PS::Comm::barrier();
    tree.calcForceAllAndWriteBack(CountNngb, system, dinfo, true);
    RearrngeParticle(system);
    WriteFile(soutput, system);
    PS::Finalize();
    return 0;
}
