#ifdef P3T_64BIT
#define CALC_EP_64bit
#define CALC_SP_64bit
#define RSQRT_NR_EPJ_X4
#define RSQRT_NR_SPJ_X4
#elif P3T_MIXBIT
#define CALC_EP_64bit
#define RSQRT_NR_EPJ_X4
#else
#define RSQRT_NR_EPJ_X2
#endif //P3T_64BIT

#ifdef ENERGY_CHECK
//#define ENERGY_DIRECT
#endif

#ifdef FORCE_CHECK
#define FORCE_DIRECT
#endif

#if defined(INTRINSIC_K) || defined(INTRINSIC_X86)
#define INTRINSIC
#endif

#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<unistd.h>
#ifdef USE_C03
#include<map>
#else //USE_C03
#include<unordered_map>
#endif //USE_C03
#include<particle_simulator.hpp>
#include"class.hpp"
#include"hard.hpp"
#include"kepler.hpp"
#include"io.hpp"
#include"profile.hpp"
#include"domain.hpp"

template<class Tsys>
void CalcOneBodyEnergy(const Tsys & system, 
		       const PS::F64 mass_sun,
		       const PS::F64vec & pos_sun,
		       const PS::F64vec & vel_sun,
		       Energy & eng){
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
	if(system[i].n_ngb <= 0){
	    PS::F64 m = system[i].mass;
	    PS::F64vec dr = system[i].pos - pos_sun;
	    PS::F64vec dv = system[i].vel - vel_sun;
	    eng.kin += 0.5*m*dv*dv;
	    eng.pot -= m*mass_sun*sqrt(1.0 / (dr*dr));
	}
    }
    eng.tot = eng.pot + eng.kin;
}

PS::F64 GetQuantizedValue(const PS::F64 & val){
    static const PS::F64 inv_log2 = 1.0 / log(2.0);
    PS::F64 log2val = log(val) * inv_log2;
    log2val = (log2val > 0.0) ? log2val : log2val - 1.0;
    PS::S32 power = (PS::S32)(log2val);
    return pow(2.0, (PS::F64)(power));

}

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

template<class Tpsys>
PS::F64 GetRootFullLenght(const Tpsys & psys, const PS::F64vec & cen){
    PS::S64 nloc = psys.getNumberOfParticleLocal();
    PS::F64 len_loc_max = 0.0;
    for(PS::S32 i=0; i<nloc; i++){
	PS::F64vec dr = psys[i].pos - cen;
	for(PS::S32 k=0; k<3; k++){
	    if(len_loc_max < dr[k]) len_loc_max = dr[k];
	}
    }
    return 2.1*fabs(PS::Comm::getMaxValue(len_loc_max));
}

template<class T>
void Print(const T str, std::ostream & fout){
#ifdef DEBUG_PRINT_PLANET
    fout<<str<<std::endl;
#endif
}

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
	//std::cout<<"eng_init.tot="<<eng_init.tot<<std::endl;
	//std::cout<<"eng_now.tot="<<eng_now.tot<<std::endl;
    }
    PS::S32 readAscii(FILE * fp){
        fscanf(fp, "%lld%lf %lf%lf%lf%lf%lf %lf%lf%lf%lf%lf\n", 
	       &n_body, &time, 
	       &eng_init.kin, &eng_init.pot, &eng_init.pot_planet, &eng_init.tot, &eng_init.disp_merge,
	       &eng_now.kin,  &eng_now.pot,  &eng_now.pot_planet,  &eng_now.tot,  &eng_now.disp_merge);
	std::cout<<"n_body="<<n_body<<" time="<<time<<std::endl;
        return n_body;
    }
    /*
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%lf\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 
		n_body, time, 
		eng_init.kin, eng_init.pot, eng_init.pot_planet, eng_init.tot, eng_init.disp,
		eng_now.kin,  eng_now.pot,  eng_now.pot_planet,  eng_now.tot,  eng_now.disp);
    }
    */
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%lf\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", 
		n_body, time, 
		eng_init.kin, eng_init.pot, eng_init.pot_planet, eng_init.tot, eng_init.disp_merge,
		eng_now.kin,  eng_now.pot,  eng_now.pot_planet,  eng_now.tot,  eng_now.disp_merge);
    }

};

template<class Tpsys>
PS::F64 GetMassMax(const Tpsys & system, const PS::S64 n){
    PS::F64 m_max_loc = -1.0;
    for(PS::S64 i=0; i<n; i++){
        if(m_max_loc < system[i].mass) m_max_loc = system[i].mass;
    }
    return PS::Comm::getMaxValue(m_max_loc);
}

template<class Tpsys>
void CalcEng(const Tpsys & psys, Energy & eng_glb, bool clear=true){
    Energy eng_loc;
    if(clear){
        eng_loc.clear();
        eng_glb.clear();
    }
    const PS::S32 n = psys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
        eng_loc.pot += 0.5 * psys[i].mass * psys[i].pot_tot;
        eng_loc.kin += 0.5 * psys[i].mass * psys[i].vel * psys[i].vel;
    }
    eng_glb.kin += PS::Comm::getSum(eng_loc.kin);
    eng_glb.pot += PS::Comm::getSum(eng_loc.pot);
    eng_glb.tot = eng_glb.kin + eng_glb.pot;
}

template<class Tpsys>
void CalcEngKepler(const Tpsys & psys, 
		   Energy & eng_glb, 
		   PS::F64 m_sun, 
		   bool clear=true){
    Energy eng_loc;
    if(clear){
        eng_loc.clear();
        eng_glb.clear();
    }
    const PS::S32 n = psys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
        eng_loc.pot_planet += 0.5 * psys[i].mass * psys[i].pot_tot;
        eng_loc.kin += 0.5 * psys[i].mass * psys[i].vel * psys[i].vel;
        eng_loc.pot -= psys[i].mass * m_sun / sqrt(psys[i].pos * psys[i].pos);
    }
    eng_loc.pot += eng_loc.pot_planet;
    eng_glb.pot_planet += PS::Comm::getSum(eng_loc.pot_planet);
    eng_glb.kin += PS::Comm::getSum(eng_loc.kin);
    eng_glb.pot += PS::Comm::getSum(eng_loc.pot);
    eng_glb.tot = eng_glb.kin + eng_glb.pot;
}

template<class Tpsys>
void CalcEngKepler(const Tpsys & psys, 
		   const PS::F64 eng_disp_loc,
                   Energy & eng_glb, 
		   PS::F64 m_sun, 
		   bool clear=true){
    PS::F64 eng_disp_cum_glb = eng_glb.disp_merge;
    Energy eng_loc;
    if(clear){
        eng_loc.clear();
        eng_glb.clear();
    }
    const PS::S32 n = psys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
        eng_loc.pot_planet += 0.5 * psys[i].mass * psys[i].pot_tot;
        eng_loc.kin += 0.5 * psys[i].mass * psys[i].vel * psys[i].vel;
        eng_loc.pot -= psys[i].mass * m_sun / sqrt(psys[i].pos * psys[i].pos);
    }
    eng_loc.pot += eng_loc.pot_planet;
    eng_loc.disp_merge = eng_disp_loc;
    eng_glb.pot_planet += PS::Comm::getSum(eng_loc.pot_planet);
    eng_glb.kin += PS::Comm::getSum(eng_loc.kin);
    eng_glb.pot += PS::Comm::getSum(eng_loc.pot);
    eng_glb.disp_merge += PS::Comm::getSum(eng_loc.disp_merge) + eng_disp_cum_glb;
    eng_glb.tot = eng_glb.kin + eng_glb.pot + eng_glb.disp_merge;
}


template<class Tpsys>
void CalcEngKepler(const Tpsys & psys, 
		   const PS::F64 eng_disp_merge_loc,
		   const PS::F64 eng_disp_aero_loc,
                   Energy & eng_glb, 
		   PS::F64 m_sun, 
		   bool clear=true){
    PS::F64 eng_disp_merge_cum_glb = eng_glb.disp_merge;
    PS::F64 eng_disp_aero_cum_glb = eng_glb.disp_aero;
    Energy eng_loc;
    if(clear){
        eng_loc.clear();
        eng_glb.clear();
    }
    const PS::S32 n = psys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
        eng_loc.pot_planet += 0.5 * psys[i].mass * psys[i].pot_tot;
        eng_loc.kin += 0.5 * psys[i].mass * psys[i].vel * psys[i].vel;
        eng_loc.pot -= psys[i].mass * m_sun / sqrt(psys[i].pos * psys[i].pos);
    }
    eng_loc.pot += eng_loc.pot_planet;
    eng_loc.disp_merge = eng_disp_merge_loc;
    eng_loc.disp_aero  = eng_disp_aero_loc;
    eng_glb.pot_planet += PS::Comm::getSum(eng_loc.pot_planet);
    eng_glb.kin += PS::Comm::getSum(eng_loc.kin);
    eng_glb.pot += PS::Comm::getSum(eng_loc.pot);
    eng_glb.disp_merge += PS::Comm::getSum(eng_loc.disp_merge) + eng_disp_merge_cum_glb;
    eng_glb.disp_aero += PS::Comm::getSum(eng_loc.disp_aero) + eng_disp_aero_cum_glb;
    eng_glb.tot = eng_glb.kin + eng_glb.pot + eng_glb.disp_merge + eng_glb.disp_aero;
}

template<class Tpsys, class Ttree, class Tforce, class Tforce_func>
void CalcForceDirectFromAllPlanet(const Tpsys & psys, 
				  const PS::DomainInfo & dinfo,
				  Ttree & tree,
				  Tforce force[],
				  Tforce_func func){

    const PS::S32 n = psys.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n; i++){
	force[i].clear();
    }
    tree.calcForceDirect(func, force, dinfo, true);
}

#ifdef ENERGY_DIRECT
template<class Tpsys, class Ttree>
void CalcEngKeplerDirect(const Tpsys & psys, 
			 const PS::DomainInfo & dinfo,
			 Ttree & tree,
			 const PS::F64 eng_disp_loc, 
			 Energy & eng_glb, 
			 PS::F64 m_sun, 
			 bool clear=true){
    PS::F64 eng_disp_cum_glb = eng_glb.disp_merge;
    Energy eng_loc;
    if(clear){
        eng_loc.clear();
        eng_glb.clear();
    }
    const PS::S32 n = psys.getNumberOfParticleLocal();
    ForceSoft * force = new ForceSoft[n+1024];
    CalcForceDirectFromAllPlanet(psys, dinfo, tree, force, CalcForceEPEPNoCutoff64bit());
    /*
    for(PS::S32 i=0; i<n+1024; i++){
	force[i].clear();
    }
    tree.calcForceDirect(CalcForceEPEPNoCutoff64bit(), force, dinfo, true);
    */
    for(PS::S32 i=0; i<n; i++){
        eng_loc.pot_planet += 0.5 * psys[i].mass * force[i].pot;
        eng_loc.kin += 0.5 * psys[i].mass * psys[i].vel * psys[i].vel;
        eng_loc.pot -= psys[i].mass * m_sun / sqrt(psys[i].pos * psys[i].pos);
    }
    eng_loc.pot += eng_loc.pot_planet;
    eng_loc.disp_merge = eng_disp_loc;
    eng_glb.pot_planet += PS::Comm::getSum(eng_loc.pot_planet);
    eng_glb.kin += PS::Comm::getSum(eng_loc.kin);
    eng_glb.pot += PS::Comm::getSum(eng_loc.pot);
    eng_glb.disp_merge += PS::Comm::getSum(eng_loc.disp_merge) + eng_disp_cum_glb;
    eng_glb.tot = eng_glb.kin + eng_glb.pot + eng_glb.disp_merge;
    delete [] force;
}
#endif

template<class Tpsys, class Ttree>
void Kick(Tpsys & system,
          const Ttree & tree,
          const PS::F64 dt){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(int i=0; i<n; i++){
	system[i].vel  += system[i].acc * dt;
    }
}

template<class Tpsys, class Ttree>
void Drift(Tpsys & system,
           const Ttree & tree,
           const PS::F64 dt){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(int i=0; i<n; i++){
        //if(tree.getForce(i).n_ngb <= 0){
	if(system[i].n_ngb <= 0){
            system[i].pos  += system[i].vel * dt;
        }
    }
}

#if 1
template<class Tpsys, class Ttree>
void DriftKepler(Tpsys & system,
                 const Ttree & tree,
                 const PS::F64 dt,
                 const PS::F64 mass_sun = 1.0){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
	if(system[i].n_ngb <= 0){
	    PS::F64vec pos0 = 0.0;
	    PS::F64vec vel0 = 0.0;
	    const PS::F64 mass1 = 0.0;
	    PS::F64vec pos1 = system[i].pos;
	    PS::F64vec vel1 = system[i].vel;
	    DriveKepler(mass_sun, mass1, pos0, pos1, vel0, vel1, dt);
	    system[i].pos = pos1;;
	    system[i].vel = vel1;
        }
    }
}

template<class Tpsys, class Ttree>
void DriftKeplerDebug(Tpsys & system,
		      const Ttree & tree,
		      const PS::F64 dt,
		      const PS::F64 mass_sun,
		      std::ostream & fout){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
	if(system[i].n_ngb <= 0){
	    fout<<"i0= "<<i<<std::endl;
	    PS::F64vec pos0 = 0.0;
	    PS::F64vec vel0 = 0.0;
	    const PS::F64 mass1 = 0.0;
	    PS::F64vec & pos1 = system[i].pos;
	    PS::F64vec & vel1 = system[i].vel;
	    fout<<"i1= "<<i<<std::endl;
	    DriveKepler(mass_sun, mass1, pos0, pos1, vel0, vel1, dt);
	    fout<<"i2= "<<i<<std::endl;
        }
    }
}
#else

template<class Tpsys, class Ttree>
void DriftKepler(Tpsys & system,
                 const Ttree & tree,
                 const PS::F64 dt,
                 const PS::F64 mass_sun = 1.0){
    const PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i=0; i<n; i++){
	static __thread PS::F64vec pos0, vel0, pos1, vel1;
	static __thread PS::F64 mass1;
        //if(tree.getForce(i).n_ngb <= 0){
	if(system[i].n_ngb <= 0){
	    pos0 = 0.0;
	    vel0 = 0.0;
	    mass1 = 0.0;
	    pos1 = system[i].pos;
	    vel1 = system[i].vel;
	    DriveKepler(mass_sun, mass1, pos0, pos1, vel0, vel1, dt);
	    system[i].pos = pos1;
	    system[i].vel = vel1;
        }
    }
}

#endif

template<class Tpsys>
void SetRoutRinRhillVdisp(const Tpsys & system_soft,
                          const PS::F64vec pos_sun,
                          const PS::F64vec vel_sun,
                          const PS::F64 mass_sun,
                          const PS::F64 ratio_r_cut,
                          const PS::F64 ratio_r_search,
                          PS::F64 & r_out,
                          PS::F64 & r_in,
                          PS::F64 & r_hill_max_glb,
                          PS::F64 & vel_disp,
			  PS::F64 & ecc_rms,
			  PS::F64 & inc_rms,
                          PS::F64 & mass_planet_max_glb,
                          PS::F64 & mass_planet_tot_glb,
			  PS::F64 & Tkep_min_glb){
    static const PS::F64 PI = 4.0*atan(1.0);
    r_out = r_in = r_hill_max_glb = vel_disp = mass_planet_max_glb = mass_planet_tot_glb = 0.0;
    const PS::S32 n_loc = system_soft.getNumberOfParticleLocal();
    const PS::S32 n_glb = system_soft.getNumberOfParticleGlobal();
    PS::F64 mass_planet_max_loc = 0.0;
    PS::F64 mass_planet_tot_loc = 0.0;
    PS::F64 r_hill_max_loc = 0.0;
    PS::F64 v_kep_max_loc = 0.0;
    PS::F64 ecc_sq_tot_loc = 0.0;
    PS::F64 inc_sq_tot_loc = 0.0;
    PS::F64 Tkep_min_loc = PS::LARGE_FLOAT;
    for(PS::S32 i=0; i<n_loc; i++){
        PS::F64 ax, ecc, inc;
        PosVel2AxEccInc(ax, ecc, inc,
                        pos_sun, system_soft[i].pos,
                        vel_sun, system_soft[i].vel,
                        mass_sun,        PS::F64(0.0));
        ecc_sq_tot_loc += ecc * ecc;
        inc_sq_tot_loc += inc * inc;
        PS::F64 v_kep  = mass_sun / ax;
        PS::F64 r_hill = ax*ax*ax*(2.0*system_soft[i].mass) / (3.0*mass_sun);
	PS::F64 Tkep  = ax*ax*ax / (system_soft[i].mass+mass_sun);
        if(r_hill_max_loc < r_hill) r_hill_max_loc = r_hill;
        if(v_kep_max_loc < v_kep) v_kep_max_loc = v_kep;
        if(Tkep < Tkep_min_loc)  Tkep_min_loc = Tkep; 
        if(mass_planet_max_loc < system_soft[i].mass) mass_planet_max_loc = system_soft[i].mass;
        mass_planet_tot_loc += system_soft[i].mass;
    }
    Tkep_min_loc = 2.0*PI*sqrt(Tkep_min_loc);
    Tkep_min_glb = PS::Comm::getMinValue(Tkep_min_loc);
    mass_planet_max_glb = PS::Comm::getMaxValue(mass_planet_max_loc);
    mass_planet_tot_glb = PS::Comm::getSum(mass_planet_tot_loc);
    r_hill_max_loc = cbrt(r_hill_max_loc);
    r_hill_max_glb = PS::Comm::getMaxValue(r_hill_max_loc);
    v_kep_max_loc = sqrt(v_kep_max_loc);
    PS::F64 v_kep_max_glb = PS::Comm::getMaxValue(v_kep_max_loc);
    ecc_rms = sqrt( PS::Comm::getSum(ecc_sq_tot_loc) / n_glb );
    inc_rms = sqrt( PS::Comm::getSum(inc_sq_tot_loc) / n_glb );
    vel_disp = v_kep_max_glb * sqrt( ecc_rms*ecc_rms + inc_rms*inc_rms);
    r_out = r_hill_max_glb * ratio_r_search;
    r_in = r_out * ratio_r_cut;
    //EPISoft::r_out = r_out;
    //EPISoft::r_in = r_in;
}

// ref Ogihara, Kobayashi and Inutsuka 2014
class GasDisk{
private:
    PS::F64 alpha;
    PS::F64 beta;
    PS::F64 Cd;
    PS::F64 f_gas;
    PS::F64 tau_gas;
    PS::F64 rho_planet;
    PS::F64 mass_sun;
    PS::F64 lum_sun;
    PS::F64vec pos_sun;
    PS::F64vec vel_sun;
public:
    void initialize(const PS::F64 _rho_planet, // g/cm^3
		    const PS::F64 _mass_sun = 1.0,
		    const PS::F64 _lum_sun = 1.0,
		    const PS::F64vec _pos_sun = 0.0,
		    const PS::F64vec _vel_sun = 0.0){
	static const PS::F64 PI = 4.0*atan(1.0);
	alpha = 11.0 / 4.0; // exponent of rho
	beta = 0.5; // exponent of temperature
	Cd = 1.0; // coeficient aero drag
	f_gas = 1.0; // scaling factor of gas density
	tau_gas = 1e6 * 2.0 * PI; // 1e6[yr] * 2PI[T]/[yr]
	rho_planet = _rho_planet * 1.49597871e13 * 1.49597871e13 * 1.49597871e13 / 1.9884e33; // [Msun/AU^3]
	mass_sun = _mass_sun;
	lum_sun  = _lum_sun;
	pos_sun  = _pos_sun;
	vel_sun  = _vel_sun;
    }
    template<class Tsys>
    void calcAeroDrag(Tsys & sys,
		      const PS::F64 time){
	// this function must be call after calcForce from planet
	static bool first = true;
	static const PS::F64 PI = 4.0*atan(1.0);
	const PS::S64 n = sys.getNumberOfParticleLocal();
	const PS::F64 coef_rho_gas = exp(-time/tau_gas) * 2400.0 * f_gas * 1.49597871e13 * 1.49597871e13 / (0.047*2 * 1.9884e33) * pow(lum_sun, -0.125) * pow(mass_sun, 0.5);
	//const PS::F64 coef_rho_gas = 2400.0 * f_gas * 1.49597871e13 * 1.49597871e13 / (0.047*2 * 1.9884e33) * pow(lum_sun, -0.125) * pow(mass_sun, 0.5);
	const PS::F64 coef_cs_vk = 1.0/29.78 * pow(lum_sun, 0.125) * pow(mass_sun, -0.5); // about 0.0018 for solar luminosity and solar mass
	// 1.0/29.78 is cs/vk at 1AU (cs=1km/sec, vkep=29.78km/sec)
	const PS::F64 coef_acc_gd = 0.5*Cd*PI;
#pragma omp parallel for
	for(PS::S64 i=0; i<n; i++){
	    PS::F64 r_sq = sys[i].pos.x*sys[i].pos.x + sys[i].pos.y*sys[i].pos.y;
	    PS::F64 inv_r = 1.0 / sqrt(r_sq);
	    PS::F64 r = r_sq * inv_r;
	    PS::F64 r_sqrt = sqrt(r);
	    PS::F64 r_4rt = sqrt(r_sqrt);
	    PS::F64vec ev(-sys[i].pos.y*inv_r, sys[i].pos.x*inv_r, 0.0); // unit vector of kepler velocity
	    PS::F64vec vkep = sqrt(mass_sun * inv_r) * ev; // kepler velocity
	    //PS::F64 cs_vk = coef_cs_vk * r_sqrt;
	    PS::F64 cs_vk = coef_cs_vk * r_4rt;
	    PS::F64 eta = 0.5*(alpha+beta)*cs_vk*cs_vk;
	    PS::F64vec vgas = (1.0 - eta)*vkep;
	    PS::F64vec u = sys[i].vel - vgas;
	    PS::F64 rplanet = cbrt(3.0*sys[i].mass/(4.0*PI*rho_planet));
	    PS::F64 rho_gas = coef_rho_gas * inv_r * inv_r * inv_r * r_4rt;
#if 1
	    sys[i].acc_gd = (-coef_acc_gd * rplanet * rplanet * rho_gas * sqrt(u*u) * u) / sys[i].mass;
	    sys[i].acc += sys[i].acc_gd;
#else
              // ref Higuchi et al.
	    PS::F64 rho_gas_tmp = 2e-6*pow(r, (-11.0/4.0)); //[kg/m^3]
	    rho_gas_tmp = rho_gas_tmp * (1.496e11*1.496e11*1.496e11/1.989e30);
	    PS::F64 eta_tmp = 0.002*sqrt(r);
	    PS::F64vec u_tmp = sys[i].vel - (1.0-eta_tmp)*vkep;
	    sys[i].acc_gd = -Cd * PI * rplanet * rplanet * rho_gas_tmp * sqrt(u_tmp*u_tmp) * u_tmp / (2.0*sys[i].mass);
	    sys[i].acc += sys[i].acc_gd;
	    if(first){
		fout_debug<<"rho_gas="<<rho_gas<<" rho_gas_tmp="<<rho_gas_tmp<<std::endl;
		fout_debug<<"eta="<<eta<<" eta_tmp="<<eta_tmp<<std::endl;
		fout_debug<<"u="<<u<<" u_tmp="<<u_tmp<<std::endl;
		first = false;
	    }
#endif
	}
    }
};

int main(int argc, char *argv[]){
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);
#ifdef CALC_HARD_ENERGY
    PS::F64 dEerr_1body_loc = 0.0;
    PS::F64 dEerr_1body_glb = 0.0;
    PS::F64 dEerr_2body_loc = 0.0;
    PS::F64 dEerr_2body_glb = 0.0;
    PS::F64 dEerr_mbody_loc = 0.0;
    PS::F64 dEerr_mbody_glb = 0.0;
#endif
    PS::F64 r_merge_factor = 1.0;
    PS::F64 ratio_r_cut = 0.1;
    PS::F64 time_sys = 0.0;
    PS::F64 theta = 0.4;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::S32 n_smp_ave = 100;
    PS::F64 time_end = 1000.0 * 6.0; // 1000[yr] * 6[T/yr]
    PS::F64 r_cut_factor = 1.0; // use 1 hill
    PS::F64 dt_soft = 1.0/16.0; // about 1/(2pi)*1/16 = 0.01 [yr]
    PS::F64 eta = 0.1;
    char dir_name[1024];
    PS::S64 n_glb = 16384;
    PS::F64 mass_sun = 1.0; //[Msun]
    PS::F64vec pos_sun = 0.0;
    PS::F64vec vel_sun = 0.0;
    PS::F64 dt_snp = 100.0; // about 100/2pi = 16 [yr]
    PS::F64 search_factor = 3.0; // 3.0*vel_disp is width of the buffer shell
    PS::F64 dens_planet = 2.0; //[g/cm^3]
    PS::F64 dt_limit_hard_factor = 4.0;
    PS::S64 snp_id = 0;
    int c;
#ifdef READ_FILE
    char sinput[2048];
    while((c=getopt(argc,argv,"i:I:o:d:D:E:t:r:T:n:s:S:l:R:r:P:X:h")) != -1){
        switch(c){
        case 'i':
            sprintf(sinput, optarg);
            std::cerr<<"sinput="<<sinput<<std::endl;
            break;
        case 'I':
            snp_id = atoi(optarg);
            std::cerr<<"snp_id="<<snp_id<<std::endl;
            break;
        case 'o':
            sprintf(dir_name, optarg);
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
        case 'P':
            dens_planet = atof(optarg);
            std::cerr<<"dens_planet="<<dens_planet<<std::endl;
            break;
        case 'X':
	    dt_limit_hard_factor = atof(optarg);
            std::cerr<<"dt_limit_hard_factor="<<dt_limit_hard_factor<<std::endl;
            break;
        case 'h':
            std::cerr<<"i: input_file"<<std::endl;
            std::cerr<<"I: snp_id"<<std::endl;
            std::cerr<<"o: dir name of output"<<std::endl;
            std::cerr<<"d: inv_dt (dafult 16 ~ 0.01yr  )"<<std::endl;
            std::cerr<<"D: dt_snp (dafult 100 ~ 16yr  )"<<std::endl;
            std::cerr<<"E: eta (dafult 0.1)"<<std::endl;
            std::cerr<<"t: theta (dafult: 0.5)"<<std::endl;
            std::cerr<<"T: time_end (dafult: 6000 ~ 1000yr)"<<std::endl;
            std::cerr<<"n: n_group_limit (dafult: 64.0)"<<std::endl;
            std::cerr<<"s: n_smp_ave (dafult: 100)"<<std::endl;
            std::cerr<<"S: search_factor (dafult: 3.0)"<<std::endl;
            std::cerr<<"l: n_leaf_limit (dafult: 8)"<<std::endl;
            std::cerr<<"r: r_merge_factor (dafult: 1.0)"<<std::endl;
            std::cerr<<"R: r_cut_factor (dafult: 3.0)"<<std::endl;
            std::cerr<<"P: dens_planet [g/cm^3](dafult: 2.0)"<<std::endl;
            std::cerr<<"X: dt_limit_hard_factor(dafult: 4.0 -> dt_limit_hard = dt_soft/4.0)"<<std::endl;
            PS::Finalize();
            return 0;
        }
    }
#else // not READ_FILE
    PS::F64 ax_in = 0.98; //[AU]
    PS::F64 ax_out = 1.02; //[AU]
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
#endif//READ_FILE
    PTCLHard::r_factor = r_merge_factor;
    PTCLHard::dens = dens_planet;

    std::ofstream fout_tcal;
    char sout_tcal[1024];
    sprintf(sout_tcal, "%s/tcal_%05d.dat", dir_name, PS::Comm::getRank());
    fout_tcal.open(sout_tcal, std::ios::app);

    std::ofstream fout_merge;
    std::ofstream fout_diag;
    std::ofstream fout_log;
#if 0
    if(PS::Comm::getRank() == 0){
        char sout_merge[1024];
        sprintf(sout_merge, "%s/merge.dat", dir_name);
        fout_merge.open(sout_merge, std::ios::app);
    }
#else
    char sout_merge[1024];
    sprintf(sout_merge, "%s/merge_%05d.dat", dir_name, PS::Comm::getRank());
    fout_merge.open(sout_merge, std::ios::app);
    fout_merge<<std::setprecision(15);

    char sout_debug[1024];
    sprintf(sout_debug, "%s/debug_%05d.dat", dir_name, PS::Comm::getRank());
    fout_debug.open(sout_debug, std::ios::app);
    fout_debug<<std::setprecision(15);

    std::ofstream fout_domain;
    char sout_domain[1024];
    sprintf(sout_domain, "%s/domain_%05d.dat", dir_name, PS::Comm::getRank());
    fout_domain.open(sout_domain, std::ios::app);
    fout_domain<<std::setprecision(15);

    if(PS::Comm::getRank() == 0){
        char sout_diag[1024];
        sprintf(sout_diag, "%s/diag.dat", dir_name);
        fout_diag.open(sout_diag, std::ios::app);
        fout_diag<<std::setprecision(15);
        char sout_log[1024];
        sprintf(sout_log, "%s/log.dat", dir_name);
        fout_log.open(sout_log, std::ios::app);
        fout_log<<std::setprecision(15);
    }
#endif
    
    PS::ParticleSystem<FPSoft> system_soft;
    system_soft.initialize();
    system_soft.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);
    PS::S32 n_loc;
    Energy eng_init, eng_now;
#ifdef READ_FILE
    const bool flag_snp_with_energy = GetFlagSnpWithEnergy(sinput);
    Print("reading file", fout_debug);
    if(flag_snp_with_energy){
	Print("with energy", fout_debug);
	FileHeaderWithEnergy file_header_read;
	system_soft.readParticleAscii(sinput, file_header_read);
	time_sys = file_header_read.time;
	eng_init = file_header_read.eng_init;
	eng_now = file_header_read.eng_now;
	PS::Comm::broadcast(&time_sys, 1, 0);
	Print(file_header_read.eng_init.disp_merge, fout_debug);
	Print(file_header_read.eng_now.disp_merge,  fout_debug);
    }
    else{
	Print("without energy", fout_debug);
	FileHeader file_header_read;
	system_soft.readParticleAscii(sinput, file_header_read);
	time_sys = file_header_read.time;
	PS::Comm::broadcast(&time_sys, 1, 0);
    }
    Print("finish reading file", fout_debug);

    Print((std::string)("time_sys= "+std::to_string(time_sys)), fout_debug);
    n_glb = system_soft.getNumberOfParticleGlobal();
    Print((std::string)("n_glb= "+std::to_string(n_glb)), fout_debug);
    n_loc = system_soft.getNumberOfParticleLocal();
    Print((std::string)("n_loc= "+std::to_string(n_loc)), fout_debug);
#else //READ_FILE
    std::cerr<<"SetParticleKeplerDisk Beg"<<std::endl;
    /*
    SetParticleKeplerDisk(system_soft, n_glb, n_loc, time_sys,
			  ax_in, ax_out, ecc_sigma_hill, inc_sigma_hill, dens, mass_sun, seed);
    */
    SetParticleKeplerDisk(system_soft, n_glb, n_loc, time_sys,
			  ax_in, ax_out, ecc_sigma_hill, inc_sigma_hill, dens, mass_sun, ax_ice, f_ice, power, seed);
    std::cerr<<"SetParticleKeplerDisk Fin"<<std::endl;
#endif
#ifdef REV_VEL
    for(PS::S32 i=0; i<n_loc; i++){
	system_soft[i].vel *= -1.0;
    }
#endif
    for(PS::S32 i=0; i<10; i++){
	std::cout<<"system_soft[i].pos="<<system_soft[i].pos<<std::endl;
	std::cout<<"system_soft[i].vel="<<system_soft[i].vel<<std::endl;
    }
    PS::Comm::barrier();


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

    PS::F64 hight_rms_loc = 0.0;
    PS::F64 r_rms_loc = 0.0;
    for(PS::S32 i=0; i<n_loc; i++){
	hight_rms_loc += system_soft[i].pos.z * system_soft[i].pos.z;
	r_rms_loc += system_soft[i].pos * system_soft[i].pos;
    }
    PS::F64 hight_rms = sqrt(PS::Comm::getSum(hight_rms_loc) / n_glb);
    PS::F64 r_rms = sqrt(PS::Comm::getSum(r_rms_loc) / n_glb);

#ifdef MERGE
    EPISoft::eps = 0.0; // eps should be zero between the sun and planets, otherwise kepler solver doesn't work well
#else //MERGE
    static const PS::F64 rho_ave = 3.0 * ( (1.49597871e13*1.49597871e13*1.49597871e13) / 1.989e33); // [Msun/AU^3]
    static const PS::F64 PI = 4.0*atan(1.0);
    static const PS::F64 C = 3.0/(4.0*PI*rho_ave);
    EPISoft::eps = cbrt(C*mass_planet_tot_glb/system_soft.getNumberOfParticleGlobal()) * r_merge_factor;
    //EPISoft::eps = cbrt(C*system_soft[0].mass) * r_merge_factor;
#endif
    std::cerr<<"EPISoft::eps="<<EPISoft::eps<<std::endl;

    EPJSoft::r_search = r_out + search_factor*vel_disp*dt_soft;
    EPISoft::r_out = r_out;
    EPISoft::r_in = r_in;
    if(r_out <= EPISoft::eps){
        EPJSoft::r_search = 0.0;
    }

    PS::F64 unit_m = 1.989e30; //[kg]
    PS::F64 unit_v = 29.7886203575; //[km/sec]
    PS::F64 unit_t = 0.15924595715; //[yr]

    if(PS::Comm::getRank() == 0){
#ifdef READ_FILE
#else //READ_FILE
        fout_log<<"ax_in= "<<ax_in<<std::endl;
        fout_log<<"ax_out= "<<ax_out<<std::endl;
        fout_log<<"inc_sigma_hill= "<<inc_sigma_hill<<std::endl;
        fout_log<<"ecc_sigma_hill= "<<ecc_sigma_hill<<std::endl;
#endif //READ_FILE
        fout_log<<"r_out= "<<r_out<<std::endl;
        fout_log<<"r_in= "<<r_in<<std::endl;
        fout_log<<"vel_disp= "<<vel_disp<<std::endl;;
        fout_log<<"r_hill_max_glb= "<<r_hill_max_glb<<std::endl;
        fout_log<<"mass_planet_max_glb= "<<mass_planet_max_glb<<std::endl;
        fout_log<<"mass_planet_tot_glb= "<<mass_planet_tot_glb<<std::endl;
        fout_log<<"ecc_rms= "<<ecc_rms<<std::endl;
        fout_log<<"inc_rms= "<<inc_rms<<std::endl;
        fout_log<<"dir_name= "<<dir_name<<std::endl;
        fout_log<<"dt_soft= "<<dt_soft<<std::endl;
        fout_log<<"dt_snp= "<<dt_snp<<std::endl;
        fout_log<<"eta= "<<eta<<std::endl;
        fout_log<<"theta= "<<theta<<std::endl;
        fout_log<<"time_end= "<<time_end<<std::endl;
        fout_log<<"n_group_limit= "<<n_group_limit<<std::endl;
        fout_log<<"n_glb= "<<n_glb<<std::endl;
        fout_log<<"n_smp_ave= "<<n_smp_ave<<std::endl;
        fout_log<<"n_leaf_limit= "<<n_leaf_limit<<std::endl;
        fout_log<<"r_cut_factor= "<<r_cut_factor<<std::endl;
        fout_log<<"PTCLHard::r_factor= "<<PTCLHard::r_factor<<std::endl;
        fout_log<<"PTCLHard::dens= "<<PTCLHard::dens<<std::endl;
        fout_log<<"EPISoft::eps= "<<EPISoft::eps<<std::endl;
        fout_log<<"EPJSoft::r_search= "<<EPJSoft::r_search<<std::endl;
        fout_log<<"unit_m= "<<unit_m<<"[kg]"<<std::endl;
        fout_log<<"unit_v= "<<unit_v<<"[km/sec]"<<std::endl;
        fout_log<<"unit_t= "<<unit_t<<"[yr]"<<std::endl;
    }
#ifdef AERO_DRAG
    GasDisk gas_disk;
    PS::F64 lum_sun = 1.0;
    gas_disk.initialize(dens_planet, mass_sun, lum_sun, pos_sun, vel_sun);
    for(PS::S32 i=0; i<n_loc; i++){
	system_soft[i].acc_gd = system_soft[i].acc_gd_prev = system_soft[i].vel_prev = 0.0;
    }
#endif //AERO_DRAG

    const PS::F32 coef_ema = 0.2;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
#ifndef DIV_FIX
    PS::S32 nx = sqrt((PS::F64)n_proc-0.000001)+1;
    while( n_proc % nx != 0) nx++;
    PS::S32 ny = n_proc / nx;
    PS::S32 nz = 1;
    if(PS::Comm::getRank() == 0){
        fout_log<<"nx,ny,nz= "<<nx<<" "<<ny<<" "<<nz<<std::endl;
    }
    dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);
    dinfo.collectSampleParticle(system_soft, true);
    dinfo.decomposeDomain();
#endif //DIV_FIX

#ifdef DIV_FIX
    if(n_proc%4==0){
	DomainDecision(dinfo, system_soft);
    }
    else{
	PS::S32 nx = sqrt((PS::F64)n_proc-0.000001)+1;
	while( n_proc % nx != 0) nx++;
	PS::S32 ny = n_proc / nx;
	PS::S32 nz = 1;
	if(PS::Comm::getRank() == 0){
	    fout_log<<"nx,ny,nz= "<<nx<<" "<<ny<<" "<<nz<<std::endl;
	}
	dinfo.setNumberOfDomainMultiDimension(nx, ny, nz);
	dinfo.collectSampleParticle(system_soft, true);
	dinfo.decomposeDomain();
    }
    PS::Comm::barrier();
    if(PS::Comm::getRank() == 0){
        for(PS::S32 i=0; i<n_proc; i++){
            std::cout<<std::setprecision(15)<<i<<"   "<<dinfo.getPosDomain(i)<<std::endl;
        }
    }
#endif // DIV_FIX

    PS::F64ort * pos_domain = new PS::F64ort[n_proc];
    for(PS::S32 i=0; i<n_proc; i++) pos_domain[i] = dinfo.getPosDomain(i);

    Print("check 1", fout_debug);

    system_soft.exchangeParticle(dinfo); 

    Print("check 2", fout_debug);

    n_loc = system_soft.getNumberOfParticleLocal();

    for(PS::S32 i=0; i<n_loc; i++) system_soft[i].adr = i;

    Print("check 3", fout_debug);
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++) system_soft[i].rank_org = PS::Comm::getRank();

#ifdef USE_QUAD
    PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::QuadrupoleWithScatterSearch tree_soft;
#else
    PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::MonopoleWithScatterSearch tree_soft;
#endif

    Print("check 4", fout_debug);

    tree_soft.initialize(n_glb, theta, n_leaf_limit, n_group_limit);

    Print("check 5", fout_debug);

    const PS::F64vec root_cen(0.0);
    PS::F64 root_len = GetRootFullLenght(system_soft, root_cen);
    tree_soft.setParticleLocalTree(system_soft);
    tree_soft.setRootCell(root_len, root_cen);
    tree_soft.mortonSortLocalTreeOnly();
    tree_soft.linkCellLocalTreeOnly();
    tree_soft.calcMomentLocalTreeOnly();
    tree_soft.exchangeLocalEssentialTree(dinfo);
    tree_soft.setLocalEssentialTreeToGlobalTree();
    tree_soft.mortonSortGlobalTreeOnly();
    tree_soft.linkCellGlobalTreeOnly();
    tree_soft.calcMomentGlobalTreeOnly();
    tree_soft.makeIPGroup();
    Print("check 6", fout_debug);

#ifdef USE_QUAD
    tree_soft.calcForce(CalcForceEPEP(), CalcForceEPSPQuad(), true);
#else
    tree_soft.calcForce(CalcForceEPEP(), CalcForceEPSPMono(), true);
#endif
    n_loc = system_soft.getNumberOfParticleLocal();


    Print("check 7", fout_debug);
    PS::F64 rout_inv = 1.0 / EPISoft::r_out;
#pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
	system_soft[i].copyFromForce(tree_soft.getForce(i));
#if 1
	system_soft[i].pot_tot += system_soft[i].mass * rout_inv; // correction for self-gravity
#endif
    }

#ifdef AERO_DRAG
    gas_disk.calcAeroDrag(system_soft, time_sys);
    for(PS::S32 i=0; i<n_loc; i++){
	system_soft[i].acc_gd_prev = system_soft[i].acc_gd;
	system_soft[i].vel_prev = system_soft[i].vel;
    }
#endif

    Print("check 8", fout_debug);

    HardSystem system_hard;
    system_hard.setSun(mass_sun, pos_sun, vel_sun);

    Print("check 9", fout_debug);

    //const PS::F64 dt_limit_hard = dt_soft / dt_limit_hard_factor;
#ifdef DEBUG2
    PS::F64 dt_limit_hard = dt_soft/dt_limit_hard_factor;
#else
    PS::F64 dt_limit_hard = std::min( GetQuantizedValue(Tkep_min/1024.0), dt_soft/dt_limit_hard_factor);
#endif
    const PS::F64 eta_s = eta * 0.01;

    Print("check 10", fout_debug);

    ///////////////
    // hard part
    system_hard.setUp(tree_soft, system_soft, system_soft.getNumberOfParticleLocal(), r_out, r_in, pos_domain, time_sys);

    Print("check 11", fout_debug);

#ifdef FORCE_CHECK
    std::ofstream fout_force;
    char sout_force[1024];
    sprintf(sout_force, "%s/force.dat", dir_name);

    ForceSoft * force_soft_direct = new ForceSoft[n_loc];
    ForceSoft * force_planet_direct = new ForceSoft[n_loc];
    PS::F64vec * force_sun = new PS::F64vec[n_loc];
    CalcForceDirectFromAllPlanet(system_soft, dinfo, tree_soft, force_soft_direct, CalcForceEPEPCutoff64bit() );
    CalcForceDirectFromAllPlanet(system_soft, dinfo, tree_soft, force_planet_direct, CalcForceEPEPNoCutoff64bit() );
    for(PS::S32 i=0; i<n_loc; i++){
	PS::F64vec r_sun = system_soft[i].pos - pos_sun;
	PS::F64 r2 = r_sun * r_sun;
	PS::F64 r = sqrt(r2);
	PS::F64 r3 = r * r2;
	force_sun[i] = -PS::F64(mass_sun) / r3 * r_sun;
    }
    PS::F64 dacc_tot_loc = 0.0;
    PS::F64 dacc_planet_loc = 0.0;
    int moment_num = 1;
#ifdef USE_QUAD
    moment_num = 4;
#endif
    for(int ip=0; ip<PS::Comm::getNumberOfProc(); ip++){
	if(PS::Comm::getRank() == ip){
	    fout_force.open(sout_force, std::ios::app);
	    fout_force<<std::setprecision(15);
	    for(int i=0; i<n_loc; i++){
		PS::F64vec dforce = system_soft[i].acc - force_soft_direct[i].acc;
		PS::F64vec force_tot = force_planet_direct[i].acc + force_sun[i];
		fout_force<<"system_soft[i].acc= "<<system_soft[i].acc
			  <<" system_soft[i].id= "<<system_soft[i].id
			  <<" force_soft_direct[i]= "<<force_soft_direct[i].acc
			  <<" force_planet_direct[i]= "<<force_planet_direct[i].acc
			  <<" force_sun[i]= "<<force_sun[i]
			  <<" dforce= "<<dforce
			  <<" moment= "<<moment_num
			  <<std::endl;
		dacc_tot_loc += dforce*dforce / (force_tot*force_tot);
		dacc_planet_loc += dforce*dforce / (force_planet_direct[i].acc*force_planet_direct[i].acc);
	    }
	    fout_force.flush();
	    fout_force.close();
	}
	PS::Comm::barrier();
    }
    PS::F64 dacc_tot_glb    = sqrt( PS::Comm::getSum(dacc_tot_loc) / system_soft.getNumberOfParticleGlobal() );
    PS::F64 dacc_planet_glb = sqrt( PS::Comm::getSum(dacc_planet_loc) / system_soft.getNumberOfParticleGlobal() );
    if(PS::Comm::getRank() == 0){
	std::ofstream fout_frcchk;
	char sout_frcchk[1024];
	sprintf(sout_frcchk, "%s/frcchk.dat", dir_name);
	fout_frcchk.open(sout_frcchk, std::ios::app);
	fout_frcchk<<std::setprecision(15);
	fout_frcchk<<"theta= "<<theta
		   <<" r_cut_factor= "<<r_cut_factor
		   <<" n_group_limit= "<<n_group_limit
		   <<" dt_soft= "<<dt_soft
		   <<" search_factor= "<<search_factor
		   <<" dacc_tot_glb= "<<dacc_tot_glb
		   <<" dacc_planet_glb= "<<dacc_planet_glb
		   <<std::endl;
	fout_frcchk.close();
    }

    delete [] force_soft_direct;
    delete [] force_planet_direct;
    delete [] force_sun;
    PS::Finalize();
    return 0;
#endif

#ifdef AERO_DRAG
    PS::F64 eng_disp_aero_loc = 0.0;
#endif
    system_hard.clearEngDisp();
#ifdef READ_FILE
    if(!flag_snp_with_energy){
#ifdef ENERGY_DIRECT
	CalcEngKeplerDirect(system_soft, dinfo, tree_soft, system_hard.eng_disp_, eng_init, PS::F64(mass_sun), true);
#else //ENERGY_DIRECT
#ifdef AERO_DRAG
	CalcEngKepler(system_soft, system_hard.eng_disp_, eng_disp_aero_loc, eng_init, PS::F64(mass_sun), true);
#else //AERO_DRAG
	CalcEngKepler(system_soft, system_hard.eng_disp_, eng_init, PS::F64(mass_sun), true);
#endif //AERO_DRAG
#endif //ENERGY_DIRECT
	eng_now = eng_init;
    }
#else //READ_FILE
#ifdef ENERGY_DIRECT
    CalcEngKeplerDirect(system_soft, dinfo, tree_soft, system_hard.eng_disp_, eng_init, PS::F64(mass_sun), true);
#else //ENERGY_DIRECT
    CalcEngKepler(system_soft, system_hard.eng_disp_, eng_init, PS::F64(mass_sun), true);
#endif //ENERGY_DIRECT
    eng_now = eng_init;
#endif //READ_FILE
    eng_init.dump(std::cerr);
    Print("check 12", fout_debug);

    //PS::S64 snp_id = 0;
    PS::S64 n_loop = 0;
    PS::S64 n_loop_offset = 0;
    Wtime::cum_offset = PS::GetWtime();
    Wtime::interval_offset = PS::GetWtime();
    Wtime::clear();
    system_hard.clear_counter();
    PS::F64 dEerr_max = 0.0;

    Print("check 13", fout_debug);


    while(time_sys < time_end){
        Print("check a", fout_debug);
        Print(n_loop, fout_debug);

        if(1){
#ifdef ENERGY_DIRECT
            CalcEngKeplerDirect(system_soft, dinfo, tree_soft, system_hard.eng_disp_, eng_now, PS::F64(mass_sun), true);
#else //ENERGY_DIRECT
#ifdef AERO_DRAG
	    CalcEngKepler(system_soft, system_hard.eng_disp_, eng_disp_aero_loc, eng_now, PS::F64(mass_sun), true);
#else //AERO_DRAG
	    CalcEngKepler(system_soft, system_hard.eng_disp_, eng_now, PS::F64(mass_sun), true);
#endif//AERO_DRAG
#endif//ENERGY_DIRECT

#ifdef READ_FILE
            if(fmod(time_sys, dt_snp) == 0.0 && n_loop != 0){
#else
            if(fmod(time_sys, dt_snp) == 0.0){
#endif
		PS::Comm::barrier();
		Wtime::take_snp_offset = PS::GetWtime();
		char file_snp[1024];
		sprintf(file_snp, "%s/snap%05d.dat", dir_name, (int)snp_id++);
		PS::Comm::broadcast(&eng_init, 1, 0);
		PS::Comm::broadcast(&eng_now, 1, 0);
		FileHeaderWithEnergy header(system_soft.getNumberOfParticleGlobal(), time_sys, eng_init, eng_now);
		system_soft.writeParticleAscii(file_snp, header);

		// new position
		if(system_hard.merge_history_.size() > 0){
		    fout_merge<<"system_hard.merge_history_.size()="<<system_hard.merge_history_.size()<<std::endl;
		    //fout_merge<<"n_loc_new="<<n_loc_new<<" n_loc_old="<<n_loc_old<<std::endl;
		    for(size_t ip=0; ip<system_hard.merge_history_.size(); ip++){
			system_hard.merge_history_[ip].dump(fout_merge);
		    }
		    fout_merge<<std::endl;
		    system_hard.merge_history_.clear(); // clear history
		}
		// new position

		PS::Comm::barrier();
		Wtime::take_snp += PS::GetWtime() - Wtime::take_snp_offset;
	    }
	    system_hard.clearEngDisp();
	    PS::F64 dEerr = (eng_now.tot-eng_init.tot)/eng_init.tot;
	    if( fabs(dEerr_max) < fabs(dEerr) ) dEerr_max = dEerr;
#ifdef ENERGY_CHECK
	    if(n_loop != 0 && fmod(time_sys, 1.0) == 0.0){
#else // no ENERGY_CHECK
#ifdef PROFILE
	    if(n_loop != 0){
#else //NOPROFILE
            if(n_loop != 0 && fmod(time_sys, 1.0) == 0.0){
#endif //PROFILE
#endif //ENERGY_CHECK
                        PS::Comm::barrier();
                        Wtime::interval = PS::GetWtime() - Wtime::interval_offset;
                        Wtime::cum = PS::GetWtime() - Wtime::cum_offset;
                        Wtime::dump(fout_tcal, time_sys, n_loop-n_loop_offset);

                        Profile<PS::DomainInfo, PS::ParticleSystem<FPSoft>, PS::TreeForForceLong<ForceSoft, EPISoft, EPJSoft>::QuadrupoleWithScatterSearch> 
                            time_profile(&dinfo, &system_soft, &tree_soft, 34, 64, 1.6e9);
                        time_profile.dump(fout_tcal, time_sys, n_loop-n_loop_offset, Wtime::interval);

                        system_hard.dump_counter(fout_tcal, n_loop-n_loop_offset);
                        fout_tcal<<"tree_soft.getNumberOfLETEPSend1stLocal()= "<<tree_soft.getNumberOfLETEPSend1stLocal()
				 <<"   "<<tree_soft.getNumberOfLETEPSend1stLocal()/(n_loop-n_loop_offset)
                                 <<" tree_soft.getNumberOfLETSPSend1stLocal()= "<<tree_soft.getNumberOfLETSPSend1stLocal()
				 <<"   "<<tree_soft.getNumberOfLETSPSend1stLocal()/(n_loop-n_loop_offset)
                                 <<" tree_soft.getNumberOfLETEPSend1stGlobal()= "<<tree_soft.getNumberOfLETEPSend1stGlobal()
				 <<"   "<<tree_soft.getNumberOfLETEPSend1stGlobal()/(n_loop-n_loop_offset)
                                 <<" tree_soft.getNumberOfLETSPSend1stGlobal()= "<<tree_soft.getNumberOfLETSPSend1stGlobal()
				 <<"   "<<tree_soft.getNumberOfLETSPSend1stGlobal()/(n_loop-n_loop_offset)
				 <<std::endl;
                        fout_tcal<<"tree_soft.getNumberOfLETEPRecv1stLocal()= "<<tree_soft.getNumberOfLETEPRecv1stLocal()
				 <<"   "<<tree_soft.getNumberOfLETEPRecv1stLocal()/(n_loop-n_loop_offset)
                                 <<" tree_soft.getNumberOfLETSPRecv1stLocal()= "<<tree_soft.getNumberOfLETSPRecv1stLocal()
				 <<"   "<<tree_soft.getNumberOfLETSPRecv1stLocal()/(n_loop-n_loop_offset)
                                 <<" tree_soft.getNumberOfLETEPRecv1stGlobal()= "<<tree_soft.getNumberOfLETEPRecv1stGlobal()
				 <<"   "<<tree_soft.getNumberOfLETEPRecv1stGlobal()/(n_loop-n_loop_offset)
                                 <<" tree_soft.getNumberOfLETSPRecv1stGlobal()= "<<tree_soft.getNumberOfLETSPRecv1stGlobal()
				 <<"   "<<tree_soft.getNumberOfLETSPRecv1stGlobal()/(n_loop-n_loop_offset)
				 <<std::endl;
                        fout_tcal<<" system_soft.getNumberOfParticleLocal()= "<<system_soft.getNumberOfParticleLocal()<<std::endl;
                        fout_tcal<<"time_sys= "<<time_sys<<" dinfo.getPosDomain(PS::Comm::getRank())="<<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
                        fout_tcal<<"-----------------------"<<std::endl;
                        fout_tcal<<std::endl;
                        fout_debug<<"time_sys= "<<time_sys<<" dinfo.getPosDomain(PS::Comm::getRank())="<<dinfo.getPosDomain(PS::Comm::getRank())<<std::endl;
                        system_hard.clear_counter();
                        
                        //CalcEngKepler(system_soft, system_hard.eng_disp_, eng, PS::F64(mass_sun), true);
                        //system_hard.clearEngDisp();
                        PS::S64 n_body_tmp = system_soft.getNumberOfParticleGlobal();
#ifdef CALC_HARD_ENERGY
                        dEerr_1body_glb = PS::Comm::getSum(dEerr_1body_loc);
                        dEerr_2body_glb = PS::Comm::getSum(dEerr_2body_loc);
                        dEerr_mbody_glb = PS::Comm::getSum(dEerr_mbody_loc);
#endif
                        if(PS::Comm::getRank() == 0){
                            fout_diag<<"time_sys= "<<time_sys
                                     <<" n_body= "<<n_body_tmp
                                     <<" mass_planet_max= "<<mass_planet_max_glb
                                     <<" r_hill_max= "<<r_hill_max_glb
                                     <<" ecc_rms= "<<ecc_rms
                                     <<" inc_rms= "<<inc_rms
                                     <<" (eng_now.tot-eng_init.tot)/eng_init.tot= "<<(eng_now.tot-eng_init.tot)/eng_init.tot
                                     <<" (eng_now.tot-eng_init.tot)/eng_init.pot_planet= "<<(eng_now.tot-eng_init.tot)/eng_init.pot_planet
                                     <<" eng_now.tot= "<<eng_now.tot
                                     <<" eng_now.pot= "<<eng_now.pot
                                     <<" eng_now.kin= "<<eng_now.kin
                                     <<" eng_now.pot_planet= "<<eng_now.pot_planet
                                     <<" eng_now.disp_merge= "<<eng_now.disp_merge
                                     <<" eng_now.disp_aero= "<<eng_now.disp_aero
                                     <<" eng_init.tot= "<<eng_init.tot
                                     <<" eng_init.pot= "<<eng_init.pot
                                     <<" eng_init.kin= "<<eng_init.kin
                                     <<" eng_init.pot_planet= "<<eng_init.pot_planet
                                     <<" eng_init.disp_merge= "<<eng_init.disp_merge
                                     <<" eng_init.disp_aero= "<<eng_init.disp_aero
                                     <<" dEerr_max= "<<dEerr_max
				     <<" hight_rms= "<<hight_rms
				     <<" r_rms= "<<r_rms
#ifdef CALC_HARD_ENERGY
                                     <<" dEerr_1body_loc= "<<dEerr_1body_loc / eng_init.tot // 39,40
                                     <<" dEerr_1body_glb= "<<dEerr_1body_glb / eng_init.tot // 41,42
                                     <<" dEerr_2body_loc= "<<dEerr_2body_loc / eng_init.tot // 43,44
                                     <<" dEerr_2body_glb= "<<dEerr_2body_glb / eng_init.tot // 45,46
                                     <<" dEerr_mbody_loc= "<<dEerr_mbody_loc / eng_init.tot // 47,48
                                     <<" dEerr_mbody_glb= "<<dEerr_mbody_glb / eng_init.tot // 49,50
#endif
                                     <<std::endl;
                        }
                        n_loop_offset = n_loop;
                        Wtime::clear();
                        Wtime::interval_offset = PS::GetWtime();
                        time_profile.clear();
                    }
                }
                Print("check b", fout_debug);

                Wtime::soft_offset = PS::GetWtime();
                /////////////
                // 1st KICK
                Kick(system_soft, tree_soft, dt_soft*0.5);
                Print("check c", fout_debug);
                system_hard.copyVelSoftToHard(system_soft);
                Print("check d", fout_debug);
                // 1st KICK
                /////////////
                Wtime::soft += PS::GetWtime() - Wtime::soft_offset;

                //////////////
                // HARD PART
                Wtime::hard_offset = PS::GetWtime();
                //////////////
                // DRIFT 1body
#ifdef CALC_HARD_ENERGY
                Energy eng_hard_1body_prev_loc;
                CalcOneBodyEnergy(system_soft, mass_sun, pos_sun, vel_sun, eng_hard_1body_prev_loc);
#endif
                Wtime::hard_1body_offset = PS::GetWtime();
		//Print(system_soft[0].pos, fout_debug);
		//Print(system_soft[0].vel, fout_debug);
                DriftKepler(system_soft, tree_soft, dt_soft);
                Wtime::hard_1body += PS::GetWtime() - Wtime::hard_1body_offset;
#ifdef CALC_HARD_ENERGY
                Energy eng_hard_1body_now_loc;
                CalcOneBodyEnergy(system_soft, mass_sun, pos_sun, vel_sun, eng_hard_1body_now_loc);
                dEerr_1body_loc += eng_hard_1body_now_loc.tot - eng_hard_1body_prev_loc.tot;
                //dEerr_1body_glb += PS::Comm::getSum(dEerr_1body_loc);
#endif
                Print("check e", fout_debug);
                // DRIFT 1body
                //////////////
                
#ifndef ONEBODY
                ////////
                // DRIFT 2body
                Wtime::hard_2body_select_system_offset = PS::GetWtime();
                system_hard.selectIsolatedSystem();
                Wtime::hard_2body_select_system += PS::GetWtime() - Wtime::hard_2body_select_system_offset;
                Print("check f", fout_debug);

                Wtime::hard_2body_offset = PS::GetWtime();
                system_hard.evolveIsolatedSystem(system_soft, r_out, r_in, mass_sun, pos_sun, vel_sun,
                                                 eta_s, eta,  dt_soft,  dt_limit_hard);
                Wtime::hard_2body += PS::GetWtime() - Wtime::hard_2body_offset;
#ifdef CALC_HARD_ENERGY
                dEerr_2body_loc += (system_hard.eng_2body_now_loc_.tot + system_hard.eng_2body_now_loc_.disp_merge) - system_hard.eng_2body_prev_loc_.tot;
                //dEerr_2body_loc += system_hard.eng_2body_now_loc_.tot - system_hard.eng_2body_prev_loc_.tot;
                //dEerr_2body_glb += PS::Comm::getSum(dEerr_2body_loc);
#endif
                Print("check g", fout_debug);
                // DRIFT 2body
                ////////

                /////////
                // GATHER HARD PTCL
                Wtime::hard_gather_data_offset = PS::GetWtime(); 
                system_hard.gatherData();
                Wtime::hard_gather_data += PS::GetWtime() - Wtime::hard_gather_data_offset;
                // GATHER HARD PTCL
                /////////

                Print("check h", fout_debug);
                ////////
                // DRIFT multi body
                Wtime::hard_multi_offset = PS::GetWtime();
                if(PS::Comm::getRank() == 0){
                    if( system_hard.ptcl_multi_glb_.size() > 0){
                        system_hard.evolveMultiSirial(dt_soft, dt_limit_hard, r_out, r_in, eta_s, eta);
                    }
                }
                Wtime::hard_multi += PS::GetWtime() - Wtime::hard_multi_offset;
#ifdef CALC_HARD_ENERGY
                dEerr_mbody_loc += (system_hard.eng_mbody_now_loc_.tot + system_hard.eng_mbody_now_loc_.disp_merge) - system_hard.eng_mbody_prev_loc_.tot;
                //dEerr_mbody_glb += PS::Comm::getSum(dEerr_mbody_loc);
#endif
                // DRIFT multi body
                ////////

                Print("check i", fout_debug);
                /////////
                // SCATTER HARD PTCL
                Wtime::hard_scatter_data_offset = PS::GetWtime(); 
                system_hard.scatterData();
                Wtime::hard_scatter_data += PS::GetWtime() - Wtime::hard_scatter_data_offset;
                // SCATTER HARD PTCL
                /////////
                Print("check j", fout_debug);

                /////////
                // COPY HARD PTCL
                Wtime::hard_copy_h2s_offset = PS::GetWtime();
                system_hard.copyPtclHardToSoft(system_soft);
                Wtime::hard_copy_h2s += PS::GetWtime() - Wtime::hard_copy_h2s_offset;
                // COPY HARD PTCL
                /////////

                Print("check k", fout_debug);
                
                /////////
                // MERGE HARD PTCL
                Wtime::hard_merge_offset = PS::GetWtime();
#ifdef MERGE
                const PS::S32 n_loc_old = system_soft.getNumberOfParticleLocal();
                PS::S32 n_loc_new = 0;
                Print("check l", fout_debug);
                for(PS::S32 ip=0; ip<n_loc_old; ip++){
                    if(system_soft[ip].mass > 0.0){ // not merged
                        system_soft[n_loc_new] = system_soft[ip];
                        n_loc_new++;
                    }
                }
                Print("check m", fout_debug);
                /*
// original
        if(system_hard.merge_history_.size() > 0){
            fout_merge<<"system_hard.merge_history_.size()="<<system_hard.merge_history_.size()<<std::endl;
            fout_merge<<"n_loc_new="<<n_loc_new<<" n_loc_old="<<n_loc_old<<std::endl;
            for(size_t ip=0; ip<system_hard.merge_history_.size(); ip++){
                system_hard.merge_history_[ip].dump(fout_merge);
            }
            fout_merge<<std::endl;
            system_hard.merge_history_.clear(); // clear history
        }
	*/
                Print("check n", fout_debug);
                PS::S32 n_glb_old = PS::Comm::getSum(n_loc_old);
                PS::S32 n_glb_new = PS::Comm::getSum(n_loc_new);
                if(n_glb_new != n_glb_old){
                    mass_planet_max_glb = GetMassMax(system_soft, n_loc_new);
                    if(PS::Comm::getRank() == 0){
                        std::cout<<time_sys<<"   "<<n_glb_new<<std::endl;
                        //fout_diag<<time_sys<<"   "<<n_glb_new<<"   "<<mass_planet_max_glb<<std::endl;
                    }
                }
                Print("check o", fout_debug);
                system_soft.setNumberOfParticleLocal(n_loc_new);
#endif // MERGE
                Wtime::hard_merge += PS::GetWtime() - Wtime::hard_merge_offset;
                system_hard.accumulate_counter();
                Wtime::hard += PS::GetWtime() - Wtime::hard_offset;
                // HARD PART
                //////////////

#endif //ONEBODY

        //////////////
        // SOFT PART
	Wtime::soft_offset = PS::GetWtime();
#ifndef DIV_FIX
	if( n_loop < 12 || n_loop % 16 == 0){
	    if(PS::Comm::getRank() == 0){
		fout_domain<<time_sys<<"   "<<tree_soft.pos_root_cell_;
		for(PS::S32 i=0; i<PS::Comm::getNumberOfProc(); i++){
		    fout_domain<<"   "<<dinfo.getPosDomain(i);
		}
		fout_domain<<std::endl;
	    }
	    dinfo.collectSampleParticle(system_soft);
	    dinfo.decomposeDomain();
	}
#endif
        //////////////////////
        // EVALUATE SOFT FORCE
	Print("check s", fout_debug);
        system_soft.exchangeParticle(dinfo);
	Print("check t", fout_debug);

	n_loc = system_soft.getNumberOfParticleLocal();
	for(PS::S32 i=0; i<n_loc; i++) system_soft[i].adr = i;

#pragma omp parallel for
	for(PS::S32 i=0; i<n_loc; i++) system_soft[i].rank_org = PS::Comm::getRank();


	if(fmod(time_sys, 100.0) == 0.0){
	    //#ifdef PROFILE
	    hight_rms_loc = 0.0;
	    r_rms_loc = 0.0;
	    for(PS::S32 i=0; i<n_loc; i++){
		hight_rms_loc += system_soft[i].pos.z * system_soft[i].pos.z;
		r_rms_loc += system_soft[i].pos * system_soft[i].pos;
	    }
	    hight_rms = sqrt(PS::Comm::getSum(hight_rms_loc) / n_glb);
	    r_rms = sqrt(PS::Comm::getSum(r_rms_loc) / n_glb);
	    //#endif //PROFILE
            SetRoutRinRhillVdisp(system_soft,           pos_sun,
                                 vel_sun,               mass_sun,
                                 ratio_r_cut,           r_cut_factor,
                                 r_out,                 r_in, 
                                 r_hill_max_glb,        vel_disp,
				 ecc_rms,               inc_rms,
                                 mass_planet_max_glb,   mass_planet_tot_glb,
				 Tkep_min);
#ifdef DEBUG2
	    dt_limit_hard = dt_soft/dt_limit_hard_factor;
#else
	    dt_limit_hard = std::min( GetQuantizedValue(Tkep_min)/1024.0, dt_soft/dt_limit_hard_factor);
#endif
            EPJSoft::r_search = r_out + search_factor*vel_disp*dt_soft;
	    EPISoft::r_out = r_out;
	    EPISoft::r_in = r_in;
            if(r_out <= EPISoft::eps){
                EPJSoft::r_search = 0.0;
            }
        }
	Print("check u", fout_debug);
	Wtime::soft_force_offset = PS::GetWtime();
#ifdef FREE_MALLOC_TREE
	tree_soft.freeMem();
	tree_soft.reallocMem();
#endif //FREE_MALLOC_TREE

	root_len = GetRootFullLenght(system_soft, root_cen);
	tree_soft.setParticleLocalTree(system_soft);
	tree_soft.setRootCell(root_len, root_cen);
	tree_soft.mortonSortLocalTreeOnly();
	tree_soft.linkCellLocalTreeOnly();
	tree_soft.calcMomentLocalTreeOnly();
	tree_soft.exchangeLocalEssentialTree(dinfo);
	tree_soft.setLocalEssentialTreeToGlobalTree();
	tree_soft.mortonSortGlobalTreeOnly();
	tree_soft.linkCellGlobalTreeOnly();
	tree_soft.calcMomentGlobalTreeOnly();
	tree_soft.makeIPGroup();
#ifdef USE_QUAD
	tree_soft.calcForce(CalcForceEPEP(), CalcForceEPSPQuad(), true); 
#else
	tree_soft.calcForce(CalcForceEPEP(), CalcForceEPSPMono(), true);
#endif
	n_loc = system_soft.getNumberOfParticleLocal();

	rout_inv = 1.0 / EPISoft::r_out;
#pragma omp parallel for
	for(PS::S32 i=0; i<n_loc; i++){
	    system_soft[i].copyFromForce(tree_soft.getForce(i));
#if 1
	    system_soft[i].pot_tot += system_soft[i].mass * rout_inv;
#endif
	}

#ifdef AERO_DRAG
	gas_disk.calcAeroDrag(system_soft, time_sys);
#endif //AERO_DRAG


	Wtime::soft_force += PS::GetWtime() - Wtime::soft_force_offset;
	Print("check v", fout_debug);
	//for(int i=0; i<system_soft.getNumberOfParticleLocal(); i++) system_soft[i].acc_pla = system_soft[i].acc;
        system_hard.setUp(tree_soft, system_soft, system_soft.getNumberOfParticleLocal(), r_out, r_in, pos_domain, time_sys);

	Print("check w", fout_debug);
        // EVALUATE SOFT FORCE
        //////////////////////

        //////////////
        // second KICK
        Kick(system_soft, tree_soft, dt_soft*0.5);
        time_sys += dt_soft;
        // second KICK
        //////////////

        //////////////
	// energy correction for aero drag
#ifdef AERO_DRAG
	eng_disp_aero_loc = 0.0;
	for(PS::S32 i=0; i<n_loc; i++){
	    eng_disp_aero_loc -= dt_soft * 0.5 * system_soft[i].mass * 
		(system_soft[i].acc_gd_prev*system_soft[i].vel_prev
		 + system_soft[i].acc_gd*system_soft[i].vel);
	    system_soft[i].acc_gd_prev = system_soft[i].acc_gd;
	    system_soft[i].vel_prev = system_soft[i].vel;
	}
#endif
	// energy correction for aero drag
        //////////////

#ifdef ONEBODY
	if(n_loop % 1000000000){
	    fout_domain<<time_sys<<"  "<<system_soft[0].pos<<"  "<<system_soft[0].vel<<std::endl;
	}
	if(system_soft[0].pos*system_soft[0].pos < 1.0) break;
#endif //ONEBODY

	Wtime::soft += PS::GetWtime() - Wtime::soft_offset;
        // SOFT PART
        //////////////
        n_loop++;
	Print("check x", fout_debug);
    }

    std::cerr<<"system_hard.merge_history_.size()="<<system_hard.merge_history_.size()<<std::endl;
    std::cerr<<"system_soft.getNumberOfParticleGlobal()="<<system_soft.getNumberOfParticleGlobal()<<std::endl;
    /*
    for(size_t ip=0; ip<system_hard.merge_history_.size(); ip++){
        system_hard.merge_history_[ip].dump();
    }
    */

#ifdef ENERGY_CHECK
    std::ofstream fout_engchk;
    char sout_engchk[1024];
    if(PS::Comm::getRank() == 0){
	sprintf(sout_engchk, "%s/engchk.dat", dir_name);
	fout_engchk.open(sout_engchk, std::ios::app);
	fout_engchk<<std::setprecision(15);
	fout_engchk<<"dt_soft= "<<dt_soft
		   <<" r_cut_factor= "<<r_cut_factor
		   <<" theta= "<<theta
		   <<" eta= "<<eta
		   <<" search_factor= "<<search_factor
		   <<" dt_limit_hard_factor= "<<dt_limit_hard_factor
		   <<" dEerr_max= "<<dEerr_max
		   <<" dt_limit_hard= "<<dt_limit_hard
		   <<" Tkep_min= "<<Tkep_min<<std::endl;

	fout_engchk.close();
    }
#endif

    PS::Finalize();
    return 0;
}

