#ifdef INTRINSIC_K
#include"phantomquad_for_p3t_k.hpp"
#endif
#ifdef INTRINSIC_X86
#include"phantomquad_for_p3t_x86.hpp"
#endif

const PS::F64 SAFTY_FACTOR_FOR_SEARCH = 1.05;
const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = SAFTY_FACTOR_FOR_SEARCH * SAFTY_FACTOR_FOR_SEARCH;

inline PS::F64 CalcK(const PS::F64 rij,
                     const PS::F64 rout,
                     const PS::F64 rin){
    PS::F64 inv_dr = 1.0 / (rout-rin);
    PS::F64 x = (rij - rin)*inv_dr;
    x = (x < 1.0) ? x : 1.0;
    x = (x > 0.0) ? x : 0.0;
    PS::F64 x2 = x*x;
    PS::F64 x4 = x2*x2;
    PS::F64 k = (((-20.0*x+70.0)*x-84.0)*x+35.0)*x4;
    std::max( std::min(k, 1.0), 0.0);
    //k = (k < 1.0) ? k : 1.0;
    //k = (k > 0.0) ? k : 0.0;
    return k;
}

class ForceSoft{
public:
    PS::F64vec acc; // soft
    PS::F64 pot; // soft
    PS::S32 n_ngb;
    void clear(){
        acc = 0.0;
        pot = 0.0;
        n_ngb = 0;
    }
};

class FPSoft{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc; // soft
#ifdef AERO_DRAG
    PS::F64vec acc_gd; // aero drag
    PS::F64vec vel_prev; // for energy check
    PS::F64vec acc_gd_prev; // aero drag
#endif
    PS::F64 pot_tot; // soft + hard
    PS::S32 rank_org;
    PS::S32 n_ngb;
    PS::S32 adr;

    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & p) { pos = p; }
    void copyFromForce(const ForceSoft & force){
        acc = force.acc;
        pot_tot = force.pot;
	n_ngb = force.n_ngb;
    }

    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %d \n", 
                this->id, this->mass, 
		this->pos.x, this->pos.y, this->pos.z,  // 3-5
		this->vel.x, this->vel.y, this->vel.z,  // 6-8
		this->acc.x, this->acc.y, this->acc.z,  // 9-11
		this->pot_tot, this->n_ngb);
    }

    void readAscii(FILE* fp) {
        fscanf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d \n", 
                &this->id, &this->mass, 
		&this->pos.x, &this->pos.y, &this->pos.z,  // 3-5
		&this->vel.x, &this->vel.y, &this->vel.z,  // 6-8
		&this->acc.x, &this->acc.y, &this->acc.z,  // 9-11
		&this->pot_tot, &this->n_ngb);
    }

    /*
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %d \n", 
                this->id, this->mass, 
		this->pos.x, this->pos.y, this->pos.z,  // 3-5
		this->vel.x, this->vel.y, this->vel.z,  // 6-8
		this->acc.x, this->acc.y, this->acc.z,  // 9-11
		this->acc_pla.x, this->acc_pla.y, this->acc_pla.z, // 12-14
		this->pot_tot, this->n_ngb);
    }
    void readAscii(FILE* fp) {
        fscanf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d \n", 
                &this->id, &this->mass, 
		&this->pos.x, &this->pos.y, &this->pos.z,  // 3-5
		&this->vel.x, &this->vel.y, &this->vel.z,  // 6-8
		&this->acc.x, &this->acc.y, &this->acc.z,  // 9-11
		&this->acc_pla.x, &this->acc_pla.y, &this->acc_pla.z, // 12-14
		&this->pot_tot, &this->n_ngb);
    }
    */

    /*
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 
                this->id, this->mass, this->pos.x, this->pos.y, this->pos.z, this->vel.x, this->vel.y, this->vel.z);
    }
    */

    void dump(std::ofstream & fout){
	fout<<"id= "<<id<<std::endl;
	fout<<"adr= "<<adr<<std::endl;
	fout<<"mass= "<<mass<<std::endl;
	fout<<"pos= "<<pos<<std::endl;
	fout<<"vel= "<<vel<<std::endl;
	fout<<"acc= "<<acc<<std::endl;
	fout<<"pot_tot= "<<pot_tot<<std::endl;
    }
};

class EPISoft{
public:
    PS::S64 id;
    PS::F64vec pos;
    static PS::F64 eps;
    static PS::F64 r_out;
    static PS::F64 r_in;
    PS::S32 rank_org;
    PS::F64vec getPos() const { return pos;}
    void copyFromFP(const FPSoft & fp){ 
        pos = fp.pos;
        id = fp.id;
	rank_org = fp.rank_org;
    }
    void dump(std::ostream & fout=std::cout) const {
	fout<<"id="<<id<<std::endl;
	fout<<"rank_org="<<rank_org<<std::endl;
	fout<<"pos="<<pos<<std::endl;
	fout<<"eps="<<eps<<std::endl;
    }
};

PS::F64 EPISoft::eps = 1.0/1024.0;
PS::F64 EPISoft::r_out = 0.0;
PS::F64 EPISoft::r_in = 0.0;



class EPJSoft{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::S32 rank_org;
    PS::S32 adr_org;
    static PS::F64 r_search;
    void copyFromFP(const FPSoft & fp){
        mass = fp.mass;
        pos = fp.pos;
        id = fp.id;
        vel = fp.vel;
	rank_org = fp.rank_org;
	adr_org = fp.adr;
    }
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec & pos_new){ pos = pos_new;}
    PS::F64 getCharge() const { return mass; }
    //PS::F64 getRSearch() const { return r_search; }
    PS::F64 getRSearch() const { return r_search * SAFTY_FACTOR_FOR_SEARCH; }
    // FORDEBUG
    void dump(std::ostream & fout=std::cout) const {
	fout<<"id="<<id<<std::endl;
	fout<<"rank_org="<<rank_org<<std::endl;
	fout<<"mass="<<mass<<std::endl;
	fout<<"pos="<<pos<<std::endl;
	fout<<"vel="<<vel<<std::endl;
    }
    void clear(){
	mass = 0.0;
	pos = vel = 0.0;
	id = rank_org = adr_org = -1;
    }
};

PS::F64 EPJSoft::r_search;

#ifdef INTRINSIC

struct CalcForceEPEP{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
#ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
#else
        static __thread PhantomGrapeQuad pg;
#endif
        pg.set_eps2(eps2);
        pg.set_r_crit2( (r_crit2+eps2)*SAFTY_FACTOR_FOR_SEARCH_SQ );
        pg.set_cutoff(EPISoft::r_out, EPISoft::r_in);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
#ifdef INTRINSIC_K
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#else
#if defined(CALC_EP_64bit) || defined(CALC_EP_MIX)
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z);
#else
	    pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#endif // CALC_EP_64bit
#endif
        }
	PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
	for(PS::S32 loop=0; loop<loop_max; loop++){
	    const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
	    const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;
	    const PS::S32 it =ih + n_jp_tmp;
	    PS::S32 i_tmp = 0;
	    for(PS::S32 i=ih; i<it; i++, i_tmp++){
		const PS::F64 m_j = ep_j[i].getCharge();
		const PS::F64vec pos_j = ep_j[i].getPos();
#ifdef INTRINSIC_K
		pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#else
#if defined(CALC_EP_64bit) || defined(CALC_EP_MIX)
		pg.set_epj_one_d(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#else
		pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#endif
#endif
	    }
#ifdef CALC_EP_64bit
	    //pg.run_epj_for_p3t_d(n_ip, n_jp_tmp);
	    pg.run_epj_for_p3t_with_linear_cutoff_d(n_ip, n_jp_tmp);
#elif CALC_EP_MIX
	    pg.run_epj_for_p3t_mix(n_ip, n_jp_tmp);
#else
            //pg.run_epj_for_p3t(n_ip, n_jp_tmp);
	    pg.run_epj_for_p3t_with_linear_cutoff(n_ip, n_jp_tmp);
#endif
	    for(PS::S32 i=0; i<n_ip; i++){
		PS::F64 * p = &(force[i].pot);
		PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
		PS::F64 n_ngb = 0; 
#ifdef CALC_EP_64bit
		pg.accum_accp_one_d(i, a[0], a[1], a[2], *p, n_ngb);
#else
		pg.accum_accp_one(i, a[0], a[1], a[2], *p, n_ngb);
#endif
		force[i].n_ngb += (PS::S32)(n_ngb*1.00001);
	    }
	}
        for(PS::S32 i=0; i<n_ip; i++){
	    force[i].n_ngb--;
	}
    }
};


struct CalcForceEPEPNoCutoff{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
#ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
#else
        static __thread PhantomGrapeQuad pg;
#endif
	//PhantomGrapeQuad pg;
        pg.set_eps2(eps2);
	pg.set_r_crit2( (r_crit2+eps2)*SAFTY_FACTOR_FOR_SEARCH_SQ );
	pg.set_cutoff(EPISoft::r_out, EPISoft::r_in);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
#ifdef INTRINSIC_K
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#else
#ifdef CALC_EP_64bit
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z);
#else
	    pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#endif // CALC_EP_64bit
#endif
        }
	PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
	for(PS::S32 loop=0; loop<loop_max; loop++){
	    const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
	    const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;
	    const PS::S32 it =ih + n_jp_tmp;
	    PS::S32 i_tmp = 0;

	    for(PS::S32 i=ih; i<it; i++, i_tmp++){
		const PS::F64 m_j = ep_j[i].getCharge();
		const PS::F64vec pos_j = ep_j[i].getPos();
#ifdef INTRINSIC_K
		pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#else
#ifdef CALC_EP_64bit
		pg.set_epj_one_d(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#else
		pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#endif
#endif
	    }
#ifdef CALC_EP_64bit
	    pg.run_epj_d(n_ip, n_jp_tmp);
#else
            pg.run_epj(n_ip, n_jp_tmp);
#endif
	    for(PS::S32 i=0; i<n_ip; i++){
		PS::F64 * p = &(force[i].pot);
		PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
		PS::F64 n_ngb = 0; 
#ifdef CALC_EP_64bit
		pg.accum_accp_one_d(i, a[0], a[1], a[2], *p, n_ngb);
#else
		pg.accum_accp_one(i, a[0], a[1], a[2], *p, n_ngb);
#endif
		force[i].n_ngb += (PS::S32)(n_ngb*1.00001);
	    }
	}
        for(PS::S32 i=0; i<n_ip; i++){
	    force[i].n_ngb--;
	}
    }
};

#ifdef INTRINSIC_X86

struct CalcForceEPEPNoCutoff64bit{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
#ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
#else
        static __thread PhantomGrapeQuad pg;
#endif
        //PhantomGrapeQuad pg;
        pg.set_eps2(eps2);
        pg.set_r_crit2( (r_crit2+eps2)*SAFTY_FACTOR_FOR_SEARCH_SQ );
        pg.set_cutoff(EPISoft::r_out, EPISoft::r_in);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z);
        }
        PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
        for(PS::S32 loop=0; loop<loop_max; loop++){
            const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
            const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;
	    const PS::S32 it =ih + n_jp_tmp;
	    PS::S32 i_tmp = 0;

	    for(PS::S32 i=ih; i<it; i++, i_tmp++){
		const PS::F64 m_j = ep_j[i].getCharge();
		const PS::F64vec pos_j = ep_j[i].getPos();
		pg.set_epj_one_d(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
	    }
	    pg.run_epj_d(n_ip, n_jp_tmp);
	    for(PS::S32 i=0; i<n_ip; i++){
		PS::F64 * p = &(force[i].pot);
		PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
		PS::F64 n_ngb = 0; 
		pg.accum_accp_one_d(i, a[0], a[1], a[2], *p, n_ngb);
		force[i].n_ngb += (PS::S32)(n_ngb*1.00001);
	    }
	}
        for(PS::S32 i=0; i<n_ip; i++){
	    force[i].n_ngb--;
	}
    }
};

#endif

//struct CalcForceEPSP{
struct CalcForceEPSPMono{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      //const PS::SPJMonopoleScatter * sp_j,
		      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
#ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
#else
        static __thread PhantomGrapeQuad pg;
#endif
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
#ifdef CALC_SP_64bit
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z);
#else
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#endif
        }
	PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
	for(PS::S32 loop=0; loop<loop_max; loop++){
	    const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
	    const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;  
	    const PS::S32 it = ih + n_jp_tmp;
	    PS::S32 i_tmp = 0;
	    for(PS::S32 i=ih; i<it; i++, i_tmp++){
		const PS::F64 m_j = sp_j[i].getCharge();
		const PS::F64vec pos_j = sp_j[i].getPos();
#ifdef CALC_SP_64bit
		pg.set_epj_one_d(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#else
		pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
#endif
		//pg.set_epj_one(i_tmp, pos_j.x, pos_j.y, pos_j.z, m_j);
	    }
#ifdef CALC_SP_64bit
	    pg.run_epj_d(n_ip, n_jp_tmp);
#else
            pg.run_epj(n_ip, n_jp_tmp);
#endif
            //pg.run_epj(n_ip, n_jp_tmp);
	    for(PS::S32 i=0; i<n_ip; i++){
		PS::F64 * p = &(force[i].pot);
		PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
#ifdef CALC_SP_64bit
		pg.accum_accp_one_d(i, a[0], a[1], a[2], *p);
#else
		pg.accum_accp_one(i, a[0], a[1], a[2], *p);
#endif
		//pg.accum_accp_one(i, a[0], a[1], a[2], *p);
	    }
	}
    }
};


// SIMD
struct CalcForceEPSPQuad{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      //const PS::SPJMonopoleScatter * sp_j,
		      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
#ifdef __HPC_ACE__
        PhantomGrapeQuad pg;
#else
        static __thread PhantomGrapeQuad pg;
#endif
        pg.set_eps2(eps2);
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec pos_i = ep_i[i].getPos();
#ifdef CALC_SP_64bit
            pg.set_xi_one_d(i, pos_i.x, pos_i.y, pos_i.z);
#else
            pg.set_xi_one(i, pos_i.x, pos_i.y, pos_i.z);
#endif
        }
	PS::S32 loop_max = (n_jp-1) / PhantomGrapeQuad::NJMAX + 1;
	for(PS::S32 loop=0; loop<loop_max; loop++){
	    const PS::S32 ih = PhantomGrapeQuad::NJMAX*loop;
	    const PS::S32 n_jp_tmp = ( (n_jp - ih) < PhantomGrapeQuad::NJMAX) ? (n_jp - ih) : PhantomGrapeQuad::NJMAX;  
	    const PS::S32 it = ih + n_jp_tmp;
	    PS::S32 i_tmp = 0;
	    for(PS::S32 i=ih; i<it; i++, i_tmp++){
		const PS::F64 m_j = sp_j[i].getCharge();
		const PS::F64vec pos_j = sp_j[i].getPos();
                const PS::F64mat q = sp_j[i].quad;
#ifdef CALC_SP_64bit
                pg.set_spj_one_d(i, pos_j.x, pos_j.y, pos_j.z, m_j,
				 q.xx, q.yy, q.zz, q.xy, q.yz, q.xz);
#else
                pg.set_spj_one(i, pos_j.x, pos_j.y, pos_j.z, m_j,
                               q.xx, q.yy, q.zz, q.xy, q.yz, q.xz);
#endif
	    }
#ifdef CALC_SP_64bit
            pg.run_spj_d(n_ip, n_jp_tmp);
#else
            pg.run_spj(n_ip, n_jp_tmp);
#endif
	    for(PS::S32 i=0; i<n_ip; i++){
		PS::F64 * p = &(force[i].pot);
		PS::F64 * a = (PS::F64 * )(&force[i].acc[0]);
#ifdef CALC_SP_64bit
		pg.accum_accp_one_d(i, a[0], a[1], a[2], *p);
#else
		pg.accum_accp_one(i, a[0], a[1], a[2], *p);
#endif
	    }
	}
    }
};


// 32 bit
struct CalcForceEPEPNOSIMD{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F32 eps2 = EPISoft::eps * EPISoft::eps;
        //const PS::F32 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
	const PS::F32 r_crit2 = EPJSoft::r_search * EPJSoft::r_search * SAFTY_FACTOR_FOR_SEARCH_SQ; 
	const PS::F32 r_out = EPISoft::r_out; 
	const PS::F32 r_in = EPISoft::r_in;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F32vec xi = ep_i[i].pos;
            PS::S64 id_i = ep_i[i].id;
            PS::F32vec ai = 0.0;
            PS::F32 poti = 0.0;
	    PS::S32 n_ngb_i = 0;
            for(PS::S32 j=0; j<n_jp; j++){
                if(id_i == ep_j[j].id) continue;
                const PS::F32vec rij = xi - ep_j[j].pos;
                const PS::F32 r2 = rij * rij;
                const PS::F32 r2_eps = rij * rij + eps2;
                if(r2 < r_crit2){
		    n_ngb_i++;
                }
                const PS::F32 r_inv = 1.0/sqrt(r2_eps);
                const PS::F32 m_r = ep_j[j].mass * r_inv;
                const PS::F32 m_r3 = m_r * r_inv * r_inv;
		const PS::F32 r_eps = r2_eps * r_inv;
		const PS::F32 k = CalcK(r_eps, r_out, r_in);
                ai -= m_r3 * rij * k;
		poti -= m_r * k;
            }
            force[i].acc += ai;
            force[i].pot += poti;
            force[i].n_ngb = n_ngb_i;
        }
    }
};
//32bit
struct CalcForceEPSPQuadNOSIMD{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
		      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F32 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 ip=0; ip<n_ip; ip++){
            PS::F32vec xi = ep_i[ip].pos;
            PS::F32vec ai = 0.0;
            PS::F32 poti = 0.0;
            for(PS::S32 jp=0; jp<n_jp; jp++){
                PS::F32 mj = sp_j[jp].mass;
                PS::F32vec xj= sp_j[jp].pos;
                PS::F32vec rij= xi - xj;
                PS::F32 r2 = rij * rij + eps2;
                PS::F32mat qj = sp_j[jp].quad;
                PS::F32 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F32 qrr = qr * rij;
                PS::F32 r_inv = 1.0f/sqrt(r2);
                PS::F32 r2_inv = r_inv * r_inv;
                PS::F32 r3_inv = r2_inv * r_inv;
                PS::F32 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F32 qrr_r5 = r5_inv * qrr;
                PS::F32 qrr_r7 = r2_inv * qrr_r5;
                PS::F32 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F32 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[ip].acc += ai;
            force[ip].pot += poti;
        }
    }
};

//64 bit
struct CalcForceEPEPCutoff64bit{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        //const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
	const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search * SAFTY_FACTOR_FOR_SEARCH_SQ;
	const PS::F64 r_out = EPISoft::r_out; 
	const PS::F64 r_in = EPISoft::r_in;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::S64 id_i = ep_i[i].id;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
	    PS::S32 n_ngb_i = 0;
            for(PS::S32 j=0; j<n_jp; j++){
                if(id_i == ep_j[j].id) continue;
                const PS::F64vec rij = xi - ep_j[j].pos;
                const PS::F64 r2 = rij * rij;
                const PS::F64 r2_eps = rij * rij + eps2;
                if(r2 < r_crit2){
		    n_ngb_i++;
                }
                const PS::F64 r_inv = 1.0/sqrt(r2_eps);
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
		const PS::F64 r_eps = r2_eps * r_inv;
		const PS::F64 k = CalcK(r_eps, r_out, r_in);
                ai -= m_r3 * rij * k;
		poti -= m_r * k;
            }
            force[i].acc += ai;
            force[i].pot += poti;
            force[i].n_ngb = n_ngb_i;
        }
    }
};

// 64 bit
struct CalcForceEPSPQuadNOSIMD64{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
		      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 ip=0; ip<n_ip; ip++){
            PS::F64vec xi = ep_i[ip].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 jp=0; jp<n_jp; jp++){
                PS::F64 mj = sp_j[jp].mass;
                PS::F64vec xj= sp_j[jp].pos;
                PS::F64vec rij= xi - xj;
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = sp_j[jp].quad;
                PS::F64 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr = qr * rij;
                PS::F64 r_inv = 1.0f/sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[ip].acc += ai;
            force[ip].pot += poti;
        }
    }
};

#else //INTRINSIC

/*
struct CalcForceEPEP{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::S64 id_i = ep_i[i].id;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S32 & n_ngb_i = force[i].n_ngb;
	    //force[i].n_ep = n_jp;
            for(PS::S32 j=0; j<n_jp; j++){
                if(id_i == ep_j[j].id) continue;
                const PS::F64vec rij = xi - ep_j[j].pos;
                const PS::F64 r2 = rij * rij;
                const PS::F64 r2_eps = rij * rij + eps2;
                if(r2 < r_crit2){
                    if(id_i != ep_j[j].id){
                        n_ngb_i++;
                    }
                    continue;
                }
                const PS::F64 r_inv = 1.0/sqrt(r2_eps);
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                ai -= m_r3 * rij;
                poti -= m_r;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
*/

#if 1
//64 bit
struct CalcForceEPEP{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search * SAFTY_FACTOR_FOR_SEARCH_SQ;
        const PS::F64 r_out = EPISoft::r_out; 
        const PS::F64 r_in = EPISoft::r_in;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::S64 id_i = ep_i[i].id;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S32 n_ngb_i = 0;
            for(PS::S32 j=0; j<n_jp; j++){
                if(id_i == ep_j[j].id){
                    n_ngb_i++;
                    continue;
                }
                const PS::F64vec rij = xi - ep_j[j].pos;
                const PS::F64 r2 = rij * rij;
                const PS::F64 r2_eps = rij * rij + eps2;
                if(r2 < r_crit2){
                    n_ngb_i++;
                }
                const PS::F64 r_inv = 1.0/sqrt(r2_eps);
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                const PS::F64 r_eps = r2_eps * r_inv;
                const PS::F64 k = CalcK(r_eps, r_out, r_in);
                ai -= m_r3 * rij * k;
                poti -= m_r * k;
            }
            force[i].acc += ai;
            force[i].pot += poti;
            force[i].n_ngb = n_ngb_i;
        }
    }
};

#else
//32 bit
struct CalcForceEPEP{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F32 eps2 = EPISoft::eps * EPISoft::eps;
        //const PS::F32 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
	const PS::F32 r_crit2 = EPJSoft::r_search * EPJSoft::r_search * SAFTY_FACTOR_FOR_SEARCH_SQ;
	const PS::F32 r_out = EPISoft::r_out; 
	const PS::F32 r_in = EPISoft::r_in;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F32vec xi = ep_i[i].pos;
            PS::S64 id_i = ep_i[i].id;
            PS::F32vec ai = 0.0;
            PS::F32 poti = 0.0;
	    PS::S32 n_ngb_i = 0;
            for(PS::S32 j=0; j<n_jp; j++){
                if(id_i == ep_j[j].id) continue;
                const PS::F32vec rij = xi - ep_j[j].pos;
                const PS::F32 r2 = rij * rij;
                const PS::F32 r2_eps = rij * rij + eps2;
                if(r2 < r_crit2){
		    n_ngb_i++;
                }
                const PS::F32 r_inv = 1.0/sqrt(r2_eps);
                const PS::F32 m_r = ep_j[j].mass * r_inv;
                const PS::F32 m_r3 = m_r * r_inv * r_inv;
		const PS::F32 r_eps = r2_eps * r_inv;
		const PS::F32 k = CalcK(r_eps, r_out, r_in);
                ai -= m_r3 * rij * k;
		poti -= m_r * k;
            }
            force[i].acc += ai;
            force[i].pot += poti;
            force[i].n_ngb = n_ngb_i;
        }
    }
};
#endif


struct CalcForceEPSP{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
		      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].getPos();
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].getCharge();
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};


#if 0
// 64 bit
struct CalcForceEPSPQuad{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
		      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 ip=0; ip<n_ip; ip++){
            PS::F64vec xi = ep_i[ip].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 jp=0; jp<n_jp; jp++){
                PS::F64 mj = sp_j[jp].mass;
                PS::F64vec xj= sp_j[jp].pos;
                PS::F64vec rij= xi - xj;
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = sp_j[jp].quad;
                PS::F64 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr = qr * rij;
                PS::F64 r_inv = 1.0f/sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[ip].acc += ai;
            force[ip].pot += poti;
        }
    }
};
#else
// 32 bit
struct CalcForceEPSPQuad{
    template<class Tsp>
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
		      const Tsp * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F32 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 ip=0; ip<n_ip; ip++){
            PS::F32vec xi = ep_i[ip].pos;
            PS::F32vec ai = 0.0;
            PS::F32 poti = 0.0;
            for(PS::S32 jp=0; jp<n_jp; jp++){
                PS::F32 mj = sp_j[jp].mass;
                PS::F32vec xj= sp_j[jp].pos;
                PS::F32vec rij= xi - xj;
                PS::F32 r2 = rij * rij + eps2;
                PS::F32mat qj = sp_j[jp].quad;
                PS::F32 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F32 qrr = qr * rij;
                PS::F32 r_inv = 1.0f/sqrt(r2);
                PS::F32 r2_inv = r_inv * r_inv;
                PS::F32 r3_inv = r2_inv * r_inv;
                PS::F32 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F32 qrr_r5 = r5_inv * qrr;
                PS::F32 qrr_r7 = r2_inv * qrr_r5;
                PS::F32 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F32 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[ip].acc += ai;
            force[ip].pot += poti;
        }
    }
};
#endif





struct CalcForceEPSP2{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopole * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].getPos();
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].getCharge();
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};


struct CalcForceEPEPSimple{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::S64 id_i = ep_i[i].id;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            //PS::S32 & n_ngb_i = force[i].n_ngb;
            for(PS::S32 j=0; j<n_jp; j++){
                if(id_i == ep_j[j].id) continue;
                const PS::F64vec rij = xi - ep_j[j].pos;
                const PS::F64 r2_eps = rij * rij + eps2;
                const PS::F64 r_inv = 1.0/sqrt(r2_eps);
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                ai -= m_r3 * rij;
                poti -= m_r;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};


struct CalcForceEPSPSimple{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopole * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                const PS::F64vec rij = xi - sp_j[j].pos;
                const PS::F64 r2_eps = rij * rij + eps2;
                const PS::F64 r_inv = 1.0/sqrt(r2_eps);
                const PS::F64 m_r = sp_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                ai -= m_r3 * rij;
                poti -= m_r;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
	/*
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].getPos();
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].getCharge();
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
	*/
    }
};

#endif //INTRINSIC

/*
struct CalcForceEPEP0{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const EPJSoft * ep_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        const PS::F64 r_crit2 = EPJSoft::r_search * EPJSoft::r_search;
        for(PS::S32 i=0; i<n_ip; i++){
            const PS::F64vec xi = ep_i[i].pos;
            PS::S64 id_i = ep_i[i].id;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            PS::S32 & n_ngb_i = force[i].n_ngb;
            //PS::S32 & id_ngb_i = force[i].id_ngb;
            for(PS::S32 j=0; j<n_jp; j++){
                if(id_i == ep_j[j].id) continue;
                const PS::F64vec rij = xi - ep_j[j].pos;
                const PS::F64 r2 = rij * rij;
                const PS::F64 r2_eps = rij * rij + eps2;
                if(r2 < r_crit2){
                    //if(r2_eps < r_crit2){
                    if(id_i != ep_j[j].id){
                        //id_ngb_i = ep_j[j].id;
                        n_ngb_i++;
                    }
                    continue;
                }
                const PS::F64 r_inv = 1.0/sqrt(r2_eps);
                const PS::F64 m_r = ep_j[j].mass * r_inv;
                const PS::F64 m_r3 = m_r * r_inv * r_inv;
                ai -= m_r3 * rij;
                poti -= m_r;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};


struct CalcForceEPSP0{
    void operator () (const EPISoft * ep_i,
                      const PS::S32 n_ip,
                      const PS::SPJMonopoleScatter * sp_j,
                      const PS::S32 n_jp,
                      ForceSoft * force){
        const PS::F64 eps2 = EPISoft::eps * EPISoft::eps;
        for(PS::S32 i=0; i<n_ip; i++){
            PS::F64vec xi = ep_i[i].pos;
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(PS::S32 j=0; j<n_jp; j++){
                PS::F64vec rij = xi - sp_j[j].getPos();
                PS::F64 r3_inv = rij * rij + eps2;
                PS::F64 r_inv = 1.0/sqrt(r3_inv);
                r3_inv = r_inv * r_inv;
                r_inv *= sp_j[j].getCharge();
                r3_inv *= r_inv;
                ai -= r3_inv * rij;
                poti -= r_inv;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};
*/

class PTCLPred{
public:
    PS::F64vec pos;
    PS::F64vec vel;
    void dump(std::ostream & fout=std::cout){
        fout<<"pos="<<pos<<std::endl;
        fout<<"vel="<<vel<<std::endl;
    }
};

class PTCLForce{
public:
    PS::F64vec acc0;
    PS::F64vec acc1;
    PS::F64vec acc0_pla;
    PS::F64vec acc1_pla;
    //PS::F64 pot;
    void clear(){
        acc0 = acc1 = 0.0;
        acc0_pla = acc1_pla = 0.0;
        //pot = 0.0;
    }
    void reduce(){
	acc0 += acc0_pla;
	acc1 += acc1_pla;
    }
    void dump(std::ostream & fout=std::cout){
        fout<<"acc0="<<acc0<<std::endl;
        fout<<"acc1="<<acc1<<std::endl;
    }
};



class PTCLHard{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64 time;
    PS::F64 dt;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc0;
    PS::F64vec acc1;
    PS::F64vec acc0_pla;
    PS::F64vec acc1_pla;
    PS::S32 n_ngb;
    PS::S32 adr_pair;
    PS::F64 pot; // debug
    PS::F64 r_merge;
    PS::S32 rank_org;
    PS::S32 adr_fp;
    static PS::F64 r_factor;
    static PS::F64 dens;

    PTCLHard(){
        //id = n_ngb = adr_pair = rank_org = adr_fp = -1;
    }

    /*
    PTCLHard(const FPSoft & fp, const PS::S32 _adr_fp){
        id = fp.id;
        mass = fp.mass;
        pos = fp.pos;
        vel = fp.vel;
	rank_org = fp.rank_org;
	adr_fp = _adr_fp;
    }
    */

    PTCLHard(const PS::S64 _id, 
	     const PS::F64 _mass, 
	     const PS::F64vec & _pos, 
	     const PS::F64vec & _vel, 
	     const PS::S32 _n_ngb,
	     const PS::S32 _rank_org,
	     const PS::S32 _adr_fp){
        id = _id;
        mass = _mass;
        pos = _pos;
        vel = _vel;
        n_ngb = _n_ngb;
        rank_org = _rank_org;
        adr_fp = _adr_fp;
    }


    void setRMerge(){
#ifdef REMOVE_BODY
	PS::F64 mass_sun = 1.0;
	PS::F64 ax = 1.0;
	PS::F64 factor = 1.1;
	r_merge = cbrt(2.0*mass/3.0*mass_sun)*ax*factor;
#else
	static const PS::F64 rho_ave = dens * ( (1.49597871e13*1.49597871e13*1.49597871e13) / 1.989e33); // [Msun/AU^3]
        static const PS::F64 PI = 4.0*atan(1.0);
        static const PS::F64 C = 3.0/(4.0*PI*rho_ave);
	r_merge = cbrt(C*mass) * r_factor; // correct
#endif
    }
    /*
    void copyFromFP(const FPSoft & fp, const PS::S32 _adr_fp){
        id = fp.id;
        mass = fp.mass;
        pos = fp.pos;
        vel = fp.vel;
	rank_org = fp.rank_org;
	adr_fp = _adr_fp;
    }
    */
    void dump(std::ostream & fout=std::cout){
/*
        fout<<"id="<<id<<std::endl;
        fout<<"n_ngb="<<n_ngb<<std::endl;
*/
        fout<<"dt="<<dt<<std::endl;
        fout<<"time="<<time<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"pos.x="<<pos.x<<std::endl;
        fout<<"mass="<<mass<<std::endl;
/*
        fout<<"vel="<<vel<<std::endl;
        fout<<"acc0="<<acc0<<std::endl;
        fout<<"acc1="<<acc1<<std::endl;
        fout<<"n_ngb="<<n_ngb<<std::endl;
        fout<<"rank_org="<<rank_org<<std::endl;
        fout<<"adr_fp="<<adr_fp<<std::endl;
*/
    }

    static PS::F64 calcDt2nd(const PS::F64vec a0,
			     const PS::F64vec a1,
                             const PS::F64 eta,
			     const PS::F64 a0_offset_sq=0.0){
        const PS::F64 s0 = a0 * a0 + a0_offset_sq;
        const PS::F64 s1 = a1 * a1;
	if(s0 == a0_offset_sq || s1 == 0.0){
            return PS::LARGE_FLOAT;
        }
        else{
	    return eta * sqrt(s0 / s1);
        }
    }

    void setDt2nd(const PTCLForce & force, 
                  const PS::F64 eta, 
                  const PS::F64 dt_limit,
		  const PS::F64 a0_offset_sq=0.0){
        //const PS::F64 dt_ref = calcDt2nd(force.acc0, force.acc1, eta, a0_offset_sq);
	//const PS::F64 dt_ref = calcDt2nd(force.acc0_pla, force.acc1_pla, eta, a0_offset_sq);
	const PS::F64 dt_ref = std::min( calcDt2nd(force.acc0_pla, force.acc1_pla, eta, a0_offset_sq), calcDt2nd(force.acc0, force.acc1, eta) );
        this->dt = dt_limit;
#if 1
        while(this->dt > dt_ref) this->dt *= 0.5;
#elif 1
	this->dt = 1.0/16384;
#else
	this->dt = dt_limit;
#endif
    }

    static PS::F64 calcDt4th(const PS::F64vec a0, const PS::F64vec a1,
                             const PS::F64vec a2, const PS::F64vec a3,
                             const PS::F64 eta,   const PS::F64 a0_offset_sq=0.0){
        const PS::F64 s0 = a0 * a0 + a0_offset_sq;
        const PS::F64 s1 = a1 * a1;
        const PS::F64 s2 = a2 * a2;
        const PS::F64 s3 = a3 * a3;
#if 1
        if(s0 == a0_offset_sq || s1 == 0.0){
            return PS::LARGE_FLOAT;
        }
        else{
            return eta * sqrt( (sqrt(s0*s2) + s1) / (sqrt(s1*s3) + s2) );
        }
#elif 1
	return 1.0/16384;
#else
	return PS::LARGE_FLOAT;
#endif
    }

    void merge(PTCLHard & ptcl_del){
#ifdef REMOVE_BODY
	// do notihing
#else
        pos = mass*pos + ptcl_del.mass*ptcl_del.pos;
        vel = mass*vel + ptcl_del.mass*ptcl_del.vel;
        mass = mass + ptcl_del.mass;
        pos /= mass;
        vel /= mass;
#endif
        setRMerge();
    }
    
    void correct(const PTCLForce & force,
                 const PS::F64 eta,
                 const PS::F64 dt_limit,
		 const PS::F64 a0_offset_sq=0.0){
        static const PS::F64 inv3 = 1.0 / 3.0;
        const PS::F64 h = 0.5 * dt;
        const PS::F64 hinv = 2.0 / dt;
        const PS::F64vec A0p = (force.acc0 + this->acc0);
        const PS::F64vec A0m = (force.acc0 - this->acc0);
        const PS::F64vec A1p = (force.acc1 + this->acc1)*h;
        const PS::F64vec A1m = (force.acc1 - this->acc1)*h;
        const PS::F64vec vel_new = this->vel + h*( A0p - inv3*A1m );
        this->pos += h*( (this->vel + vel_new) + h*(-inv3*A0m));
        this->vel = vel_new;
        this->acc0 = force.acc0;
        this->acc1 = force.acc1;
#ifdef FORDEBUG
        this->pot = force.pot;
#endif
        this->time += dt;
        const PS::F64vec acc3 = (1.5*hinv*hinv*hinv) * (A1p - A0m);
        const PS::F64vec acc2 = (0.5*hinv*hinv) * A1m + h*acc3;
        //const PS::F64 dt_ref = calcDt4th(this->acc0, this->acc1, acc2, acc3, eta, a0_offset_sq);
        const PS::F64vec A0m_pla = (force.acc0_pla - this->acc0_pla);
        const PS::F64vec A1p_pla = (force.acc1_pla + this->acc1_pla)*h;
        const PS::F64vec A1m_pla = (force.acc1_pla - this->acc1_pla)*h;
        const PS::F64vec acc3_pla = (1.5*hinv*hinv*hinv) * (A1p_pla - A0m_pla);
        const PS::F64vec acc2_pla = (0.5*hinv*hinv) * A1m_pla + h*acc3_pla;
	//const PS::F64 dt_ref = calcDt4th(this->acc0_pla, this->acc1_pla, acc2_pla, acc3_pla, eta, a0_offset_sq);
	const PS::F64 dt_ref = std::min( calcDt4th(this->acc0_pla, this->acc1_pla, acc2_pla, acc3_pla, eta, a0_offset_sq), calcDt4th(this->acc0, this->acc1, acc2, acc3, eta) ); 
        this->acc0_pla = force.acc0_pla;
        this->acc1_pla = force.acc1_pla;
        const PS::F64 dt_old = this->dt;
        assert(dt_old != 0.0);
        this->dt = dt_limit;
        while(this->dt > dt_ref) this->dt *= 0.5;
        this->dt = dt_old*2 < this->dt ?  dt_old*2 : this->dt;
    }
    bool isDead(){
	return (mass > 0.0) ? false : true;
    }
};

PS::F64 PTCLHard::r_factor = 1.0;
PS::F64 PTCLHard::dens = 2.0;

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

class PTCLHardComm{
public:
    PS::S64 id;
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::S32 n_ngb;
    PS::S32 rank_org;
    PS::S32 adr_fp;
    PTCLHardComm(){}
    PTCLHardComm(const PTCLHard & ph): 
	id(ph.id), 
	mass(ph.mass),
	pos(ph.pos),
	vel(ph.vel),
	n_ngb(ph.n_ngb),
	rank_org(ph.rank_org),
	adr_fp(ph.adr_fp) { }
    void copyFromFP(const FPSoft & fp){
	id = fp.id;
	mass = fp.mass;
	pos = fp.pos;
	vel = fp.vel;
	n_ngb = fp.n_ngb;
	rank_org = fp.rank_org;
	adr_fp = fp.adr;
    }
    void copyFromPH(const PTCLHard & ph){
	id = ph.id;
	mass = ph.mass;
	pos = ph.pos;
	vel = ph.vel;
	n_ngb = ph.n_ngb;
	rank_org = ph.rank_org;
	adr_fp = ph.adr_fp;
    }
    void copyFromPHC(const PTCLHardComm & ph){
	id = ph.id;
	mass = ph.mass;
	pos = ph.pos;
	vel = ph.vel;
	n_ngb = ph.n_ngb;
	rank_org = ph.rank_org;
	adr_fp = ph.adr_fp;
    }
    /*
    PTCLHardComm(const PS::S64 _id, 
		 const PS::F64 _mass, 
		 const PS::F64vec & _pos, 
		 const PS::F64vec & _vel, 
		 const PS::S32 _n_ngb,
		 const PS::S32 _rank_org,
		 const PS::S32 _adr_fp){
	id = _id;
	mass = _mass;
	pos = _pos;
	vel = _vel;
	n_ngb = _n_ngb;
	rank_org = _rank_org;
	adr_fp = _adr_fp;
    }
    */
    /*
    PTCLHardComm(const PTCLHardComm & phc){
	id = phc.id;
	mass = phc.mass;
	pos = phc.pos;
	vel = phc.vel;
	n_ngb = phc.n_ngb;
	rank_org = phc.rank_org;
	adr_fp = phc.adr_fp;
    }
    */
    /*
    const PTCLHardComm & operator = (const PTCLHardComm & phc){
	id = phc.id;
	mass = phc.mass;
	pos = phc.pos;
	vel = phc.vel;
	n_ngb = phc.n_ngb;
	rank_org = phc.rank_org;
	adr_fp = phc.adr_fp;
	return (*this)
    }
    */
};


class MergeLog{
public:
    PS::F64 time;
    PTCLHard ptcl_merged;
    PTCLHard ptcl_dead;
    MergeLog(){
	time = -1.0;
    }
    MergeLog(const PS::F64 t, const PTCLHard & p_m, const PTCLHard & p_d){
	time = t;
	ptcl_merged = p_m;
	ptcl_dead = p_d;
    }
    void dump(std::ostream & fout=std::cout){
        fout<<"time= "<<time<<std::endl;
        fout<<"ptcl_merged.id= "      <<ptcl_merged.id<<std::endl;
        fout<<"ptcl_merged.mass= "    <<ptcl_merged.mass<<std::endl;
        fout<<"ptcl_merged.pos= "     <<ptcl_merged.pos<<std::endl;
        fout<<"ptcl_merged.vel= "     <<ptcl_merged.vel<<std::endl;
        fout<<"ptcl_merged.rank_org= "<<ptcl_merged.rank_org<<std::endl;
        fout<<"ptcl_dead.id= "      <<ptcl_dead.id<<std::endl;
        fout<<"ptcl_dead.mass= "    <<ptcl_dead.mass<<std::endl;
        fout<<"ptcl_dead.pos= "     <<ptcl_dead.pos<<std::endl;
        fout<<"ptcl_dead.vel= "     <<ptcl_dead.vel<<std::endl;
        fout<<"ptcl_dead.rank_org= "<<ptcl_dead.rank_org<<std::endl;
    }
};
