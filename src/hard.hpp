/*
system_soft[i].n_epj_: without particle itself
pair_hadr_loc_: without particle itself
epj_ngb_arra_; without particle itself
dead particle has zero-mass
*/

#define PARALLEL_2BODY
//#define ORIGINAL

#ifdef USE_INTRINSIC_FOR_X86
#include<immintrin.h>
#endif

#include"kepler.hpp"

std::ofstream fout_debug;

template<class T>
void Print(const T str, std::ostream & fout);

class SortAdr{
public:
    static PS::F64 * time;
    bool operator() (const PS::S32 & left, const PS::S32 & right) const {
        return time[left] < time[right];
    }
};

PS::F64 * SortAdr::time;

// in the case of merge, dt_limit must be changed. 
inline PS::F64 CalcDtLimit(const PS::F64 time_sys, 
                           const PS::F64 time_sync, 
                           const PS::F64 dt_limit_org){
    PS::F64 dt_limit_ret = dt_limit_org;
    PS::F64 s = time_sys / dt_limit_ret;
    while(s != PS::F64(PS::S64(s))){
        s *= 2.0;
        dt_limit_ret *= 0.5;
    }
    return dt_limit_ret;
}

#include"force.hpp"
class HardSystem{
private:
    void CalcTwoBodyEnergy(Energy & eng,
			   const PS::F64 r_in,
			   const PS::F64 r_out){
	eng.clear();
        const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;
	PS::S32 n_pair = pair_adr_2body_loc_.size();
        for(PS::S32 np=0; np<n_pair; np++){
            const PS::S64 adr_i = pair_adr_2body_loc_[np].first;
            const PS::S64 adr_j = pair_adr_2body_loc_[np].second;
            if(adr_i > adr_j) continue;
            const PTCLHard & ptcl_i = ptcl_2body_loc_[adr_i];
            const PTCLHard & ptcl_j = ptcl_2body_loc_[adr_j];
	    if(ptcl_i.mass > 0.0){
		eng.kin += 0.5*ptcl_i.mass*ptcl_i.vel*ptcl_i.vel;
		const PS::F64vec ri = ptcl_i.pos - pos_sun_;
		eng.pot -= ptcl_i.mass*mass_sun_*sqrt(1.0/(ri*ri));
	    }
	    if(ptcl_j.mass > 0.0){
		eng.kin += 0.5*ptcl_j.mass*ptcl_j.vel*ptcl_j.vel;
		const PS::F64vec rj = ptcl_j.pos - pos_sun_;
		eng.pot -= ptcl_j.mass*mass_sun_*sqrt(1.0/(rj*rj));
	    }
	    if(ptcl_i.mass > 0.0 && ptcl_j.mass > 0.0){
		const PS::F64vec rij = ptcl_i.pos - ptcl_j.pos;
		const PS::F64 q = r_in / r_out;
		const PS::F64 y = sqrt(rij*rij+eps_sq) / r_out;
		eng.pot -= (1.0-CalcW(y,q)) * (ptcl_i.mass*ptcl_j.mass*sqrt(1.0/(rij*rij+eps_sq)));
	    }
	}
	eng.tot = eng.kin + eng.pot;
    }

    void CalcMultiBodyEnergy(Energy & eng,
			     const PS::F64 r_in,
			     const PS::F64 r_out){
	eng.clear();
        const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;
        const PS::S32 ni = ptcl_multi_glb_.size();
	for(PS::S32 i=0; i<ni; i++){
	    PTCLHard & ptcl_i = ptcl_multi_glb_[i];
	    if(ptcl_i.mass <= 0.0) continue;
	    PS::F64vec ri = ptcl_i.pos - pos_sun_;
	    eng.kin += 0.5*ptcl_i.mass*ptcl_i.vel*ptcl_i.vel;
	    eng.pot -= mass_sun_*ptcl_i.mass*sqrt(1.0/(ri*ri));
            const PS::S32 nj = ptcl_i.n_ngb;
            const PS::S32 adr_pair = ptcl_i.adr_pair;
            for(PS::S32 j=0; j<nj; j++){
                const PS::S32 adr_j = pair_adr_multi_glb_[adr_pair+j].second;
		PTCLHard & ptcl_j = ptcl_multi_glb_[adr_j];
		if(ptcl_j.mass <= 0.0) continue;
		PS::F64vec rij = ptcl_i.pos - ptcl_j.pos;
		PS::F64 q = r_in / r_out;
		PS::F64 y = sqrt(rij*rij+eps_sq) / r_out;
		eng.pot -= (1.0-CalcW(y,q)) * (0.5*ptcl_i.mass*ptcl_j.mass*sqrt(1.0/(rij*rij+eps_sq)));
	    }
	}
	eng.tot = eng.kin + eng.pot;
    }

    size_t n_interaction_2body_loc_;
    size_t n_interaction_2body_glb_;
    size_t n_interaction_multi_loc_;
    size_t n_interaction_multi_glb_;

    void makeDictionaryGlobalImpl(std::vector<PTCLHard> & ptcl,
				  const std::vector< std::pair<PS::S64, PS::S64> > & pair_id,
				  std::vector< std::pair<PS::S64, PS::S64> > & pair_adr){
	static std::unordered_map<PS::S64, PS::S64> idx_to_adr;
        idx_to_adr.clear();
        pair_adr.clear();
        const PS::S32 n_ptcl = ptcl.size();
        for(PS::S32 i=0; i<n_ptcl; i++){
            idx_to_adr.insert( std::pair<PS::S64, PS::S64>(ptcl[i].id, i) );
        }
        const PS::S32 n_interaction = pair_id.size();
        PS::S32 first_adr_old = -1;
        for(PS::S32 i=0; i<n_interaction; i++){
            //std::unordered_map<PS::S64, PS::S64>::iterator itr = idx_to_adr.find(pair_id[i].second);
            const PS::S64 first_adr = idx_to_adr[pair_id[i].first];
            if(first_adr_old != first_adr){
                ptcl[first_adr].adr_pair = i;
                first_adr_old = first_adr;
            }
#if 1
	    assert( idx_to_adr.find(pair_id[i].second) != idx_to_adr.end());
	    pair_adr.push_back( std::pair<PS::S64, PS::S64>(first_adr,
							    idx_to_adr[pair_id[i].second]));
#else
            if(itr == idx_to_adr.end()){
		// no correspondig particle
		std::cout<<"global"<<std::endl;
		std::cout<<"i="<<i<<std::endl;
		std::cout<<"n_interaction="<<n_interaction<<std::endl;
		std::cout<<"first_adr="<<first_adr<<std::endl;
		std::cout<<"pair_id[i].first="<<pair_id[i].first<<std::endl;
		std::cout<<"pair_id[i].second="<<pair_id[i].second<<std::endl;
                pair_adr.push_back( std::pair<PS::S64, PS::S64>(first_adr, -1));
            }
            else{
                pair_adr.push_back( std::pair<PS::S64, PS::S64>(first_adr,
                                                                idx_to_adr[pair_id[i].second]));
            }
#endif
        }
    }

#if 1
    // construct pair_adr from pair_id
    // not "make Dicitionay"
    void makeDictionaryLocalImpl(std::vector<PTCLHard> & ptcl,
				 const std::vector< std::pair<PS::S64, PS::S64> > & pair_id,
				 const std::vector< std::pair<PS::S32, PS::S32> > & pair_rank,
				 std::vector< std::pair<PS::S64, PS::S64> > & pair_adr){
	static std::unordered_map<PS::S64, PS::S64> idx_to_adr;
	const PS::S32 my_rank = PS::Comm::getRank();
        idx_to_adr.clear();
        pair_adr.clear();
        const PS::S32 n_ptcl = ptcl.size();
        for(PS::S32 i=0; i<n_ptcl; i++){
            idx_to_adr.insert( std::pair<PS::S64, PS::S64>(ptcl[i].id, i) );
        }
        const PS::S32 n_interaction = pair_id.size();
        PS::S32 first_adr_old = -1;
        for(PS::S32 i=0; i<n_interaction; i++){
            //std::unordered_map<PS::S64, PS::S64>::iterator itr = idx_to_adr.find(pair_id[i].second);
            const PS::S64 first_adr = idx_to_adr[pair_id[i].first];
            if(first_adr_old != first_adr){
                ptcl[first_adr].adr_pair = i;
                first_adr_old = first_adr;
            }
#if 1 
	    if(pair_rank[i].second == my_rank){
                pair_adr.push_back( std::pair<PS::S64, PS::S64>(first_adr,
                                                                idx_to_adr[pair_id[i].second]));
	    }
	    else{
                pair_adr.push_back( std::pair<PS::S64, PS::S64>(first_adr, -1));
	    }
#else
            if(itr == idx_to_adr.end()){
		// no correspondig particle
		/*
		std::cout<<"local"<<std::endl;
		std::cout<<"i="<<i<<std::endl;
		std::cout<<"n_interaction="<<n_interaction<<std::endl;
		std::cout<<"first_adr="<<first_adr<<std::endl;
		std::cout<<"pair_id[i].first="<<pair_id[i].first<<std::endl;
		std::cout<<"pair_id[i].second="<<pair_id[i].second<<std::endl;
		*/
                pair_adr.push_back( std::pair<PS::S64, PS::S64>(first_adr, -1));
            }
            else{
                pair_adr.push_back( std::pair<PS::S64, PS::S64>(first_adr,
                                                                idx_to_adr[pair_id[i].second]));
            }
#endif
        }
    }
#endif

    PS::F64 mass_sun_;
    PS::F64vec pos_sun_;
    PS::F64vec vel_sun_;

#if 1
    template<class Tptcl>
    void merge2body(Tptcl & ptcl0,
		    Tptcl & ptcl1, 
		    PS::F64 & eng_disp, 
		    std::vector<MergeLog> & merge_history){
        const PS::F64 m0 = ptcl0.mass;
        const PS::F64 m1 = ptcl1.mass;
        const PS::S64 id0 = ptcl0.id;
        const PS::S64 id1 = ptcl1.id;

        Tptcl & ptcl_merge = (m0>m1) ? ptcl0 : ( (m1>m0) ? ptcl1 : ((id0<id1) ? ptcl0 : ptcl1) );
        Tptcl & ptcl_dead = (&ptcl_merge == &ptcl0) ? ptcl1 : ptcl0;

        const PS::F64vec r01 = ptcl0.pos - ptcl1.pos;
        const PS::F64vec v01 = ptcl0.vel - ptcl1.vel;
        const PS::F64vec r0 = ptcl0.pos - pos_sun_;
        const PS::F64vec r1 = ptcl1.pos - pos_sun_;
        const PS::F64 pot0 = -mass_sun_ * m0 / sqrt(r0*r0);
        const PS::F64 pot1 = -mass_sun_ * m1 / sqrt(r1*r1);
        const PS::F64 pot_old = pot0 + pot1 - m0*m1/sqrt(r01*r01);
        const PS::F64 m_red = (m0*m1) / (m0+m1);

        const PS::F64 merge_time = ptcl0.time + time_origin_;
        merge_history.push_back( MergeLog(merge_time, ptcl_merge, ptcl_dead) );
        ptcl_merge.merge(ptcl_dead);
        ptcl_dead.mass = 0.0;

        const PS::F64vec r_new = ptcl_merge.pos - pos_sun_;
        const PS::F64 pot_new = -mass_sun_ * ptcl_merge.mass / sqrt(r_new*r_new);
        eng_disp += 0.5*m_red*v01*v01 - (pot_new - pot_old);
    }
#endif

#if 0
    template<class Tptcl>
    void merge2bodyWithHardEnergy(Tptcl & ptcl0, 
				  Tptcl & ptcl1, 
				  PS::F64 & eng_disp, 
				  std::vector<MergeLog> & merge_history,
				  const PS::F64 r_in,
				  const PS::F64 r_out,
				  PS::F64 & eng_hard_disp){

        const PS::F64 m0 = ptcl0.mass;
        const PS::F64 m1 = ptcl1.mass;
        const PS::S64 id0 = ptcl0.id;
        const PS::S64 id1 = ptcl1.id;

        Tptcl & ptcl_merge = (m0>m1) ? ptcl0 : ( (m1>m0) ? ptcl1 : ((id0<id1) ? ptcl0 : ptcl1) );
        Tptcl & ptcl_dead = (&ptcl_merge == &ptcl0) ? ptcl1 : ptcl0;

        const PS::F64vec r01 = ptcl0.pos - ptcl1.pos;
        const PS::F64vec v01 = ptcl0.vel - ptcl1.vel;
        const PS::F64vec r0 = ptcl0.pos - pos_sun_;
        const PS::F64vec r1 = ptcl1.pos - pos_sun_;
        const PS::F64 pot0 = -mass_sun_ * m0 / sqrt(r0*r0);
        const PS::F64 pot1 = -mass_sun_ * m1 / sqrt(r1*r1);
        const PS::F64 pot_old = pot0 + pot1 - m0*m1/sqrt(r01*r01);
        const PS::F64 m_red = (m0*m1) / (m0+m1);

        const PS::F64 merge_time = ptcl0.time + time_origin_;
        merge_history.push_back( MergeLog(merge_time, ptcl_merge, ptcl_dead) );
        ptcl_merge.merge(ptcl_dead);
        ptcl_dead.mass = 0.0;

        const PS::F64vec r_new = ptcl_merge.pos - pos_sun_;
        const PS::F64 pot_new = -mass_sun_ * ptcl_merge.mass / sqrt(r_new*r_new);
        eng_disp += 0.5*m_red*v01*v01 - (pot_new - pot_old);

	const PS::F64 q = r_in/r_out;
	const PS::F64 y = sqrt(r01*r01)/r_out;
        const PS::F64 pot_hard_old = pot0 + pot1 - (1-CalcW(y,q))*m0*m1/sqrt(r01*r01);
        eng_hard_disp += 0.5*m_red*v01*v01 - (pot_new - pot_hard_old);
    }
#endif

#if 1
    template<class Tptcl>
    void merge2body(std::vector<Tptcl> & ptcl, 
                    const std::vector< std::pair<PS::S64, PS::S64> > pair_adr,
                    const PS::S64 adr0, 
		    const PS::S64 adr1, 
                    PS::F64 & eng_disp, 
		    std::vector<MergeLog> & merge_history){
        Tptcl & ptcl0 = ptcl[adr0];
        Tptcl & ptcl1 = ptcl[adr1];
        const PS::F64 m0 = ptcl0.mass;
        const PS::F64 m1 = ptcl1.mass;
        const PS::S64 id0 = ptcl0.id;
        const PS::S64 id1 = ptcl1.id;
        Tptcl & ptcl_merge = (m0>m1) ? ptcl0 : ( (m1>m0) ? ptcl1 : ((id0<id1) ? ptcl0 : ptcl1) );
        Tptcl & ptcl_dead = (&ptcl_merge == &ptcl0) ? ptcl1 : ptcl0;

        std::set<PS::S64> adr_jp;
        adr_jp.clear();

        const PS::S32 n0 = ptcl0.n_ngb;
        const PS::S32 adr_pair_ngb_0 = ptcl0.adr_pair;
        for(PS::S32 j=0; j<n0; j++){
            PS::S32 adr = pair_adr[adr_pair_ngb_0+j].second;
            const Tptcl & ptclj = ptcl[adr];
            if( ptcl1.id == ptclj.id) continue;
            adr_jp.insert(adr);
        }
        const PS::S32 n1 = ptcl1.n_ngb;
        const PS::S32 adr_pair_ngb_1 = ptcl1.adr_pair;
        for(PS::S32 j=0; j<n1; j++){
            PS::S32 adr = pair_adr[adr_pair_ngb_1+j].second;
            const Tptcl & ptclj = ptcl[adr];
            if( ptcl0.id == ptclj.id) continue;
            adr_jp.insert(adr);
        }
        PS::F64 pot_old = 0.0;
        std::set<PS::S64>::iterator adr_itr = adr_jp.begin();
#if 0
        const size_t n_ngb = adr_jp.size();
        for(size_t j=0; j<n_ngb; j++, adr_itr++){
            const Tptcl & ptclj = ptcl[*adr_itr];
            //std::cerr<<"---------------"<<std::endl;
            //std::cerr<<"---------------"<<std::endl;
            if(ptcl0.id != ptclj.id && ptcl1.id != ptclj.id){
                //std::cerr<<"ptcl0.id="<<ptcl0.id<<" ptcl1.id="<<ptcl1.id<<" ptclj.id="<<ptclj.id<<std::endl;
                //std::cerr<<"ptcl0.pos="<<ptcl0.pos<<" ptcl1.pos="<<ptcl1.pos<<" ptclj.pos="<<ptclj.pos<<std::endl;
                const PS::F64vec dr0 = ptcl0.pos - ptclj.pos;
                pot_old -= (ptcl0.mass * ptclj.mass) / sqrt(dr0*dr0);
                const PS::F64vec dr1 = ptcl1.pos - ptclj.pos;
                pot_old -= (ptcl1.mass * ptclj.mass) / sqrt(dr1*dr1);
                //std::cerr<<"dr0="<<dr0<<" dr1="<<dr1<<std::endl;
            }
            //std::cerr<<"---------------"<<std::endl;
            //std::cerr<<"---------------"<<std::endl;
        }
#endif

        const PS::F64vec r0 = ptcl0.pos - pos_sun_;
        pot_old -= mass_sun_ * m0 / sqrt(r0*r0);
        const PS::F64vec r1 = ptcl1.pos - pos_sun_;
        pot_old -= mass_sun_ * m1 / sqrt(r1*r1);
        const PS::F64vec r01 = ptcl0.pos - ptcl1.pos;
        pot_old -= m0*m1/sqrt(r01*r01);
        const PS::F64vec v01 = ptcl0.vel - ptcl1.vel;
        const PS::F64 m_red = (m0*m1) / (m0+m1);
        PS::F64 kin = 0.5*m_red*v01*v01;

        const PS::F64 merge_time = ptcl0.time + time_origin_;
        merge_history.push_back( MergeLog(merge_time, ptcl_merge, ptcl_dead) );
        ptcl_merge.merge(ptcl_dead);
        ptcl_dead.mass = 0.0;

        const PS::F64vec r_new = ptcl_merge.pos - pos_sun_;
        PS::F64 pot_new = -mass_sun_ * ptcl_merge.mass / sqrt(r_new*r_new); 
        adr_itr = adr_jp.begin();
#if 0
        //std::cerr<<"---------------"<<std::endl;
        //std::cerr<<"---------------"<<std::endl;
        for(size_t j=0; j<n_ngb; j++, adr_itr++){
            Tptcl & ptclj = ptcl[*adr_itr];
            if(ptcl0.id != ptclj.id && ptcl1.id != ptclj.id){
                const PS::F64vec dr = ptcl_merge.pos - ptclj.pos;
                pot_new -= (ptcl_merge.mass * ptclj.mass) / sqrt(dr*dr);
                //std::cerr<<"dr="<<dr<<std::endl;
            }
        }
#endif
        //std::cerr<<"---------------"<<std::endl;
        //std::cerr<<"---------------"<<std::endl;
        eng_disp += kin - (pot_new - pot_old);
    }
#endif

    template<class Tptcl>
    void merge2bodyWithHardEnergy(std::vector<Tptcl> & ptcl, 
				  const std::vector< std::pair<PS::S64, PS::S64> > pair_adr,
				  const PS::S64 adr0, 
				  const PS::S64 adr1, 
				  PS::F64 & eng_disp, 
				  std::vector<MergeLog> & merge_history,
				  const PS::F64 r_in,
				  const PS::F64 r_out,
				  PS::F64 & eng_hard_disp){
	const PS::F64 q = r_in / r_out;
        Tptcl & ptcl0 = ptcl[adr0];
        Tptcl & ptcl1 = ptcl[adr1];
        const PS::F64 m0 = ptcl0.mass;
        const PS::F64 m1 = ptcl1.mass;
        const PS::S64 id0 = ptcl0.id;
        const PS::S64 id1 = ptcl1.id;
        Tptcl & ptcl_merge = (m0>m1) ? ptcl0 : ( (m1>m0) ? ptcl1 : ((id0<id1) ? ptcl0 : ptcl1) );
        Tptcl & ptcl_dead = (&ptcl_merge == &ptcl0) ? ptcl1 : ptcl0;

        std::set<PS::S64> adr_jp;
        adr_jp.clear();

        const PS::S32 n0 = ptcl0.n_ngb;
        const PS::S32 adr_pair_ngb_0 = ptcl0.adr_pair;
        for(PS::S32 j=0; j<n0; j++){
            PS::S32 adr = pair_adr[adr_pair_ngb_0+j].second;
            const Tptcl & ptclj = ptcl[adr];
            if( ptcl1.id == ptclj.id) continue;
            adr_jp.insert(adr);
        }
        const PS::S32 n1 = ptcl1.n_ngb;
        const PS::S32 adr_pair_ngb_1 = ptcl1.adr_pair;
        for(PS::S32 j=0; j<n1; j++){
            PS::S32 adr = pair_adr[adr_pair_ngb_1+j].second;
            const Tptcl & ptclj = ptcl[adr];
            if( ptcl0.id == ptclj.id) continue;
            adr_jp.insert(adr);
        }
        std::set<PS::S64>::iterator adr_itr = adr_jp.begin();

	PS::F64 pot_old = 0.0;
	PS::F64 pot_hard_old = 0.0;

        const size_t n_ngb = adr_jp.size();
        for(size_t j=0; j<n_ngb; j++, adr_itr++){
            const Tptcl & ptclj = ptcl[*adr_itr];
            if(ptcl0.id != ptclj.id && ptcl1.id != ptclj.id){
                const PS::F64vec dr0 = ptcl0.pos - ptclj.pos;
                const PS::F64vec dr1 = ptcl1.pos - ptclj.pos;
                pot_old -= (ptcl0.mass * ptclj.mass) / sqrt(dr0*dr0);
                pot_old -= (ptcl1.mass * ptclj.mass) / sqrt(dr1*dr1);

		const PS::F64 y0 = sqrt(dr0*dr0) / r_out;
		const PS::F64 y1 = sqrt(dr1*dr1) / r_out;
                pot_hard_old -= (1.0-CalcW(y0, q))*(ptcl0.mass * ptclj.mass) / sqrt(dr0*dr0);
                pot_hard_old -= (1.0-CalcW(y1, q))*(ptcl1.mass * ptclj.mass) / sqrt(dr1*dr1);
            }
        }


        const PS::F64vec r0 = ptcl0.pos - pos_sun_;
        const PS::F64vec r1 = ptcl1.pos - pos_sun_;
        const PS::F64vec r01 = ptcl0.pos - ptcl1.pos;
	const PS::F64 pot0 = -mass_sun_ * m0 / sqrt(r0*r0);
	const PS::F64 pot1 = -mass_sun_ * m1 / sqrt(r1*r1);
	const PS::F64 pot01 = -m0 * m1 / sqrt(r01*r01);
        pot_old += pot0 + pot1 + pot01;
	const PS::F64 y = sqrt(r01*r01) / r_out;
        pot_hard_old += pot0 + pot1 + (1.0-CalcW(y,q))*pot01;
        const PS::F64vec v01 = ptcl0.vel - ptcl1.vel;
        const PS::F64 m_red = (m0*m1) / (m0+m1);
        PS::F64 kin = 0.5*m_red*v01*v01;

        const PS::F64 merge_time = ptcl0.time + time_origin_;
        merge_history.push_back( MergeLog(merge_time, ptcl_merge, ptcl_dead) );
        ptcl_merge.merge(ptcl_dead);
        ptcl_dead.mass = 0.0;

        const PS::F64vec r_new = ptcl_merge.pos - pos_sun_;
        PS::F64 pot_new = -mass_sun_ * ptcl_merge.mass / sqrt(r_new*r_new); 
        PS::F64 pot_hard_new = -mass_sun_ * ptcl_merge.mass / sqrt(r_new*r_new); 
        adr_itr = adr_jp.begin();

        for(size_t j=0; j<n_ngb; j++, adr_itr++){
            Tptcl & ptclj = ptcl[*adr_itr];
            if(ptcl0.id != ptclj.id && ptcl1.id != ptclj.id){
                const PS::F64vec dr = ptcl_merge.pos - ptclj.pos;
                pot_new -= (ptcl_merge.mass * ptclj.mass) / sqrt(dr*dr);
		const PS::F64 y = sqrt(dr*dr) / r_out;
                pot_hard_new -= (1.0-CalcW(y,q)) * (ptcl_merge.mass * ptclj.mass) / sqrt(dr*dr);
            }
        }

        eng_disp += kin - (pot_new - pot_old);
        eng_hard_disp += kin - (pot_hard_new - pot_hard_old);
    }

    template<class Tptcl>
    void merge2body(Tptcl & ptcl0, Tptcl & ptcl1, std::vector<MergeLog> & merge_history){
        const PS::F64 m0 = ptcl0.mass;
        const PS::F64 m1 = ptcl1.mass;
        const PS::S64 id0 = ptcl0.id;
        const PS::S64 id1 = ptcl1.id;
        Tptcl & ptcl_merge = (m0>m1) ? ptcl0 : ( (m1>m0) ? ptcl1 : ((id0<id1) ? ptcl0 : ptcl1) );
        Tptcl & ptcl_dead = (&ptcl_merge == &ptcl0) ? ptcl1 : ptcl0;
        const PS::F64 merge_time = ptcl0.time + time_origin_;
        merge_history.push_back( MergeLog(merge_time, ptcl_merge, ptcl_dead) );
        ptcl_merge.merge(ptcl_dead);
        ptcl_dead.mass = 0.0;
    }

    template<class Tptcl>
    void merge2body(Tptcl & ptcl0, Tptcl & ptcl1){
        const PS::F64 m0 = ptcl0.mass;
        const PS::F64 m1 = ptcl1.mass;
        const PS::S64 id0 = ptcl0.id;
        const PS::S64 id1 = ptcl1.id;
        Tptcl & ptcl_merge = (m0>m1) ? ptcl0 : ( (m1>m0) ? ptcl1 : ((id0<id1) ? ptcl0 : ptcl1) );
        Tptcl & ptcl_dead = (&ptcl_merge == &ptcl0) ? ptcl1 : ptcl0;
        ptcl_merge.merge(ptcl_dead);
        ptcl_dead.mass = 0.0;
    }

    void calcAcc0AndAcc1FromSun(const PS::F64vec & pos,
                                const PS::F64vec & vel,
                                PS::F64vec & acc0,
                                PS::F64vec & acc1){
        const PS::F64vec rij = pos - pos_sun_;
        const PS::F64vec vij = vel - vel_sun_;
        const PS::F64 r_inv = 1.0 / sqrt(rij * rij);
        const PS::F64 r2_inv = r_inv * r_inv;
        const PS::F64 r3_inv = r2_inv * r_inv;
        const PS::F64 m_r3 = mass_sun_ * r3_inv;
        const PS::F64vec F0 = -m_r3*rij;
        const PS::F64vec F1 = -m_r3*vij - 3.0*rij*vij*r2_inv*F0;
        acc0 += F0;
        acc1 += F1;
    }

public:
#ifdef CALC_HARD_ENERGY
    Energy eng_2body_prev_loc_;
    Energy eng_2body_now_loc_;
    Energy eng_mbody_prev_loc_;
    Energy eng_mbody_now_loc_;
#endif
    class Wtime{
    public:
        PS::F64 predict_;
        PS::F64 force_;
        PS::F64 force_kepler_;
        PS::F64 correct_;
        PS::F64 sort_;
        void clear(){
            predict_ = force_ = force_kepler_ = correct_ = sort_ = 0.0;
        }
        void dump(std::ostream & fout){
            fout<<"predict_= "<<predict_<<" force_= "<< force_<<" force_kepler_= "<<force_kepler_
                <<" correct_= "<<correct_<<" sort_= "<<sort_<<std::endl;
        }
    } wtime_profile_;

    void setSun(const PS::F64 m, const PS::F64vec & p, const PS::F64vec & v){
        mass_sun_ = m;
        pos_sun_ = p;
        vel_sun_ = v;
    }

    void clearEngDisp(){
        eng_disp_ = 0.0;
    }

    std::vector<PS::S32> adr_ptcl_loc_;
    std::vector< std::pair<PS::S64, PS::S64> > pair_id_loc_;
    std::vector< std::pair<PS::S32, PS::S32> > pair_adr_loc_org_; // M.I.
    std::vector< std::pair<PS::S32, PS::S32> > pair_rank_loc_; // M.I.
    std::vector<EPJSoft> epj_ngb_array_;
    std::vector<PTCLHard> ptcl_loc_;
    std::vector<PS::S32> rank_neighbor_;
    PS::S64 n_loop_;
    PS::F64 a0_offset_sq_;
    PS::F64 time_origin_;

    std::vector<PS::S32> rank_org_j_; // new
    template<class Ttree, class Tsystem>
    void setUp(Ttree & tree_soft,
               Tsystem & system_soft,
               const PS::S32 n_tot_loc, 
               const PS::F64 r_out,
               const PS::F64 r_in,
	       const PS::F64ort pos_domain[],
               const PS::F64 time_origin=0.0){
        time_origin_ = time_origin;
        const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;
        adr_ptcl_loc_.clear();
	epj_ngb_array_.clear();
	pair_id_loc_.clear();
	pair_adr_loc_org_.clear();
	pair_rank_loc_.clear();
        ptcl_loc_.clear();
	rank_org_j_.clear();
	ptcl_loc_org_.resize(n_tot_loc);
	for(PS::S32 i=0; i<n_tot_loc; i++){
	    PS::S32 n_ngb_init = system_soft[i].n_ngb;
            if( n_ngb_init > 0){
		EPJSoft * nbl = NULL;
		PS::S32 n_ngb = tree_soft.getNeighborListOneParticle(system_soft[i], nbl);
		/*
		if(n_ngb <= 0){
		    std::cout<<"rank, n_ngb, i, pos, n_ngb="<<PS::Comm::getRank()<<"  "<<n_ngb<<"  "<<i<<"  "<<system_soft[i].pos<<"  "<<system_soft[i].n_ngb<<std::endl;
		}
		*/
		assert(n_ngb>0);
		system_soft[i].n_ngb = n_ngb-1;
		if( n_ngb == 1) continue; // neighbor is itself
		for(PS::S32 j=0; j<n_ngb; j++){
		    if( system_soft[i].id == (nbl+j)->id) continue;
		    const PS::F64vec pos_i = system_soft[i].pos;
		    PS::F64 & pot_tot_i = system_soft[i].pot_tot;
		    PS::F64vec & acc_pla_i = system_soft[i].acc;
		    CalcAccPotShortWithLinearCutoff
			(pos_i,            acc_pla_i,		 pot_tot_i, 
			 (nbl+j)->pos,     (nbl+j)->mass,
			 eps_sq,           r_out, 
			 r_in);
		    const PS::F64vec dr = system_soft[i].pos - (nbl+j)->pos;
		    const PS::F64 reps_sq = dr*dr + eps_sq;
		    if(reps_sq > EPJSoft::r_search*EPJSoft::r_search){
			system_soft[i].n_ngb--;
			continue;
		    }
		    epj_ngb_array_.push_back(nbl[j]);
		    pair_id_loc_.push_back (std::pair<PS::S64, PS::S64>    ( (PS::S64)system_soft[i].id,  (PS::S64)(nbl+j)->id) );
		    pair_adr_loc_org_.push_back(std::pair<PS::S32, PS::S32>( (PS::S32)system_soft[i].adr, (PS::S32)(nbl+j)->adr_org));
		    pair_rank_loc_.push_back (std::pair<PS::S64, PS::S64>  ( (PS::S64)system_soft[i].rank_org,  (PS::S64)(nbl+j)->rank_org) );
		    rank_org_j_.push_back((nbl+j)->rank_org); // new
		}
		if(system_soft[i].n_ngb > 0){
		    adr_ptcl_loc_.push_back(i);
		    ptcl_loc_.push_back( PTCLHard(system_soft[i].id,    system_soft[i].mass,
						  system_soft[i].pos,   system_soft[i].vel,
						  system_soft[i].n_ngb, system_soft[i].rank_org,
						  system_soft[i].adr) );
		}
	    }
	} // end of for over i

	/*
	for(PS::S32 i=0; i<n_tot_loc; i++){
	    PS::S32 n_ngb_init = system_soft[i].n_ngb;
            if( n_ngb_init > 0){
		EPJSoft * nbl = NULL;
		PS::S32 n_ngb = tree_soft.getNeighborListOneParticle(system_soft[i], nbl);
		assert(n_ngb>0);
		system_soft[i].n_ngb = n_ngb-1;
		if( n_ngb == 1) continue; // neighbor is itself
		for(PS::S32 j=0; j<n_ngb; j++){
		    const PS::F64vec dr = system_soft[i].pos - (nbl+j)->pos;
		    const PS::F64 reps_sq = dr*dr + eps_sq;
		    if( system_soft[i].id == (nbl+j)->id) continue;
		    if(reps_sq > EPJSoft::r_search*EPJSoft::r_search){
			system_soft[i].n_ngb--;
			continue;
		    }
		    epj_ngb_array_.push_back(nbl[j]);
		    pair_id_loc_.push_back (std::pair<PS::S64, PS::S64>    ( (PS::S64)system_soft[i].id,  (PS::S64)(nbl+j)->id) );
		    pair_adr_loc_org_.push_back(std::pair<PS::S32, PS::S32>( (PS::S32)system_soft[i].adr, (PS::S32)(nbl+j)->adr_org));
		    pair_rank_loc_.push_back (std::pair<PS::S64, PS::S64>  ( (PS::S64)system_soft[i].rank_org,  (PS::S64)(nbl+j)->rank_org) );
		    rank_org_j_.push_back((nbl+j)->rank_org); // new
		    const PS::F64vec pos_i = system_soft[i].pos;
		    PS::F64 & pot_tot_i = system_soft[i].pot_tot;
#if 1
		    PS::F64vec & acc_tot_i = system_soft[i].acc;
		    CalcAccPotShortWithLinearCutoff
			(pos_i,            acc_tot_i,		 pot_tot_i, 
			 (nbl+j)->pos,     (nbl+j)->mass,
			 eps_sq,           r_out, 
			 r_in);
#else
		    CalcPotShort(pos_i,            pot_tot_i, 
				 (nbl+j)->pos,     (nbl+j)->mass,
				 eps_sq,           r_out, 
				 r_in);
#endif
		}
		if(system_soft[i].n_ngb > 0){
		    adr_ptcl_loc_.push_back(i);
		    ptcl_loc_.push_back( PTCLHard(system_soft[i].id,    system_soft[i].mass,
						  system_soft[i].pos,   system_soft[i].vel,
						  system_soft[i].n_ngb, system_soft[i].rank_org,
						  system_soft[i].adr) );
		}
	    }
	} // end of for over i
	*/

	/// set a0 offset
	PS::F64 m_min_loc = PS::LARGE_FLOAT;
	const PS::S32 n = ptcl_loc_.size();
	for(PS::S32 i=0; i<n; i++){
	    if(m_min_loc > ptcl_loc_[i].mass && ptcl_loc_[i].mass > 0.0) m_min_loc = ptcl_loc_[i].mass;
	}
	PS::F64 a0_offset_sq_loc = 0.1 * m_min_loc / (r_out*r_out);
	a0_offset_sq_ = PS::Comm::getMinValue(a0_offset_sq_loc);
	////////////
	// add M.I.
	const PS::S32 n_proc = PS::Comm::getNumberOfProc();
	const PS::S32 my_rank = PS::Comm::getRank();
	const PS::F64 r_crit_sq = (EPJSoft::r_search * EPJSoft::r_search + eps_sq)*SAFTY_FACTOR_FOR_SEARCH_SQ*1.01;
	const PS::F64ort pos_my_domain = pos_domain[my_rank];
	rank_neighbor_.clear();
	for(int i=0; i<n_proc; i++){
	    if(i==my_rank) continue;
	    else if( r_crit_sq >= pos_my_domain.getDistanceMinSQ(pos_domain[i]) ){
		rank_neighbor_.push_back(i);
		//fout_debug<<"rank_neighbor_.back()"<<rank_neighbor_.back()<<std::endl;
	    }
	}
    }

    std::vector<PTCLHardComm> ptcl_loc_org_;
    template<class Tsystem>
    void copyVelSoftToHard(const Tsystem & system){
	const PS::S32 n = adr_ptcl_loc_.size();
#pragma omp parallel for
	for(PS::S32 i=0; i<n; i++){
	    const PS::S32 adr = adr_ptcl_loc_[i];
	    ptcl_loc_[i].vel = system[adr].vel;
	}
	const PS::S32 n_loc_tot = system.getNumberOfParticleLocal();
#pragma omp parallel for 
	for(PS::S32 i=0; i<n_loc_tot; i++){
	    ptcl_loc_org_[i].copyFromFP(system[i]);
	}
    }
    
#if 1
    std::vector< std::pair<PS::S64, PS::S64> > pair_adr_loc_;
    void makeDictionaryLocal(){
	//makeDictionaryLocalImpl(ptcl_loc_, pair_id_loc_, pair_adr_loc_);
	makeDictionaryLocalImpl(ptcl_loc_, pair_id_loc_, pair_rank_loc_, pair_adr_loc_);
    }
#endif


    std::vector<PS::S32> n_ptcl_array_;
    std::vector< std::pair<PS::S32, PS::S32> > n_int_n_ptcl_array_;
    std::vector<PS::S32> n_ptcl_disp_;
    std::vector<PTCLHard> ptcl_multi_glb_;
    std::vector< std::pair<PS::S64, PS::S64> > pair_id_multi_glb_;
    void gatherData(){
	static std::vector<PS::S32> n_interaction_array;
	static std::vector<PS::S32> n_interaction_disp;
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        n_int_n_ptcl_array_.reserve(n_proc);
        std::pair<PS::S32, PS::S32> n_int_n_ptcl(pair_id_multi_loc_.size(), ptcl_multi_loc_.size());
        n_interaction_array.reserve(n_proc);
        n_ptcl_array_.reserve(n_proc);
        n_interaction_disp.reserve(n_proc+1);
        n_ptcl_disp_.reserve(n_proc+1);
#if 1
	PS::Comm::gather(&n_int_n_ptcl, 1, &n_int_n_ptcl_array_[0]);
        n_interaction_disp[0] = n_ptcl_disp_[0] = 0;
	if(PS::Comm::getRank() == 0){
	    for(PS::S32 i=0; i<n_proc; i++){
		n_interaction_array[i] = n_int_n_ptcl_array_[i].first;
		n_ptcl_array_[i] = n_int_n_ptcl_array_[i].second;
		n_interaction_disp[i+1] = n_interaction_disp[i] + n_interaction_array[i];
		n_ptcl_disp_[i+1] = n_ptcl_disp_[i] + n_ptcl_array_[i];
	    }
	    pair_id_multi_glb_.resize(n_interaction_disp[n_proc]);
	    ptcl_multi_glb_.resize(n_ptcl_disp_[n_proc]);
	}
#else
        PS::Comm::allGather(&n_int_n_ptcl, 1, &n_int_n_ptcl_array_[0]);
        n_interaction_disp[0] = n_ptcl_disp[0] = 0;
	for(PS::S32 i=0; i<n_proc; i++){
	    n_interaction_array[i] = n_int_n_ptcl_array_[i].first;
	    n_ptcl_array_[i] = n_int_n_ptcl_array_[i].second;
	    n_interaction_disp[i+1] = n_interaction_disp[i] + n_interaction_array[i];
	    n_ptcl_disp_[i+1] = n_ptcl_disp_[i] + n_ptcl_array_[i];
	}
        pair_id_multi_glb_.resize(n_interaction_disp[n_proc]);
        ptcl_multi_glb_.resize(n_ptcl_disp_[n_proc]);
#endif

        PS::Comm::gatherV(&pair_id_multi_loc_[0], pair_id_multi_loc_.size(), &pair_id_multi_glb_[0],
                          &n_interaction_array[0], &n_interaction_disp[0]);
        PS::Comm::gatherV(&ptcl_multi_loc_[0], ptcl_multi_loc_.size(),  &ptcl_multi_glb_[0],
                          &n_ptcl_array_[0], &n_ptcl_disp_[0]);
    }

    //std::unordered_map<PS::S64, PS::S64> idx_to_adr_multi_glb_;
    std::vector< std::pair<PS::S64, PS::S64> > pair_adr_multi_glb_;
    void makeDictionaryGlobal(){
        //makeDictionary(ptcl_multi_glb_, pair_id_multi_glb_, idx_to_adr_multi_glb_, pair_adr_multi_glb_);
	makeDictionaryGlobalImpl(ptcl_multi_glb_, pair_id_multi_glb_, pair_adr_multi_glb_);
    }

    std::vector<PTCLPred> ptcl_pred_;
    std::vector<PTCLForce> ptcl_force_;
    std::vector<PS::F64> time_next_;
    std::vector<PS::S32> adr_sorted_;
    PS::F64 time_sys_;
    PS::F64 time_end_;
    PS::F64 dt_limit_;
    PS::F64 time_sync_;
    PS::S32 n_active_;

    std::vector<bool> merge_flag_glb_;
    bool merge_state_;
    void setUpGlb(const PS::F64 time_end, 
                  const PS::F64 dt_limit){
        const PS::S32 n = ptcl_multi_glb_.size();
        n_active_ = n;
        time_sys_ = 0.0;
        time_end_ = time_end;
        dt_limit_ = dt_limit;
        time_sync_ = time_sys_ + dt_limit_;
        adr_sorted_.clear();
        adr_sorted_.resize(n);
        ptcl_pred_.resize(n);
        ptcl_force_.resize(n);
        time_next_.resize(n);
        merge_flag_glb_.clear();
        merge_flag_glb_.resize(n);
        for(PS::S32 ip=0; ip<n; ip++){
            adr_sorted_[ip] = ip;
            ptcl_pred_[ip].pos = ptcl_multi_glb_[ip].pos;
            ptcl_pred_[ip].vel = ptcl_multi_glb_[ip].vel;
            ptcl_force_[ip].acc0 = ptcl_force_[ip].acc1 = 0.0;
            ptcl_multi_glb_[ip].time = time_next_[ip] = 0.0;
            ptcl_multi_glb_[ip].setRMerge();
            merge_flag_glb_[ip] = false;
        }
        merge_state_ = false;
        n_loop_ = 0;
        wtime_profile_.clear();
    }

    std::vector< std::pair<PS::S64, PS::S64> > pair_merge_adr_;
    void calcAcc0AndAcc1(const PS::F64 r_out,
                         const PS::F64 r_in){
        const PS::S32 n_thread_max = PS::Comm::getNumberOfThread();
        static std::vector< std::pair<PS::S64, PS::S64> > * pair_merge_adr_omp;
        static bool first = true;
        if(first){
            pair_merge_adr_omp = new std::vector< std::pair<PS::S64, PS::S64> >[n_thread_max];
            first = false;
        }
        for(PS::S32 i=0; i<n_thread_max; i++){
            pair_merge_adr_omp[i].clear();
        }
        pair_merge_adr_.clear();
        const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;
        const PS::S32 ni = n_active_;
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            //const PS::S32 ith = PS::Comm::getThreadNum();
            const PS::S32 adr_i = adr_sorted_[ip];
            ptcl_force_[adr_i].clear();
            const PS::S32 nj = ptcl_multi_glb_[adr_i].n_ngb;
            const PS::S32 adr_pair = ptcl_multi_glb_[adr_i].adr_pair;
            for(PS::S32 jp=0; jp<nj; jp++){
                //PS::F64 r2 = 0.0;
                const PS::S32 adr_j = pair_adr_multi_glb_[adr_pair+jp].second;
                if(   ptcl_multi_glb_[adr_i].id != pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair].first
                      || ptcl_multi_glb_[adr_j].id != pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second){
		    std::cerr<<"normal"<<std::endl;
		    std::cerr<<"ip="<<ip<<" jp="<<jp<<" nj="<<nj<<" adr_pair="<<adr_pair<<std::endl;
                    std::cerr<<"adr_i="<<adr_i<<" adr_j="<<adr_j<<std::endl;
                    std::cerr<<"ptcl_multi_glb_[adr_i].id="<<ptcl_multi_glb_[adr_i].id<<" ptcl_multi_glb_[adr_j].id="<<ptcl_multi_glb_[adr_j].id<<std::endl;
                    std::cerr<<"pair_multi_id_glb_[ptcl_multi_glb_[adr_i].adr_pair].first="<<pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair].first
                             <<" pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second="<<pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second<<std::endl;
                    std::cerr<<std::endl;
                }
                assert(ptcl_multi_glb_[adr_i].id == pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair].first);
                assert(ptcl_multi_glb_[adr_j].id == pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second);
#ifdef FORDEBUG
                CalcAcc0AndAcc1Cutoff(ptcl_pred_[adr_i].pos,   ptcl_pred_[adr_i].vel,
                                      ptcl_force_[adr_i].acc0, ptcl_force_[adr_i].acc1,
                                      ptcl_force_[adr_i].pot,  ptcl_pred_[adr_j].pos,
                                      ptcl_pred_[adr_j].vel,   ptcl_multi_glb_[adr_j].mass,
                                      eps_sq,                  r_out, r_in);		
#else // FORDEBUG
#ifdef MERGE
		PS::F64 r2 = 0.0;
                CalcAcc0Acc1AndR2Cutoff(ptcl_pred_[adr_i].pos,   ptcl_pred_[adr_i].vel,
                                        ptcl_force_[adr_i].acc0_pla, ptcl_force_[adr_i].acc1_pla,
                                        r2, 
                                        ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                        ptcl_multi_glb_[adr_j].mass,   eps_sq,
                                        r_out, r_in);
                const PS::F64 r_merge = (ptcl_multi_glb_[adr_i].r_merge + ptcl_multi_glb_[adr_j].r_merge);
                const PS::F64 r_merge_2 = r_merge * r_merge;
                if(r2 < r_merge_2){
		    const PS::S32 ith = PS::Comm::getThreadNum();
                    //pair_merge_adr_.push_back( std::pair<PS::S64, PS::S64>(adr_i, adr_j) );
                    pair_merge_adr_omp[ith].push_back( std::pair<PS::S64, PS::S64>(adr_i, adr_j) );
                }
#else // MERGE
                CalcAcc0AndAcc1Cutoff(ptcl_pred_[adr_i].pos,   ptcl_pred_[adr_i].vel,
                                      ptcl_force_[adr_i].acc0_pla, ptcl_force_[adr_i].acc1_pla,
                                      ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                      ptcl_multi_glb_[adr_j].mass,   eps_sq,
                                      r_out, r_in);
#endif // MERGE
#endif // FORDEBUG
            }
        }
#ifdef MERGE
        for(PS::S32 i=0; i<n_thread_max; i++){
            const size_t j_max = pair_merge_adr_omp[i].size();
            for(size_t j=0; j<j_max; j++){
                pair_merge_adr_.push_back(pair_merge_adr_omp[i][j]);
            }
        }
        const PS::S32 list_len = pair_merge_adr_.size();
        bool sort_again = false;
        for(PS::S32 i=0; i<list_len; i++){
            const PS::S32 adr_i = pair_merge_adr_[i].first;
            const PS::S32 adr_j = pair_merge_adr_[i].second;
            if( time_next_[adr_i] != time_next_[adr_j] && !ptcl_multi_glb_[adr_j].isDead()){
		// if the merging two particles are not integraed at the same time,
		// the particle non-integrated particle is forced to be integraed.
#if 0
                std::cout<<"adr_i="<<adr_i
                         <<" adr_j="<<adr_j<<std::endl;
                std::cout<<"time_next_[adr_i]="<<time_next_[adr_i]
                         <<" time_next_[adr_j]="<<time_next_[adr_j]<<std::endl;
                std::cout<<"ptcl_multi_glb_[adr_i].n_ngb="<<ptcl_multi_glb_[adr_i].n_ngb<<std::endl;
#endif
                time_next_[adr_j] = time_next_[adr_i];
                ptcl_multi_glb_[adr_j].dt = time_next_[adr_j] - ptcl_multi_glb_[adr_j].time;
                ptcl_force_[adr_j].clear();
                const PS::S32 nk = ptcl_multi_glb_[adr_j].n_ngb;
                const PS::S32 adr_pair = ptcl_multi_glb_[adr_j].adr_pair;
                for(PS::S32 kp=0; kp<nk; kp++){
                    // Here, not care about second (j-th and k-th) merger.
                    const PS::S32 adr_k = pair_adr_multi_glb_[adr_pair+kp].second;
                    CalcAcc0AndAcc1Cutoff(ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                          ptcl_force_[adr_j].acc0_pla, ptcl_force_[adr_j].acc1_pla,
                                          ptcl_pred_[adr_k].pos,   ptcl_pred_[adr_k].vel,
                                          ptcl_multi_glb_[adr_k].mass,   eps_sq,
                                          r_out, r_in);
                }
                n_active_++;
                sort_again = true;
            }
        }
        if(sort_again){
            const PS::S32 n = ptcl_multi_glb_.size();
            SortAdr::time = &time_next_[0];
            std::sort(&adr_sorted_[0], &adr_sorted_[0]+n, SortAdr());
        }
#endif
    }




    void calcAcc0AndAcc1_tmp(const PS::F64 r_out,
                         const PS::F64 r_in){
        const PS::S32 n_thread_max = PS::Comm::getNumberOfThread();
        static std::vector< std::pair<PS::S64, PS::S64> > * pair_merge_adr_omp;
        static bool first = true;
        if(first){
            pair_merge_adr_omp = new std::vector< std::pair<PS::S64, PS::S64> >[n_thread_max];
            first = false;
        }
        for(PS::S32 i=0; i<n_thread_max; i++){
            pair_merge_adr_omp[i].clear();
        }
        pair_merge_adr_.clear();
        const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;
        const PS::S32 ni = n_active_;
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            //const PS::S32 ith = PS::Comm::getThreadNum();
            const PS::S32 adr_i = adr_sorted_[ip];
            ptcl_force_[adr_i].clear();
            const PS::S32 nj = ptcl_multi_glb_[adr_i].n_ngb;
            const PS::S32 adr_pair = ptcl_multi_glb_[adr_i].adr_pair;
            for(PS::S32 jp=0; jp<nj; jp++){
                //PS::F64 r2 = 0.0;
                const PS::S32 adr_j = pair_adr_multi_glb_[adr_pair+jp].second;
                if( ptcl_multi_glb_[adr_i].id != pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair].first
		    || ptcl_multi_glb_[adr_j].id != pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second ){
		    std::cerr<<"tmp"<<std::endl;
		    std::cerr<<"ip="<<ip<<" jp="<<jp<<" nj="<<nj<<" adr_pair="<<adr_pair<<std::endl;
                    std::cerr<<"adr_i="<<adr_i<<" adr_j="<<adr_j<<std::endl;
                    std::cerr<<"ptcl_multi_glb_[adr_i].id="<<ptcl_multi_glb_[adr_i].id<<" ptcl_multi_glb_[adr_j].id="<<ptcl_multi_glb_[adr_j].id<<std::endl;
                    std::cerr<<"pair_multi_id_glb_[ptcl_multi_glb_[adr_i].adr_pair].first="<<pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair].first
                             <<" pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second="<<pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second<<std::endl;
                    std::cerr<<std::endl;
                }
                assert(ptcl_multi_glb_[adr_i].id == pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair].first);
                assert(ptcl_multi_glb_[adr_j].id == pair_id_multi_glb_[ptcl_multi_glb_[adr_i].adr_pair+jp].second);
#ifdef FORDEBUG
                CalcAcc0AndAcc1Cutoff(ptcl_pred_[adr_i].pos,   ptcl_pred_[adr_i].vel,
                                      ptcl_force_[adr_i].acc0, ptcl_force_[adr_i].acc1,
                                      ptcl_force_[adr_i].pot,  ptcl_pred_[adr_j].pos,
                                      ptcl_pred_[adr_j].vel,   ptcl_multi_glb_[adr_j].mass,
                                      eps_sq,                  r_out, r_in);		
#else // FORDEBUG
#ifdef MERGE
		PS::F64 r2 = 0.0;
                CalcAcc0Acc1AndR2Cutoff(ptcl_pred_[adr_i].pos,   ptcl_pred_[adr_i].vel,
                                        ptcl_force_[adr_i].acc0_pla, ptcl_force_[adr_i].acc1_pla,
                                        r2, 
                                        ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                        ptcl_multi_glb_[adr_j].mass,   eps_sq,
                                        r_out, r_in);
                const PS::F64 r_merge = (ptcl_multi_glb_[adr_i].r_merge + ptcl_multi_glb_[adr_j].r_merge);
                const PS::F64 r_merge_2 = r_merge * r_merge;
                if(r2 < r_merge_2){
		    const PS::S32 ith = PS::Comm::getThreadNum();
                    //pair_merge_adr_.push_back( std::pair<PS::S64, PS::S64>(adr_i, adr_j) );
                    pair_merge_adr_omp[ith].push_back( std::pair<PS::S64, PS::S64>(adr_i, adr_j) );
                }
#else // MERGE
                CalcAcc0AndAcc1Cutoff(ptcl_pred_[adr_i].pos,   ptcl_pred_[adr_i].vel,
                                      ptcl_force_[adr_i].acc0_pla, ptcl_force_[adr_i].acc1_pla,
                                      ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                      ptcl_multi_glb_[adr_j].mass,   eps_sq,
                                      r_out, r_in);
#endif // MERGE
#endif // FORDEBUG
            }
        }
#ifdef MERGE
        for(PS::S32 i=0; i<n_thread_max; i++){
            const size_t j_max = pair_merge_adr_omp[i].size();
            for(size_t j=0; j<j_max; j++){
                pair_merge_adr_.push_back(pair_merge_adr_omp[i][j]);
            }
        }
        const PS::S32 list_len = pair_merge_adr_.size();
        bool sort_again = false;
        for(PS::S32 i=0; i<list_len; i++){
            const PS::S32 adr_i = pair_merge_adr_[i].first;
            const PS::S32 adr_j = pair_merge_adr_[i].second;
            if( time_next_[adr_i] != time_next_[adr_j] && !ptcl_multi_glb_[adr_j].isDead()){
		// if the merging two particles are not integraed at the same time,
		// the particle non-integrated particle is forced to be integraed.
#if 0
                std::cout<<"adr_i="<<adr_i
                         <<" adr_j="<<adr_j<<std::endl;
                std::cout<<"time_next_[adr_i]="<<time_next_[adr_i]
                         <<" time_next_[adr_j]="<<time_next_[adr_j]<<std::endl;
                std::cout<<"ptcl_multi_glb_[adr_i].n_ngb="<<ptcl_multi_glb_[adr_i].n_ngb<<std::endl;
#endif
                time_next_[adr_j] = time_next_[adr_i];
                ptcl_multi_glb_[adr_j].dt = time_next_[adr_j] - ptcl_multi_glb_[adr_j].time;
                ptcl_force_[adr_j].clear();
                const PS::S32 nk = ptcl_multi_glb_[adr_j].n_ngb;
                const PS::S32 adr_pair = ptcl_multi_glb_[adr_j].adr_pair;
                for(PS::S32 kp=0; kp<nk; kp++){
                    // Here, not care about second (j-th and k-th) merger.
                    const PS::S32 adr_k = pair_adr_multi_glb_[adr_pair+kp].second;
                    CalcAcc0AndAcc1Cutoff(ptcl_pred_[adr_j].pos,   ptcl_pred_[adr_j].vel,
                                          ptcl_force_[adr_j].acc0_pla, ptcl_force_[adr_j].acc1_pla,
                                          ptcl_pred_[adr_k].pos,   ptcl_pred_[adr_k].vel,
                                          ptcl_multi_glb_[adr_k].mass,   eps_sq,
                                          r_out, r_in);
                }
                n_active_++;
                sort_again = true;
            }
        }
        if(sort_again){
            const PS::S32 n = ptcl_multi_glb_.size();
            SortAdr::time = &time_next_[0];
            std::sort(&adr_sorted_[0], &adr_sorted_[0]+n, SortAdr());
        }
#endif
    }


    void calcAcc0AndAcc1Kepler(const PS::F64 mass_sun,
                               const PS::F64vec & pos_sun,
                               const PS::F64vec & vel_sun){
        const PS::S32 ni = n_active_;
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            PTCLPred & pred = ptcl_pred_[adr];
            PTCLForce & force = ptcl_force_[adr];
            const PS::F64vec rij = pred.pos - pos_sun;
            const PS::F64vec vij = pred.vel - vel_sun;
            const PS::F64 r_inv = 1.0 / sqrt(rij * rij);
            const PS::F64 r2_inv = r_inv * r_inv;
            const PS::F64 r3_inv = r2_inv * r_inv;
            const PS::F64 m_r3 = mass_sun * r3_inv;
            const PS::F64vec F0 = -m_r3*rij;
            const PS::F64vec F1 = -m_r3*vij - 3.0*rij*vij*r2_inv*F0;
            force.acc0 += F0;
            force.acc1 += F1;
        }
    }
    
    
    void setInitDt(const PS::F64 eta){
        const PS::S32 ni = n_active_;
        const PS::F64 dt_limit_new = CalcDtLimit(time_sys_, time_sync_, dt_limit_);
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            ptcl_multi_glb_[adr].setDt2nd(ptcl_force_[adr], eta, dt_limit_new, a0_offset_sq_);
        }
    }

    void copyForceToPtcl(){
        const PS::S32 ni = n_active_;
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            ptcl_multi_glb_[adr].acc0 = ptcl_force_[adr].acc0;
            ptcl_multi_glb_[adr].acc1 = ptcl_force_[adr].acc1;
            ptcl_multi_glb_[adr].acc0_pla = ptcl_force_[adr].acc0_pla;
            ptcl_multi_glb_[adr].acc1_pla = ptcl_force_[adr].acc1_pla;
#ifdef FORDEBUG
            ptcl_multi_glb_[adr].pot = ptcl_force_[adr].pot;
#endif
        }
    }

    void sortAndSelectIp(){
        const PS::S32 ni_old = n_active_;
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni_old; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            time_next_[adr] += ptcl_multi_glb_[adr].dt;
        }
        SortAdr::time = &time_next_[0];
        if(merge_state_){
            PS::S32 n_tot = ptcl_multi_glb_.size();
            std::sort(&adr_sorted_[0], &adr_sorted_[0]+n_tot, SortAdr());
        }
        else{
            std::sort(&adr_sorted_[0], &adr_sorted_[0]+ni_old, SortAdr());
        }
        const PS::F64 t_n = time_next_[adr_sorted_[0]];
        const PS::S32 n = ptcl_multi_glb_.size();
        for(n_active_=1; n_active_<n; n_active_++){
            if(t_n < time_next_[adr_sorted_[n_active_]]) {
                break;
            }
        }
    }

#if 1
    void predictAll(){
        static const PS::F64 inv3 = 1.0 / 3.0;
        const PS::S32 n = ptcl_multi_glb_.size();
        const PS::F64 time_next = time_next_[adr_sorted_[0]];
#pragma omp parallel for
        for(PS::S32 i=0; i<n/2; i++){
            const PTCLHard & p_0 = ptcl_multi_glb_[i*2];
            const PTCLHard & p_1 = ptcl_multi_glb_[i*2+1];
            const PS::F64 dt_0 = time_next - p_0.time;
            const PS::F64 dt_1 = time_next - p_1.time;
            ptcl_pred_[i*2].pos = p_0.pos + dt_0*(p_0.vel  + 0.5*dt_0*(p_0.acc0 + inv3*dt_0*p_0.acc1));
            ptcl_pred_[i*2].vel = p_0.vel + dt_0*(p_0.acc0 + 0.5*dt_0*p_0.acc1);
            ptcl_pred_[i*2+1].pos = p_1.pos + dt_1*(p_1.vel  + 0.5*dt_1*(p_1.acc0 + inv3*dt_1*p_1.acc1));
            ptcl_pred_[i*2+1].vel = p_1.vel + dt_1*(p_1.acc0 + 0.5*dt_1*p_1.acc1);
#if 0
            const PS::F64 c1_0 = 0.5 * dt_0 * dt_0;
            const PS::F64 c2_0 = c1_0 * inv3 * dt_0;
            const PS::F64 c1_1 = 0.5 * dt_1 * dt_1;
            const PS::F64 c2_1 = c1_1 * inv3 * dt_1;
            ptcl_pred_[i*2].pos = p_0.pos + dt_0*p_0.vel  + c1_0*p_0.acc0 + c2_0*p_0.acc1;
            ptcl_pred_[i*2].vel = p_0.vel + dt_0*p_0.acc0 + c1_0*p_0.acc1;
            ptcl_pred_[i*2+1].pos = p_1.pos + dt_1*p_1.vel  + c1_1*p_1.acc0 + c2_1*p_1.acc1;
            ptcl_pred_[i*2+1].vel = p_1.vel + dt_1*p_1.acc0 + c1_1*p_1.acc1;
#endif
        }
        if(n%2 != 0){
            const PTCLHard & p_0 = ptcl_multi_glb_[n-1];
            const PS::F64 dt_0 = time_next - p_0.time;
            ptcl_pred_[n-1].pos = p_0.pos + dt_0*(p_0.vel  + 0.5*dt_0*(p_0.acc0 + inv3*dt_0*p_0.acc1));
            ptcl_pred_[n-1].vel = p_0.vel + dt_0*(p_0.acc0 + 0.5*dt_0*p_0.acc1);
        }
        /*
        //#pragma omp parallel for
        for(PS::S32 i=0; i<n; i++){
        const PTCLHard & p = ptcl_multi_glb_[i];
        const PS::F64 dt = time_next - p.time;
	    ptcl_pred_[i].pos = p.pos + dt*(p.vel  + 0.5*dt*(p.acc0 + inv3*dt*p.acc1));
        ptcl_pred_[i].vel = p.vel + dt*(p.acc0 + 0.5*dt*p.acc1);
        }
        */
    }
#else
#ifdef INTRINSIC_X86
    void predictAll(){
        asm("###predictAll");
        typedef double v2df __attribute__((vector_size(16)));
        typedef double v4df __attribute__((vector_size(32)));
        const PS::F64 tn = time_next_[adr_sorted_[0]];
        const PS::S32 n = ptcl_multi_glb_.size();
        for(PS::S32 i=0; i<n; i++){
            const PTCLHard & p_0 = ptcl_multi_glb_[i];
            const v4df pos{p_0.pos.x, p_0.pos.y, p_0.pos.z, 0.0};
            const v4df vel{p_0.vel.x, p_0.vel.y, p_0.vel.z, 0.0};
            const v4df acc0{p_0.acc0.x, p_0.acc0.y, p_0.acc0.z, 0.0};
            const v4df acc1{p_0.acc1.x, p_0.acc1.y, p_0.acc1.z, 0.0};
            const PS::F64 dt = tn - p_0.time;
            const v4df c0{dt, dt, dt, 0.0};
            const v4df c1 = c0 * (v4df){0.5, 0.5, 0.5, 0.0};
            const v4df c2 = c0 * (v4df){1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0};
            const v4df pos_pre = pos + c0*(vel  + c1*(acc0 + c2*acc1));
            const v4df vel_pre = vel + c0*(acc0 + c1*acc1);
            PS::F64vec & pos_pred = ptcl_pred_[i].pos;
            PS::F64vec & vel_pred = ptcl_pred_[i].vel;
            const v2df pos_pre_l = __builtin_ia32_vextractf128_pd256(pos_pre, 0);
            const v2df pos_pre_h = __builtin_ia32_vextractf128_pd256(pos_pre, 1);
            const v2df vel_pre_l = __builtin_ia32_vextractf128_pd256(vel_pre, 0);
            const v2df vel_pre_h = __builtin_ia32_vextractf128_pd256(vel_pre, 1);
            pos_pred.x = __builtin_ia32_vec_ext_v2df( pos_pre_l, 0 );
            pos_pred.y = __builtin_ia32_vec_ext_v2df( pos_pre_l, 1 );
            pos_pred.z = __builtin_ia32_vec_ext_v2df( pos_pre_h, 0 );
            vel_pred.x = __builtin_ia32_vec_ext_v2df( vel_pre_l, 0 );
            vel_pred.y = __builtin_ia32_vec_ext_v2df( vel_pre_l, 1 );
            vel_pred.z = __builtin_ia32_vec_ext_v2df( vel_pre_h, 0 );
        }
    }
#endif

#endif

    std::vector<MergeLog> merge_history_;
    PS::F64 eng_disp_;
    PS::F64 eng_hard_multi_disp_;
    HardSystem(){
        merge_history_.clear();
        merge_history_.reserve(1000);
        eng_disp_ = 0.0;
        eng_hard_multi_disp_ = 0.0;
	/*
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
	n_2body_cand_send_ = new PS::S32[n_proc];
	n_2body_cand_send_disp_ = new PS::S32[n_proc+1];
	n_2body_cand_recv_ = new PS::S32[n_proc];
	n_2body_cand_recv_disp_ = new PS::S32[n_proc+1];
	*/
    }
    
    void correctIp(const PS::F64 eta,
                   const PS::F64 dt_limit_org,
                   const PS::F64 r_out,
                   const PS::F64 r_in){
        const PS::S32 ni = n_active_;
        time_sys_ = time_next_[adr_sorted_[0]];
        if(time_sys_ == time_sync_) time_sync_ += dt_limit_org;
        const PS::F64 dt_limit_new = CalcDtLimit(time_sys_, time_sync_, dt_limit_org);
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            ptcl_multi_glb_[adr].correct(ptcl_force_[adr], eta, dt_limit_new, a0_offset_sq_);
        }
#ifdef MERGE
        merge_state_ = false;
        const PS::S32 n_int = pair_merge_adr_.size();
        for(PS::S32 i=0; i<n_int; i++){
            const PS::S32 adr_i = pair_merge_adr_[i].first;
            const PS::S32 adr_j = pair_merge_adr_[i].second;
            if(merge_flag_glb_[adr_i] == true || merge_flag_glb_[adr_j] == true) continue;
#ifdef CALC_HARD_ENERGY
            merge2bodyWithHardEnergy(ptcl_multi_glb_, pair_adr_multi_glb_, adr_i, adr_j, eng_disp_, merge_history_, r_in, r_out, eng_hard_multi_disp_);
#else
            merge2body(ptcl_multi_glb_, pair_adr_multi_glb_, adr_i, adr_j, eng_disp_, merge_history_);
#endif
            PS::S64 adr_dead = (ptcl_multi_glb_[adr_i].mass == 0.0) ? adr_i : adr_j;
            merge_flag_glb_[adr_dead] = true;
            merge_state_ = true;
        }

        if(merge_state_){
            const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;
            const PS::S32 ni = n_active_;
            for(PS::S32 ip=0; ip<ni; ip++){
                const PS::S32 adr_i = adr_sorted_[ip];
                if(ptcl_multi_glb_[adr_i].isDead()){
                    ptcl_multi_glb_[adr_i].dt = PS::LARGE_FLOAT;
                }
                else{
                    ptcl_force_[adr_i].clear();
                    const PS::S32 nj = ptcl_multi_glb_[adr_i].n_ngb;
                    const PS::S32 adr_pair = ptcl_multi_glb_[adr_i].adr_pair;
                    for(PS::S32 jp=0; jp<nj; jp++){
                        const PS::S32 adr_j = pair_adr_multi_glb_[adr_pair+jp].second;
                        if(ptcl_multi_glb_[adr_j].isDead()) continue;
                        CalcAcc0AndAcc1Cutoff(ptcl_multi_glb_[adr_i].pos,    ptcl_multi_glb_[adr_i].vel,
                                              ptcl_force_[adr_i].acc0, ptcl_force_[adr_i].acc1,
                                              ptcl_multi_glb_[adr_j].pos,   ptcl_multi_glb_[adr_j].vel,
                                              ptcl_multi_glb_[adr_j].mass,   eps_sq,
                                              r_out, r_in);
                    }
                    calcAcc0AndAcc1FromSun(ptcl_multi_glb_[adr_i].pos, 
                                           ptcl_multi_glb_[adr_i].vel,
                                           ptcl_force_[adr_i].acc0,      
                                           ptcl_force_[adr_i].acc1);
                    ptcl_multi_glb_[adr_i].acc0 = ptcl_force_[adr_i].acc0;
                    ptcl_multi_glb_[adr_i].acc1 = ptcl_force_[adr_i].acc1;
                    ptcl_multi_glb_[adr_i].setDt2nd(ptcl_force_[adr_i], 0.01*eta, dt_limit_new, a0_offset_sq_);
                }
            }
        }
#endif
    }

    bool evolve(const PS::F64 r_out, const PS::F64 r_in,
                const PS::F64 eta,   const PS::F64 dt_limit_org){
        PS::F64 offset = PS::GetWtimeNoBarrier();
        predictAll();
        wtime_profile_.predict_ += PS::GetWtimeNoBarrier() - offset;

        offset = PS::GetWtimeNoBarrier();
        calcAcc0AndAcc1(r_out, r_in);
        wtime_profile_.force_ += PS::GetWtimeNoBarrier() - offset;

        offset = PS::GetWtimeNoBarrier();
        correctIp(eta, dt_limit_org, r_out, r_in);
        wtime_profile_.correct_ += PS::GetWtimeNoBarrier() - offset;
	
#ifdef FORDEBUG	
        if(n_active_ == ptcl_multi_glb_.size()){
            calc_eng();
            std::cout<<"time_sys_="<<time_sys_
		     <<" time_end_="<<time_end_
		     <<" loop="<< loop
		     <<" eng_tot_init_="<<eng_tot_init_
		     <<" eng_tot_="<<eng_tot_
		     <<" (eng_tot_init_ - eng_tot_)/eng_tot_init_="
                     <<(eng_tot_init_ - eng_tot_)/eng_tot_init_<<std::endl;
            std::cout<<"ptcl_multi_glb_[0].pos="<<ptcl_multi_glb_[0].pos<<std::endl;
            std::cout<<std::endl;
        }
#endif
        offset = PS::GetWtimeNoBarrier();
        sortAndSelectIp();
        wtime_profile_.sort_ += PS::GetWtimeNoBarrier() - offset;
	
        n_loop_++;
        if(time_sys_ == time_end_) return true;
        else return false;
    }

    bool evolveKepler(const PS::F64 r_out,        const PS::F64 r_in,
		      const PS::F64 eta,          const PS::F64 dt_limit_org,
		      const PS::F64 mass_sun=1.0, const PS::F64vec pos_sun=0.0,
		      const PS::F64vec vel_sun=0.0){

        PS::F64 offset = PS::GetWtimeNoBarrier();
        predictAll();
        wtime_profile_.predict_ += PS::GetWtimeNoBarrier() - offset;

        offset = PS::GetWtimeNoBarrier();
        calcAcc0AndAcc1_tmp(r_out, r_in);
        wtime_profile_.force_ += PS::GetWtimeNoBarrier() - offset;

        offset = PS::GetWtimeNoBarrier();
        calcAcc0AndAcc1Kepler(mass_sun, pos_sun, vel_sun);
        wtime_profile_.force_kepler_ += PS::GetWtimeNoBarrier() - offset;

	reduce_force();

        offset = PS::GetWtimeNoBarrier();
        correctIp(eta, dt_limit_org, r_out, r_in);
        wtime_profile_.correct_ += PS::GetWtimeNoBarrier() - offset;

        offset = PS::GetWtimeNoBarrier();
        sortAndSelectIp();
        wtime_profile_.sort_ += PS::GetWtimeNoBarrier() - offset;

        //n_loop_++;
        if(time_sys_ == time_end_) return true;
        else return false;
    }

    void scatterData(){
        const PS::S32 n_loc = ptcl_multi_loc_.size();
        PS::Comm::scatterV(&ptcl_multi_glb_[0], &n_ptcl_array_[0], &n_ptcl_disp_[0],
                           &ptcl_multi_loc_[0], n_loc);

	
#if 0
        if(PS::Comm::getRank() == 1){
            for(PS::S32 i=0; i<n_loc; i++){
                ptcl_hard_loc[i].dump();
            }
        }
#endif
    }

    template<class Tsoft>
    void copyPtclHardToSoft(Tsoft & system_soft){
        const PS::S32 n = ptcl_multi_loc_.size();
#pragma omp parallel for
        for(PS::S32 i=0; i<n; i++){
            //const PS::S32 adr = adr_ptcl_multi_loc_[i]; // by M.I.
            const PS::S32 adr = ptcl_multi_loc_[i].adr_fp; // by M.I.
            assert(system_soft[adr].id == ptcl_multi_loc_[i].id);
            system_soft[adr].mass = ptcl_multi_loc_[i].mass;
            system_soft[adr].pos = ptcl_multi_loc_[i].pos;
            system_soft[adr].vel = ptcl_multi_loc_[i].vel;
        }
    }

    PS::F64 eng_tot_init_;
    PS::F64 eng_kin_init_;
    PS::F64 eng_pot_init_;

    PS::F64 eng_tot_;
    PS::F64 eng_kin_;
    PS::F64 eng_pot_;

    static void calc_eng_impl(const std::vector<PTCLHard> ptcl, PS::F64 & et, PS::F64 & ek, PS::F64 & ep){
        const PS::S32 n = ptcl.size();
        et = ek = ep = 0.0;
        for(PS::S32 i=0; i<n; i++){
            const PTCLHard & p = ptcl[i];
            ek += 0.5 * p.mass * p.vel * p.vel;
            ep += 0.5 * p.mass * p.pot;
        }
        et = ek + ep;
    }

    void calc_eng_init(){
        calc_eng_impl(ptcl_multi_glb_, eng_tot_init_, eng_kin_init_, eng_pot_init_);
    }

    void calc_eng(){
        calc_eng_impl(ptcl_multi_glb_, eng_tot_, eng_kin_, eng_pot_);
    }

    void predict2body(PTCLPred & pred_i,
		      PTCLPred & pred_j,
                      const PTCLHard & ptcl_i, 
		      const PTCLHard & ptcl_j, 
		      const PS::F64 dt){
        static const PS::F64 inv3 = 1.0 / 3.0;
        pred_i.pos = ptcl_i.pos + dt*(ptcl_i.vel  + 0.5*dt*(ptcl_i.acc0 + inv3*dt*ptcl_i.acc1));
        pred_i.vel = ptcl_i.vel + dt*(ptcl_i.acc0 + 0.5*dt* ptcl_i.acc1);
        pred_j.pos = ptcl_j.pos + dt*(ptcl_j.vel  + 0.5*dt*(ptcl_j.acc0 + inv3*dt*ptcl_j.acc1));
        pred_j.vel = ptcl_j.vel + dt*(ptcl_j.acc0 + 0.5*dt* ptcl_j.acc1);
    }

    std::vector< std::pair<PS::S64, PS::S64> > pair_id_2body_loc_;
    std::vector< std::pair<PS::S64, PS::S64> > pair_id_multi_loc_;
    std::vector< std::pair<PS::S64, PS::S64> > pair_adr_2body_loc_; // point to ptcl_loc_
    std::vector<PTCLHard> ptcl_2body_loc_;
    std::vector<PTCLHard> ptcl_multi_loc_;
    //std::vector<PS::S64> adr_ptcl_2body_loc_;
    //std::vector<PS::S64> adr_ptcl_multi_loc_;



#ifndef ORIGINAL
    // new 
    void selectIsolatedSystem(){
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
	//static std::vector<PTCLHard> ptcl_2body_cand_send;
	//static std::vector<PTCLHard> ptcl_2body_cand_recv;
	static std::vector<PTCLHardComm> ptcl_2body_cand_send;
	static std::vector<PTCLHardComm> ptcl_2body_cand_recv;
	static std::vector<PS::S32>  rank_j_2body_cand;
	static std::vector<PS::S32>  adr_i_2body_cand_buf;
	static std::vector<PS::S32>  adr_i_2body_cand_sort;
	//static std::vector<PS::S64> id_j_2body_cand_buf; // not needed
	//static std::vector<PS::S64> id_j_2body_cand_send;
	//static std::vector<PS::S64> id_j_2body_cand_recv;
	static std::vector<PS::S32> adr_j_2body_cand_buf;
	static std::vector<PS::S32> adr_j_2body_cand_send;
	static std::vector<PS::S32> adr_j_2body_cand_recv;
	static MPI_Request * req_send;
	static MPI_Request * req_recv;
	static MPI_Status * stat_send;
	static MPI_Status * stat_recv;
	static PS::S32 * n_2body_cand_send;
	static PS::S32 * n_2body_cand_send_disp;
	static PS::S32 * n_2body_cand_recv;
	static PS::S32 * n_2body_cand_recv_disp;
        static std::unordered_map<PS::S64, PS::S64> idx_to_adr_2body;
	static bool first = true;
	if(first){
	    req_send  = new MPI_Request[n_proc+10];
	    req_recv  = new MPI_Request[n_proc+10];
	    stat_send = new MPI_Status[n_proc+10];
	    stat_recv = new MPI_Status[n_proc+10];
	    n_2body_cand_send = new PS::S32[n_proc+10];
	    n_2body_cand_send_disp = new PS::S32[n_proc+10+1];
	    n_2body_cand_recv = new PS::S32[n_proc+10];
	    n_2body_cand_recv_disp = new PS::S32[n_proc+10+1];
	    first = false;
	}

#pragma omp parallel for
	for(PS::S32 i=0; i<n_proc; i++){
	    n_2body_cand_send[i] = n_2body_cand_recv[i] = 0;
	}
        makeDictionaryLocal();
        pair_id_multi_loc_.clear();
        pair_id_2body_loc_.clear();
        ptcl_2body_loc_.clear();
        ptcl_multi_loc_.clear();
        idx_to_adr_2body.clear();

	ptcl_2body_cand_send.clear();
	//id_j_2body_cand_buf.clear();
	//id_j_2body_cand_send.clear();
	adr_j_2body_cand_buf.clear();
	adr_j_2body_cand_send.clear();
	rank_j_2body_cand.clear();
	adr_i_2body_cand_buf.clear();
	adr_i_2body_cand_sort.clear();
        PS::S64 adr_i_old = -1;
        const size_t n_pair = pair_id_loc_.size();
        PS::S64 n_cnt_2body = 0;

        for(size_t i=0; i<n_pair; i++){
            const PS::S64 adr_i = pair_adr_loc_[i].first;
            const PS::S64 adr_j = pair_adr_loc_[i].second;
            if( ptcl_loc_[adr_i].n_ngb == 1 && adr_j >= 0 && ptcl_loc_[adr_j].n_ngb == 1){
		// both particles are in this domain.
                pair_id_2body_loc_.push_back( pair_id_loc_[i] );
                ptcl_2body_loc_.push_back( ptcl_loc_[adr_i] );
                idx_to_adr_2body.insert(std::pair<PS::S64, PS::S64>(ptcl_loc_[adr_i].id, n_cnt_2body));
                n_cnt_2body++;
            }
	    else if(ptcl_loc_[adr_i].n_ngb == 1 && adr_j < 0){
		//neighbor is in another domain.
		PS::S32 rank_j = rank_org_j_[i];
		PS::S32 flag = 0;
		for(PS::S32 i0=0; i0<static_cast<PS::S32>(rank_neighbor_.size()); i0++){
		    if(rank_j == rank_neighbor_[i0]) flag++;
		}
		assert(flag==1);
		//id_j_2body_cand_buf.push_back(pair_id_loc_[i].second);
		adr_j_2body_cand_buf.push_back(pair_adr_loc_org_[i].second);
		rank_j_2body_cand.push_back(rank_j);
		n_2body_cand_send[ rank_j ]++;
		adr_i_2body_cand_buf.push_back(adr_i);
	    }
            else{
		// the particle has more than one neighbors.
                pair_id_multi_loc_.push_back(pair_id_loc_[i]);
                if(adr_i != adr_i_old){
		    // push back ptcl_loc_[adr_i] once.
                    ptcl_multi_loc_.push_back( ptcl_loc_[adr_i] );
                    adr_i_old = adr_i;
                }
            }
        }
	/////////
	//// add
	PS::S32 n_proc_send = 0;
	PS::S32 n_proc_recv = 0;
	const PS::S32 n_proc_neighbor = rank_neighbor_.size();
	for(PS::S32 i=0; i<n_proc_neighbor; i++){
	    PS::S32 tag = 0;
	    PS::S32 rank = rank_neighbor_[i];
	    MPI_Isend(&n_2body_cand_send[rank], 1, PS::GetDataType<PS::S32>(), rank, tag, MPI_COMM_WORLD, req_send+n_proc_send);
	    n_proc_send++;
	    MPI_Irecv(&n_2body_cand_recv[rank], 1, PS::GetDataType<PS::S32>(), rank, tag, MPI_COMM_WORLD, req_recv+n_proc_recv);
	    n_proc_recv++;
	}
	Print(n_proc_neighbor, fout_debug);
	MPI_Waitall(n_proc_send, req_send, stat_send);
	MPI_Waitall(n_proc_recv, req_recv, stat_recv);

	// pack index of j-particle
	// NOTE: n_2body_cand_XXX with id of which j-process is not neighbor is undefined.
	const PS::S32 rank_head = rank_neighbor_[0];

	n_2body_cand_send_disp[rank_head] = n_2body_cand_recv_disp[rank_head] = 0;
	for(PS::S32 i=1; i<n_proc_neighbor; i++){
	    PS::S32 rank_prev = rank_neighbor_[i-1];
	    PS::S32 rank      = rank_neighbor_[i];
	    n_2body_cand_send_disp[rank]  = n_2body_cand_send_disp[rank_prev] + n_2body_cand_send[rank_prev];
	    n_2body_cand_recv_disp[rank]  = n_2body_cand_recv_disp[rank_prev] + n_2body_cand_recv[rank_prev];
	    n_2body_cand_send[rank_prev] = 0;
	}

	const PS::S32 rank_tail = rank_neighbor_[n_proc_neighbor-1];
	PS::S32 n_send_tot = n_2body_cand_send_disp[rank_tail] + n_2body_cand_send[rank_tail];
	PS::S32 n_recv_tot = n_2body_cand_recv_disp[rank_tail] + n_2body_cand_recv[rank_tail];
	n_2body_cand_send[rank_tail] = 0;

	Print("debug e5", fout_debug);

	//id_j_2body_cand_send.resize(n_send_tot);
	adr_j_2body_cand_send.resize(n_send_tot);
	adr_i_2body_cand_sort.resize(n_send_tot);
	for(PS::S32 i=0; i<n_send_tot; i++){
	    PS::S32 rank = rank_j_2body_cand[i];
	    PS::S32 adr = n_2body_cand_send_disp[rank] + n_2body_cand_send[rank];
	    //id_j_2body_cand_send[adr]   = id_j_2body_cand_buf[i];
	    adr_j_2body_cand_send[adr]  = adr_j_2body_cand_buf[i];
	    adr_i_2body_cand_sort[adr]  = adr_i_2body_cand_buf[i];
	    n_2body_cand_send[rank]++;
	}
	// exchange index of j-particle
	//id_j_2body_cand_recv.resize(n_recv_tot);
	adr_j_2body_cand_recv.resize(n_recv_tot);
	n_proc_send = n_proc_recv = 0;
	for(PS::S32 i=0; i<n_proc_neighbor; i++){
	    PS::S32 rank = rank_neighbor_[i];
	    PS::S32 tag = 0;
	    if(n_2body_cand_send[rank] > 0 ){
		//MPI_Isend(&id_j_2body_cand_send[n_2body_cand_send_disp[rank]], n_2body_cand_send[rank], PS::GetDataType<PS::S64>(), rank, 10, MPI_COMM_WORLD, req_send+n_proc_send);
		MPI_Isend(&adr_j_2body_cand_send[n_2body_cand_send_disp[rank]], n_2body_cand_send[rank], PS::GetDataType<PS::S32>(), rank, tag, MPI_COMM_WORLD, req_send+n_proc_send);
		n_proc_send++;
	    }
	    if(n_2body_cand_recv[rank] > 0 ){
		//MPI_Irecv(&id_j_2body_cand_recv[n_2body_cand_recv_disp[rank]], n_2body_cand_recv[rank], PS::GetDataType<PS::S64>(), rank, 10, MPI_COMM_WORLD, req_recv+n_proc_recv);
		MPI_Irecv(&adr_j_2body_cand_recv[n_2body_cand_recv_disp[rank]], n_2body_cand_recv[rank], PS::GetDataType<PS::S32>(), rank, tag, MPI_COMM_WORLD, req_recv+n_proc_recv);
		n_proc_recv++;
	    }
	}
	MPI_Waitall(n_proc_send, req_send, stat_send);
	MPI_Waitall(n_proc_recv, req_recv, stat_recv);

	// set j-particles to send to i-process
	ptcl_2body_cand_send.clear();
	ptcl_2body_cand_send.resize(n_recv_tot);
#pragma omp parallel for
	for(PS::S32 i=0; i<n_recv_tot; i++){
	    //const PS::S32 id = id_j_2body_cand_recv[i];
	    //const PS::S32 adr_b = idx_to_adr_loc_[id];
	    //assert(adr_b >= 0);
	    //assert(ptcl_loc_[adr_b].n_ngb > 0);
	    //assert( id == ptcl_loc_[adr_b].id);
	    //ptcl_2body_cand_send[i].copyFromPH(ptcl_loc_[adr_b]);


	    const PS::S32 adr_org = adr_j_2body_cand_recv[i];
	    assert(adr_org >= 0);
	    ptcl_2body_cand_send[i].copyFromPHC(ptcl_loc_org_[adr_org]);

	    /*
	    fout_debug<<"ptcl_loc_[adr_b].id="<<ptcl_loc_[adr_b].id<<" ptcl_loc_org_[adr_org].id="<<ptcl_loc_org_[adr_org].id<<std::endl;
	    fout_debug<<"ptcl_loc_[adr_b].mass="<<ptcl_loc_[adr_b].mass<<" ptcl_loc_org_[adr_org].mass="<<ptcl_loc_org_[adr_org].mass<<std::endl;
	    fout_debug<<"ptcl_loc_[adr_b].pos="<<ptcl_loc_[adr_b].pos<<" ptcl_loc_org_[adr_org].pos="<<ptcl_loc_org_[adr_org].pos<<std::endl;
	    fout_debug<<"ptcl_loc_[adr_b].vel="<<ptcl_loc_[adr_b].vel<<" ptcl_loc_org_[adr_org].vel="<<ptcl_loc_org_[adr_org].vel<<std::endl;
	    fout_debug<<"ptcl_loc_[adr_b].n_ngb="<<ptcl_loc_[adr_b].n_ngb<<" ptcl_loc_org_[adr_org].n_ngb="<<ptcl_loc_org_[adr_org].n_ngb<<std::endl;
	    fout_debug<<"ptcl_loc_[adr_b].rank_org="<<ptcl_loc_[adr_b].rank_org<<" ptcl_loc_org_[adr_org].rank_org="<<ptcl_loc_org_[adr_org].rank_org<<std::endl;
	    fout_debug<<"ptcl_loc_[adr_b].adr_fp="<<ptcl_loc_[adr_b].adr_fp<<" ptcl_loc_org_[adr_org].adr_fp="<<ptcl_loc_org_[adr_org].adr_fp<<std::endl;
	    assert(ptcl_loc_[adr_b].id==ptcl_loc_org_[adr_org].id);
	    assert(ptcl_loc_[adr_b].n_ngb==ptcl_loc_org_[adr_org].n_ngb);
	    assert(ptcl_loc_[adr_b].mass==ptcl_loc_org_[adr_org].mass);
	    assert(ptcl_loc_[adr_b].pos.x==ptcl_loc_org_[adr_org].pos.x);
	    assert(ptcl_loc_[adr_b].vel.x==ptcl_loc_org_[adr_org].vel.x);
	    assert(ptcl_loc_[adr_b].rank_org==ptcl_loc_org_[adr_org].rank_org);
	    assert(ptcl_loc_[adr_b].adr_fp==ptcl_loc_org_[adr_org].adr_fp);
	    */
	}

	Print("debug e6", fout_debug);

	///// send back particle with the index of which the particles received
	ptcl_2body_cand_recv.resize(n_send_tot);
	n_proc_send = n_proc_recv = 0;
	for(PS::S32 i=0; i<n_proc_neighbor; i++){
	    PS::S32 rank = rank_neighbor_[i];
	    PS::S32 tag = 0;
	    if(n_2body_cand_recv[rank] > 0 ){
		MPI_Isend(&ptcl_2body_cand_send[n_2body_cand_recv_disp[rank]], n_2body_cand_recv[rank], PS::GetDataType<PTCLHardComm>(), rank, tag, MPI_COMM_WORLD, req_send+n_proc_send);
		n_proc_send++;
	    }
	    if(n_2body_cand_send[rank] > 0 ){
		MPI_Irecv(&ptcl_2body_cand_recv[n_2body_cand_send_disp[rank]], n_2body_cand_send[rank], PS::GetDataType<PTCLHardComm>(), rank, tag, MPI_COMM_WORLD, req_recv+n_proc_recv);
		n_proc_recv++;
	    }
	}
	MPI_Waitall(n_proc_send, req_send, stat_send);
	MPI_Waitall(n_proc_recv, req_recv, stat_recv);
#ifdef PARALLEL_2BODY
	PS::S32 corr = 0; // for debug
	for(PS::S32 i=0; i<n_send_tot; i++){
	    assert(ptcl_2body_cand_recv[i].n_ngb >= 1);
	    if(ptcl_2body_cand_recv[i].n_ngb > 1){
		// the particle has more than one neighbors.
		const PS::S64 id_j = ptcl_2body_cand_recv[i].id; //ptcl_2body_cand_recv is j particle
		const PS::S32 adr_i = adr_i_2body_cand_sort[i];
                pair_id_multi_loc_.push_back( std::pair<PS::S64, PS::S64>(ptcl_loc_[adr_i].id, id_j) );
		ptcl_multi_loc_.push_back( ptcl_loc_[adr_i] );
	    }
	    else if(ptcl_2body_cand_recv[i].n_ngb == 1){
		// this system can be isolated
		//const PTCLHard & ptcl_j = ptcl_2body_cand_recv[i];
		const PTCLHardComm & ptcl_j = ptcl_2body_cand_recv[i];
		const PS::S64 id_j = ptcl_j.id; //ptcl_2body_cand_recv is j particle
		const PS::S32 adr_i = adr_i_2body_cand_sort[i];
		const PTCLHard & ptcl_i = ptcl_loc_[adr_i];
		const PS::S64 id_i = ptcl_i.id;
		const PS::F64 mass_i = ptcl_i.mass;
		const PS::F64 mass_j = ptcl_j.mass;
		assert(id_i != id_j);
		if( (mass_i > mass_j) || ( (mass_i == mass_j) && (id_i < id_j) ) ){
		    pair_id_2body_loc_.push_back( std::pair<PS::S64, PS::S64>(id_i, id_j) );
		    ptcl_2body_loc_.push_back( ptcl_i );
		    idx_to_adr_2body.insert(std::pair<PS::S64, PS::S64>(id_i, n_cnt_2body));
		    n_cnt_2body++;

		    pair_id_2body_loc_.push_back( std::pair<PS::S64, PS::S64>(id_j, id_i) );
		    //ptcl_2body_loc_.push_back( ptcl_j );
		    ptcl_2body_loc_.push_back( PTCLHard(ptcl_j.id,    ptcl_j.mass,
							ptcl_j.pos,   ptcl_j.vel,
							ptcl_j.n_ngb, ptcl_j.rank_org,
							ptcl_j.adr_fp) );
		    idx_to_adr_2body.insert(std::pair<PS::S64, PS::S64>(id_j, n_cnt_2body));
		    n_cnt_2body++;
		    corr++;
		}
		else corr--;
	    }
	}
	if(pair_id_loc_.size()+corr != pair_id_multi_loc_.size() + pair_id_2body_loc_.size()){
	    std::cout<<"pair_id_loc_.size()="<<pair_id_loc_.size()
		     <<" pair_id_multi_loc_.size()="<<pair_id_multi_loc_.size()
		     <<" pair_id_2body_loc_.size()="<<pair_id_2body_loc_.size()<<std::endl;
	}
	if( ptcl_loc_.size()+corr != ptcl_multi_loc_.size() + ptcl_2body_loc_.size() ){
	    std::cout<<"ptcl_loc_.size()="<<ptcl_loc_.size()
		     <<" ptcl_multi_loc_.size()="<<ptcl_multi_loc_.size()
		     <<" ptcl_2body_loc_.size()="<<ptcl_2body_loc_.size()<<std::endl;
	}
        assert( pair_id_loc_.size()+corr == pair_id_multi_loc_.size() + pair_id_2body_loc_.size() );
        assert( ptcl_loc_.size()+corr    == ptcl_multi_loc_.size() + ptcl_2body_loc_.size() );
#else //PARALLEL_2BODY
	for(PS::S32 i=0; i<n_send_tot; i++){
	    assert(ptcl_2body_cand_recv[i].n_ngb >= 1);
	    if(1){ // for debug
		//const PTCLHard & ptcl_j = ptcl_2body_cand_recv[i];
		const PTCLHardComm & ptcl_j = ptcl_2body_cand_recv[i];
		const PS::S64 id_j = ptcl_j.id; //ptcl_2body_cand_recv is j particle
		const PS::S32 adr_i = adr_i_2body_cand_sort[i];
		const PTCLHard & ptcl_i = ptcl_loc_[adr_i];
		const PS::S64 id_i = ptcl_i.id;
                pair_id_multi_loc_.push_back( std::pair<PS::S64, PS::S64>(id_i, id_j) );
		ptcl_multi_loc_.push_back( ptcl_i );
	    }
	}
	if(pair_id_loc_.size() != pair_id_multi_loc_.size() + pair_id_2body_loc_.size()){
	    std::cout<<"pair_id_loc_.size()="<<pair_id_loc_.size()
		     <<" pair_id_multi_loc_.size()="<<pair_id_multi_loc_.size()
		     <<" pair_id_2body_loc_.size()="<<pair_id_2body_loc_.size()<<std::endl;
	}
	if( ptcl_loc_.size() != ptcl_multi_loc_.size() + ptcl_2body_loc_.size() ){
	    std::cout<<"ptcl_loc_.size()="<<ptcl_loc_.size()
		     <<" ptcl_multi_loc_.size()="<<ptcl_multi_loc_.size()
		     <<" ptcl_2body_loc_.size()="<<ptcl_2body_loc_.size()<<std::endl;
	}
        assert( pair_id_loc_.size() == pair_id_multi_loc_.size() + pair_id_2body_loc_.size() );
        assert( ptcl_loc_.size() == ptcl_multi_loc_.size() + ptcl_2body_loc_.size() );
#endif //PARALLEL_2BODY
	////// check if the particles is isolated or not.
	//// add
	/////////
        pair_adr_2body_loc_.clear();
        const size_t n_pair_2body = pair_id_2body_loc_.size();
        for(size_t i=0; i<n_pair_2body; i++){
            pair_adr_2body_loc_.push_back(std::pair<PS::S64, PS::S64>
                                          (idx_to_adr_2body[pair_id_2body_loc_[i].first],
                                           idx_to_adr_2body[pair_id_2body_loc_[i].second]));
        }
    }

#else //ORIGINAL
    // original
    void selectIsolatedSystem(){
        std::unordered_map<PS::S64, PS::S64> idx_to_adr_2body;
        makeDictionaryLocal();
        pair_id_multi_loc_.clear();
        pair_id_2body_loc_.clear();
        ptcl_2body_loc_.clear();
        ptcl_multi_loc_.clear();
        idx_to_adr_2body.clear();
        PS::S64 adr_i_old = -1;
        const size_t n_pair = pair_id_loc_.size();
        PS::S64 n_cnt_2body = 0;
        for(size_t i=0; i<n_pair; i++){
            const PS::S64 adr_i = pair_adr_loc_[i].first;
            const PS::S64 adr_j = pair_adr_loc_[i].second;
            if( ptcl_loc_[adr_i].n_ngb == 1 && adr_j >= 0 && ptcl_loc_[adr_j].n_ngb == 1){
                pair_id_2body_loc_.push_back( pair_id_loc_[i] );
                ptcl_2body_loc_.push_back( ptcl_loc_[adr_i] );
                idx_to_adr_2body.insert(std::pair<PS::S64, PS::S64>(ptcl_loc_[adr_i].id, n_cnt_2body));
                n_cnt_2body++;
            }
            else{
                pair_id_multi_loc_.push_back(pair_id_loc_[i]);
                if(adr_i != adr_i_old){
                    ptcl_multi_loc_.push_back( ptcl_loc_[adr_i] );
                    adr_i_old = adr_i;
                }
            }
        }
	if(pair_id_loc_.size() != pair_id_multi_loc_.size() + pair_id_2body_loc_.size()){
	    std::cerr<<"pair_id_loc_.size()="<<pair_id_loc_.size()
		     <<" pair_id_multi_loc_.size()="<<pair_id_multi_loc_.size()
		     <<" pair_id_2body_loc_.size()="<<pair_id_2body_loc_.size()<<std::endl;
	}
	if( ptcl_loc_.size() != ptcl_multi_loc_.size() + ptcl_2body_loc_.size() ){
	    std::cerr<<"ptcl_loc_.size()="<<ptcl_loc_.size()
		     <<" ptcl_multi_loc_.size()="<<ptcl_multi_loc_.size()
		     <<" ptcl_2body_loc_.size()="<<ptcl_2body_loc_.size()<<std::endl;
	}
        assert( pair_id_loc_.size() == pair_id_multi_loc_.size() + pair_id_2body_loc_.size() );
        assert( ptcl_loc_.size() == ptcl_multi_loc_.size() + ptcl_2body_loc_.size() );
        pair_adr_2body_loc_.clear();
        const size_t n_pair_2body = pair_id_2body_loc_.size();
        for(size_t i=0; i<n_pair_2body; i++){
            pair_adr_2body_loc_.push_back(std::pair<PS::S64, PS::S64>
                                          (idx_to_adr_2body[pair_id_2body_loc_[i].first],
                                           idx_to_adr_2body[pair_id_2body_loc_[i].second]));
        }
    }
#endif


    
    std::vector<PTCLPred> ptcl_pred_2body_loc_;
    std::vector<PTCLForce> ptcl_force_2body_loc_;
    std::vector<bool> merge_flag_2body_loc_;
#ifndef ORIGINAL
    template<class Tsoft>
    void evolveIsolatedSystem(Tsoft & system,
                              const PS::F64 r_out, 
			      const PS::F64 r_in,
                              const PS::F64 mass_sun, 
			      const PS::F64vec pos_sun,
                              const PS::F64vec vel_sun,
                              const PS::F64 eta_s, 
			      const PS::F64 eta,
                              const PS::F64 time_end, 
			      const PS::F64 dt_limit_org){
        const PS::S32 n_thread_max = PS::Comm::getNumberOfThread();
        static std::vector<MergeLog> * merge_log;
        static PS::F64 * eng_disp_omp;
#ifdef CALC_HARD_ENERGY
        PS::F64 eng_hard_disp = 0.0;
        static PS::F64 * eng_hard_disp_omp;
#endif
        static bool first = true;
        if(first){
            merge_log = new std::vector<MergeLog>[n_thread_max];
            eng_disp_omp = new PS::F64[n_thread_max];
#ifdef CALC_HARD_ENERGY
            eng_hard_disp_omp = new PS::F64[n_thread_max];
#endif
            first = false;
        }
        for(PS::S32 i=0; i<n_thread_max; i++){
            merge_log[i].clear();
            eng_disp_omp[i] = 0.0;
#ifdef CALC_HARD_ENERGY
            eng_hard_disp_omp[i] = 0.0;
#endif
        }
        const size_t n_2body = ptcl_2body_loc_.size();
        ptcl_pred_2body_loc_.clear();
        ptcl_pred_2body_loc_.resize(n_2body);
        ptcl_force_2body_loc_.clear();
        ptcl_force_2body_loc_.resize(n_2body);
        merge_flag_2body_loc_.clear();
        merge_flag_2body_loc_.resize(n_2body);
#pragma omp parallel for
        for(size_t ip=0; ip<n_2body; ip++){
            ptcl_pred_2body_loc_[ip].pos = ptcl_2body_loc_[ip].pos;
            ptcl_pred_2body_loc_[ip].vel = ptcl_2body_loc_[ip].vel;
            ptcl_2body_loc_[ip].time = 0.0;
            ptcl_2body_loc_[ip].setRMerge();
            merge_flag_2body_loc_[ip] = false;
        }
        const PS::S32 n_pair = pair_adr_2body_loc_.size();
        const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;

#ifdef CALC_HARD_ENERGY
        CalcTwoBodyEnergy(eng_2body_prev_loc_, r_in, r_out);
#endif

#pragma omp parallel for
        for(PS::S32 np=0; np<n_pair; np++){
            //PS::S32 ith = PS::Comm::getThreadNum();
            PS::F64 time_sys_tmp = 0.0;
            PS::F64 time_sync_tmp = time_sys_tmp + dt_limit_org;
            PS::F64 time_end_tmp = time_end;

            const PS::S64 adr_i = pair_adr_2body_loc_[np].first;
            const PS::S64 adr_j = pair_adr_2body_loc_[np].second;
            assert(adr_i >= 0);
            assert(adr_j >= 0);
            if(adr_i > adr_j) continue;
            PTCLHard & ptcl_i = ptcl_2body_loc_[adr_i];
            PTCLHard & ptcl_j = ptcl_2body_loc_[adr_j];
/*
            if(ptcl_i.rank_org != PS::Comm::getRank() || ptcl_j.rank_org != PS::Comm::getRank()){
                std::cerr<<"PS::Comm::getRank()="<<PS::Comm::getRank()<<" ptcl_i.rank_org="<<ptcl_i.rank_org<<" ptcl_j.rank_org="<<ptcl_j.rank_org<<std::endl;
                std::cerr<<"PS::Comm::getRank()="<<PS::Comm::getRank()<<" ptcl_i.pos="<<ptcl_i.pos<<" ptcl_j.pos="<<ptcl_j.pos<<std::endl;
            }
            assert(ptcl_i.rank_org == PS::Comm::getRank());
            assert(ptcl_j.rank_org == PS::Comm::getRank());
*/
/*
            if(ptcl_i.rank_org != PS::Comm::getRank() || ptcl_j.rank_org != PS::Comm::getRank()){
                ptcl_i.dump(fout_debug);
                ptcl_j.dump(fout_debug);
            }
*/

            PTCLForce & force_i = ptcl_force_2body_loc_[adr_i];
            PTCLForce & force_j = ptcl_force_2body_loc_[adr_j];
            PTCLPred & pred_i = ptcl_pred_2body_loc_[adr_i];
            PTCLPred & pred_j = ptcl_pred_2body_loc_[adr_j];
            PS::F64 r2 = 0.0;
            force_i.clear();
            force_j.clear();
            CalcAcc0Acc1AndR2CutoffPair(ptcl_i.mass, pred_i.pos, pred_i.vel, force_i.acc0_pla, force_i.acc1_pla,
                                        ptcl_j.mass, pred_j.pos, pred_j.vel, force_j.acc0_pla, force_j.acc1_pla,
                                        r2,   eps_sq, r_out, r_in);

            CalcAcc0Acc1AndR2(pred_i.pos,    pred_i.vel,
                              force_i.acc0,  force_i.acc1,
                              pos_sun, vel_sun, mass_sun, 0.0);
	    
            CalcAcc0Acc1AndR2(pred_j.pos,     pred_j.vel,
                              force_j.acc0,   force_j.acc1,
                              pos_sun, vel_sun, mass_sun, 0.0);
            force_i.reduce();
            force_j.reduce();

            PS::F64 dt_limit_new = CalcDtLimit(time_sys_tmp, time_sync_tmp, dt_limit_org);
            ptcl_i.setDt2nd(force_i, eta_s, dt_limit_new, a0_offset_sq_);
            ptcl_j.setDt2nd(force_j, eta_s, dt_limit_new, a0_offset_sq_);
            PS::F64 dt = (ptcl_i.dt < ptcl_j.dt) ? ptcl_i.dt : ptcl_j.dt;
            ptcl_i.dt = ptcl_j.dt = dt;
            ptcl_i.acc0 = force_i.acc0;
            ptcl_i.acc1 = force_i.acc1;
            ptcl_j.acc0 = force_j.acc0;
            ptcl_j.acc1 = force_j.acc1;
            ptcl_i.acc0_pla = force_i.acc0_pla;
            ptcl_i.acc1_pla = force_i.acc1_pla;
            ptcl_j.acc0_pla = force_j.acc0_pla;
            ptcl_j.acc1_pla = force_j.acc1_pla;
            while(time_sys_tmp != time_end_tmp){
                n_loop_++;
                n_interaction_2body_loc_ += 2; // plus 2 means using force symmetry.
                predict2body(pred_i, pred_j,
                             ptcl_i, ptcl_j,
                             dt);
                force_i.clear();
                force_j.clear();
                CalcAcc0Acc1AndR2CutoffPair(ptcl_i.mass, pred_i.pos, pred_i.vel, force_i.acc0_pla, force_i.acc1_pla,
                                            ptcl_j.mass, pred_j.pos, pred_j.vel, force_j.acc0_pla, force_j.acc1_pla,
                                            r2,   eps_sq, r_out, r_in);
                CalcAcc0Acc1AndR2(pred_i.pos,	  pred_i.vel,
                                  force_i.acc0,   force_i.acc1,
                                  pos_sun, vel_sun, mass_sun, 0.0);
                CalcAcc0Acc1AndR2(pred_j.pos,     pred_j.vel,
                                  force_j.acc0,   force_j.acc1,
                                  pos_sun, vel_sun, mass_sun, 0.0);
                force_i.reduce();
                force_j.reduce();
#ifdef MERGE
                const PS::F64 r_merge = ptcl_i.r_merge + ptcl_j.r_merge;
                const PS::F64 r_merge_2 = r_merge * r_merge;
                bool merge = false;
                if(r2 < r_merge_2){
                    //pair_isolated_merge_adr_.push_back( std::pair<PS::S64, PS::S64>(adr_i, adr_j) );
                    merge = true;
                }
#endif //MERGE
                time_sys_tmp += dt;
                if(time_sys_tmp == time_sync_tmp) time_sync_ += dt_limit_org;
                dt_limit_new = CalcDtLimit(time_sys_tmp, time_sync_tmp, dt_limit_org);
                ptcl_i.correct(force_i, eta, dt_limit_new, a0_offset_sq_);
                ptcl_j.correct(force_j, eta, dt_limit_new, a0_offset_sq_);
                dt = (ptcl_2body_loc_[adr_i].dt < ptcl_2body_loc_[adr_j].dt) ? ptcl_2body_loc_[adr_i].dt : ptcl_2body_loc_[adr_j].dt;
                ptcl_i.dt = ptcl_j.dt = dt;
#ifdef MERGE
                if(merge){
		    PS::S32 ith = PS::Comm::getThreadNum();
#ifdef CALC_HARD_ENERGY
                    merge2bodyWithHardEnergy(ptcl_i, ptcl_j, eng_disp_omp[ith], merge_log[ith], r_in, r_out, eng_hard_disp_omp[ith]);
#else
                    merge2body(ptcl_i, ptcl_j, eng_disp_omp[ith], merge_log[ith]);
#endif
                    PTCLHard & ptcl_merge = (ptcl_i.mass == 0.0) ? ptcl_j : ptcl_i;
                    const PS::F64 dt_rem = time_end_tmp - time_sys_tmp;

                    DriveKepler(mass_sun, ptcl_merge.mass,
                                pos_sun,  ptcl_merge.pos,
                                vel_sun,  ptcl_merge.vel, dt_rem);

                    ptcl_merge.time = time_end_tmp;
                    time_sys_tmp = time_end_tmp;
                }
#endif //MERGE
            }
        }
        for(PS::S32 i=0; i<n_thread_max; i++){
            const size_t j_max = merge_log[i].size();
            for(size_t j=0; j<j_max; j++){
                merge_history_.push_back(merge_log[i][j]);
            }
            eng_disp_ += eng_disp_omp[i];
#ifdef CALC_HARD_ENERGY
            eng_hard_disp += eng_hard_disp_omp[i];
#endif
        }

#ifdef PARALLEL_2BODY
        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        static bool first2 = true;
        static PS::S32 * n_send;
        static PS::S32 * n_recv;
        static PS::S32 * n_send_disp;
        static PS::S32 * n_recv_disp;
        static MPI_Request * req_send;
        static MPI_Request * req_recv;
        static MPI_Status * stat_send;
        static MPI_Status * stat_recv;
        static std::vector<PTCLHardComm> ptcl_buf;
        static std::vector<PTCLHardComm> ptcl_send;
        static std::vector<PTCLHardComm> ptcl_recv;
        static std::vector<PS::S32> rank_j;
        if(first2){
            n_send = new PS::S32[n_proc];
            n_recv = new PS::S32[n_proc];
            n_send_disp = new PS::S32[n_proc+1];
            n_recv_disp = new PS::S32[n_proc+1];
            req_send  = new MPI_Request[n_proc];
            req_recv  = new MPI_Request[n_proc];
            stat_send = new MPI_Status[n_proc];
            stat_recv = new MPI_Status[n_proc];
            first2 = false;
        }
        ptcl_buf.clear();
        rank_j.clear();
#pragma omp parallel for
        for(PS::S32 i=0; i<n_proc; i++){
            n_send[i] = n_recv[i] = 0;
        }
        for(size_t i=0; i<n_2body; i++){
            if(ptcl_2body_loc_[i].rank_org == PS::Comm::getRank()){
                const PS::S64 adr = ptcl_2body_loc_[i].adr_fp; // M.I.
                assert(ptcl_2body_loc_[i].id == system[adr].id);
                system[adr].mass = ptcl_2body_loc_[i].mass;
                system[adr].pos = ptcl_2body_loc_[i].pos;
                system[adr].vel = ptcl_2body_loc_[i].vel;
            }
            else{
                n_send[ ptcl_2body_loc_[i].rank_org ]++;
                ptcl_buf.push_back(ptcl_2body_loc_[i]);
                rank_j.push_back(ptcl_2body_loc_[i].rank_org);
            }
        }
        PS::S32 n_proc_send = 0;
        PS::S32 n_proc_recv = 0;
        const PS::S32 n_proc_neighbor = rank_neighbor_.size();
        for(PS::S32 i=0; i<n_proc_neighbor; i++){
            PS::S32 tag = 0;
            PS::S32 rank = rank_neighbor_[i];
            MPI_Isend(&n_send[rank], 1, PS::GetDataType<PS::S32>(), rank, tag, MPI_COMM_WORLD, req_send+n_proc_send);
            n_proc_send++;
            MPI_Irecv(&n_recv[rank], 1, PS::GetDataType<PS::S32>(), rank, tag, MPI_COMM_WORLD, req_recv+n_proc_recv);
            n_proc_recv++;
        }
        MPI_Waitall(n_proc_send, req_send, stat_send);
        MPI_Waitall(n_proc_recv, req_recv, stat_recv);

        const PS::S32 rank_head = rank_neighbor_[0];
        n_send_disp[rank_head] = n_recv_disp[rank_head] = 0;
        for(PS::S32 i=1; i<n_proc_neighbor; i++){
            PS::S32 rank_prev = rank_neighbor_[i-1];
            PS::S32 rank      = rank_neighbor_[i];
            n_send_disp[rank]  = n_send_disp[rank_prev] + n_send[rank_prev];
            n_recv_disp[rank]  = n_recv_disp[rank_prev] + n_recv[rank_prev];
            n_send[rank_prev] = 0;
        }
        const PS::S32 rank_tail = rank_neighbor_[n_proc_neighbor-1];
        PS::S32 n_send_tot = n_send_disp[rank_tail] + n_send[rank_tail];
        PS::S32 n_recv_tot = n_recv_disp[rank_tail] + n_recv[rank_tail];
        n_send[rank_tail] = 0;
        ptcl_send.resize(n_send_tot);
        for(PS::S32 i=0; i<n_send_tot; i++){
            PS::S32 rank = rank_j[i];
            PS::S32 adr = n_send_disp[rank] + n_send[rank];
            ptcl_send[adr] = ptcl_buf[i];
            n_send[rank]++;
        }
        ptcl_recv.resize(n_recv_tot);
        n_proc_send = n_proc_recv = 0;
        for(PS::S32 i=0; i<n_proc_neighbor; i++){
            PS::S32 tag = 0;
            PS::S32 rank = rank_neighbor_[i];
            if(n_send[rank] > 0 ){
                MPI_Isend(&ptcl_send[n_send_disp[rank]], n_send[rank], PS::GetDataType<PTCLHardComm>(), rank, tag, MPI_COMM_WORLD, req_send+n_proc_send);
                n_proc_send++;
            }
            if(n_recv[rank] > 0 ){
                MPI_Irecv(&ptcl_recv[n_recv_disp[rank]], n_recv[rank], PS::GetDataType<PTCLHardComm>(), rank, tag, MPI_COMM_WORLD, req_recv+n_proc_recv);
                n_proc_recv++;
            }
        }
        MPI_Waitall(n_proc_send, req_send, stat_send);
        MPI_Waitall(n_proc_recv, req_recv, stat_recv);
#pragma omp parallel for
        for(PS::S32 i=0; i<n_recv_tot; i++){
            const PS::S64 adr = ptcl_recv[i].adr_fp; // M.I.
            system[adr].mass  = ptcl_recv[i].mass;
            system[adr].pos   = ptcl_recv[i].pos;
            system[adr].vel   = ptcl_recv[i].vel;
        }
#else //PARALLEL_2BODY
        // original
        for(size_t i=0; i<n_2body; i++){
            const PS::S64 adr = ptcl_2body_loc_[i].adr_fp; // M.I.
            assert(ptcl_2body_loc_[i].id == system[adr].id);
            system[adr].mass = ptcl_2body_loc_[i].mass;
            system[adr].pos = ptcl_2body_loc_[i].pos;
            system[adr].vel = ptcl_2body_loc_[i].vel;
        }
#endif //PARALLEL_2BODY

#ifdef CALC_HARD_ENERGY
        CalcTwoBodyEnergy(eng_2body_now_loc_, r_in, r_out);
        eng_2body_now_loc_.disp_merge = eng_hard_disp;
#endif
    }

#else //ORIGINAL
    // original
    template<class Tsoft>
    void evolveIsolatedSystem(Tsoft & system,
                              const PS::F64 r_out, 
			      const PS::F64 r_in,
                              const PS::F64 mass_sun, 
			      const PS::F64vec pos_sun,
                              const PS::F64vec vel_sun,
                              const PS::F64 eta_s, 
			      const PS::F64 eta,
                              const PS::F64 time_end, 
			      const PS::F64 dt_limit_org){
        const PS::S32 n_thread_max = PS::Comm::getNumberOfThread();
        static std::vector<MergeLog> * merge_log;
        static PS::F64 * eng_disp_omp;
        static bool first = true;
        if(first){
            merge_log = new std::vector<MergeLog>[n_thread_max];
            eng_disp_omp = new PS::F64[n_thread_max];
            first = false;
        }
        for(PS::S32 i=0; i<n_thread_max; i++){
            merge_log[i].clear();
            eng_disp_omp[i] = 0.0;
        }
        const size_t n_2body = ptcl_2body_loc_.size();
        ptcl_pred_2body_loc_.clear();
        ptcl_pred_2body_loc_.resize(n_2body);
        ptcl_force_2body_loc_.clear();
        ptcl_force_2body_loc_.resize(n_2body);
        merge_flag_2body_loc_.clear();
        merge_flag_2body_loc_.resize(n_2body);
#pragma omp parallel for
        for(size_t ip=0; ip<n_2body; ip++){
            ptcl_pred_2body_loc_[ip].pos = ptcl_2body_loc_[ip].pos;
            ptcl_pred_2body_loc_[ip].vel = ptcl_2body_loc_[ip].vel;
            ptcl_2body_loc_[ip].time = 0.0;
            ptcl_2body_loc_[ip].setRMerge();
            merge_flag_2body_loc_[ip] = false;
        }
        const PS::S32 n_pair = pair_adr_2body_loc_.size();
        const PS::F64 eps_sq = EPISoft::eps * EPISoft::eps;
#pragma omp parallel for
        for(PS::S32 np=0; np<n_pair; np++){
            //PS::S32 ith = PS::Comm::getThreadNum();
            PS::F64 time_sys_tmp = 0.0;
            PS::F64 time_sync_tmp = time_sys_tmp + dt_limit_org;
            PS::F64 time_end_tmp = time_end;

            const PS::S64 adr_i = pair_adr_2body_loc_[np].first;
            const PS::S64 adr_j = pair_adr_2body_loc_[np].second;
            if(adr_i > adr_j) continue;

            PTCLHard & ptcl_i = ptcl_2body_loc_[adr_i];
            PTCLHard & ptcl_j = ptcl_2body_loc_[adr_j];
            PTCLForce & force_i = ptcl_force_2body_loc_[adr_i];
            PTCLForce & force_j = ptcl_force_2body_loc_[adr_j];
            PTCLPred & pred_i = ptcl_pred_2body_loc_[adr_i];
            PTCLPred & pred_j = ptcl_pred_2body_loc_[adr_j];
            PS::F64 r2 = 0.0;
            force_i.clear();
            force_j.clear();
            CalcAcc0Acc1AndR2CutoffPair(ptcl_i.mass, pred_i.pos, pred_i.vel, force_i.acc0_pla, force_i.acc1_pla,
                                        ptcl_j.mass, pred_j.pos, pred_j.vel, force_j.acc0_pla, force_j.acc1_pla,
                                        r2,   eps_sq, r_out, r_in);

            CalcAcc0Acc1AndR2(pred_i.pos,    pred_i.vel,
                              force_i.acc0,  force_i.acc1,
                              pos_sun, vel_sun, mass_sun, 0.0);
	    
            CalcAcc0Acc1AndR2(pred_j.pos,     pred_j.vel,
                              force_j.acc0,   force_j.acc1,
                              pos_sun, vel_sun, mass_sun, 0.0);
	    force_i.reduce();
	    force_j.reduce();
            PS::F64 dt_limit_new = CalcDtLimit(time_sys_tmp, time_sync_tmp, dt_limit_org);
            ptcl_i.setDt2nd(force_i, eta_s, dt_limit_new, a0_offset_sq_);
            ptcl_j.setDt2nd(force_j, eta_s, dt_limit_new, a0_offset_sq_);
            PS::F64 dt = (ptcl_i.dt < ptcl_j.dt) ? ptcl_i.dt : ptcl_j.dt;
            ptcl_i.dt = ptcl_j.dt = dt;
            ptcl_i.acc0 = force_i.acc0;
            ptcl_i.acc1 = force_i.acc1;
            ptcl_j.acc0 = force_j.acc0;
            ptcl_j.acc1 = force_j.acc1;
            ptcl_i.acc0_pla = force_i.acc0_pla;
            ptcl_i.acc1_pla = force_i.acc1_pla;
            ptcl_j.acc0_pla = force_j.acc0_pla;
            ptcl_j.acc1_pla = force_j.acc1_pla;
            while(time_sys_tmp != time_end_tmp){
		n_loop_++;
		n_interaction_2body_loc_ += 2; // plus 2 means using force symmetry.
                predict2body(pred_i, pred_j,
                             ptcl_i, ptcl_j,
                             dt);
                force_i.clear();
                force_j.clear();
                CalcAcc0Acc1AndR2CutoffPair(ptcl_i.mass, pred_i.pos, pred_i.vel, force_i.acc0_pla, force_i.acc1_pla,
                                            ptcl_j.mass, pred_j.pos, pred_j.vel, force_j.acc0_pla, force_j.acc1_pla,
                                            r2,   eps_sq, r_out, r_in);
                CalcAcc0Acc1AndR2(pred_i.pos,	  pred_i.vel,
                                  force_i.acc0,   force_i.acc1,
                                  pos_sun, vel_sun, mass_sun, 0.0);
                CalcAcc0Acc1AndR2(pred_j.pos,     pred_j.vel,
                                  force_j.acc0,   force_j.acc1,
                                  pos_sun, vel_sun, mass_sun, 0.0);
		force_i.reduce();
		force_j.reduce();
#ifdef MERGE
                const PS::F64 r_merge = ptcl_i.r_merge + ptcl_j.r_merge;
                const PS::F64 r_merge_2 = r_merge * r_merge;
                bool merge = false;
                if(r2 < r_merge_2){
                    //pair_isolated_merge_adr_.push_back( std::pair<PS::S64, PS::S64>(adr_i, adr_j) );
                    merge = true;
                }
#endif //MERGE
                time_sys_tmp += dt;
                if(time_sys_tmp == time_sync_tmp) time_sync_ += dt_limit_org;
                dt_limit_new = CalcDtLimit(time_sys_tmp, time_sync_tmp, dt_limit_org);
                ptcl_i.correct(force_i, eta, dt_limit_new, a0_offset_sq_);
                ptcl_j.correct(force_j, eta, dt_limit_new, a0_offset_sq_);
                dt = (ptcl_2body_loc_[adr_i].dt < ptcl_2body_loc_[adr_j].dt) ? ptcl_2body_loc_[adr_i].dt : ptcl_2body_loc_[adr_j].dt;
                ptcl_i.dt = ptcl_j.dt = dt;
#ifdef MERGE
                if(merge){
		    PS::S32 ith = PS::Comm::getThreadNum();
                    //merge2body(ptcl_i, ptcl_j, merge_log[ith]);
                    merge2body(ptcl_i, ptcl_j, eng_disp_omp[ith], merge_log[ith]);
                    //merge2body(ptcl_2body_loc_, adr_i, adr_j, eng_disp_omp[ith], merge_log[ith]);
                    PTCLHard & ptcl_merge = (ptcl_i.mass == 0.0) ? ptcl_j : ptcl_i;
                    const PS::F64 dt_rem = time_end_tmp - time_sys_tmp;
                    DriveKepler(mass_sun, ptcl_merge.mass,
                                pos_sun, ptcl_merge.pos,
                                vel_sun, ptcl_merge.vel, dt_rem);
                    ptcl_merge.time = time_end_tmp;
                    time_sys_tmp = time_end_tmp;
                }
#endif //MERGE
            }
        }
        for(PS::S32 i=0; i<n_thread_max; i++){
            const size_t j_max = merge_log[i].size();
            for(size_t j=0; j<j_max; j++){
                merge_history_.push_back(merge_log[i][j]);
            }
            eng_disp_ += eng_disp_omp[i];
        }
        for(size_t i=0; i<n_2body; i++){
            //const PS::S64 adr = adr_ptcl_2body_loc_[i]; //M.I.
            const PS::S64 adr = ptcl_2body_loc_[i].adr_fp; //M.I.
            assert(ptcl_2body_loc_[i].id == system[adr].id);
            system[adr].mass = ptcl_2body_loc_[i].mass;
            system[adr].pos = ptcl_2body_loc_[i].pos;
            system[adr].vel = ptcl_2body_loc_[i].vel;
        }
    }
#endif



    void reduce_force(){
        const PS::S32 ni = n_active_;
#pragma omp parallel for
        for(PS::S32 ip=0; ip<ni; ip++){
            const PS::S32 adr = adr_sorted_[ip];
            ptcl_force_[adr].reduce();
	}
    }

    void evolveMultiSirial(const PS::F64 t_end,
			   const PS::F64 dt_limit_hard,
			   const PS::F64 r_out,
			   const PS::F64 r_in,
			   const PS::F64 eta_s,
			   const PS::F64 eta){

	makeDictionaryGlobal();
#ifdef CALC_HARD_ENERGY
	// NOTE: the function is call after setting of ptcl_multi_glb_[i].adr_pair.
	eng_hard_multi_disp_ = 0.0;
	CalcMultiBodyEnergy(eng_mbody_prev_loc_, r_in, r_out);
#endif
	setUpGlb(t_end, dt_limit_hard);
	calcAcc0AndAcc1(r_out, r_in); // BUGs are
	calcAcc0AndAcc1Kepler(mass_sun_, pos_sun_, vel_sun_);
	reduce_force();
	setInitDt(eta_s);
	copyForceToPtcl();
	sortAndSelectIp();
	bool end_hard_part = false;
	while( !end_hard_part ){
	    n_loop_++;
	    n_interaction_multi_loc_ += n_active_; 
	    end_hard_part = evolveKepler(r_out, r_in, eta, dt_limit_hard);
	}
#ifdef CALC_HARD_ENERGY
	CalcMultiBodyEnergy(eng_mbody_now_loc_, r_in, r_out);
	eng_mbody_now_loc_.disp_merge = eng_hard_multi_disp_;
#endif
    }

    PS::S64 n_loop_cum_;
    size_t n_ptcl_2body_loc_;
    size_t n_pair_2body_loc_;
    size_t n_ptcl_2body_glb_;
    size_t n_pair_2body_glb_;
    size_t n_ptcl_multi_;
    size_t n_pair_multi_;
    void accumulate_counter(){
	n_loop_cum_ += n_loop_;
	n_ptcl_2body_loc_ += ptcl_2body_loc_.size();
	n_pair_2body_loc_ += pair_id_2body_loc_.size();
	n_ptcl_multi_ += ptcl_multi_glb_.size();
	n_pair_multi_ += pair_id_multi_glb_.size();
    }

    void dump_counter(std::ofstream & fout, const PS::S64 n_interval){
	n_ptcl_2body_glb_ = PS::Comm::getSum(n_ptcl_2body_loc_);
	n_pair_2body_glb_ = PS::Comm::getSum(n_pair_2body_loc_);
	n_interaction_2body_glb_ = PS::Comm::getSum(n_interaction_2body_loc_);
	n_interaction_multi_glb_ = PS::Comm::getSum(n_interaction_multi_loc_);
	fout<<"n_loop_cum= "        <<n_loop_cum_      <<" "<< (PS::F64)n_loop_cum_ / n_interval
	    <<" n_ptcl_2body_loc_= "<<n_ptcl_2body_loc_<<" "<< (PS::F64)n_ptcl_2body_loc_ / n_interval
	    <<" n_pair_2body_loc_= "<<n_pair_2body_loc_<<" "<< (PS::F64)n_pair_2body_loc_ / n_interval
	    <<" n_ptcl_2body_glb_= "<<n_ptcl_2body_glb_<<" "<< (PS::F64)n_ptcl_2body_glb_ / n_interval
	    <<" n_pair_2body_glb_= "<<n_pair_2body_glb_<<" "<< (PS::F64)n_pair_2body_glb_ / n_interval
	    <<" n_ptcl_multi_= "    <<n_ptcl_multi_    <<" "<< (PS::F64)n_ptcl_multi_ / n_interval
	    <<" n_pair_multi_= "    <<n_pair_multi_    <<" "<< (PS::F64)n_pair_multi_ / n_interval
	    <<" n_interaction_2body_loc_= "    <<n_interaction_2body_loc_    <<" "<< (PS::F64)n_interaction_2body_loc_ / n_interval
	    <<" n_interaction_2body_glb_= "    <<n_interaction_2body_glb_    <<" "<< (PS::F64)n_interaction_2body_glb_ / n_interval
	    <<" n_interaction_multi_loc_= "    <<n_interaction_multi_loc_    <<" "<< (PS::F64)n_interaction_multi_loc_ / n_interval
	    <<" n_interaction_multi_glb_= "    <<n_interaction_multi_glb_    <<" "<< (PS::F64)n_interaction_multi_glb_ / n_interval<<std::endl;
    }

    void clear_counter(){
	n_loop_cum_ = 0;
	n_ptcl_2body_loc_ = n_pair_2body_loc_ 
	    = n_ptcl_2body_glb_ = n_pair_2body_glb_
	    = n_ptcl_multi_ =  n_pair_multi_ = 0;
	n_interaction_2body_loc_ = n_interaction_2body_glb_ = 0;
	n_interaction_multi_loc_ = n_interaction_multi_glb_ = 0;
    }

#if 0
    template<class Tptcl>
    void calcEnergyHard(const Tptcl ptcl[],
			const PS::S32 n,
			const PS::F64 r_out,
			const PS::F64 r_in,
			const PS::F64 m_sun,
			const PS::F64vec pos_sun,
			const PS::F64vec vel_sun){
	const PS::F64 q = r_in / r_out;
	for(PS::S32 i=0; i<n; i++){
	    PS::F64vec rij;
	    PS::F64 dr = sqrt(rij*rij + eps_sq);
	    y = dr / r_out;
	    PS::F64 w = CalcW(y, q);

	    pot_hard -= (ptcl[j].mass / dr)*w;

	    pot_hard -= (m_sun / dr); // from sun
	}
    }
#endif

};
