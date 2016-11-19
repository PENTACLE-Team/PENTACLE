inline PS::F64 cutoff_poly_3rd(const PS::F64 rij,
                               const PS::F64 rout,
                               const PS::F64 rin){
    PS::F64 inv_dr = 1.0 / (rout-rin);
    PS::F64 x = (rij - rin)*inv_dr;
    x = (x < 1.0) ? x : 1.0;
    x = (x > 0.0) ? x : 0.0;
    PS::F64 x2 = x*x;
    PS::F64 x4 = x2*x2;
    PS::F64 k = (((-20.0*x+70.0)*x-84.0)*x+35.0)*x4;
    return k;
}


inline PS::F64 cutoff_poly_3rd(const PS::F64 rij,
                               const PS::F64 rout,
                               const PS::F64 rin,
                               const PS::F64 inv_dr){
    PS::F64 x = (rij - rin)*inv_dr;
    x = (x < 1.0) ? x : 1.0;
    x = (x > 0.0) ? x : 0.0;
    PS::F64 x2 = x*x;
    PS::F64 x4 = x2*x2;
    PS::F64 k = (((-20.0*x+70.0)*x-84.0)*x+35.0)*x4;
    return k;
}

/*
inline void calc_acc_split(const PS::F64vec & posi,
                           PS::F64vec & acci_long,
                           PS::F64vec & acci_short,
                           PS::F64 & poti_tot,
                           const PS::F64vec & posj,
                           const PS::F64 massj,
                           const PS::F64 eps2,
                           const PS::F64 rcut_out,
                           const PS::F64 rcut_in,
                           const PS::F64 rsearch2){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2 = rij * rij;
    const PS::F64 r2_eps = r2 + eps2;
    //const PS::F64 r2_eps = rij*rij + eps2;
    const PS::F64 inv_dr = 1.0 / (rcut_out-rcut_in);
    //if(r2_eps <= rsearch2){
    if(r2 <= rsearch2){ // to be consisten with CalcForceEPEP()
        PS::F64 r_eps = sqrt(r2_eps);
        PS::F64 R = 1.0/r_eps;
        PS::F64 R2 = R*R;
        PS::F64 R3 = R2*R;
#ifdef FORDEBUG
        PS::F64 K = 0.0; // for debug
#else
        //PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in, inv_dr);
	PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
#endif
        PS::F64vec F0 = -massj * R3 * rij;
        acci_short += F0 * (1.0-K);
        acci_long += F0 * K;
        poti_tot -= massj * R;
    }
}
*/

/*
inline void CalcAcc0Short(const PS::F64vec & posi,
			  PS::F64vec & acci_short,
			  PS::F64 & poti_tot,
			  const PS::F64vec & posj,
			  const PS::F64 massj,
			  const PS::F64 eps2,
			  const PS::F64 rcut_out,
			  const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2 = rij * rij;
    const PS::F64 r2_eps = r2 + eps2;
    const PS::F64 inv_dr = 1.0 / (rcut_out-rcut_in);
    PS::F64 r_eps = sqrt(r2_eps);
    PS::F64 R = 1.0/r_eps;
    PS::F64 R2 = R*R;
    PS::F64 R3 = R2*R;
    //PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in, inv_dr);
    PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
    PS::F64vec F0 = -massj * R3 * rij;
    acci_short += F0 * (1.0-K);
    poti_tot -= massj * R;
}
*/

inline void CalcPotShort(const PS::F64vec & posi,
			 PS::F64 & poti_tot,
			 const PS::F64vec & posj,
			 const PS::F64 massj,
			 const PS::F64 eps2,
			 const PS::F64 rcut_out,
			 const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2 = rij * rij;
    const PS::F64 r2_eps = r2 + eps2;
    //const PS::F64 inv_dr = 1.0 / (rcut_out-rcut_in);
    PS::F64 R = 1.0/sqrt(r2_eps);
    PS::F64 r_eps = R * r2_eps;
    //PS::F64 k = cutoff_poly_3rd(r_eps, rcut_out, rcut_in, inv_dr);
    PS::F64 k = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
    poti_tot -= massj * R * (1.0-k);
}

inline void CalcAccPotShortWithLinearCutoff(const PS::F64vec & posi,
					    PS::F64vec & acci_pla,
					    PS::F64 & poti_tot,
					    const PS::F64vec & posj,
					    const PS::F64 massj,
					    const PS::F64 eps2,
					    const PS::F64 rcut_out,
					    const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2 = rij * rij;
    const PS::F64 r2_eps = r2 + eps2;
    const PS::F64 rcut2_out = rcut_out * rcut_out;
    const PS::F64 R = 1.0/sqrt(r2_eps);
    const PS::F64 Rm = massj * R;
    const PS::F64 R2 = R * R;
    const PS::F64 Rm3 = Rm * R2;
    const PS::F64 r_eps = R * r2_eps;
    const PS::F64 k = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
    const PS::F64 r2_max = (r2_eps > rcut2_out) ? r2_eps : rcut2_out;
    const PS::F64 R_max = 1.0/sqrt(r2_max);
    const PS::F64 Rm_max = massj * R_max;
    const PS::F64 R2_max = R_max * R_max;
    const PS::F64 Rm3_max = Rm_max * R2_max;
    poti_tot -= (Rm - Rm_max);
    acci_pla -= (Rm3*k - Rm3_max)*rij;
}


inline PS::F64 cutoff_poly_3rd_dot(const PS::F64 &rij,
                                   const PS::F64 &rijvij,
                                   const PS::F64 &_rout,
                                   const PS::F64 &_rin){
    PS::F64 rout = _rout;
    PS::F64 rin = _rin;
    PS::F64 inv_dr = 1.0/(rout-rin);
    PS::F64 x = (rij - rin)*inv_dr;
    PS::F64 xdot = rijvij/rij*inv_dr;
    PS::F64 Kdot = 0.0;
    if(x <= 0.0)
        Kdot = 0.0;
    else if(1.0 <= x)
        Kdot = 0.0;
    else{
        PS::F64 x2 = x*x;
        PS::F64 x3 = x2*x;
        PS::F64 x4 = x2*x2;
        PS::F64 x5 = x4*x;
        PS::F64 x6 = x4*x2;
        Kdot = (-140.0*x6 + 420.0*x5 - 420.0*x4 + 140.0*x3) * xdot;
    }
    return Kdot;
}

#if 1
inline void CalcAcc0AndAcc1Cutoff(const PS::F64vec posi,
                                  const PS::F64vec veli,
                                  PS::F64vec & acci,
                                  PS::F64vec & jrki,
                                  const PS::F64vec posj, 
                                  const PS::F64vec velj, 
                                  const PS::F64 massj, 
                                  const PS::F64 eps2, 
                                  const PS::F64 rcut_out,
                                  const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2_eps = rij*rij + eps2;
    if(r2_eps <= rcut_out*rcut_out){
        const PS::F64vec vij = veli - velj;
        const PS::F64 rijvij = rij * vij;
        const PS::F64 r_eps = sqrt(r2_eps);
        const PS::F64 R = 1.0/r_eps;
	//const PS::F64 R = 1.0 / sqrt(r2_eps);
        const PS::F64 R2 = R*R;
        const PS::F64 R3 = R2*R;
        const PS::F64 A = (rijvij)*R2;
#ifdef FORDEBUG
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -massj*R3*rij*(1.0-K);
        const PS::F64vec F1 = -massj*R3*vij*(1.0-K) - 3.0*A*F0 + massj*R3*rij*Kdot;
        acci += F0;
        jrki += F1;
    }
}
#endif

inline void CalcAcc0Acc1AndR2Cutoff(const PS::F64vec posi,
				    const PS::F64vec veli,
				    PS::F64vec & acci,
				    PS::F64vec & jrki,
				    PS::F64 & r2,
				    const PS::F64vec posj, 
				    const PS::F64vec velj,
				    const PS::F64 massj,
				    const PS::F64 eps2,
				    const PS::F64 rcut_out,
				    const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    r2 = rij*rij;
    const PS::F64 r2_eps = r2 + eps2;
    if(r2_eps <= rcut_out*rcut_out){
        const PS::F64vec vij = veli - velj;
        const PS::F64 rijvij = rij * vij;
        const PS::F64 r_eps = sqrt(r2_eps);
        const PS::F64 R = 1.0/r_eps;
	//const PS::F64 R = 1.0 / sqrt(r2_eps);
        const PS::F64 R2 = R*R;
        const PS::F64 R3 = R2*R;
        const PS::F64 A = (rijvij)*R2;
#ifdef FORDEBUG
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -massj*R3*rij*(1.0-K);
        const PS::F64vec F1 = -massj*R3*vij*(1.0-K) - 3.0*A*F0 + massj*R3*rij*Kdot;
        acci += F0;
        jrki += F1;
    }
}

inline void CalcAcc0Acc1AndR2CutoffPair(const PS::F64 massi,
					const PS::F64vec posi,
					const PS::F64vec veli,
					PS::F64vec & acci,
					PS::F64vec & jrki,
					const PS::F64 massj,
					const PS::F64vec posj,
					const PS::F64vec velj,
					PS::F64vec & accj,
					PS::F64vec & jrkj,
					PS::F64 & r2,
					const PS::F64 eps2,
					const PS::F64 rcut_out,
					const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    r2 = rij*rij;
    const PS::F64 r2_eps = r2 + eps2;
    if(r2_eps <= rcut_out*rcut_out){
        const PS::F64vec vij = veli - velj;
        const PS::F64 rijvij = rij * vij;
        const PS::F64 r_eps = sqrt(r2_eps);
        const PS::F64 R = 1.0/r_eps;
	//const PS::F64 R = 1.0 / sqrt(r2_eps);
        const PS::F64 R2 = R*R;
        const PS::F64 R3 = R2*R;
        const PS::F64 A = (rijvij)*R2;
	//#ifdef FORDEBUG
#if 0
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -R3*rij*(1.0-K);
        const PS::F64vec F1 = -R3*vij*(1.0-K) - 3.0*A*F0 + R3*rij*Kdot;
        acci += massj*F0;
        jrki += massj*F1;
	accj -= massi*F0;
        jrkj -= massi*F1;
    }
}

inline void CalcAcc0Acc1AndPotCutoff(const PS::F64vec posi,
                                  const PS::F64vec veli,
                                  PS::F64vec & acci,
                                  PS::F64vec & jrki,
                                  PS::F64 & poti,
                                  const PS::F64vec posj, 
                                  const PS::F64vec velj, 
                                  const PS::F64 massj,
                                  const PS::F64 eps2,
                                  const PS::F64 rcut_out,
                                  const PS::F64 rcut_in){
    const PS::F64vec rij = posi - posj;
    const PS::F64 r2_eps = rij*rij + eps2;
#ifdef FORDEBUG
    if(1){
#else
    if(r2_eps <= rcut_out*rcut_out){
#endif
        const PS::F64vec vij = veli - velj;
        const PS::F64 rijvij = rij * vij;
        const PS::F64 r_eps = sqrt(r2_eps);
        const PS::F64 R = 1.0/r_eps;
	//const PS::F64 R = 1.0 / sqrt(r2_eps);
        const PS::F64 R2 = R * R;
        const PS::F64 R3 = R2 * R;
        const PS::F64 A = rijvij * R2;
#ifdef FORDEBUG
        PS::F64 K = 0.0; // for debug
        PS::F64 Kdot = 0.0; // for debug
        poti -= massj * R;
#else
        const PS::F64 K = cutoff_poly_3rd(r_eps, rcut_out, rcut_in);
        const PS::F64 Kdot = cutoff_poly_3rd_dot(r_eps, rijvij, rcut_out, rcut_in);
#endif
        const PS::F64vec F0 = -massj*R3*rij*(1.0-K);
        const PS::F64vec F1 = -massj*R3*vij*(1.0-K) - 3.0*A*F0 + massj*R3*rij*Kdot;
        acci += F0;
        jrki += F1;
    }
}

inline void CalcAcc0Acc1AndR2(const PS::F64vec posi,
			      const PS::F64vec veli,
			      PS::F64vec & acci,
			      PS::F64vec & jrki,
			      const PS::F64vec posj,
			      const PS::F64vec velj,
			      const PS::F64 massj,
			      const PS::F64 eps2){
    const PS::F64vec rij = posi - posj;
    const PS::F64vec vij = veli - velj;
    //const PS::F64 r2_eps = rij*rij + eps2;
    const PS::F64 r_inv = 1.0 / sqrt(rij * rij);
    const PS::F64 r2_inv = r_inv * r_inv;
    const PS::F64 r3_inv = r2_inv * r_inv;
    const PS::F64 m_r3 = massj * r3_inv;
    const PS::F64vec F0 = -m_r3*rij;
    const PS::F64vec F1 = -m_r3*vij - 3.0*rij*vij*r2_inv*F0;
    acci += F0;
    jrki += F1;
}

/* 
class CalcW{
private:
    PS::F64 A7, A6, A5, A4, A3, A2, A1, A0, B1, A1_dash;
    PS::F64 q; // rin/rout
public:
    void setParam(const PS::F64 r_in, const PS::F64 r_out){
    }
};
*/

//#ifdef CALC_HARD_ENERGY
// y: reps/rout [reps: sqrt(r^2+eps^2)], q: rin/rout
inline PS::F64 CalcW(const PS::F64 y, const PS::F64 q=0.1){
     PS::F64 q2 = q*q;
     PS::F64 q3 = q2*q;
     PS::F64 q4 = q2*q2;
     PS::F64 q5 = q3*q2;
     PS::F64 q6 = q3*q3;
     PS::F64 q7 = q4*q3;
     PS::F64 denominator = (q-1.0)*(q-1.0)*(q-1.0)*(q-1.0)*(q-1.0)*(q-1.0)*(q-1.0);
     PS::F64 A7 = 20.0/denominator/-6;
     PS::F64 A6 = (-70.0*q - 70.0)/denominator/-5;
     PS::F64 A5 = (84.0*q2 + 252.0*q + 84.0)/denominator/-4;
     PS::F64 A4 = (-35.0*q3 - 315.0*q2 - 315.0*q - 35.0)/denominator/-3;
     PS::F64 A3 = (140.0*q3 + 420.0*q2 + 140.0*q)/denominator/-2;
     PS::F64 A2 = (-210*q3 - 210.0*q2)/denominator/-1;
     PS::F64 A1 = (140*q3)/denominator*-1;
     PS::F64 A0 = (-35.0*q4 + 21.0*q5 - 7.0*q6 + q7)/denominator;
     PS::F64 x = 1.0; // x=rout/rout
     PS::F64 B1 = 1.0 - ( (((((((A7*x + A6)*x + A5)*x + A4)*x + A3)*x + A2)*x + A1*log(x))*x) + A0 ); // to W(r>rout) = 1.0
     PS::F64 A1_dash = -7*(60*q3*log(q) - q6 + 9.0*q5 - 45.0*q4 + 45.0*q2 - 9.0*q + 1.0)/(3.0*denominator);
     if(y <= q) return A1_dash*y;
     else if(y >= 1.0) return 1.0;
     else return (((((((A7*y + A6)*y + A5)*y + A4)*y + A3)*y + A2)*y + A1*log(y) + B1)*y) + A0;
}
//#endif

#if 0
inline PS::F64 calc_G(const PS::F64 y, const PS::F64 gamma){
    const PS::F64 y2 = y * y;
    const PS::F64 y3 = y2 * y;
    const PS::F64 y4 = y2 * y2;
    const PS::F64 y5 = y3 * y2;
    const PS::F64 y6 = y3 * y3;
    const PS::F64 y7 = y3 * y3;
    const PS::F64 gamma2 = gamma * gamma;
    const PS::F64 gamma3 = gamma2 * gamma;
    const PS::F64 c = 1.0 - gamma;
    const PS::F64 c7 = 1.0 / (c*c*c*c*c*c*c);
    PS::F64 G =
	(-10.0/3.0*y7
	 + 14.0*(gamma+1.0)*y6
	 - 21.0*(gamma2+3.0*gamma+1.0)*y5
	 + 35.0/3.0*(gamma3+9.0*gamma2+9.0*gamma+1.0)*y4
	 - 70.0*(gamma3+3.0*gamma2+gamma)*y3
	 + 210.0*(gamma3+gamma2)*y2
	 - 140.0*gamma3*y*log(y)
	 + (gamma7-7.0*gamma6+21.0*gamma5-35.0*gamma4)) * c7;
    return G;
}

inline PS::F64 calc_W(const PS::F64 y, const PS::F64 gamma){
    if(1.0 < y){
	return 1.0;
    }
    else if(gamma < y){
	calc_G
    }
    else{
	PS::F64 ret = 7.0*(gamma6 - 9.0*gamma5 + 45.0*gamma4)
    }
    const PS::F64 y2 = y * y;
    const PS::F64 y3 = y2 * y;
    const PS::F64 y4 = y2 * y2;
    const PS::F64 y5 = y3 * y2;
    const PS::F64 y6 = y3 * y3;
    const PS::F64 y7 = y3 * y3;
    const PS::F64 gamma2 = gamma * gamma;
    const PS::F64 gamma3 = gamma2 * gamma;
    const PS::F64 c = 1.0 - gamma;
    const PS::F64 c7 = 1.0 / (c*c*c*c*c*c*c);
    PS::F64 G =
	(-10.0/3.0*y7
	 + 14.0*(gamma+1.0)*y6
	 - 21.0*(gamma2+3.0*gamma+1.0)*y5
	 + 35.0/3.0*(gamma3+9.0*gamma2+9.0*gamma+1.0)*y4
	 - 70.0*(gamma3+3.0*gamma2+gamma)*y3
	 + 210.0*(gamma3+gamma2)*y2
	 - 140.0*gamma3*y*log(y)
	 + (gamma7-7.0*gamma6+21.0*gamma5-35.0*gamma4)) * c7;
    return G;
}


 
#endif
