// updateDerivs.cpp
// TWS, December 2019
// Updated 2024 for LVRV model
/***************************************************
Definition of variables
	y[0]: a1 - overall elongation
	y[1]: a2 - overall contraction
	y[2]: a3 - torsion
LVRV model
	y[3]: a4 - RV expansion
	y[4]: a5 - septal expansion
LV model
	y[3]: a4 - RV diameter
	y[4] not used
Both models
	y[5]: a6 - LA diameter
	y[6]: a7 - RA diameter
	y[7]: q_pv - pulmonary vein flow
	y[8]: q_sv - vena cava flow
	y[9]: q_lv - aortic valve flow
	y[10]: q_la - mitral valve flow
	y[11]: q_rv - pulmonary valve flow
	y[12]: q_ra - tricuspid valve flow
	y[13]: Psp - systemic arterial pressure
	y[14]: Psv - systemic venous pressure
	y[15]: Ppp - pulmonary arterial pressure
	y[16]: Ppv - pulmonary venous presure
	y[17]: Zaov - aortic valve state
	y[18]: Zmiv - mitral valve state
	y[19]: Zpuv - pulmonary valve state
	y[20]: Ztcv - tricuspid valve state
****************************************************/
#include <math.h>
#include <stdio.h>
#include <string>
#include "nrutil.h"

using namespace std;

void gaussj(double **a, int n, double *b);

int updateDerivs(double t) {
	extern int Nmat, Nvars, Npar;
	extern double** J, ** J1;
	extern double* f, * X, * dy;
	extern double Plv, Psp, Psv, Pra, Prv, Ppp, Ppv, Pla, Psa, Ppa;
	extern double q_la, q_lv, q_ra, q_rv, q_LVAD;
	extern double q_sp, q_sv, q_pp, q_pv;
	extern double Csa, Csv, Cpa, Cpv;
	extern double** alpha_mat_aij, * eta_lv_ai, * eta_rv_ai, * kappa_ai;
	extern double eta_la, kappa_la, chi_la;
	extern double eta_rv, kappa_rv, chi_rv;
	extern double eta_ra, kappa_ra, chi_ra;
	extern double B_lv, B_la, B_rv, B_ra, L_lv, L_la, L_rv, L_ra, L_pv, L_sv;
	extern double Rsp, Rsv, Rpa, Rpp, Rpv, Rra, f_rv;
	extern double Vlv, Vrv, Vla, Vra;
	extern double VA_aovo, VA_aovc, VA_mivo, VA_mivc, VA_puvo, VA_puvc, VA_tcvo, VA_tcvc, VL_aov, VL_miv, VL_puv, VL_tcv, rho;
	extern double Zmiv, Zaov, Ztcv, Zpuv, Kaovo, Kaovc, Kmivo, Kmivc, Kpuvo, Kpuvc, Ktcvo, Ktcvc;
	extern double Z0start, Z0end, P_ao_0, q_lv_0, P_echo, Q_trans;
	extern double V0_peri, s_peri, Pperi, Tc;
	extern double q0_LVAD, Ps_LVAD, t0_LVAD;
	extern string odeSystem, septalDisp, LVAD, pericardium;

	int i, j, ip, k, Ndiff = Nmat - Npar;
	double linerror, linerror2;
	double err_la, err_lv, err_sp, err_sv, err_ra, err_rv, err_pp, err_pv;
	double VA_aov, VA_miv, VA_puv, VA_tcv;
	double decayrate = 100.;

	//system for LVRV model (Npar = 5, Ndiff = 2) and LV model (Npar = 3, Ndiff = 1) - Eq. (S94)
	//for LVRV model without septal displacement Npar = 4, Ndiff = 2
	for (i = 1; i <= Npar; i++) {
		for (j = 1; j <= Npar; j++) J[i][j] = alpha_mat_aij[i][j];
		f[i] = -kappa_ai[i];
		J[i][Npar + 1] = -eta_lv_ai[i];
		J[Npar + 1][i] = -eta_lv_ai[i];
		if (Ndiff == 2) {
			J[i][Npar + 2] = -eta_rv_ai[i];
			J[Npar + 2][i] = -eta_rv_ai[i];
		}
	}
	for (i = 1; i <= Ndiff; i++) for (j = 1; j <= Ndiff; j++) J[Npar + i][Npar + j] = 0.;
	f[Npar + 1] = -q_la + q_lv + q_LVAD;
	if (Ndiff == 2) f[Npar + 2] = -q_ra + q_rv;

	for (i = 1; i <= Nmat; i++) {
		X[i] = f[i];
		for (j = 1; j <= Nmat; j++) J1[i][j] = J[i][j];
	}
	gaussj(J1, Nmat, X);	//using Gaussian elimination to solve system of (S94,S125-128)
	for (ip = 1; ip <= Npar; ip++) dy[ip - 1] = X[ip]; // dai_dt, LVRV model to a5, LV model to a3

	if (pericardium == "on") Pperi = exp(s_peri * ((Vlv + Vrv + Vla + Vra)/V0_peri - 1.));	//(S91)
	else Pperi = 0.;

	Plv = Pperi + X[Npar + 1];
	if (Ndiff == 2) Prv = Pperi + X[Npar + 2];
	if (septalDisp == "off") dy[4] = 0.;

	if (Plv < -1.333)
		printf("*** Warning: Plv = %f kPa < -10 mmHg\n", Plv);	//pressures can be negative but should not go to -10 mmHg
	linerror2 = 0.;
	for (i = 1; i <= Nmat; i++) {	//check that result is correct
		linerror = f[i];
		for (k = 1; k <= Nmat; k++) linerror -= J[i][k] * X[k];
		linerror2 += DSQR(linerror);
	}

	if (linerror2 > 1.e-12 || linerror2 != linerror2) {
		printf("*** Error in linear solver: LV\n");
		return 1;
	}
	if(odeSystem == "LVmodel2024") {   //system for RV diameter
		dy[3] = (q_ra - q_rv) / (eta_rv * f_rv);	//eq. (S114), 2023 paper
		Prv = (kappa_rv + chi_rv * dy[3]) / eta_rv;	//eq. (S124), 2023 paper
		dy[4] = 0.;									//not used
	}
	//system for LA diameter
	dy[5] = (q_pv - q_la) / eta_la;				//Eq. (S127)
	Pla = Pperi + (kappa_la + chi_la * dy[5]) / eta_la;	//Eq. (S136)
	//system for RA diameter
	dy[6] = (q_sv - q_ra) / eta_ra;				//Eq. (S128)
	Pra = Pperi + (kappa_ra + chi_ra * dy[6]) / eta_ra;	//Eq. (S137)
	//q_pv - pulmonary vein flow, inertia
	dy[7] = (Ppv - Pla - Rpv * q_pv) / L_pv;	//Eq. (S133)	
	//q_sv - vena cava flow, inertia
	dy[8] = (Psv - Pra - Rsv * q_sv) / L_sv;	//Eq. (S134)
	//valve equations
	//Valve areas
	VA_aov = Zaov * VA_aovo + (1. - Zaov) * VA_aovc;	//Eq. (S120)
	VA_miv = Zmiv * VA_mivo + (1. - Zmiv) * VA_mivc;
	VA_puv = Zpuv * VA_puvo + (1. - Zpuv) * VA_puvc;
	VA_tcv = Ztcv * VA_tcvo + (1. - Ztcv) * VA_tcvc;
	// kinetic energy coefficients (convert dyn/cm2 to kPa)
	B_lv = 1.e-4 * rho / (2. * DSQR(VA_aov));			//Eq. (S119)
	B_la = 1.e-4 * rho / (2. * DSQR(VA_miv));
	B_rv = 1.e-4 * rho / (2. * DSQR(VA_puv));
	B_ra = 1.e-4 * rho / (2. * DSQR(VA_tcv));
	// inertia coefficients (convert dyn/cm2 to kPa)
	L_lv = 1.e-4 * rho * VL_aov / VA_aov;				//Eq. (S118)
	L_la = 1.e-4 * rho * VL_miv / VA_miv;
	L_rv = 1.e-4 * rho * VL_puv / VA_puv;
	L_ra = 1.e-4 * rho * VL_tcv / VA_tcv;
	//pressures downstream of valves
	//q_lv_0 includes LVAD contribution
	Psa = P_ao_0 + P_echo + Z0start * 1.e-4 * (q_lv + q_LVAD - q_lv_0);	//subtract DC component. Eq. (S115)
	Ppa = Ppp + Rpa * q_rv;																//Eq. (S138)
	//flow derivatives
	// d / dt[q_lv] - aortic valve flow
	if (VA_aov > 1e-3) dy[9] = 1. / L_lv * (Plv - Psa - B_lv * fabs(q_lv) * q_lv);		//Eq. (S116)
	else dy[9] = -decayrate * q_lv;
	// d / dt[q_la] - mitral valve flow
	if (VA_miv > 1e-3) dy[10] = 1. / L_la * (Pla - Plv - B_la * fabs(q_la) * q_la);		//Eq. (S116)
	else dy[10] = -decayrate * q_la;
	// d / dt[q_rv] - pulmonary valve flow
	if (VA_puv > 1e-3) dy[11] = 1. / L_rv * (Prv - Ppa - B_rv * fabs(q_rv) * q_rv);		//Eq. (S117)
	else dy[11] = -decayrate * q_rv;
	// d / dt[q_ra] - tricuspid valve flow
	if (VA_tcv > 1e-3) dy[12] = 1. / L_ra * (Pra - Prv - B_ra * fabs(q_ra) * q_ra);		//Eq. (S117)
	else dy[12] = -decayrate * q_ra;
	if (LVAD == "on" && t > t0_LVAD * Tc) {		// turn on LVAD after t0_LVAD cycles
		if (Psa - Plv < Ps_LVAD * 0.1333) q_LVAD = q0_LVAD * sqrt(1. + (Plv - Psa) / (Ps_LVAD * 0.1333));
	}
	else q_LVAD = 0.;
	//compliance pressures
	// d / dt[Psp] - systemic arterial pressure
	dy[13] = 1. / Csa * (Q_trans - (Psp - Psv) / Rsp);	//Eq. (S129)
	// d / dt[Psv] - systemic venous pressure
	dy[14] = 1. / Csv * ((Psp - Psv) / Rsp - q_sv);							//Eq. (S130)
	// d / dt[Ppp] - pulmonary arterial pressure
	dy[15] = 1. / Cpa * (q_rv - (Ppp - Ppv) / Rpp);							//Eq. (S131)
	// d / dt[Ppv] - pulmonary venous presure
	dy[16] = 1. / Cpv * ((Ppp - Ppv) / Rpp - q_pv);							//Eq. (S132)
	q_sp = (Psp - Psv) / Rsp;
	q_pp = (Ppp - Ppv) / Rpp;

	//derivatives of valve states
	//Zaov - aortic valve state
	if (Plv > Psa) dy[17] = (1. - Zaov) * Kaovo * (Plv - Psa);				//Eq. (S121)
	else dy[17] = Zaov * Kaovc * (Plv - Psa);
	//Zmiv - mitral valve state
	if (Pla > Plv) dy[18] = (1. - Zmiv) * Kmivo * (Pla - Plv);				//Eq. (S122)
	else dy[18] = Zmiv * Kmivc * (Pla - Plv);
	//Zpuv - pulmonary valve state
	if (Prv > Ppa) dy[19] = (1. - Zpuv) * Kpuvo * (Prv - Ppa);				//Eq. (S123)
	else dy[19] = Zpuv * Kpuvc * (Prv - Ppa);
	//Ztcv - tricuspid valve state
	if (Pra > Prv) dy[20] = (1. - Ztcv) * Ktcvo * (Pra - Prv);				//Eq. (S124)
	else dy[20] = Ztcv * Ktcvc * (Pra - Prv);

	//useful for debugging
	err_lv = q_la - q_lv;
	for (i = 1; i <= Npar; i++) err_lv -= eta_lv_ai[i] * dy[i - 1];
	err_rv = q_ra - q_rv;
	if(odeSystem == "LVmodel2024") err_rv -= f_rv * eta_rv * dy[3];
	else for (i = 1; i <= Npar; i++) err_rv -= eta_rv_ai[i] * dy[i - 1];
	err_la = q_pv - q_la - eta_la * dy[5];
	err_ra = q_sv - q_ra - eta_ra * dy[6];
	err_sp = q_lv - q_sp - Csa * dy[13];
	err_sv = q_sp - q_sv - Csv * dy[14];
	err_pp = q_rv - q_pp - Cpa * dy[15];
	err_pv = q_pp - q_pv - Cpv * dy[16];

	return 0;
}
