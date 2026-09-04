// odeSolver.cpp
// edited by TWS, November 2019
#define _CRT_SECURE_NO_WARNINGS	//TWS
#include <stdio.h>
#include "nrutil.h"
#include <math.h>
#include <string>
#include <sstream>
#include <chrono>

using namespace std;

int updateModel(double t);

int ODEsolver() {
	extern int Nvars, Nmat, Nr, Ncycles, Nt, Nta, iRecord, Nmu, Nmu1, Nnu, Nstore, Npar, Nphi, Ncycles_show;
	extern double Tc, dtRecord, dtMax;
	extern double Vwall_rv, Vwall_lv, Vwall_se, Vlv, Vrv, Vla, Vra;
	extern double Csa, Csv, Cpa, Cpv;
	extern double Plv, Pla, Ppa, Psa, Psp, Psv, Pra, Ppv, Ppp, Prv, Pperi;
	extern double Rsp, Rsv, Rra, Rpa, Rpp, Rpv, Zmiv, Zaov, Ztcv, Zpuv;
	extern double a0, a1, a2, a3, a4, a5, a6, a7, At, At_la;
	extern double q_ra, q_rv, q_la, q_lv, q_LVAD;
	extern double q_sp, q_sv, q_pp, q_pv;
	extern double *eta_lv_ai, *eta_rv_ai, nu_up;
	extern double eta_la, eta_rv, f_rv, eta_ra;
	extern double VA_aovo, VA_aovc, VA_mivo, VA_mivc, VA_puvo, VA_puvc, VA_tcvo, VA_tcvc;
	extern double Pext, Psp0, Psv0, Ppp0, Ppv0;
	extern double LV_SW, RV_SW;

	extern double *ytest, *dy, *r;
	extern double *y, *y_init, *tvec, *a1vec, *a2vec, *a3vec, *a4vec, *a5vec, *a6vec, *a7vec;
	extern double **store, **Y, ***chmu, *cnu;
	extern double *echo, *trans, *q_lv_hist, *p_ao_hist, P_ao_0, q_lv_0, P_echo, Q_trans, Vaorta, Z0end;
	extern string odeSystem;
	extern double ***ell, ***ell0, ***shmu;
	extern double ***ch_store, ***sh_store, ***sh_store2;

	int cycleNum = 1;
	int i, k, ip, error;
	int kk = 0, j, Nphi2 = Nphi / 2;
	int Nmu0 = Nmu / 2 - 1, Nnu0 = Nnu / 2 - 1;
	int Nnu_mid = (PI / 2. - nu_up) / (PI - nu_up) * (Nnu - 1) + 0.5;		//find index on equator

	double tRecord, t = 0., tstart = Tc * (Ncycles - 1), tstart0 = tstart;
	double Amin = 0.25;	//echo doppler sampling area (cm2) in calculation of blood velocity
	double a, Vtot = 0., dVtot = 0., dVaorta;
	double SV_lv, EF_lv, SV_rv, EF_rv, dVlv, dVrv;
	FILE *ofp;
	char fname[80];
	
	//initial values
	//y = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 10., 0.55, 2., 1., 0., 1., 0., 1.]
	for (k = 0; k < Nvars; k++) y[k] = 0.;
	y[13] = Psp0;
	y[14] = Psv0;
	y[15] = Ppp0;
	y[16] = Ppv0;
	y[18] = 1.;
	y[20] = 1.;

	P_ao_0 = 0.;
	q_lv_0 = 0.;
	P_echo = 0.;

	LV_SW = 0.;
	RV_SW = 0.;
	   
	for (iRecord = 0; iRecord < Nt; iRecord++) {
		for (i = 0; i < Nstore; i++) store[i][iRecord] = 0.;
		tvec[iRecord] = 0.;
		a1vec[iRecord] = 0.;
		a2vec[iRecord] = 0.;
		a3vec[iRecord] = 0.;
		a4vec[iRecord] = 0.;
		a5vec[iRecord] = 0.;
		a6vec[iRecord] = 0.;
		a7vec[iRecord] = 0.;
	}
	tRecord = (Ncycles - Ncycles_show) * Tc;
	iRecord = 0;

	for (i = 0; i < Nta; i++) {			// don't store history from previous run
		q_lv_hist[i] = 0.;
		p_ao_hist[i] = 0.;
	}
	//////////////////// start of time integration loop //////////////////
	while (t < Tc * Ncycles) {

		error = updateModel(t);

		Vtot = Vla + Vlv + Vra + Vrv + Csa * Psp + Csv * Psv + Cpa * Ppp + Cpv * Ppv;
		dVtot = eta_la * dy[6] + eta_ra * dy[7] + Csa * dy[13] + Csv * dy[14] + Cpa * dy[15] + Cpv * dy[16];
		if (odeSystem == "LVmodel2024") dVtot += f_rv * eta_rv * dy[3];
		else for (ip = 1; ip <= Npar; ip++) dVtot += (eta_lv_ai[ip] + eta_rv_ai[ip]) * dy[ip - 1];
				 
		//update q_lv history vector, calculate reflected pressure
		for (i = Nta - 1; i > 0; i--) {
			q_lv_hist[i] = q_lv_hist[i - 1];	//shift contents by one position
			p_ao_hist[i] = p_ao_hist[i - 1];
		}
		//q_lv_hist includes LVAD contribution
		q_lv_hist[0] = q_lv + q_LVAD;
		p_ao_hist[0] = Psp;
		P_echo = 0.;
		Q_trans = 0.;
		for (i = 0; i < Nta; i++) {
			P_echo += echo[i] * q_lv_hist[i] * dtMax;	//convolution, eq. (S113)
			Q_trans += trans[i] * q_lv_hist[i] * dtMax;	//convolution, eq. (S114)
		}
		P_ao_0 += (p_ao_hist[0] - p_ao_hist[Nta - 1]) / Nta;	//update pressure average over one cycle at downstream end of aorta
		//q_lv_0 includes LVAD contribution
		q_lv_0 += (q_lv_hist[0] - q_lv_hist[Nta - 1]) / Nta;	//update LV outflow average over one cycle
		//dVaorta = q_lv - P_trans * 1.e4 / Z0end;				// d/dt aorta volume
		dVaorta = q_lv + q_LVAD - Q_trans;				// d/dt aorta volume
		Vaorta += dVaorta * dtMax;				// update aorta volume
		
		Vtot += Vaorta;
		dVtot += dVaorta;
		//if (t > 0.1 && fabs(dVtot) > 0.01) printf("*** Error: Volume not conserved %f\n", dVtot);
		//with aorta model, volume conservation is not exact. Check graph.
		// Store results, include only the last Ncycles_show
		if (t >= tRecord) {			// Store results
			if (iRecord > Nt) {
				printf("***Error: storage arrays not big enough ****\n");
				return 1;
			}

			for (k = 0; k < Nvars; k++)	Y[k][iRecord] = y[k];
			tvec[iRecord] = t;
			if (a1 != a1) return 0;
			a = a0 + a1;
			if (abs(q_lv) < 1e-6) q_lv = 0.;		//avoid large negative exponents
			if (abs(q_rv) < 1e-6) q_rv = 0.;
			if (abs(q_la) < 1e-6) q_la = 0.;
			if (abs(q_ra) < 1e-6) q_ra = 0.;
			a1vec[iRecord] = a1;
			a2vec[iRecord] = a2;
			a3vec[iRecord] = a3;
			a4vec[iRecord] = a4;
			a5vec[iRecord] = a5;
			a6vec[iRecord] = a6;
			a7vec[iRecord] = a7;
			store[0][iRecord] = Vtot;
			store[1][iRecord] = Vlv;
			store[2][iRecord] = Plv + Pext;
			store[3][iRecord] = Vla;
			store[4][iRecord] = Pla + Pext;
			store[5][iRecord] = Vrv;
			store[6][iRecord] = Prv + Pext;
			store[7][iRecord] = Vra;
			store[8][iRecord] = Pra + Pext;
			store[9][iRecord] = Pperi;
			store[10][iRecord] = Psa + Pext;
			store[11][iRecord] = Psp + Pext;
			store[12][iRecord] = Psv + Pext;
			store[13][iRecord] = Ppa + Pext;
			store[14][iRecord] = Ppp + Pext;
			store[15][iRecord] = Ppv + Pext;
			store[16][iRecord] = q_la;
			store[17][iRecord] = q_lv;
			store[18][iRecord] = q_ra;
			store[19][iRecord] = q_rv;
			store[20][iRecord] = q_LVAD;
			store[21][iRecord] = q_sp;
			store[22][iRecord] = q_sv;
			store[23][iRecord] = 0.;
			store[24][iRecord] = q_pp;
			store[25][iRecord] = q_pv;
			store[26][iRecord] = Zmiv;
			store[27][iRecord] = Zaov;
			store[28][iRecord] = At;
			store[29][iRecord] = At_la;
			store[33][iRecord] = Ztcv;
			store[34][iRecord] = Zpuv;
			//base to apex length and time derivative
			store[35][iRecord] = a * (chmu[0][0][0] * cnu[0] + chmu[0][Nnu - 1][0]);
			if (iRecord == 0) store[36][iRecord] = 0.;
			else store[36][iRecord] = -(store[35][iRecord] - store[35][iRecord - 1]) / (tvec[iRecord] - tvec[iRecord - 1]);
			store[37][iRecord] = P_echo + Pext;
			store[38][iRecord] = Q_trans;
			store[39][iRecord] = P_ao_0 + Pext;
			store[40][iRecord] = q_lv_0;

			//blood flow velocities in valves
				//assume fully open if q > 0, fully closed if q < 0
				//use of actual valve areas causes artifacts due to time lag of flow changes
				//set minimum value of area to Amin
			if (q_lv > 0.) store[41][iRecord] = q_lv / DMAX(Amin, VA_aovo);
			else store[41][iRecord] = q_lv / DMAX(Amin, VA_aovc);
			if (q_la > 0.) store[42][iRecord] = q_la / DMAX(Amin, VA_mivo);
			else store[42][iRecord] = q_la / DMAX(Amin, VA_mivc);
			if (q_rv > 0.) store[43][iRecord] = q_rv / DMAX(Amin, VA_puvo);
			else store[43][iRecord] = q_rv / DMAX(Amin, VA_puvc);
			if (q_ra > 0.) store[44][iRecord] = q_ra / DMAX(Amin, VA_tcvo);
			else store[44][iRecord] = q_ra / DMAX(Amin, VA_tcvc);
			//internal diameter of LV at equator
			store[45][iRecord] = 2. * a * shmu[0][Nnu_mid][0];  //corrected August 2025
			//external diameter of LV at equator
			store[46][iRecord] = 2. * a * shmu[0][Nnu_mid][Nmu-1];  //corrected August 2025
			store[47][iRecord] = Vaorta;
			tRecord += dtRecord;		// Reset record variables
			iRecord++;
		}
		//integration using RK2 scheme //////////////////
		for (k = 0; k < Nvars; k++) {
			ytest[k] = y[k];
			y[k] += dtMax * dy[k] / 2.;
		}
		for (k = Nvars - 4; k < Nvars; k++) {	//avoid overshoot of valve states
			if (y[k] > 1.) y[k] = 1.;
			if (y[k] < 0.) y[k] = 0.;
		}
		error = updateModel(t + dtMax / 2.);
		for (k = 0; k < Nvars; k++) y[k] = ytest[k] + dtMax * dy[k];

		///////////////////////////////////////////
		for (k = Nvars - 4; k < Nvars; k++) {	//avoid overshoot of valve states
			if (y[k] > 1.) y[k] = 1.;
			if (y[k] < 0.) y[k] = 0.;
		}
		t += dtMax;
		if (t > tstart) {		//in last cycle of run, compute stroke work		
			dVlv = q_la - q_lv - q_LVAD;			// d/dt LV volume
			LV_SW -= Plv * dVlv * dtMax;			// update LV stroke work
			dVrv = q_ra - q_rv;						// d/dt RV volume
			RV_SW -= Prv * dVrv * dtMax;			// update RV stroke work
		}
		if (t > tstart0) {		//in last cycle of run, save data for short-axis plots
			kk++;
			tstart0 += Tc / 4.;
			if (kk <= 4) {
				a = a0 + a1;
				for (i = 0; i < Nphi; i++) {
					sh_store2[i][kk][1] = a * shmu[i][Nnu_mid][0];
					sh_store2[i][kk][2] = a * shmu[i][Nnu_mid][Nmu1 - 1];
					sh_store2[i][kk][3] = a * shmu[i][Nnu_mid][Nmu1];
					sh_store2[i][kk][4] = a * shmu[i][Nnu_mid][Nmu - 1];
				}
				for (j = 0; j < Nnu; j++) {
					ch_store[j][kk][1] = a * chmu[0][j][0];
					ch_store[j][kk][2] = a * chmu[0][j][Nmu1 - 1];
					ch_store[j][kk][3] = a * chmu[0][j][Nmu1];
					ch_store[j][kk][4] = a * chmu[0][j][Nmu - 1];
					ch_store[j][kk][5] = a * chmu[Nphi2][j][0];
					ch_store[j][kk][6] = a * chmu[Nphi2][j][Nmu1 - 1];
					ch_store[j][kk][7] = a * chmu[Nphi2][j][Nmu1];
					ch_store[j][kk][8] = a * chmu[Nphi2][j][Nmu - 1];
					sh_store[j][kk][1] = a * shmu[0][j][0];
					sh_store[j][kk][2] = a * shmu[0][j][Nmu1 - 1];
					sh_store[j][kk][3] = a * shmu[0][j][Nmu1];
					sh_store[j][kk][4] = a * shmu[0][j][Nmu - 1];
					sh_store[j][kk][5] = a * shmu[Nphi2][j][0];
					sh_store[j][kk][6] = a * shmu[Nphi2][j][Nmu1 - 1];
					sh_store[j][kk][7] = a * shmu[Nphi2][j][Nmu1];
					sh_store[j][kk][8] = a * shmu[Nphi2][j][Nmu - 1];
				}
			}
		}
	}
}
