                                                                                                       // equilibriumIntegrals.cpp
// edited by TWS, November 2019, June 2025
#include <stdio.h>
#include <math.h>
#include "nrutil.h"

double Integral3(double ***f);

void updateIntegrals(double t) {
	extern int Nmu, Nnu, Nphi, Npar;
	extern double** alpha_mat_aij, * kappa_ai;
	extern double*** Se11, *** Se12, *** Se13, *** Se22, *** Se23, *** Se33;
	extern double**** dE11_dai, **** dE12_dai, **** dE13_dai, **** dE22_dai, **** dE23_dai, **** dE33_dai;
	extern double*** Sf11_const, *** Sf12_const, *** Sf13_const, *** Sf22_const, *** Sf23_const, *** Sf33_const;
	extern double**** Sf11_ai, **** Sf12_ai, **** Sf13_ai, **** Sf22_ai, **** Sf23_ai, **** Sf33_ai;
	extern double**** Sv11_ai, **** Sv12_ai, **** Sv13_ai, **** Sv22_ai, **** Sv23_ai, **** Sv33_ai;
	extern double*** Frefinv11, *** Frefinv12, *** Frefinv13, *** Frefinv22, *** Frefinv33;
	extern double*** integrand, *** integratingFactor;
	extern double*** S0_11, *** S0_12, *** S0_13, *** S0_22, *** S0_23, *** S0_33;
	extern double**** S0_11_ai, **** S0_12_ai, **** S0_13_ai, **** S0_22_ai, **** S0_23_ai, **** S0_33_ai;

	int k, j, i, ip, jp;
	double integ;
	double S1_11, S1_12, S1_13, S1_22, S1_23, S1_33;
	double S1_11_ai, S1_12_ai, S1_13_ai, S1_22_ai, S1_23_ai, S1_33_ai;
	double e11, e12, e13, e21, e22, e23, e31, e32, e33;

	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
		// Elastic stresses in reference configuration
		S1_11 = Se11[i][j][k] + Sf11_const[i][j][k];
		S1_12 = Se12[i][j][k] + Sf12_const[i][j][k];
		S1_13 = Se13[i][j][k] + Sf13_const[i][j][k];
		S1_22 = Se22[i][j][k] + Sf22_const[i][j][k];
		S1_23 = Se23[i][j][k] + Sf23_const[i][j][k];
		S1_33 = Se33[i][j][k] + Sf33_const[i][j][k];
		// Pullback to axisymmetric configuration S0 = F0^(-1) S F0^(-T), eq. (S24.4)
		e11 = S1_11 * Frefinv11[i][j][k] + S1_12 * Frefinv12[i][j][k] + S1_13 * Frefinv13[i][j][k];
		e12 = S1_12 * Frefinv22[i][j][k];
		e13 = S1_13 * Frefinv33[i][j][k];
		e21 = S1_12 * Frefinv11[i][j][k] + S1_22 * Frefinv12[i][j][k] + S1_23 * Frefinv13[i][j][k];
		e22 = S1_22 * Frefinv22[i][j][k];
		e23 = S1_23 * Frefinv33[i][j][k];
		e31 = S1_13 * Frefinv11[i][j][k] + S1_23 * Frefinv12[i][j][k] + S1_33 * Frefinv13[i][j][k];
		e32 = S1_23 * Frefinv22[i][j][k];
		e33 = S1_33 * Frefinv33[i][j][k];
		S0_11[i][j][k] = Frefinv11[i][j][k] * e11 + Frefinv12[i][j][k] * e21 + Frefinv13[i][j][k] * e31;
		S0_12[i][j][k] = Frefinv11[i][j][k] * e12 + Frefinv12[i][j][k] * e22 + Frefinv13[i][j][k] * e32;
		S0_13[i][j][k] = Frefinv11[i][j][k] * e13 + Frefinv12[i][j][k] * e23 + Frefinv13[i][j][k] * e33;
		S0_22[i][j][k] = Frefinv22[i][j][k] * e22;
		S0_23[i][j][k] = Frefinv22[i][j][k] * e23;
		S0_33[i][j][k] = Frefinv33[i][j][k] * e33;
		for (ip = 1; ip <= Npar; ip++) {
			//Viscous stress derivatives in reference frame
			S1_11_ai = Sv11_ai[i][j][k][ip] + Sf11_ai[i][j][k][ip];
			S1_12_ai = Sv12_ai[i][j][k][ip] + Sf12_ai[i][j][k][ip];
			S1_13_ai = Sv13_ai[i][j][k][ip] + Sf13_ai[i][j][k][ip];
			S1_22_ai = Sv22_ai[i][j][k][ip] + Sf22_ai[i][j][k][ip];
			S1_23_ai = Sv23_ai[i][j][k][ip] + Sf23_ai[i][j][k][ip];
			S1_33_ai = Sv33_ai[i][j][k][ip] + Sf33_ai[i][j][k][ip];
			// Pullback to axisymmetric configuration S0 = F0^(-1) S F0^(-T), eq. (S24.4)
			e11 = S1_11_ai * Frefinv11[i][j][k] + S1_12_ai * Frefinv12[i][j][k] + S1_13_ai * Frefinv13[i][j][k];
			e12 = S1_12_ai * Frefinv22[i][j][k];
			e13 = S1_13_ai * Frefinv33[i][j][k];
			e21 = S1_12_ai * Frefinv11[i][j][k] + S1_22_ai * Frefinv12[i][j][k] + S1_23_ai * Frefinv13[i][j][k];
			e22 = S1_22_ai * Frefinv22[i][j][k];
			e23 = S1_23_ai * Frefinv33[i][j][k];
			e31 = S1_13_ai * Frefinv11[i][j][k] + S1_23_ai * Frefinv12[i][j][k] + S1_33_ai * Frefinv13[i][j][k];
			e32 = S1_23_ai * Frefinv22[i][j][k];
			e33 = S1_33_ai * Frefinv33[i][j][k];
			S0_11_ai[i][j][k][ip] = Frefinv11[i][j][k] * e11 + Frefinv12[i][j][k] * e21 + Frefinv13[i][j][k] * e31;
			S0_12_ai[i][j][k][ip] = Frefinv11[i][j][k] * e12 + Frefinv12[i][j][k] * e22 + Frefinv13[i][j][k] * e32;
			S0_13_ai[i][j][k][ip] = Frefinv11[i][j][k] * e13 + Frefinv12[i][j][k] * e23 + Frefinv13[i][j][k] * e33;
			S0_22_ai[i][j][k][ip] = Frefinv22[i][j][k] * e22;
			S0_23_ai[i][j][k][ip] = Frefinv22[i][j][k] * e23;
			S0_33_ai[i][j][k][ip] = Frefinv33[i][j][k] * e33;
		}
	}
	// KAPPA INTEGRALS
	for (ip = 1; ip <= Npar; ip++) {
		for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
			integ = dE11_dai[i][j][k][ip] * S0_11[i][j][k] + dE22_dai[i][j][k][ip] * S0_22[i][j][k]
				+ dE33_dai[i][j][k][ip] * S0_33[i][j][k] + 2. * (dE12_dai[i][j][k][ip] * S0_12[i][j][k]
				+ dE13_dai[i][j][k][ip] * S0_13[i][j][k] + dE23_dai[i][j][k][ip] * S0_23[i][j][k]);		//Eq. (S93)
			integrand[i][j][k] = integ * integratingFactor[i][j][k];
		}
		kappa_ai[ip] = Integral3(integrand);
	}
	// ALPHA_MAT INTEGRALS - use symmetry wrt ip, jp
	for (ip = 1; ip <= Npar; ip++) for (jp = ip; jp <= Npar; jp++) {
		if (ip == jp || ip <= 3) {
			for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
				integ = dE11_dai[i][j][k][ip] * S0_11_ai[i][j][k][jp] + dE22_dai[i][j][k][ip] * S0_22_ai[i][j][k][jp]
					+ dE33_dai[i][j][k][ip] * S0_33_ai[i][j][k][jp] + 2. * (dE12_dai[i][j][k][ip] * S0_12_ai[i][j][k][jp]
					+ dE13_dai[i][j][k][ip] * S0_13_ai[i][j][k][jp] + dE23_dai[i][j][k][ip] * S0_23_ai[i][j][k][jp]);	//Eq. (S93)
				integrand[i][j][k] = integ * integratingFactor[i][j][k];
			}
			alpha_mat_aij[ip][jp] = Integral3(integrand);
		}
		else alpha_mat_aij[ip][jp] = 0.;
	}
	for (ip = 2; ip <= Npar; ip++) for (jp = 1; jp < ip; jp++) alpha_mat_aij[ip][jp] = alpha_mat_aij[jp][ip];
}
