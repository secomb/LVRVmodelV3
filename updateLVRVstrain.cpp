// updateLVRVstrain.cpp
// TWS, February 2024, June 2025
// equation numbers refer LVRV model supplementary material (S)
#include <math.h>
#include <stdio.h>
#include "nrutil.h"

// Deformation gradient tensor F0 or F1, relative to axisymmetric configuration
void updateDefGradTensor(double t) {
	extern int Nmu, Nnu, Nphi, Npar;
	extern double a0, a1;
	extern double*** gmu0, *** gphi0, * snu;
	extern double*** gmu, *** gphi, **** dgmu_dai, **** dgphi_dai;
	extern double*** shmu, *** chmu, *** sum_sq, *** root_sum_sq;
	extern double*** dmu_dmu0, *** dmu_dphi0, **** dmu_dai, *** dmu_dnu, * dphi_dnu;
	extern double**** ddmu_dmu0dai, **** ddmu_dphi0dai, **** ddmu_dnudai;
	extern double*** F11, *** F12, *** F13, *** F22, *** F32, *** F33;
	extern double**** dF11_dai, **** dF12_dai, **** dF13_dai, **** dF22_dai, **** dF32_dai, **** dF33_dai;

	int k, j, i, ip;
	double a = a0 + a1;
	double gmu0inv, gphi0inv;

	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
		gmu[i][j][k] = a * root_sum_sq[i][j][k];							//Eq. after (S1)		
		gphi[i][j][k] = a * shmu[i][j][k] * snu[j];							//Eq. after (S1)
		gmu0inv = 1. / gmu0[i][j][k];
		gphi0inv = 1. / gphi0[i][j][k];
		//deformation gradient tensor F0 or F1
		F11[i][j][k] = gmu[i][j][k] * gmu0inv * dmu_dmu0[i][j][k];			//Eq. (S13)		
		F12[i][j][k] = gmu[i][j][k] * gmu0inv * dmu_dnu[i][j][k];
		F13[i][j][k] = gmu[i][j][k] * gphi0inv * dmu_dphi0[i][j][k];
		F22[i][j][k] = gmu[i][j][k] * gmu0inv;
		F32[i][j][k] = gphi[i][j][k] * gmu0inv * dphi_dnu[j];
		F33[i][j][k] = gphi[i][j][k] * gphi0inv;
		for (ip = 1; ip <= Npar; ip++) {
			dF11_dai[i][j][k][ip] = (dgmu_dai[i][j][k][ip] * dmu_dmu0[i][j][k] + gmu[i][j][k] * ddmu_dmu0dai[i][j][k][ip]) * gmu0inv;	//Eq. (S38,52)
			dF12_dai[i][j][k][ip] = (dgmu_dai[i][j][k][ip] * dmu_dnu[i][j][k] + gmu[i][j][k] * ddmu_dnudai[i][j][k][ip]) * gmu0inv;		//Eq. (S39,53)
			dF13_dai[i][j][k][ip] = (dgmu_dai[i][j][k][ip] * dmu_dphi0[i][j][k] + gmu[i][j][k] * ddmu_dphi0dai[i][j][k][ip]) * gphi0inv;//Eq. (S40,54)
			dF22_dai[i][j][k][ip] = dgmu_dai[i][j][k][ip] * gmu0inv;																	//Eq. (S41,55)
			if (ip == 3) dF32_dai[i][j][k][ip] = -snu[j] * gphi[i][j][k] * gmu0inv;														//Eq. (S44)
			else dF32_dai[i][j][k][ip] = dgphi_dai[i][j][k][ip] * dphi_dnu[j] * gmu0inv;												//Eq. (S42,56)
			dF33_dai[i][j][k][ip] = dgphi_dai[i][j][k][ip] * gphi0inv;																	//Eq. (S43,57)
		}
	}
}

// Right Cauchy-Green tensor C0 = F0^T F0 or C1 = F1^T F1, relative to axisymmetric configuration
void updateCauchyTensor() {
	extern int Nmu, Nnu, Nphi;
	extern double*** F11, *** F12, *** F13, *** F22, *** F32, *** F33;
	extern double*** C11, *** C12, *** C13, *** C22, *** C23, *** C33;
	extern double*** Cinv11, *** Cinv12, *** Cinv13, *** Cinv22, *** Cinv23, *** Cinv33;

	int k, j, i;
	double detC = 0.;
	double byDetC = 0.;

	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
		C11[i][j][k] = DSQR(F11[i][j][k]);
		C12[i][j][k] = F11[i][j][k] * F12[i][j][k];
		C13[i][j][k] = F11[i][j][k] * F13[i][j][k];
		C22[i][j][k] = DSQR(F12[i][j][k]) + DSQR(F22[i][j][k]) + DSQR(F32[i][j][k]);
		C23[i][j][k] = F12[i][j][k] * F13[i][j][k] + F32[i][j][k] * F33[i][j][k];
		C33[i][j][k] = DSQR(F13[i][j][k]) + DSQR(F33[i][j][k]);
		detC = DSQR(F11[i][j][k] * F22[i][j][k] * F33[i][j][k]);
		byDetC = 1. / detC;
		/* Inverse of a symmetric matrix = 1/determinant *
		[- c23^2 + c22*c33, c13*c23 - c12*c33, c12*c23 - c13*c22]
		[c13*c23 - c12*c33, - c13^2 + c11*c33, c12*c13 - c11*c23]
		[c12*c23 - c13*c22, c12*c13 - c11*c23, - c12^2 + c11*c22]
		*/
		Cinv11[i][j][k] = byDetC * (C22[i][j][k] * C33[i][j][k] - C23[i][j][k] * C23[i][j][k]);
		Cinv12[i][j][k] = byDetC * (C13[i][j][k] * C23[i][j][k] - C12[i][j][k] * C33[i][j][k]);
		Cinv13[i][j][k] = byDetC * (C12[i][j][k] * C23[i][j][k] - C13[i][j][k] * C22[i][j][k]);
		Cinv22[i][j][k] = byDetC * (C11[i][j][k] * C33[i][j][k] - C13[i][j][k] * C13[i][j][k]);
		Cinv23[i][j][k] = byDetC * (C12[i][j][k] * C13[i][j][k] - C11[i][j][k] * C23[i][j][k]);
		Cinv33[i][j][k] = byDetC * (C11[i][j][k] * C22[i][j][k] - C12[i][j][k] * C12[i][j][k]);
	}
}

// Compute derivatives of Green-Lagrange strain tensor E0 or E1 with respect to ai.
// Also compute E2 and rotate to fiber coordinates.
// Explicit calculation of Greens strain and pullback to axisymmetric configuration, Eq. (S24.3), are
// not needed because only derivatives of E wrt a_i are used in computation.
void updateGreenStrain(double t) {
	extern int Nmu, Nnu, Nphi, Npar;
	extern double*** F11, *** F12, *** F13, *** F22, *** F32, *** F33;
	extern double**** dF11_dai, **** dF12_dai, **** dF13_dai, **** dF22_dai, **** dF32_dai, **** dF33_dai;
	extern double**** dE11_dai, **** dE12_dai, **** dE13_dai, **** dE22_dai, **** dE23_dai, **** dE33_dai;
	extern double*** Ess, *** Esn, *** Enn, *** Esf, *** Eff, *** Enf;
	extern double*** Qs1, *** Qn1, *** Qf1, *** Qs2, *** Qn2, *** Qf2, *** Qs3, *** Qn3, *** Qf3;
	extern double*** Frefinv11, *** Frefinv12, *** Frefinv13, *** Frefinv22, *** Frefinv33;
	
	int k, j, i, ip;
	double e1s, e1n, e1f, e2s, e2n, e2f, e3s, e3n, e3f;
	double F2_11, F2_12, F2_13, F2_22, F2_32, F2_33;
	double C2_11, C2_12, C2_13, C2_22, C2_23, C2_33;
	double E2_11, E2_12, E2_13, E2_22, E2_23, E2_33;

	if (t >= 0.) {
		for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {		
			//Derivatives of Greens strain
			for (ip = 1; ip <= Npar; ip++) {			//Eqs. (S57.1) - (S57.6)
				dE11_dai[i][j][k][ip] = F11[i][j][k] * dF11_dai[i][j][k][ip];
				dE12_dai[i][j][k][ip] = 0.5*(F12[i][j][k] * dF11_dai[i][j][k][ip] + F11[i][j][k] * dF12_dai[i][j][k][ip]);
				dE13_dai[i][j][k][ip] = 0.5*(F13[i][j][k] * dF11_dai[i][j][k][ip] + F11[i][j][k] * dF13_dai[i][j][k][ip]);
				dE22_dai[i][j][k][ip] = F12[i][j][k] * dF12_dai[i][j][k][ip] + F22[i][j][k] * dF22_dai[i][j][k][ip] + F32[i][j][k] * dF32_dai[i][j][k][ip];
				dE23_dai[i][j][k][ip] = 0.5*(F32[i][j][k] * dF33_dai[i][j][k][ip] + F33[i][j][k] * dF32_dai[i][j][k][ip]
					 + F12[i][j][k] * dF13_dai[i][j][k][ip] + F13[i][j][k] * dF12_dai[i][j][k][ip]);
				dE33_dai[i][j][k][ip] = F13[i][j][k] * dF13_dai[i][j][k][ip] + F33[i][j][k] * dF33_dai[i][j][k][ip];
			}
			//deformation gradient tensor F2 = F1 F0^(-1). Eq. after (S24.1)	
			F2_11 = F11[i][j][k] * Frefinv11[i][j][k];
			F2_12 = F11[i][j][k] * Frefinv12[i][j][k] + F12[i][j][k] * Frefinv22[i][j][k];
			F2_13 = F11[i][j][k] * Frefinv13[i][j][k] + F13[i][j][k] * Frefinv33[i][j][k];
			F2_22 = F22[i][j][k] * Frefinv22[i][j][k];
			F2_32 = F32[i][j][k] * Frefinv22[i][j][k];
			F2_33 = F33[i][j][k] * Frefinv33[i][j][k];
			// Right Cauchy-Green tensor C2 = F2^T F2, relative to reference configuration
			C2_11 = DSQR(F2_11);
			C2_12 = F2_11 * F2_12;
			C2_13 = F2_11 * F2_13;
			C2_22 = DSQR(F2_12) + DSQR(F2_22) + DSQR(F2_32);
			C2_23 = F2_12 * F2_13 + F2_32 * F2_33;
			C2_33 = DSQR(F2_13) + DSQR(F2_33);
			// Green-Lagrange strain tensor E2 = (1/2) (C2 - 1), relative to reference configuration
			E2_11 = 0.5 * (C2_11 - 1);
			E2_12 = 0.5 * C2_12;
			E2_13 = 0.5 * C2_13;
			E2_22 = 0.5 * (C2_22 - 1);
			E2_23 = 0.5 * C2_23;
			E2_33 = 0.5 * (C2_33 - 1);
			// Rotate E2_prolate to E2_fiber for use in elastic stress computation.
			// E2(snf) = Q E2(mu nu phi) Q^T, eq. (S66)
			e1s = E2_11 * Qs1[i][j][k] + E2_12 * Qs2[i][j][k] + E2_13 * Qs3[i][j][k];
			e1n = E2_11 * Qn1[i][j][k] + E2_12 * Qn2[i][j][k] + E2_13 * Qn3[i][j][k];
			e1f = E2_11 * Qf1[i][j][k] + E2_12 * Qf2[i][j][k] + E2_13 * Qf3[i][j][k];
			e2s = E2_12 * Qs1[i][j][k] + E2_22 * Qs2[i][j][k] + E2_23 * Qs3[i][j][k];
			e2n = E2_12 * Qn1[i][j][k] + E2_22 * Qn2[i][j][k] + E2_23 * Qn3[i][j][k];
			e2f = E2_12 * Qf1[i][j][k] + E2_22 * Qf2[i][j][k] + E2_23 * Qf3[i][j][k];
			e3s = E2_13 * Qs1[i][j][k] + E2_23 * Qs2[i][j][k] + E2_33 * Qs3[i][j][k];
			e3n = E2_13 * Qn1[i][j][k] + E2_23 * Qn2[i][j][k] + E2_33 * Qn3[i][j][k];
			e3f = E2_13 * Qf1[i][j][k] + E2_23 * Qf2[i][j][k] + E2_33 * Qf3[i][j][k];
			Ess[i][j][k] = Qs1[i][j][k] * e1s + Qs2[i][j][k] * e2s + Qs3[i][j][k] * e3s;
			Esn[i][j][k] = Qs1[i][j][k] * e1n + Qs2[i][j][k] * e2n + Qs3[i][j][k] * e3n;
			Esf[i][j][k] = Qs1[i][j][k] * e1f + Qs2[i][j][k] * e2f + Qs3[i][j][k] * e3f;
			Enn[i][j][k] = Qn1[i][j][k] * e1n + Qn2[i][j][k] * e2n + Qn3[i][j][k] * e3n;
			Enf[i][j][k] = Qn1[i][j][k] * e1f + Qn2[i][j][k] * e2f + Qn3[i][j][k] * e3f;
			Eff[i][j][k] = Qf1[i][j][k] * e1f + Qf2[i][j][k] * e2f + Qf3[i][j][k] * e3f;
		}
	}
}

void updateFiberStrain(double t) {
	extern int Nmu, Nnu, Nphi, Npar;
	extern double a0, a1, a2, a3, Tc, Ta, dtMax;
	extern double*** shmu, *** chmu, *** sh_sq, * snu, * snu_sq, *** sum_sq, *** ell, *** ell_ref, *** ell0, *** eps_f, *** eps_fed, *** dmu_dnu;
	extern double**** dell_dai, *** omega;
	extern double**** dmu_dai, **** ddmu_dnudai, **** ddmu_dphi0dai;
	extern double*** gmu, *** gphi, **** dgmu_dai, **** dgphi_dai;
	extern double*** dmu_dphi0, * dphi_dnu, *** gnu0, *** gphi0;
	extern double*** root_sum_sq;

	int k, j, i, ip;
	double a = a0 + a1;
	double dphi0_dnu0, denom, ddenom_dai;
	double dmu_dnu_t, dphi_dnu_t, ddmu_dnudai_t, ddphi_dnudai_t, truncT = t - (int)(t / Tc) * Tc; ;

	//fiber stretch in reference and deformed configuration
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
		dphi0_dnu0 = -omega[i][j][k] * snu[j];										//Eq. (S58)
		dmu_dnu_t = dmu_dnu[i][j][k] + dmu_dphi0[i][j][k] * dphi0_dnu0;				//Total derivative, Eq. (S64.2)
		dphi_dnu_t = dphi_dnu[j] + dphi0_dnu0;										//Eq. (S64.2)
		denom = DSQR(gmu[i][j][k]) * (1 + DSQR(dmu_dnu_t)) + DSQR(gphi[i][j][k] * dphi_dnu_t);	//Eq. (S64.1)
		ell[i][j][k] = sqrt(denom);	//Eq. (S64.1)
		if (t >= 0.) {
			// Store end-diastolic strain at appropriate time
			eps_f[i][j][k] = 0.5 * (DSQR(ell[i][j][k]) / DSQR(ell_ref[i][j][k]) - 1.);		//fiber strain, after Eq. (S85)
			if (fabs(truncT - Tc + Ta) < dtMax) eps_fed[i][j][k] = eps_f[i][j][k];
			//Derivatives of eq. (S64.1) wrt ai
			for (ip = 1; ip <= Npar; ip++) {
				ddmu_dnudai_t = ddmu_dnudai[i][j][k][ip] + ddmu_dphi0dai[i][j][k][ip] * dphi0_dnu0;//Eq. (S64.2)
				if (ip == 3) ddphi_dnudai_t = -snu[j];
				else ddphi_dnudai_t = 0.;
				ddenom_dai = 2. * (gmu[i][j][k] * dgmu_dai[i][j][k][ip] * (1 + DSQR(dmu_dnu_t)) 
					+ DSQR(gmu[i][j][k]) * dmu_dnu_t * ddmu_dnudai_t 
					+ gphi[i][j][k] * dgphi_dai[i][j][k][ip] * DSQR(dphi_dnu_t)
					+ DSQR(gphi[i][j][k]) * dphi_dnu_t * ddphi_dnudai_t);
				dell_dai[i][j][k][ip] = ddenom_dai / ell[i][j][k] / 2.;
			}
		}
	}
}
