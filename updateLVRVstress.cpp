// updateLVstress.cpp
// edited by TWS, November 2019, March 2024, June 2025
#include <math.h>
#include <stdio.h>
#include "nrutil.h"

void computeElastic() {												//Eq. (S76)
	extern int Nmu, Nnu, Nphi, Nmu1;
	extern double c1, bff, bxx, bfx;
	extern double c1_rv, br_rv, bperp_rv;
	extern double bff_rv, bxx_rv, bfx_rv;
	extern double* phi0, phi_s;
	extern double*** Ess, *** Esn, *** Enn, *** Esf, *** Eff, *** Enf;
	extern double*** Se11, *** Se12, *** Se13, *** Se22, *** Se23, *** Se33;
	extern double*** Qs1, *** Qn1, *** Qf1, *** Qs2, *** Qn2, *** Qf2, *** Qs3, *** Qn3, *** Qf3;

	double c10, bff0, bxx0, bfx0, exponent;
	double es1, es2, es3, en1, en2, en3, ef1, ef2, ef3;
	double Sess, Sesn, Senn, Sesf, Seff, Senf, eW;

	int k, j, i;
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
		if (phi0[i] >= phi_s && phi0[i] <= 2. * PI - phi_s && k >= Nmu1) {	//in RV free wall
			c10 = c1_rv;
			bff0 = bff_rv;
			bxx0 = bxx_rv;
			bfx0 = bfx_rv;
		}
		else {
			c10 = c1;
			bff0 = bff;
			bxx0 = bxx;
			bfx0 = bfx;
		}
		exponent = bff0 * DSQR(Eff[i][j][k]) + bxx0 * (DSQR(Enn[i][j][k]) + DSQR(Ess[i][j][k])
			+ 2. * DSQR(Esn[i][j][k])) + 2. * bfx0*(DSQR(Esf[i][j][k]) + DSQR(Enf[i][j][k]));	//Eq. (S75)
		eW = c10 * exp(exponent);																//Eq. (S76)
		Sess = eW * bxx0 * Ess[i][j][k];
		Senn = eW * bxx0 * Enn[i][j][k];
		Seff = eW * bff0 * Eff[i][j][k];
		Sesn = eW * bxx0 * Esn[i][j][k];
		Sesf = eW * bfx0 * Esf[i][j][k];
		Senf = eW * bfx0 * Enf[i][j][k];
	// Rotate from fiber to prolate spheroidal coordinates
	// Se(mu nu phi) = Q^T Se(snf) Q. Note Se is symmetric. Eq. (S76.1)
		es1 = Sess * Qs1[i][j][k] + Sesn * Qn1[i][j][k] + Sesf * Qf1[i][j][k];
		es2 = Sess * Qs2[i][j][k] + Sesn * Qn2[i][j][k] + Sesf * Qf2[i][j][k];
		es3 = Sess * Qs3[i][j][k] + Sesn * Qn3[i][j][k] + Sesf * Qf3[i][j][k];
		en1 = Sesn * Qs1[i][j][k] + Senn * Qn1[i][j][k] + Senf * Qf1[i][j][k];
		en2 = Sesn * Qs2[i][j][k] + Senn * Qn2[i][j][k] + Senf * Qf2[i][j][k];
		en3 = Sesn * Qs3[i][j][k] + Senn * Qn3[i][j][k] + Senf * Qf3[i][j][k];
		ef1 = Sesf * Qs1[i][j][k] + Senf * Qn1[i][j][k] + Seff * Qf1[i][j][k];
		ef2 = Sesf * Qs2[i][j][k] + Senf * Qn2[i][j][k] + Seff * Qf2[i][j][k];
		ef3 = Sesf * Qs3[i][j][k] + Senf * Qn3[i][j][k] + Seff * Qf3[i][j][k];
		Se11[i][j][k] = Qs1[i][j][k] * es1 + Qn1[i][j][k] * en1 + Qf1[i][j][k] * ef1;
		Se12[i][j][k] = Qs1[i][j][k] * es2 + Qn1[i][j][k] * en2 + Qf1[i][j][k] * ef2;
		Se13[i][j][k] = Qs1[i][j][k] * es3 + Qn1[i][j][k] * en3 + Qf1[i][j][k] * ef3;
		Se22[i][j][k] = Qs2[i][j][k] * es2 + Qn2[i][j][k] * en2 + Qf2[i][j][k] * ef2;
		Se23[i][j][k] = Qs2[i][j][k] * es3 + Qn2[i][j][k] * en3 + Qf2[i][j][k] * ef3;
		Se33[i][j][k] = Qs3[i][j][k] * es3 + Qn3[i][j][k] * en3 + Qf3[i][j][k] * ef3;
	}
}

void computeViscous() {
	//Eq. (S84)
	extern int Nmu, Nnu, Nphi, Npar;
	extern double kv;
	extern double ***Cinv11, ***Cinv12, ***Cinv13, ***Cinv22, ***Cinv23, ***Cinv33;
	extern double ****dE11_dai, ****dE12_dai, ****dE13_dai, ****dE22_dai, ****dE23_dai, ****dE33_dai;
	extern double ****Sv11_ai, ****Sv12_ai, ****Sv13_ai, ****Sv22_ai, ****Sv23_ai, ****Sv33_ai;

	int k, j, i, ip; 
	double kvtwo = 2. * kv;
	double p11, p12, p13, p21, p22, p23, p31, p32, p33;

	// Compute Sv with ai derivatives
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) for (ip = 1; ip <= Npar; ip++) {
		p11 = Cinv11[i][j][k] * dE11_dai[i][j][k][ip] + Cinv12[i][j][k] * dE12_dai[i][j][k][ip] + Cinv13[i][j][k] * dE13_dai[i][j][k][ip];
		p12 = Cinv11[i][j][k] * dE12_dai[i][j][k][ip] + Cinv12[i][j][k] * dE22_dai[i][j][k][ip] + Cinv13[i][j][k] * dE23_dai[i][j][k][ip];
		p13 = Cinv11[i][j][k] * dE13_dai[i][j][k][ip] + Cinv12[i][j][k] * dE23_dai[i][j][k][ip] + Cinv13[i][j][k] * dE33_dai[i][j][k][ip];
		p21 = Cinv12[i][j][k] * dE11_dai[i][j][k][ip] + Cinv22[i][j][k] * dE12_dai[i][j][k][ip] + Cinv23[i][j][k] * dE13_dai[i][j][k][ip];
		p22 = Cinv12[i][j][k] * dE12_dai[i][j][k][ip] + Cinv22[i][j][k] * dE22_dai[i][j][k][ip] + Cinv23[i][j][k] * dE23_dai[i][j][k][ip];
		p23 = Cinv12[i][j][k] * dE13_dai[i][j][k][ip] + Cinv22[i][j][k] * dE23_dai[i][j][k][ip] + Cinv23[i][j][k] * dE33_dai[i][j][k][ip];
		p31 = Cinv13[i][j][k] * dE11_dai[i][j][k][ip] + Cinv23[i][j][k] * dE12_dai[i][j][k][ip] + Cinv33[i][j][k] * dE13_dai[i][j][k][ip];
		p32 = Cinv13[i][j][k] * dE12_dai[i][j][k][ip] + Cinv23[i][j][k] * dE22_dai[i][j][k][ip] + Cinv33[i][j][k] * dE23_dai[i][j][k][ip];
		p33 = Cinv13[i][j][k] * dE13_dai[i][j][k][ip] + Cinv23[i][j][k] * dE23_dai[i][j][k][ip] + Cinv33[i][j][k] * dE33_dai[i][j][k][ip];

		Sv11_ai[i][j][k][ip] = kvtwo * (Cinv11[i][j][k] * p11 + Cinv12[i][j][k] * p12 + Cinv13[i][j][k] * p13);
		Sv12_ai[i][j][k][ip] = kvtwo * (Cinv12[i][j][k] * p11 + Cinv22[i][j][k] * p12 + Cinv23[i][j][k] * p13);
		Sv13_ai[i][j][k][ip] = kvtwo * (Cinv13[i][j][k] * p11 + Cinv23[i][j][k] * p12 + Cinv33[i][j][k] * p13);
		Sv22_ai[i][j][k][ip] = kvtwo * (Cinv12[i][j][k] * p21 + Cinv22[i][j][k] * p22 + Cinv23[i][j][k] * p23);
		Sv23_ai[i][j][k][ip] = kvtwo * (Cinv13[i][j][k] * p21 + Cinv23[i][j][k] * p22 + Cinv33[i][j][k] * p23);
		Sv33_ai[i][j][k][ip] = kvtwo * (Cinv13[i][j][k] * p31 + Cinv23[i][j][k] * p32 + Cinv33[i][j][k] * p33);
	}
}

void computeActive() {
	extern int Nmu1, Nmu2, Nmu, Nnu, Nphi, Npar;
	extern double dnu, dphi, phi_s;
	extern double Ls0, Lsmax, Lsw, km, kav, km_rv, kav_rv;
	extern double*** ell, *** ell_ref, **** dell_dai, *** integrand, *** integratingFactor;
	extern double*** Sf_ff_const, **** Sf_ff_ai;
	extern double*** Sf11_const, *** Sf12_const, *** Sf13_const, *** Sf22_const, *** Sf23_const, *** Sf33_const;
	extern double**** Sf11_ai, **** Sf12_ai, **** Sf13_ai, **** Sf22_ai, **** Sf23_ai, **** Sf33_ai;
	extern double *** SigmaF, **** SigmaF_ai;
	extern double*** sps, *** cps, *** sps_sq, *** cps_sq;
	extern double*** mu0, * nu, * phi0;
	extern double*** Sf_temp;
	extern double*** Qs1, *** Qn1, *** Qf1, *** Qs2, *** Qn2, *** Qf2, *** Qs3, *** Qn3, *** Qf3;
	extern double At, At_la;

	int k, j, i, ip;
	double ellRat, ellRatLess, activationkm, activationkv, G, LS;
	double twoLsw_sq = 2. * DSQR(Lsw);

	// Compute active stress in fiber direction
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
		LS = Ls0 * ell[i][j][k] / ell_ref[i][j][k];		//Eq. after (S86)
		G = exp(-DSQR(LS - Lsmax) / twoLsw_sq);			//Eq. (S86)

		if (phi0[i] >= phi_s && phi0[i] <= 2. * PI - phi_s && k >= Nmu1) {
			activationkm = At * G * km_rv;
			activationkv = At * G * kav_rv;
		}
		else {
			activationkm = At * G * km;
			activationkv = At * G * kav;
		}
		ellRat = DSQR(ell_ref[i][j][k] / ell[i][j][k]);
		ellRatLess = ell[i][j][k] / DSQR(ell_ref[i][j][k]);

		// Compute active stress (Cauchy and PK2)
		SigmaF[i][j][k] = activationkm;													//Eq. (S85) - Cauchy, config. 2
		Sf_ff_const[i][j][k] = activationkm * ellRat;									//Eq. (S89) - PK2, config. 1
		for (ip = 1; ip <= Npar; ip++) {
			SigmaF_ai[i][j][k][ip] = activationkv * dell_dai[i][j][k][ip] * ellRatLess;	//Eq. (S85)
			Sf_ff_ai[i][j][k][ip] = SigmaF_ai[i][j][k][ip] * ellRat;					//Eq. (S90)
		}
	}	
	// Rotate to prolate spheroidal coordinates, eq. (S76.1)
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
		Sf11_const[i][j][k] = Qf1[i][j][k] * Sf_ff_const[i][j][k] * Qf1[i][j][k];
		Sf12_const[i][j][k] = Qf1[i][j][k] * Sf_ff_const[i][j][k] * Qf2[i][j][k];
		Sf13_const[i][j][k] = Qf1[i][j][k] * Sf_ff_const[i][j][k] * Qf3[i][j][k];
		Sf22_const[i][j][k] = Qf2[i][j][k] * Sf_ff_const[i][j][k] * Qf2[i][j][k];
		Sf23_const[i][j][k] = Qf2[i][j][k] * Sf_ff_const[i][j][k] * Qf3[i][j][k];
		Sf33_const[i][j][k] = Qf3[i][j][k] * Sf_ff_const[i][j][k] * Qf3[i][j][k];
		for (ip = 1; ip <= Npar; ip++) {
			Sf11_ai[i][j][k][ip] = Qf1[i][j][k] * Sf_ff_ai[i][j][k][ip] * Qf1[i][j][k];
			Sf12_ai[i][j][k][ip] = Qf1[i][j][k] * Sf_ff_ai[i][j][k][ip] * Qf2[i][j][k];
			Sf13_ai[i][j][k][ip] = Qf1[i][j][k] * Sf_ff_ai[i][j][k][ip] * Qf3[i][j][k];
			Sf22_ai[i][j][k][ip] = Qf2[i][j][k] * Sf_ff_ai[i][j][k][ip] * Qf2[i][j][k];
			Sf23_ai[i][j][k][ip] = Qf2[i][j][k] * Sf_ff_ai[i][j][k][ip] * Qf3[i][j][k];
			Sf33_ai[i][j][k][ip] = Qf3[i][j][k] * Sf_ff_ai[i][j][k][ip] * Qf3[i][j][k];
		}
	}
}
