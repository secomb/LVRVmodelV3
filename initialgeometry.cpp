// initialgeometry.cpp
// TWS, May 2024
// Note: phi_s = beginning of transition in fiber directions, beginning of phi discretization
// phi_t = begining of RV free wall, septum
// This is configuration 0, "axisymmetric" (except when RV is hypertrophic)

#include <math.h>
#include <iostream>
#include <string>
#include <iomanip>
#include "nrutil.h"
#include <complex>

using namespace std;
typedef std::complex<double> CX;

double tripleIntegral(double*** f, double **dmu, double dnu, double dphi, int Nmus, int Nmue, int Nnu, int Nphis, int Nphie);

void initialgeometry(int Nmu1, int Nmu2, int Nmu, int Nnu, int Nphi, int Nr) {
	extern double nu_up, muin0, mumid0, muout0, a0, phi_t, phi_s, phi_d, alpha, psi_in_b0, psi_out_b0, rin0, rout0, a60;
	extern double rin0_la, rout0_la, rin0_rv, rout0_rv, rin0_ra, rout0_ra;
	extern double dnu, dphi;
	extern double *phi0, *ggphi0, *ggphi0p, ***omega;
	extern double ***mu0, *nu, *jjnu, *jjnup;
	extern double **dmu1, **dmu2, **dmu_ones;
	extern double ***cps, ***sps, ***cps_sq, ***sps_sq; 
	extern double *snu, *cnu, ***shmu0, ***chmu0;
	extern double *snu_sq, *cnu_sq, ***sum_sq0, ***root_sum_sq0, **muout;
	extern double ***ell0, ***gmu0, ***gnu0, ***gphi0, ***integratingFactor;
	extern double *r0_la, *r0_rv, *r0_ra, ***chmu;
	extern double Vwall_lv, Vwall_se, Vwall_rv;
	extern string RVhypertrophy;

	int k, j, i;
	int Nnu_mid = (PI / 2. - nu_up) / (PI - nu_up) * (Nnu - 1) + 0.5;		//find index on equator
	int Nphi_st = (Nphi - 2) * (phi_t - phi_s) / (PI - phi_s) / 2. + 0.5;	//for integration strictly over RV and septum
	double denom, phi_r, mu_val, phi_val, spl_fac, lam, lam_lv, psi_eq, psi_eq_lv;
	double a0sqr = DSQR(a0), a0cube = DCUB(a0);
	double third = 1. / 3., r3 = pow(2., 1. / 3.);
	double chmuout0 = cosh(muout0), chmuout;
	double C0, RHS, b, S;
	CX Sc, root3;

	// Calculate nu and phi0 vectors, j(nu) and g(phi0) functions
	dnu = (PI - nu_up) / (Nnu - 1);
	dphi = 2. * (PI - phi_s) / (Nphi - 2);
	for (j = 0; j < Nnu; j++) {
		nu[j] = nu_up + j * dnu * 0.999;		//to avoid singulatity at pi
		if (alpha * nu[j] <= PI) {
			jjnu[j] = DSQR(sin(alpha * nu[j]));		//Eq. (S11)
			jjnup[j] = 2. * sin(alpha * nu[j]) * cos(alpha * nu[j]) * alpha;		//Eq. (S11)
		}
		else {
			jjnu[j] = 0.;
			jjnup[j] = 0.;
		}
	}

	phi0[0] = 0.;		//represents LV
	ggphi0[0] = 0.;
	ggphi0p[0] = 0.;
	for (i = 1; i < Nphi; i++) {
		phi0[i] = phi_s + (i - 1) * dphi;	// phi ranges from phi_s to 2*PI - phi_s
		phi_val = phi0[i];
		if (phi_val > PI) phi_val = 2. * PI - phi_val;
		if (phi_val <= phi_t) {
			ggphi0[i] = 0.;
			ggphi0p[i] = 0.;
		}
		else if (phi_val < phi_t + phi_d) {
			phi_r = PI / 2. * (phi_val - phi_t) / phi_d;
			ggphi0[i] = DSQR(sin(phi_r));	//Eq. (S12)
			ggphi0p[i] = PI / phi_d * sin(phi_r) * cos(phi_r);
			if (phi0[i] > PI) ggphi0p[i] = -ggphi0p[i];
		}
		else {
			ggphi0[i] = 1.;
			ggphi0p[i] = 0.;
		}
	}

	// Set up initial spherical chamber values
	for (k = 0; k < Nr; k++) {
		r0_la[k] = rin0_la + k * (rout0_la - rin0_la) / (Nr - 1);
		r0_rv[k] = rin0_rv + k * (rout0_rv - rin0_rv) / (Nr - 1);
		r0_ra[k] = rin0_ra + k * (rout0_ra - rin0_ra) / (Nr - 1);
	}
	// Set up functions of nu0
	for (j = 0; j < Nnu; j++) {
		snu[j] = sin(nu[j]);
		cnu[j] = cos(nu[j]);
		snu_sq[j] = DSQR(snu[j]);
		cnu_sq[j] = DSQR(cnu[j]);
	}

	// Set up muout for RV hypertrophy
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) {
		if (RVhypertrophy == "on") {
			C0 = third - cnu_sq[j];
			b = -cnu_sq[j];
			RHS = a0cube * (chmuout0 * (third * DSQR(chmuout0) - cnu_sq[j]) - C0) + a60 * jjnu[j] * ggphi0[i];
			S = -(RHS / a0cube + C0);
			Sc = S;
			root3 = pow(sqrt(4. * DCUB(b) + 9. * Sc * Sc) - 3. * Sc, third);
			chmuout = real(root3 / r3 - r3 * b / root3);
			muout[i][j] = acosh(chmuout);
		}
		else muout[i][j] = muout0;
	};

	// Calculate mu0 vector for each i, j
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) {
		dmu1[i][j] = (mumid0 - muin0) / (Nmu1 - 1);
		dmu2[i][j] = (muout[i][j] - mumid0) / (Nmu2 - 1);
		dmu_ones[i][j] = 1.;
		for (k = 0; k < Nmu1; k++) mu0[i][j][k] = muin0 + k * dmu1[i][j];
		for (k = 0; k < Nmu2; k++) mu0[i][j][Nmu1 + k] = mumid0 + k * dmu2[i][j];
	}		
	
	// Set up functions in axisymmetric configuration
	for (i = 0; i < Nphi; i++) {
		phi_val = phi0[i];
		if (phi_val > PI) phi_val = 2. * PI - phi_val;
		if (i == 0) spl_fac = 0;	// factor for transition from single gradient to double gradient
		else {
			spl_fac = (phi_val - phi_s) / (phi_t - phi_s);
			spl_fac = DMIN(spl_fac, 1.);
		}
		for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
			// Calculate fiber angles at nu = pi/2 with linear dependence on mu values in LV free wall
			mu_val = mu0[i][Nnu_mid][k];
			//modified for clarity, August 2025. See Eq. (62)
			if (k < Nmu1) lam = (mu_val - muin0) / (mumid0 - muin0);//position in inner layer
			else lam = (mu_val - mumid0) / (muout[i][j] - mumid0);	//position in outer layer
			psi_eq = psi_in_b0 * (1. - lam) + psi_out_b0 * lam;		//based on RV or septum
			lam_lv = (mu_val - muin0) / (muout[i][j] - muin0);		//position in combined layers
			psi_eq_lv = psi_in_b0 * (1. - lam_lv) + psi_out_b0 * lam_lv;	//based on LV
			psi_eq = spl_fac * psi_eq + (1. - spl_fac) * psi_eq_lv;	// including transition
			//with luck, psi_eq is never exactly zero
			omega[i][j][k] = 1 / (tan(psi_eq) * tanh(mu0[i][j][k]));//Eq. (S63);
			shmu0[i][j][k] = sinh(mu0[i][j][k]);
			chmu0[i][j][k] = cosh(mu0[i][j][k]);
			sum_sq0[i][j][k] = DSQR(shmu0[i][j][k]) + snu_sq[j];
			root_sum_sq0[i][j][k] = sqrt(sum_sq0[i][j][k]);
			gmu0[i][j][k] = a0 * root_sum_sq0[i][j][k];			//Eq. (S13.1)
			gnu0[i][j][k] = gmu0[i][j][k];						//Eq. (S13.1)
			gphi0[i][j][k] = a0 * shmu0[i][j][k] * snu[j];			//Eq. (S13.1)
			integratingFactor[i][j][k] = gmu0[i][j][k] * gnu0[i][j][k] * gphi0[i][j][k];
			denom = sqrt(sum_sq0[i][j][k] + DSQR(shmu0[i][j][k] * snu_sq[j] * omega[i][j][k]));
			cps[i][j][k] = shmu0[i][j][k] * snu_sq[j] * omega[i][j][k] / denom;	//Eq. (S60)
			sps[i][j][k] = root_sum_sq0[i][j][k] / denom;	//Eq. (S60.1)
			cps_sq[i][j][k] = DSQR(cps[i][j][k]);
			sps_sq[i][j][k] = DSQR(sps[i][j][k]);
			ell0[i][j][k] = a0 * denom;				//Eq. (S64)
			chmu[i][j][k] = chmu0[i][j][k];		//initialize cosh(mu)
		}
	}

	//calculate wall volumes
	Vwall_lv = 2. * tripleIntegral(integratingFactor, dmu1, dnu, phi_s, 0, Nmu1, Nnu, 0, 1);					// LV free wall
	Vwall_lv += 2. * tripleIntegral(integratingFactor, dmu2, dnu, phi_s, Nmu1, Nmu, Nnu, 0, 1);					// LV free wall
	Vwall_lv += tripleIntegral(integratingFactor, dmu1, dnu, dphi, 0, Nmu1, Nnu, 1, Nphi_st + 1);				// LV free wall
	Vwall_lv += tripleIntegral(integratingFactor, dmu1, dnu, dphi, 0, Nmu1, Nnu, Nphi - Nphi_st, Nphi);			// LV free wall
	Vwall_lv += tripleIntegral(integratingFactor, dmu2, dnu, dphi, Nmu1, Nmu, Nnu, 1, Nphi_st + 1);				// LV free wall
	Vwall_lv += tripleIntegral(integratingFactor, dmu2, dnu, dphi, Nmu1, Nmu, Nnu, Nphi - Nphi_st, Nphi);		// LV free wall
	Vwall_se = tripleIntegral(integratingFactor, dmu1, dnu, dphi, 0, Nmu1, Nnu, Nphi_st + 1, Nphi - Nphi_st);	// septum
	Vwall_rv = tripleIntegral(integratingFactor, dmu2, dnu, dphi, Nmu1, Nmu, Nnu, Nphi_st + 1, Nphi - Nphi_st);	// RV free wall
	printf("Vwall_lv = %8.4f, Vwall_se = %8.4f, Vwall_rv = %8.4f,", Vwall_lv, Vwall_se, Vwall_rv);
}
