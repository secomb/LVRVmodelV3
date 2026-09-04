// updateMapping.cpp
// edited by TWS, November 2019; Feb. 2024
// Mapping is relative to configuration 0, "axisymmetric"
#include <math.h>
#include <complex>
#include "nrutil.h"

using namespace std;
typedef std::complex<double> CX;

double simpsons(double* f, double dx, int Npoints);
double doubleIntegral(double** f, double dnu, double dphi, int Nnu, int Nphis, int Nphie);
double tripleIntegral(double*** f, double **dmu, double dnu, double dphi, int Nmus, int Nmue, int Nnu, int Nphis, int Nphie);

//double calculate_mu(double mu0, double mu0in, double nu0, double a0, double a, double a2){
//	//this is less efficient than the version in updateMapping but useful as a utility
//	double mu, cnu_sq = DSQR(cos(nu0)), chin0 = cosh(mu0in), chmu0 = cosh(mu0);
//	double C0, F, RHS, b;
//	double r3 = pow(2., 1. / 3.);
//	CX S, root3, muCX;
//
//	C0 = 1./3. - cnu_sq;
//	F = a2 * DSQR(a0) * (chin0 * (1./3. * DSQR(chin0) - cnu_sq) - C0);
//	RHS = DCUB(a0) * (chmu0 * (1./3. * DSQR(chmu0) - cnu_sq) - C0) - F;
//	b = -cnu_sq;
//	S = -(RHS / DCUB(a) + C0);	// Intermediate root calculations are complex
//	root3 = pow(sqrt(4. * DCUB(b) + 9. * S * S) - 3. * S, 1./3.);
//	muCX = acosh(root3 / r3 - r3 * b / root3);
//	mu = real(muCX);		// The imaginary part is 0, but still need to take only real part
//	return mu;
//}

void updateMapping(double t) {
	extern int Nmu1, Nmu2, Nmu, Nnu, Nphi, Npar;
	extern double **dmu1, **dmu2, **dmu_ones, dnu, dphi, phi_s;
	extern double a0, a1, a2, a3, a40, a4, a5, muin0;
	extern double Vwall_rv, Vwall_lv, Vwall_se, Vlv, Vrv, Vla, Vra;;
	extern double *nu, *vol_integrand, *dphi_dnu;
	extern double *snu, *cnu, *snu_sq, *cnu_sq, ***shmu0, ***chmu0, ***sum_sq0;
	extern double ***shmu, ***chmu, ***sh_sq, ***ch_sq, ***sum_sq, ***root_sum_sq, ***dmu_dmu0, ***dmu_dphi0, ***dmu_dnu, ****dmu_dai;
	extern double ****ddmu_dnudai, ****ddmu_dmu0dai, ****ddmu_dphi0dai;
	extern double ***gmu, ***gphi, ****dgmu_dai, ****dgphi_dai;
	extern double **vIntegrand_dai, ***integrand, **integrand2;
	extern double *eta_lv_ai, *eta_rv_ai, *dff2_dai, *dff3_dai;
	extern double *jjnu, *jjnup, *ggphi0, *ggphi0p, *hhmu0;
	extern string odeSystem, septalDisp;

	int k, j, i, n, ip, n_volint = Nmu, flag = 0;
	double a, asqr, a0sqr, acube, a0cube;
	double chin0 = cosh(muin0);
	double ff2, ff3, ff3inv, ff32inv, ff1, ffac, ff4;
	double third = 1. / 3., r3 = pow(2., 1. / 3.);
	double C0, RHS, RHS0, b, LHS, Xmu, S, err;
	double c_in, c_mid1, c_mid2, coshmu, dcoshmu;
	CX Sc, root3, muCX;

	a = a0 + a1;
	asqr = DSQR(a);
	a0sqr = DSQR(a0);
	acube = DCUB(a);
	a0cube = DCUB(a0);

	//update h(mu0)
	if (odeSystem == "LVRVmodel2024") {
		for (k = 0; k < Nmu1; k++) hhmu0[k] = a5;
		for (k = 0; k < Nmu2; k++) hhmu0[Nmu1 + k] = a40 + a4;
	}
	else for (k = 0; k < Nmu; k++) hhmu0[k] = 0.;

	// Update mu values using explicit solution for root of cubic, eq. (S6)
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
		C0 = third - cnu_sq[j];
		b = -cnu_sq[j];
		RHS0 = a0cube *(chmu0[i][j][k] * (third * DSQR(chmu0[i][j][k]) - cnu_sq[j]) - C0);
		RHS0 -= a2 * a0sqr * (chin0 * (third * DSQR(chin0) - cnu_sq[j]) - C0);		//Eq. (S7) 
		RHS = RHS0 + hhmu0[k] * jjnu[j] * ggphi0[i];								//Eq. (S9)
		S = -(RHS / acube + C0);
		Xmu = chmu[i][j][k];
		LHS = Xmu * (DSQR(Xmu) * third + b) + S;
		err = abs(LHS);
		if (err > 0.01) {		//use exact solution, intermediate root calculations are complex
			Sc = S;
			root3 = pow(sqrt(4. * DCUB(b) + 9. * Sc * Sc) - 3. * Sc, third);
			Xmu = real(root3 / r3 - r3 * b / root3); // The imaginary part is 0, but still need to take only real part
		}
		else {				//use Newton iteration if possible
			n = 0;
			do {
				Xmu -= LHS / (DSQR(Xmu) + b);
				LHS = Xmu * (DSQR(Xmu) * third + b) + S;
				err = abs(LHS);
				n++;
			} while (err > 1.e-8 && n < 100);
		}
		if (Xmu < 1.) {		//Newton failed, use exact solution
			Sc = S;
			root3 = pow(sqrt(4. * DCUB(b) + 9. * Sc * Sc) - 3. * Sc, third);
			Xmu = real(root3 / r3 - r3 * b / root3);
		}			
		if (Xmu < 1.) {
			Xmu = 1.;	//just in case
			if (flag == 0) printf("*** Error: could not evaluate mu in updateMapping, t = %f\n", t);
			flag = 1;
		}
		chmu[i][j][k] = Xmu;
		shmu[i][j][k] = sqrt(DSQR(Xmu) - 1.);
	}

	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
		dphi_dnu[j] = -a3 * snu[j];														//Eq. (S14)		
		ff1 = a0cube * shmu0[i][j][k] * sum_sq0[i][j][k];								//Eq. (S16)	 
		sh_sq[i][j][k] = DSQR(shmu[i][j][k]);
		ch_sq[i][j][k] = DSQR(chmu[i][j][k]);
		sum_sq[i][j][k] = sh_sq[i][j][k] + snu_sq[j];
		root_sum_sq[i][j][k] = sqrt(sum_sq[i][j][k]);
		ffac = 2. * cnu[j] * snu[j];
		ff2 = ffac * (acube * (1. - chmu[i][j][k]) + a0cube * (chmu0[i][j][k] - 1.)
			- a2 * a0sqr * (chin0 - 1.));												//Eq. (S17)
		ff3 = acube * shmu[i][j][k] * sum_sq[i][j][k];									//Eq. (S18)
		ff4 = acube * chmu[i][j][k] * (3. * sh_sq[i][j][k] + snu_sq[j]);				//Eq. (S29.1)
		if (ff3 > 0.) {
			ff3inv = 1. / ff3;
			ff32inv = 1. / DSQR(ff3);
			dmu_dmu0[i][j][k] = ff1 * ff3inv;											//Eq. (S14)
			dmu_dnu[i][j][k] = (ff2 + hhmu0[k] * jjnup[j] * ggphi0[i]) * ff3inv;		//Eq. (S15)
			dmu_dphi0[i][j][k] = hhmu0[k] * jjnu[j] * ggphi0p[i] * ff3inv;				//Eq. (S15)
			dmu_dai[i][j][k][1] = asqr * (1. - DCUB(chmu[i][j][k])
				+ 3. * (chmu[i][j][k] - 1.) * cnu_sq[j]) * ff3inv;						//Eq. (S25)
			dmu_dai[i][j][k][2] = a0sqr * (third * (1. - DCUB(chin0))
				+ (chin0 - 1.) * cnu_sq[j]) * ff3inv;									//Eq. (S26)
			dmu_dai[i][j][k][3] = 0.;
			dmu_dai[i][j][k][4] = 0.;
			dmu_dai[i][j][k][5] = 0.;
			if (odeSystem == "LVRVmodel2024") {				
				if (k < Nmu1) {															//Eq. (S45), septum
					if (septalDisp == "on") dmu_dai[i][j][k][5] = jjnu[j] * ggphi0[i] * ff3inv;	
				}
				else dmu_dai[i][j][k][4] = jjnu[j] * ggphi0[i] * ff3inv;				//Eq. (S45), RV free wall
			}
			for (ip = 1; ip <= Npar; ip++) {
				dff2_dai[ip] = -acube * ffac * shmu[i][j][k] * dmu_dai[i][j][k][ip];	//Eq. (S27,28,46);
				dff3_dai[ip] = ff4 * dmu_dai[i][j][k][ip];								//Eq. (S30,47)
				dgmu_dai[i][j][k][ip] = a * chmu[i][j][k] * shmu[i][j][k] * dmu_dai[i][j][k][ip] / root_sum_sq[i][j][k];	//Eq. (S35,51)
				dgphi_dai[i][j][k][ip] = a * chmu[i][j][k] * snu[j] * dmu_dai[i][j][k][ip];		//Eq. (S37)
			}
			dff2_dai[1] += ffac * 3. * asqr * (1. - chmu[i][j][k]);						//Eq. (S27)
			dff2_dai[2] -= ffac * a0sqr * (chin0 - 1.);									//Eq. (S28)		
			dff3_dai[1] += 3. * asqr * shmu[i][j][k] * sum_sq[i][j][k];					//Eq. (S29)
			for (ip = 1; ip <= Npar; ip++) {
				ddmu_dmu0dai[i][j][k][ip] = -ff1 * ff32inv * dff3_dai[ip];				//Eq. (S31,48)
				ddmu_dnudai[i][j][k][ip] = dff2_dai[ip] * ff3inv
					- (ff2 + hhmu0[k] * jjnup[j] * ggphi0[i]) * dff3_dai[ip] * ff32inv;	//Eq. (S32,49)
				ddmu_dphi0dai[i][j][k][ip] = -hhmu0[k] * jjnu[j] * ggphi0p[i] * dff3_dai[ip] * ff32inv;	//Eq. (S33,S50)
			}
			if (odeSystem == "LVRVmodel2024") {
				if (k < Nmu1) {
					if (septalDisp == "on") ddmu_dnudai[i][j][k][5] += jjnup[j] * ggphi0[i] * ff3inv;	//Eq. (S49)
				}					
				else ddmu_dnudai[i][j][k][4] += jjnup[j] * ggphi0[i] * ff3inv;			//Eq. (S49)
			}
			dgmu_dai[i][j][k][1] += root_sum_sq[i][j][k];								//Eq. (S34)
			dgphi_dai[i][j][k][1] += snu[j] * shmu[i][j][k];							//Eq. (S36)
		}		
	}

	//calculate cavity volume - LV
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) {
		c_in = chmu[i][j][0];
		dcoshmu = (c_in - 1.) / (n_volint - 1);
		for (k = 0; k < n_volint; k++) {
			coshmu = 1. + k * dcoshmu;
			integrand[i][j][k] = dcoshmu * acube * snu[j] * (DSQR(coshmu) - cnu_sq[j]);
		}
	}
	Vlv = 2. * tripleIntegral(integrand, dmu_ones, dnu, phi_s, 0, n_volint, Nnu, 0, 1);
	Vlv += tripleIntegral(integrand, dmu_ones, dnu, dphi, 0, n_volint, Nnu, 1, Nphi);

	//calculate cavity volume - RV, eq. (S67)
	for (j = 0; j < Nnu; j++) for (i = 1; i < Nphi; i++) {
		c_mid1 = chmu[i][j][Nmu1 - 1];
		c_mid2 = chmu[i][j][Nmu1];
		dcoshmu = (c_mid2 - c_mid1) / (n_volint - 1);
		for (k = 0; k < n_volint; k++) {
			coshmu = c_mid1 + k * dcoshmu;
			integrand[i][j][k] = dcoshmu * acube * snu[j] * (DSQR(coshmu) - cnu_sq[j]);
		}
	}
	Vrv = tripleIntegral(integrand, dmu_ones, dnu, dphi, 0, n_volint, Nnu, 1, Nphi);		// RV
	if (t < 0.) printf(" Vref_lv = %f, Vref_rv = %f\n", Vlv, Vrv);
	
	//derivatives of LV cavity volume
	eta_lv_ai[1] = 0.;
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++)
		integrand2[i][j] = a0sqr * ((1. - DCUB(chin0))/3. + (chin0 - 1.) * cnu_sq[j]) * snu[j];	//Eq. (S71)
	eta_lv_ai[2] = 2. * doubleIntegral(integrand2, dnu, phi_s, Nnu, 0, 1);
	eta_lv_ai[2] += doubleIntegral(integrand2, dnu, dphi, Nnu, 1, Nphi);
	eta_lv_ai[3] = 0.;
	if (odeSystem == "LVRVmodel2024") {
		eta_lv_ai[4] = 0.;
		for (j = 0; j < Nnu; j++) for (i = 1; i < Nphi; i++) integrand2[i][j] = jjnu[j] * ggphi0[i] * snu[j];	//Eq. (S73)
		eta_rv_ai[4] = doubleIntegral(integrand2, dnu, dphi, Nnu, 1, Nphi);
	}

	//derivatives of RV cavity volume
	eta_rv_ai[1] = 0.;
	eta_rv_ai[2] = 0.;
	eta_rv_ai[3] = 0.;
	if (odeSystem == "LVRVmodel2024" && septalDisp == "on") {
		eta_lv_ai[5] = eta_rv_ai[4];
		eta_rv_ai[5] = -eta_lv_ai[5];
	}
}
