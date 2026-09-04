//mathFunctions.cpp
//edited by TWS, November 2019, February 2024
#include <math.h>
#include <stdio.h>
#include "nrutil.h"

double simpsons(double* f, double dx, int Npoints) {
	double out = 0.;
	int k;
	if (Npoints >= 9) {	// Use simpsons rule if there are sufficient points
		out = 17 * f[0] + 59 * f[1] + 43 * f[2] + 49 * f[3]
			+ 49 * f[Npoints - 4] + 43 * f[Npoints - 3] + 59 * f[Npoints - 2] + 17 * f[Npoints - 1];
		for (k = 4; k < Npoints - 4; k++) out = out + 48 * f[k];
		out = out * dx / 48;
	}
	else {// Use trapezoid rule if there are fewer points, works even if Npoints = 1
		for (k = 1; k < Npoints - 1; k++) out = out + f[k];
		out = out * 2 + f[0] + f[Npoints - 1];
		out = out * dx / 2;
	}
	return out;
}

double tripleIntegral(double*** f, double **dmu, double dnu, double dphi, int Nmus, int Nmue, int Nnu, int Nphis, int Nphie) {
	extern double **integrand2, *integrand3, *integrand1;

	int k, j, i, Nmud = Nmue - Nmus, Nphid = Nphie - Nphis;
	double value;

	for (i = 0; i < Nphid; i++) {
		for (j = 0; j < Nnu; j++) {
			for (k = 0; k < Nmud; k++) integrand1[k] = f[Nphis + i][j][Nmus + k];
			integrand2[i][j] = simpsons(integrand1, dmu[i][j], Nmud);
		}
		integrand3[i] = simpsons(integrand2[i], dnu, Nnu);
	}
	value = simpsons(integrand3, dphi, Nphid);
	return value;
}

double doubleIntegral(double** f, double dnu, double dphi, int Nnu, int Nphis, int Nphie) {
	extern double *integrand3;

	int i, Nphid = Nphie - Nphis;
	double value;

	for (i = 0; i < Nphid; i++) integrand3[i] = simpsons(f[Nphis + i], dnu, Nnu);
	value = simpsons(integrand3, dphi, Nphid);
	return value;
}

double Integral3(double ***f) {			//integrate over entire wall domain
	extern int Nmu1, Nmu, Nnu, Nphi;
	extern double **dmu1, **dmu2, dnu, dphi, phi_s;
	double value;
	
	value = tripleIntegral(f, dmu1, dnu, dphi, 0, Nmu1, Nnu, 1, Nphi);
	value += 2. * tripleIntegral(f, dmu1, dnu, phi_s, 0, Nmu1, Nnu, 0, 1);
	value += tripleIntegral(f, dmu2, dnu, dphi, Nmu1, Nmu, Nnu, 1, Nphi);
	value += 2. * tripleIntegral(f, dmu2, dnu, phi_s, Nmu1, Nmu, Nnu, 0, 1);
	return value;
}
