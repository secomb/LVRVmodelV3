// referenceConfig.cpp
// TWS, March 2024
// This is configuration 1, "reference" including RV cavity
#include <math.h>
#include <stdio.h>
#include <string>
#include "nrutil.h"

using namespace std;

#include <chrono>

void updateMapping(double t);
void updateDefGradTensor(double t);
void updateCauchyTensor();
void updateGreenStrain(double t);
void updateFiberStrain(double t);

void referenceConfig() {
	extern int Nmu, Nnu, Nphi;
	extern double a1, a2, a3, a4, a5;
	extern double ***Ess, ***Enn, ***Eff, ***Esn, ***Enf, ***Esf;
	extern double ***ell, ***ell_ref;
	extern double*** sps, *** cps, *** cps_sq, *** sps_sq;
	extern double*** F11, *** F12, *** F13, *** F22, *** F32, *** F33;
	extern double*** Frefinv11, *** Frefinv12, *** Frefinv13, *** Frefinv22, *** Frefinv33;
	extern double*** Qs1, *** Qn1, *** Qf1,*** Qs2, *** Qn2, *** Qf2, *** Qs3, *** Qn3, *** Qf3;

	int i, j, k;
	double es1, en1, ef1, es2, en2, ef2, es3, en3, ef3, esm, enm, efm;

	a1 = 0.;
	a2 = 0.;
	a3 = 0.;
	a4 = 0.;
	a5 = 0.;
	updateMapping(-1.);		//gives configuration 0 to 1 mapping, i.e. axiysmmetric to reference
	updateDefGradTensor(-1.);
	// store values of inverse DefGradTensor in reference configuration
	// no torsion in reference configuration, so F32 = 0
	// inverse of | a b c | is | 1/a -b/ad -c/af|
	//            | 0 d 0 |    | 0    1/d    0  |
	//            | 0 0 f |    | 0    0     1/f |
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
		Frefinv11[i][j][k] = 1. / F11[i][j][k];
		Frefinv12[i][j][k] = -F12[i][j][k] / (F11[i][j][k] * F22[i][j][k]);
		Frefinv13[i][j][k] = -F13[i][j][k] / (F11[i][j][k] * F33[i][j][k]);
		Frefinv22[i][j][k] = 1. / F22[i][j][k];
		Frefinv33[i][j][k] = 1. / F33[i][j][k];
	}
	updateCauchyTensor();

	//calculate rotation matrix from (mu,nu,phi) to (s,n,f) in reference configuration 1
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) {
		//f direction, eq. (S65)
		ef1 = F12[i][j][k] * sps[i][j][k] - F13[i][j][k] * cps[i][j][k];
		ef2 = F22[i][j][k] * sps[i][j][k];
		ef3 = -F33[i][j][k] * cps[i][j][k];
		efm = sqrt(SQR(ef1) + SQR(ef2) + SQR(ef3));
		//s direction
		es1 = 0.;
		es2 = -ef3;
		es3 = ef2;
		esm = sqrt(SQR(es2) + SQR(es3));
		//n direction
		en1 = SQR(ef2) + SQR(ef3);
		en2 = -ef2 * ef1;
		en3 = -ef3 * ef1;
		enm = sqrt(SQR(en1) + SQR(en2) + SQR(en3));
		//Q matrix, eqn. (65.1)
		Qs1[i][j][k] = es1 / esm;
		Qs2[i][j][k] = es2 / esm;
		Qs3[i][j][k] = es3 / esm;
		Qn1[i][j][k] = en1 / enm;
		Qn2[i][j][k] = en2 / enm;
		Qn3[i][j][k] = en3 / enm;
		Qf1[i][j][k] = ef1 / efm;
		Qf2[i][j][k] = ef2 / efm;
		Qf3[i][j][k] = ef3 / efm;
	}
	updateGreenStrain(-1.);		//t < 0, signal that this is reference state
	updateFiberStrain(-1.);		//t < 0, signal that this is reference state

	//store ell values in configuration 1 (reference) 
	for (i = 0; i < Nphi; i++) for (j = 0; j < Nnu; j++) for (k = 0; k < Nmu; k++) ell_ref[i][j][k] = ell[i][j][k];
}

