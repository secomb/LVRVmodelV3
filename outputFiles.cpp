// outputFiles.cpp
// TWS, November 2019
#define _CRT_SECURE_NO_WARNINGS	//TWS
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

void outputFiles(int run) {		//TWS August 2025
	extern int Nt, Nnu, Nr, Nstore, Nparam, Nopt, Ntarget;
	extern int* ParamVar;
	extern double muin0, mumid0, muout0, a0, a40, nu_up, Ls0, Lsmax, Lsw, Tc, Psv0, km, km_rv;
	extern double rin0_la, rout0_la, rin0_rv, rout0_rv, rin0_ra, rout0_ra;
	extern double* tvec, * a1vec, * a2vec, * a3vec, * a4vec, * a5vec, * a6vec, * a7vec;
	extern double* targetVector, * targetWeight, * outputVector, * chisq;
	extern double* xx, * xxlb, * xxub;
	extern double* Param0, * Param;
	extern double** store, ** Y;
	extern std::string ParamName[128];	//note: zero-based arrays, include extra element
	extern std::string TargetName[49];
	extern string Optimize, Normalize, odeSystem;

	int i, k, Nt1 = Nt - 1;
	double result, target, residual, res_sum, weight;
	char fname[80];
	FILE* ofp;
	fstream file1;

	// Write Matlab output file
	sprintf(fname, "outputFileA%03i.m", run);
	fopen_s(&ofp, fname, "w");
	fprintf(ofp, "t = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", tvec[i]);
		else fprintf(ofp, "%g ]; \n", tvec[i]);
	}
	fprintf(ofp, "a1 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a1vec[i]);
		else fprintf(ofp, "%g ]; \n", a1vec[i]);
	}
	fprintf(ofp, "a2 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a2vec[i]);
		else fprintf(ofp, "%g ]; \n", a2vec[i]);
	}
	fprintf(ofp, "a3 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a3vec[i]);
		else fprintf(ofp, "%g ]; \n", a3vec[i]);
	}
	fprintf(ofp, "a4 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a4vec[i]);
		else fprintf(ofp, "%g ]; \n", a4vec[i]);
	}
	fprintf(ofp, "a5 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a5vec[i]);
		else fprintf(ofp, "%g ]; \n", a5vec[i]);
	}
	fprintf(ofp, "a6 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a6vec[i]);
		else fprintf(ofp, "%g ]; \n", a6vec[i]);
	}
	fprintf(ofp, "a7 = [ ");
	for (i = 0; i < Nt1; i++) {
		if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
		if (i < Nt1 - 1) fprintf(ofp, "%g, ", a7vec[i]);
		else fprintf(ofp, "%g ]; \n", a7vec[i]);
	}
	fclose(ofp);

	sprintf(fname, "outputFileStore%03i.m", run);
	fopen_s(&ofp, fname, "w");
	fprintf(ofp, "store = [ ");
	for (k = 0; k < Nstore; k++) {
		for (i = 0; i < Nt1; i++) {
			if (i % 10 == 0 && i > 0) fprintf(ofp, "...\n");
			if (i < Nt1 - 1) fprintf(ofp, "%g, ", store[k][i]);
			else if (k < Nstore - 1) fprintf(ofp, "%g; ...\n", store[k][i]);
			else fprintf(ofp, "%g ]; \n", store[k][i]);
		}
	}
	fclose(ofp);

	//Write combined text file - May 2022
	sprintf(fname, "outputFile_All%03i.txt", run);
	fopen_s(&ofp, fname, "w");
	fprintf(ofp, "t a1 a2 a3 a4 a5 a6 a7");
	fprintf(ofp, " Vtot Vlv Plv Vla Pla Vrv Prv Vra Pra Pperi Psa Psp Psv Ppa Ppp Ppv q_la q_lv q_ra q_rv q_LVAD q_sp q_sv null q_pp q_pv Zmiv Zaov");
	fprintf(ofp, " At    At_la null null null Ztcv Zpuv b-a_length dt(length) P_echo Q_trans P_ao_0 Q_lv_0");
	fprintf(ofp, " vel_lv vel_la vel_rv vel_ra LV_diam FV_fiberstress Vaorta\n");
	for (i = 0; i < Nt1; i++) {
		fprintf(ofp, "%g %g %g %g %g %g %g %g", tvec[i], a1vec[i], a2vec[i], a3vec[i], a4vec[i], a5vec[i], a6vec[i], a7vec[i]);
		for (k = 0; k < Nstore; k++) fprintf(ofp, " %g", store[k][i]);
		fprintf(ofp, "\n");
	}
	fclose(ofp);

	sprintf(fname, "outputList%03i.txt", run);
	file1.open(fname, ios_base::out);
	if (Optimize == "on") {
		if (Optimize == "on") file1 << "**** Outputs normalized to target values in optimization ****" << endl;
		file1 << "chi2 = " << chisq[0] << endl;
		file1 << "**** Optimized multipliers ****" << endl;
		for (i = 1; i <= Nopt; i++) file1 << " xx[" << i << "] = " << xx[i];
		file1 << endl;
		for (i = 1; i <= Nopt; i++) file1 << " lb[" << i << "] = " << xxlb[i];
		file1 << endl;
		for (i = 1; i <= Nopt; i++) file1 << " ub[" << i << "] = " << xxub[i];
		file1 << endl;
	}
	file1 << "**** Input parameters ****" << endl;
	for (i = 1; i <= Nparam; i++) {
		file1 << ParamName[i] << " = " << Param[i];
		if (ParamVar[i] && Optimize == "on") file1 << " optimized by xx[" << ParamVar[i] << "]";
		file1 << endl;
	}
	file1 << "**** Output variables ****" << endl;
	for (i = 1; i <= Ntarget; i++) {
		target = targetVector[i];
		if (i == 23 || i == 24) result = outputVector[i] / (1. - outputVector[i]);
		else result = outputVector[i];
		file1 << TargetName[i] << " = " << result;
		if (targetWeight[i] > 0. && Optimize == "on") {
			if (i == 23 || i == 24) result = targetVector[i] / (1. - targetVector[i]);
			else result = targetVector[i];
			file1 << ", target = " << result;
			if (Normalize == "on") file1 << ", weight = " << targetWeight[i] * DSQR(targetVector[i]);
			else file1 << ", weight = " << targetWeight[i];
		}
		file1 << endl;
	}
	file1.close();

	//if (Optimize == "on") {
		sprintf(fname, "residuals%03i.txt", run);
		ofp = fopen(fname, "w");
		res_sum = 0.;
		fprintf(ofp, "     Variable   Target    Fitted     Residual   Weight\n");
		for (i = 1; i <= Ntarget; i++) if (targetWeight[i] > 0.) {
			target = targetVector[i];
			result = outputVector[i];
			residual = result - target;
			res_sum += DSQR(residual) * targetWeight[i];
			weight = targetWeight[i];
			if (Normalize == "on") {
				residual /= target;
				weight *= DSQR(target);
			}
			fprintf(ofp, "%13s %9.4f %9.4f %9.4f %9.4f\n",
				TargetName[i], target, result, residual, weight);
		}
		fprintf(ofp, "Objective function = %9.4f\n", res_sum);
		fclose(ofp);
	//}

		sprintf(fname, "stroke_work%03i.txt", run);
		ofp = fopen(fname, "w");
		fprintf(ofp, "%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",
			km, km_rv, Psv0, outputVector[1], outputVector[45], outputVector[5], outputVector[46]);
		fclose(ofp);
	
	printf("PAOMEAN = %g\n", outputVector[36]);
	printf("PSVMEAN = %g\n", outputVector[48]);
	printf("PVMEAN = %g\n", outputVector[28]);
}
