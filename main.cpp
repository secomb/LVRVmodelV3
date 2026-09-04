// Main program
// edited by TWS, November 2019
// LVRV version, TWS 2024, with parameter estimation, August 2025
// This version allows variable muout, to represent RV hypertrophy
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <iomanip>
#include "nrutil.h"

using namespace std;

void readParameters(int run);
void readTargets(int run);
void LevMarq(double* targetVector, double* targetWeight, double* outputVector, int Ntarget,
	double* xx, double* xxlb, double* xxub, int Nopt, double * chisq, double tol, int Niter, int *ia, int run);
void setuparrays(int Nmu1, int Nmu2, int Nmu, int Nnu, int Nphi, int Npar, int Nr, int Nmat, int Nvars, int Nt, int Nstore);
void freearrays(int Nmu1, int Nmu2, int Nmu, int Nnu, int Nphi, int Npar, int Nr, int Nmat, int Nvars, int Nt, int Nstore);
int aorta(double* Param0, int* ParamVar, int Nparam, double* xx);
void runModel();
void makeGraphs(int run);
void outputFiles(int run);
void updateParams(double* Param, int Nparam);
void updateTargets(double* outputVector, int Ntarget);

int Nvars, Nmat, iRecord, Nstore, Nr, Ncycles, Nt, Nta, Ncycles_show;
double Tc, dtRecord, dtMax;
double Ta, Tca, Tca_shift, f_rv;
double *y_init, *y, *ytest, *dy;
//storage of results
double *tvec, *a1vec, *a2vec, *a3vec, *a4vec, *a5vec, *a6vec, *a7vec, **Y, **store;
double ***ch_store, ***sh_store, ***sh_store2;

// LVRV parameters
int Nmu1, Nmu2, Nmu, Nnu, Nphi, Npar;
double dnu, dphi;
double nu_up, a0, a40, a60, muin0, mumid0, muout0, phi_s, phi_t, phi_d, alpha;
double psi_in_b0, psi_out_b0;
double Ls0, Lsmax, Lsw;
double km, kav, kdfac;
double kv, c1, bff, bfx, bxx;

// Other chambers
double km_la, kav_la, kv_la, br_la, c1_la, bperp_la, rin0_la, rout0_la;
double km_rv, kav_rv, kv_rv, br_rv, c1_rv, bperp_rv, rin0_rv, rout0_rv;
double bff_rv, bxx_rv, bfx_rv;
double km_ra, kav_ra, kv_ra, br_ra, c1_ra, bperp_ra, rin0_ra, rout0_ra;
double eta_la, kappa_la, chi_la;
double eta_rv, kappa_rv, chi_rv;
double eta_ra, kappa_ra, chi_ra;

double *phi0, *ggphi0, *ggphi0p, ***mu0, *nu, *jjnu, *jjnup, *hhmu0;
double **dmu1, **dmu2, **dmu_ones;
double ***omega;
double ***cps, ***sps, ***cps_sq, ***sps_sq; 
double *snu, *cnu, ***shmu0, ***chmu0, *dphi_dnu, **muout;
double *snu_sq, *cnu_sq, ***sum_sq0, ***root_sum_sq0;
double ***ell0, ***gmu0, ***gnu0, ***gphi0, ***integratingFactor;
double ***gmu, ***gphi, ****dgmu_dai, ****dgphi_dai;
double *r0_la, *r0_ra, *r0_rv;

double *X, *f;
double **J, **J1;
double a1, a2, a3, a4, a5, a6, a7;

// Inputs
double **alpha_mat_aij, **alpha_temp, *eta_lv_ai, *eta_rv_ai, *kappa_ai;

// Lumped model variables
double Vwall_rv, Vwall_lv, Vwall_se, Vlv, Vrv, Vla, Vra;
double Csa, Csv, Cpa, Cpv;
double Rsp, Rsv, Rpa, Rpp, Rpv, Rra;
double Ppa, Ppp, Psp, Pra, Ppv, Psa, Prv, Psv, Plv, Pla;
double q_la, q_lv, q_ra, q_rv, q_LVAD = 0.;
double q_sp, q_sv, q_pp, q_pv;

double ***shmu, ***chmu;
double ***sh_sq, ***ch_sq, ***sum_sq, ***root_sum_sq;
double ***dmu_dmu0, ***dmu_dphi0, ***dmu_dnu, ****dmu_dai;
double ****ddmu_dnudai, ****ddmu_dmu0dai, ****ddmu_dphi0dai;
double *dff2_dai, *dff3_dai;

// Deformation gradient tensor (asymmetric)
double*** F11, *** F12, *** F13, *** F22, *** F32, *** F33;
double*** Frefinv11, *** Frefinv12, *** Frefinv13, *** Frefinv22, *** Frefinv33;
double**** dF11_dai, **** dF12_dai, **** dF13_dai, **** dF22_dai, **** dF32_dai, **** dF33_dai;

// Cauchy deformation tensor (symmetric)
double ***C11, ***C12, ***C13, ***C22, ***C23, ***C33;
double ***Cinv11, ***Cinv12, ***Cinv13, ***Cinv22, ***Cinv23, ***Cinv33;

// Green strain (symmetric)
double ****dE11_dai, ****dE12_dai, ****dE13_dai, ****dE22_dai, ****dE23_dai, ****dE33_dai;
double ***Ess, ***Enn, ***Eff, ***Esn, ***Enf, ***Esf;

// Fiber rotation matrix in reference configuration
double*** Qs1, *** Qn1, *** Qf1, *** Qs2, *** Qn2, *** Qf2, *** Qs3, *** Qn3, *** Qf3;

// Fiber strain
double ***ell, ***ell_ref, ***eps_f, ***eps_fed;
double ****dell_dai;

// 2nd Piola-Kirchhoff viscous stress (symmetric)
double ****Sv11_ai, ****Sv12_ai, ****Sv13_ai, ****Sv22_ai, ****Sv23_ai, ****Sv33_ai;

// 2nd Piola-Kirchhoff elastic stress (symmetric)
double ***Se11, ***Se12, ***Se13, ***Se22, ***Se23, ***Se33;

// Stress pullback to axisymmetric configuration
double*** S0_11, *** S0_12, *** S0_13, *** S0_22, *** S0_23, *** S0_33;
double**** S0_11_ai, **** S0_12_ai, **** S0_13_ai, **** S0_22_ai, **** S0_23_ai, **** S0_33_ai;

// Activation function
double At_la, At;

// Active fiber stress
double ***Sf_ff_const, ***Sf_temp;
double ****Sf_ff_ai;
double ***Sf11_const, ***Sf12_const, ***Sf13_const, ***Sf22_const, ***Sf23_const, ***Sf33_const;
double ****Sf11_ai, ****Sf12_ai, ****Sf13_ai, ****Sf22_ai, ****Sf23_ai, ****Sf33_ai;

// Active stress variables
double ***SigmaF, ****SigmaF_ai;

// Spherical chamber variables
double *r, *kappaj_integrand, *chij_integrand;

// Integration
double ***integrand, **integrand2, *integrand3, *integrand1;

// Volume computation 
double *vol_integrand, **vIntegrand_dai;

// Valve variables
double B_lv, B_la, B_rv, B_ra, L_lv, L_la, L_rv, L_ra, L_pv, L_sv;
double Zmiv, Zaov, Ztcv, Zpuv;
double Kaovo, Kaovc, Kmivo, Kmivc, Kpuvo, Kpuvc, Ktcvo, Ktcvc;
double VA_aovo, VA_aovc, VA_mivo, VA_mivc, VA_puvo, VA_puvc, VA_tcvo, VA_tcvc, VL_aov, VL_miv, VL_puv, VL_tcv, rho;

//aorta simulation
double dx, tpulse, Gaorta0, Gaorta1, Z0start, Z0end, Laorta, A00, A01;
double P_echo, Q_trans, P_ao_0, q_lv_0, Vaorta, impedancefac, fac;
double *echo, *trans, *q_lv_hist, *p_ao_hist, *tveca, VELAORTAWAVE, LV_SW, RV_SW;
int Nn;

//External tissue pressure, initial pressures
double Pext, Psp0, Psv0, Ppp0, Ppv0;

//pericardium parameters
double V0_peri, s_peri, Pperi;

//LVAD parameters
double q0_LVAD, Ps_LVAD, t0_LVAD;

//Variables used in optimization
int Nparam, Nopt, Ntarget, NtargetVar, Neval = 0, aortaVar, Niter;
int* ParamVar, * ia;
double tol;
double* Param, * xx, * xxlb, * xxub, * outputVector, * outputVector0, * chisq;
double* Param0, * targetVector, * targetWeight;
double** derivs;
std::string ParamName[128];	//note: zero-based arrays, include extra element
std::string TargetName[49];
string Optimize, Normalize, Skip;
string odeSystem, septalDisp, LVAD, RVhypertrophy, pericardium;

int main(int argc, char* argv[]) {
	int iout, Nfit = 0, j;
	int run = 0, runflag = 0;
	char fname[80];
	FILE *ofp;

	do {
		sprintf(fname, "parameters%03i.txt", run);
		//test whether file exists
		ofp = fopen(fname, "r");
		if (ofp) {
			runflag = 1;
			fclose(ofp);
		}
		else {
			runflag = 0;
			printf("*************************\n");
			printf("No more input files found\n");
			printf("*************************\n");
		}
		if(runflag) {
			printf("*******************\n");
			printf("****** %03i ********\n", run);
			printf("*******************\n");	

			readParameters(run);

			if (Skip != "on") {

				updateParams(Param0, Nparam);

				readTargets(run);

				Param = dvector(1, Nparam);		//use for model parameters after adjustment
				outputVector = dvector(1, Ntarget);
				if (Optimize == "on") {
					outputVector0 = dvector(1, Ntarget);
					derivs = dmatrix(1, Ntarget, 1, Nopt);
					chisq = dvector(0, 1); //for first iteration, need to return two values
				}
				Nmu = Nmu1 + Nmu2;
				Nta = aorta(Param0, ParamVar, Nparam, xx);
				Nt = 1 + (int)(Tc * Ncycles_show / dtRecord);	// only last cycle if Ncycles_show = 1
				Nstore = 48;
				Nvars = 21;	//no. of variables in system of ODEs
				if (odeSystem == "LVRVmodel2024") {
					if (septalDisp == "on") {
						Npar = 5;	//no. of ai(t) variables for LVRV
						Nmat = 7;	//no. of ai(t) variables for LVRV plus Plv and Prv
					}
					else {
						Npar = 4;	//no. of ai(t) variables for LVRV
						Nmat = 6;	//no. of ai(t) variables for LVRV plus Plv and Prv
					}
				}
				else {
					Npar = 3;   //no. of ai(t) variables for LV
					Nmat = 4;	//no. of ai(t) variables for LV plus Plv
				}
				setuparrays(Nmu1, Nmu2, Nmu, Nnu, Nphi, Npar, Nr, Nmat, Nvars, Nt, Nstore);

				if (Optimize == "on") {
					if (Nopt == 0) printf("*** Error: No optimization variables, cannot optimize\n");
					if (NtargetVar == 0) printf("*** Error: No target variables, cannot optimize\n");
				}
				if (Nopt == 0 || NtargetVar == 0) Optimize = "off";

				if (Optimize == "on") {
					for (j = 1; j <= Nopt; j++) if (ia[j]) Nfit++;
					printf("Optimizing, %i variable parameters, %i target variables\n", Nfit, NtargetVar);
					for (iout = 1; iout <= Ntarget; iout++) outputVector[iout] = 0.;

					LevMarq(targetVector, targetWeight, outputVector, Ntarget,
						xx, xxlb, xxub, Nopt, chisq, tol, Niter, ia, run);
				}
				else {
					readTargets(run);

					runModel();

					updateTargets(outputVector, Ntarget);
				}
				makeGraphs(run);

				outputFiles(run);

				freearrays(Nmu1, Nmu2, Nmu, Nnu, Nphi, Npar, Nr, Nmat, Nvars, Nt, Nstore);
				free_dvector(echo, 0, Nta - 1);
				free_dvector(trans, 0, Nta - 1);
				free_dvector(q_lv_hist, 0, Nta - 1);
				free_dvector(p_ao_hist, 0, Nta - 1);
				free_dvector(tveca, 0, Nta - 1);

				free_dvector(Param, 1, Nparam);
				free_dvector(Param0, 1, Nparam);
				free_ivector(ParamVar, 1, Nparam);
				free_dvector(outputVector, 1, Ntarget);

				if (Optimize == "on") {
					free_dvector(xx, 1, Nopt);
					free_dvector(xxlb, 1, Nopt);
					free_dvector(xxub, 1, Nopt);
					free_dvector(outputVector0, 1, Ntarget);
					free_dvector(targetVector, 1, Ntarget);
					free_dvector(targetWeight, 1, Ntarget);
					free_dvector(chisq, 0, 1);
					free_ivector(ia, 1, Nopt);
				}
			}
			else printf("*** File skipped **\n");
		}
		run++;
	} while (runflag);
	return(0);
}
