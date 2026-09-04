// setuparrays.cpp
// TWS, March 2024

#include <stdio.h>
#include <math.h>
#include "nrutil.h"

void setuparrays(int Nmu1, int Nmu2, int Nmu, int Nnu, int Nphi, int Npar, int Nr, int Nmat, int Nvars, int Nt, int Nstore) {
	extern double *phi0, *ggphi0, *ggphi0p, ***mu0, *nu, *jjnu, *jjnup, *hhmu0;
	extern double **dmu1, **dmu2, **dmu_ones;
	extern double **alpha_mat_aij, *eta_lv_ai, *eta_rv_ai, *kappa_ai, *X, *f;
	extern double **J, **J1, **muout;
	extern double ***omega;
	extern double ***cps, ***sps, ***cps_sq, ***sps_sq; 
	extern double *snu, *cnu, ***shmu0, ***chmu0, *dphi_dnu;
	extern double *snu_sq, *cnu_sq, ***sum_sq0, ***root_sum_sq0;
	extern double ***ell0, ***gmu0, ***gnu0, ***gphi0, ***integratingFactor;
	extern double *r0_la, *r0_rv, *r0_ra;
	extern double *y, *ytest;

	//storage of results
	extern double *tvec, *a1vec, *a2vec, *a3vec, *a4vec, *a5vec, *a6vec, *a7vec, **Y, **store;
	extern double ***ch_store, ***sh_store, ***sh_store2;

	extern double ***shmu, ***chmu;
	extern double ***sh_sq, ***ch_sq, ***sum_sq, ***root_sum_sq;
	extern double ***dmu_dmu0, ***dmu_dphi0, ***dmu_dnu, ****dmu_dai;
	extern double ****ddmu_dnudai, ****ddmu_dmu0dai, ****ddmu_dphi0dai;
	extern double *dff2_dai, *dff3_dai;
	extern double ***gmu, ***gphi, ****dgmu_dai, ****dgphi_dai;

	// Deformation gradient tensor (symmetric)
	extern double*** F11, *** F12, *** F13, *** F22, *** F32, *** F33;
	extern double*** Frefinv11, *** Frefinv12, *** Frefinv13, *** Frefinv22, *** Frefinv33;
	extern double ****dF11_dai, ****dF12_dai, ****dF13_dai, ****dF22_dai, ****dF32_dai, ****dF33_dai;

	// Cauchy deformation tensor (symmetric)
	extern double ***C11, ***C12, ***C13, ***C22, ***C23, ***C33;
	extern double ***Cinv11, ***Cinv12, ***Cinv13, ***Cinv22, ***Cinv23, ***Cinv33;

	// Green strain (symmetric)
	extern double ****dE11_dai, ****dE12_dai, ****dE13_dai, ****dE22_dai, ****dE23_dai, ****dE33_dai;
	extern double ***Ess, ***Enn, ***Eff, ***Esn, ***Enf, ***Esf;

	// Fiber strain
	extern double ***ell, ***ell_ref, ***eps_f, ***eps_fed;
	extern double ****dell_dai;

	// Fiber rotation matrix in reference configuration
	extern double*** Qs1, *** Qn1, *** Qf1, *** Qs2, *** Qn2, *** Qf2, *** Qs3, *** Qn3, *** Qf3;

	// 2nd Piola-Kirchhoff viscous stress (symmetric)
	extern double ****Sv11_ai, ****Sv12_ai, ****Sv13_ai, ****Sv22_ai, ****Sv23_ai, ****Sv33_ai;

	// 2nd Piola-Kirchhoff elastic stress (symmetric)
	extern double ***Se11, ***Se12, ***Se13, ***Se22, ***Se23, ***Se33;

	// Active fiber stress
	extern double ***Sf_ff_const, ***Sf_temp;;
	extern double ****Sf_ff_ai;
	extern double ***Sf11_const, ***Sf12_const, ***Sf13_const, ***Sf22_const, ***Sf23_const, ***Sf33_const;
	extern double ****Sf11_ai, ****Sf12_ai, ****Sf13_ai, ****Sf22_ai, ****Sf23_ai, ****Sf33_ai;

	// Total work variables
	extern double ****Se, ****Sv_prev, ****Sv, ****Sf_prev, ****Sf;

	// Active stress variables
	extern double ***SigmaF, ****SigmaF_ai;

	// Spherical chamber variables
	extern double *r, *kappaj_integrand, *chij_integrand;

	// System variables
	extern double *dy;

	// Integration
	extern double ***integrand, **integrand2, *integrand3, *integrand1;

	// Stress pullback to axisymmetric configuration
	extern double*** S0_11, *** S0_12, *** S0_13, *** S0_22, *** S0_23, *** S0_33;
	extern double**** S0_11_ai, **** S0_12_ai, **** S0_13_ai, **** S0_22_ai, **** S0_23_ai, **** S0_33_ai;

	// Volume computation 
	extern double *vol_integrand, **vIntegrand_dai;
/////////////////////////////////////////////////////////////////////////////////////////
	dy = dvector(0, Nvars - 1);
	y = dvector(0, Nvars - 1);
	ytest = dvector(0, Nvars - 1);
	
	phi0 = dvector(0, Nphi - 1);
	ggphi0 = dvector(0, Nphi - 1);
	ggphi0p = dvector(0, Nphi - 1);

	alpha_mat_aij = dmatrix(1, Npar, 1, Npar);
	kappa_ai = dvector(1, Npar);
	eta_lv_ai = dvector(1, Npar);
	eta_rv_ai = dvector(1, Npar);
	dff2_dai = dvector(1, Npar);
	dff3_dai = dvector(1, Npar);

	//Linear system
	X = dvector(1, Nmat);
	f = dvector(1, Nmat);
	J = dmatrix(1, Nmat, 1, Nmat);
	J1 = dmatrix(1, Nmat, 1, Nmat);

	omega = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	cps = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	sps = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	cps_sq = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	sps_sq = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	mu0 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	dmu1 = dmatrix(0, Nphi - 1, 0, Nnu - 1); 
	dmu2 = dmatrix(0, Nphi - 1, 0, Nnu - 1);
	dmu_ones = dmatrix(0, Nphi - 1, 0, Nnu - 1);
	muout = dmatrix(0, Nphi - 1, 0, Nnu - 1);
	hhmu0 = dvector(0, Nmu - 1);
	nu = dvector(0,Nnu-1);
	snu = dvector(0, Nnu-1);
	cnu = dvector(0, Nnu - 1);
	shmu0 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	chmu0 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	snu_sq = dvector(0, Nnu - 1);
	cnu_sq = dvector(0, Nnu - 1);
	jjnu = dvector(0, Nnu - 1);
	jjnup = dvector(0, Nnu - 1);
	dphi_dnu = dvector(0, Nnu - 1);
	sum_sq0 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	root_sum_sq0 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	ell0 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	ell_ref = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	gmu0 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	gnu0 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	gphi0 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	integratingFactor = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);

	// Initialize subelements
	shmu = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	chmu = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	sh_sq = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	ch_sq = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	sum_sq = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	root_sum_sq = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	gmu = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	gphi = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	dmu_dmu0 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	dmu_dphi0 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	dmu_dnu = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	dmu_dai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	ddmu_dnudai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	ddmu_dmu0dai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	ddmu_dphi0dai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	dgmu_dai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	dgphi_dai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);

	// deformation gradient tensor
	F11 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	F12 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	F13 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	F22 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	F32 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	F33 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);

	Frefinv11 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Frefinv12 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Frefinv13 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Frefinv22 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Frefinv33 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);

	dF11_dai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	dF12_dai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	dF13_dai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	dF22_dai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	dF32_dai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	dF33_dai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);

	// Cauchy deformation tensor
	C11 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	C12 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	C13 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	C22 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	C23 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	C33 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Cinv11 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Cinv12 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Cinv13 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Cinv22 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Cinv23 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Cinv33 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	
	// Green strain tensor derivatives
	dE11_dai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	dE12_dai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	dE13_dai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	dE22_dai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	dE23_dai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	dE33_dai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	Ess = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Enn = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Eff = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Esn = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Enf = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Esf = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);

	// Fiber rotation matrix in reference configuration
	Qs1 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Qn1 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Qf1 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Qs2 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Qn2 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Qf2 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Qs3 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Qn3 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	Qf3 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	
	// Fiber strain
	ell = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	dell_dai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	eps_f = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	eps_fed = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);

	// 2nd Piola-Kirchhoff viscous stress
	Sv11_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	Sv12_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	Sv13_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	Sv22_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	Sv23_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	Sv33_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);

	// 2nd Piola-Kirchhoff elastic stress 
	Se11 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Se12 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Se13 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Se22 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Se23 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Se33 = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);

	// active fiber strain vectors
	Sf_ff_const = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Sf_temp = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Sf_ff_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);

	Sf11_const = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Sf12_const = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Sf13_const = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Sf22_const = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Sf23_const = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	Sf33_const = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);

	Sf11_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	Sf12_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	Sf13_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	Sf22_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	Sf23_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);
	Sf33_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);

	// Stress pullback to axisymmetric configuration
	S0_11 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	S0_12 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	S0_13 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	S0_22 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	S0_23 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	S0_33 = d3tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);

	S0_11_ai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	S0_12_ai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	S0_13_ai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	S0_22_ai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	S0_23_ai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	S0_33_ai = d4tensor(0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);

	// Active stress variables
	SigmaF = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	SigmaF_ai = d4tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);

	// Initialize equilibrium integrands
	integrand = d3tensor(0,Nphi-1,0,Nnu-1,0,Nmu-1);
	integrand2 = dmatrix(0,Nphi-1,0,Nnu-1);
	integrand3 = dvector(0,Nphi-1);
	integrand1 = dvector(0,Nmu-1);

	// Initialize volume derivative
	vIntegrand_dai = dmatrix(0,Nnu-1,1, Npar);
	vol_integrand = dvector(0,Nnu-1);

	// Initialize spherical chamber variables
	r = dvector(0, Nr - 1);
	r0_la = dvector(0, Nr - 1);
	r0_rv = dvector(0, Nr - 1);
	r0_ra = dvector(0, Nr - 1);
	kappaj_integrand = dvector(0,Nr-1);
	chij_integrand = dvector(0,Nr-1);

	// Set up solution storage vector
	Y = dmatrix(0, Nvars - 1, 0, Nt - 1);
	// Set up storage for other variables
	store = dmatrix(0,Nstore-1,0,Nt-1);

	// Set up storage for a1,a2,a3,a4,a5,a6,a7
	a1vec = dvector(0, Nt - 1);
	a2vec = dvector(0, Nt - 1);
	a3vec = dvector(0, Nt - 1);
	a4vec = dvector(0, Nt - 1);
	a5vec = dvector(0, Nt - 1);
	a6vec = dvector(0, Nt - 1);
	a7vec = dvector(0, Nt - 1);

	// Set up storage for the time variable
	tvec = dvector(0, Nt - 1);

	//storage for short-axis and long-axis plots
	ch_store = d3tensor(0, Nnu - 1, 1, 4, 1, 8);
	sh_store = d3tensor(0, Nnu - 1, 1, 4, 1, 8);
	sh_store2 = d3tensor(0, Nphi - 1, 1, 4, 1, 4);
}
