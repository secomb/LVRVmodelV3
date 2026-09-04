// freearrays.cpp
// TWS, March 2024

#include <stdio.h>
#include <math.h>
#include "nrutil.h"

void freearrays(int Nmu1, int Nmu2, int Nmu, int Nnu, int Nphi, int Npar, int Nr, int Nmat, int Nvars, int Nt, int Nstore) {
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
	extern double ***F11, ***F12, ***F13, ***F22, ***F32, ***F33;
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

	// 2nd Piola-Kirchhoff viscous stress (symmetric)
	extern double ****Sv11_ai, ****Sv12_ai, ****Sv13_ai, ****Sv22_ai, ****Sv23_ai, ****Sv33_ai;

	// 2nd Piola-Kirchhoff elastic stress (symmetric)
	extern double ***Se11, ***Se12, ***Se13, ***Se22, ***Se23, ***Se33;

	// Active fiber stress
	extern double ***Sf_ff_const, ***Sf_temp;;
	extern double ****Sf_ff_ai;
	extern double ***Sf11_const, ***Sf12_const, ***Sf13_const, ***Sf22_const, ***Sf23_const, ***Sf33_const;
	extern double ****Sf11_ai, ****Sf12_ai, ****Sf13_ai, ****Sf22_ai, ****Sf23_ai, ****Sf33_ai;

	// Stress pullback to axisymmetric configuration
	extern double*** S0_11, *** S0_12, *** S0_13, *** S0_22, *** S0_23, *** S0_33;
	extern double**** S0_11_ai, **** S0_12_ai, **** S0_13_ai, **** S0_22_ai, **** S0_23_ai, **** S0_33_ai;

	// Active stress variables
	extern double ***SigmaF, ****SigmaF_ai;

	// Spherical chamber variables
	extern double *r, *kappaj_integrand, *chij_integrand;

	// System variables
	extern double *dy;

	// Integration
	extern double ***integrand, **integrand2, *integrand3, *integrand1;


	// Volume computation 
	extern double *vol_integrand, **vIntegrand_dai;
/////////////////////////////////////////////////////////////////////////////////////////
	free_dvector(dy, 0, Nvars - 1);
	free_dvector(y, 0, Nvars - 1);
	free_dvector(ytest, 0, Nvars - 1);
	
	free_dvector(phi0, 0, Nphi - 1);
	free_dvector(ggphi0, 0, Nphi - 1);
	free_dvector(ggphi0p, 0, Nphi - 1);

	free_dmatrix(alpha_mat_aij, 1, Npar, 1, Npar);
	free_dvector(kappa_ai, 1, Npar);
	free_dvector(eta_lv_ai, 1, Npar);
	free_dvector(eta_rv_ai, 1, Npar);
	free_dvector(dff2_dai, 1, Npar);
	free_dvector(dff3_dai, 1, Npar);

	//Linear system	
	free_dvector(X, 1, Nmat);
	free_dvector(f, 1, Nmat);
	free_dmatrix(J, 1, Nmat, 1, Nmat);
	free_dmatrix(J1, 1, Nmat, 1, Nmat);

	free_d3tensor(omega, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(cps, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(sps, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(cps_sq, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(sps_sq, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(mu0, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_dmatrix(dmu1, 0, Nphi - 1, 0, Nnu - 1);
	free_dmatrix(dmu2, 0, Nphi - 1, 0, Nnu - 1);
	free_dmatrix(dmu_ones, 0, Nphi - 1, 0, Nnu - 1);
	free_dmatrix(muout, 0, Nphi - 1, 0, Nnu - 1);
	free_dvector(hhmu0, 0, Nmu - 1);
	free_dvector(nu, 0, Nnu - 1);
	free_dvector(snu, 0, Nnu - 1);
	free_dvector(cnu, 0, Nnu - 1);
	free_d3tensor(shmu0, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(chmu0, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_dvector(snu_sq, 0, Nnu - 1);
	free_dvector(cnu_sq, 0, Nnu - 1);
	free_dvector(jjnu, 0, Nnu - 1);
	free_dvector(jjnup, 0, Nnu - 1);
	free_dvector(dphi_dnu, 0, Nnu - 1);
	free_d3tensor(sum_sq0, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(root_sum_sq0, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(ell0, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(ell_ref, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(gmu0, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(gnu0, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(gphi0, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(integratingFactor, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);

	// Initialize subelements
	free_d3tensor(shmu, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(chmu, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(sh_sq, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(ch_sq, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(sum_sq, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(root_sum_sq, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(gmu, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(gphi, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(dmu_dmu0, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(dmu_dphi0, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(dmu_dnu, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d4tensor(dmu_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(ddmu_dnudai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(ddmu_dmu0dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(ddmu_dphi0dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dgmu_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dgphi_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);

	// deformation gradient tensor
	free_d3tensor(F11, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(F12, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(F13, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(F22, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(F32, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(F33, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d4tensor(dF11_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dF12_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dF13_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dF22_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dF32_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dF33_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);

	// Cauchy deformation tensor
	free_d3tensor(C11, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(C12, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(C13, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(C22, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(C23, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(C33, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Cinv11, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Cinv12, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Cinv13, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Cinv22, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Cinv23, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Cinv33, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);

	// Green strain tensor
	free_d4tensor(dE11_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dE12_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dE13_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dE22_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dE23_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(dE33_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d3tensor(Ess, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Enn, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Eff, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Esn, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Enf, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Esf, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);

	// Fiber strain
	free_d3tensor(ell, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d4tensor(dell_dai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d3tensor(eps_f, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(eps_fed, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);

	// 2nd Piola-Kirchhoff viscous stress
	free_d4tensor(Sv11_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(Sv12_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(Sv13_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(Sv22_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(Sv23_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(Sv33_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);

	// 2nd Piola-Kirchhoff elastic stress 
	free_d3tensor(Se11, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Se12, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Se13, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Se22, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Se23, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Se33, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);

	// Stress pullback to axisymmetric configuration
	free_d3tensor(S0_11, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(S0_12, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(S0_13, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(S0_22, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(S0_23, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(S0_33, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d4tensor(S0_11_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar); 
	free_d4tensor(S0_12_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar); 
	free_d4tensor(S0_13_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar); 
	free_d4tensor(S0_22_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar); 
	free_d4tensor(S0_23_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar); 
	free_d4tensor(S0_33_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);

	// active fiber strain vectors
	free_d3tensor(Sf_ff_const, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Sf_temp, 0,Nphi-1,0,Nnu-1,0,Nmu-1);
	free_d4tensor(Sf_ff_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);

	free_d3tensor(Sf11_const, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Sf12_const, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Sf13_const, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Sf22_const, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Sf23_const, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_d3tensor(Sf33_const, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);

	free_d4tensor(Sf11_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(Sf12_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(Sf13_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(Sf22_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(Sf23_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	free_d4tensor(Sf33_ai, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1, 1, Npar);
	// Active stress variables
	free_d3tensor(SigmaF, 0,Nphi-1,0,Nnu-1,0,Nmu-1);
	free_d4tensor(SigmaF_ai, 0,Nphi-1,0,Nnu-1,0,Nmu-1, 1, Npar);

	// Initialize equilibrium integrands
	free_d3tensor(integrand, 0, Nphi - 1, 0, Nnu - 1, 0, Nmu - 1);
	free_dmatrix(integrand2, 0,Nphi-1,0,Nnu-1);
	free_dvector(integrand3, 0,Nphi-1);
	free_dvector(integrand1, 0,Nmu-1);

	// Initialize volume derivative
	free_dmatrix(vIntegrand_dai, 0, Nnu - 1, 1, Npar);
	free_dvector(vol_integrand, 0, Nnu - 1);

	// Initialize spherical chamber variables
	free_dvector(r, 0, Nr - 1);
	free_dvector(r0_la, 0, Nr - 1);
	free_dvector(r0_rv, 0, Nr - 1);
	free_dvector(r0_ra, 0, Nr - 1);
	free_dvector(kappaj_integrand, 0, Nr - 1);
	free_dvector(chij_integrand, 0, Nr - 1);

	// Set up solution storage vector
	free_dmatrix(Y, 0, Nvars - 1, 0, Nt - 1);
	// Set up storage for other variables
	free_dmatrix(store, 0, Nstore - 1, 0, Nt - 1);

	// Set up storage for a1,a2,a3,a4,a5,a6,a7
	free_dvector(a1vec, 0, Nt - 1);
	free_dvector(a2vec, 0, Nt - 1);
	free_dvector(a3vec, 0, Nt - 1);
	free_dvector(a4vec, 0, Nt - 1);
	free_dvector(a5vec, 0, Nt - 1);
	free_dvector(a6vec, 0, Nt - 1);
	free_dvector(a7vec, 0, Nt - 1);

	// Set up storage for the time variable
	free_dvector(tvec, 0, Nt - 1);

	//storage for short-axis and long-axis plots
	free_d3tensor(ch_store, 0, Nnu - 1, 1, 4, 1, 8);
	free_d3tensor(sh_store, 0, Nphi - 1, 1, 4, 1, 4);
}
