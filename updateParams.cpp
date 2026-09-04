// updateParams - using values stored in Param array
// TWS, July 2025

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <complex>
#include "nrutil.h"

typedef std::complex<double> CX;
using namespace std;

void updateParams(double* Param, int Nparam)
{
	extern int Nvars, Nmat, iRecord, Nstore, Nr, Ncycles, Nt, Nta, Ncycles_show, aortaVar;
	extern double Tc, dtRecord, dtMax;
	extern double Ta, Tca, Tca_shift, f_rv;

	// LV parameters
	extern int Nmu1, Nmu2, Nphi, Nnu;
	extern double a0, a40, a60, nu_up, muin0, mumid0, muout0;
	extern double psi_in_b0, psi_out_b0, phi_d, phi_s, phi_t, alpha;
	extern double Ls0, Lsmax, Lsw;
	extern double km, kav, kdfac;
	extern double kv, c1, bff, bfx, bxx;

	// Other chamber parameters
	extern double km_la, kav_la, kv_la, br_la, c1_la, bperp_la, rin0_la, rout0_la;
	extern double km_rv, kav_rv, kv_rv, br_rv, c1_rv, bperp_rv, rin0_rv, rout0_rv;
	extern double bff_rv, bxx_rv, bfx_rv;
	extern double km_ra, kav_ra, kv_ra, br_ra, c1_ra, bperp_ra, rin0_ra, rout0_ra;
	extern double At_la;

	// Lumped model parameters
	extern double Vla, Vra, Vrv, Vlv;
	extern double Csa, Csv, Cpa, Cpv;
	extern double Rsp, Rsv, Rpa, Rpp, Rpv, Rra;
	extern double q_la, q_lv, q_ra, q_rv;
	extern double q_sp, q_sv, q_pp, q_pv;

	// Valve/vein parameters
	extern double B_lv, B_la, B_rv, B_ra, L_lv, L_la, L_rv, L_ra, L_pv, L_sv;
	extern double Zmiv, Zaov, Ztcv, Zpuv;
	extern double Kaovo, Kaovc, Kmivo, Kmivc, Kpuvo, Kpuvc, Ktcvo, Ktcvc;
	extern double VA_aovo, VA_aovc, VA_mivo, VA_mivc, VA_puvo, VA_puvc, VA_tcvo, VA_tcvc, VL_aov, VL_miv, VL_puv, VL_tcv, rho;

	// Initial pressures, external tissue pressure
	extern double Psp0, Psv0, Ppp0, Ppv0, Pext;

	//pericardium parameters
	extern double V0_peri, s_peri, Pperi;

	//LVAD parameters
	extern double q0_LVAD, Ps_LVAD, t0_LVAD;

	extern string odeSystem, septalDisp, LVAD, RVhypertrophy, pericardium;

	double RefVolLV, RefVolLA, RefVolRV, RefVolRA, WallVolLV, WallVolLA, WallVolRV, WallVolRA, WallVolSep;
	double Vmuout0 = 0., Vmumid0, CSS, LAtoSA, RVfract, third = 1. / 3.;
	double VA_pv, VA_sv, VL_pv, VL_sv;

	// used for cubic solving
	double c0, pp, qq, acube, chout0;
	CX sdel, root;

	double kp;	//not used

	Ncycles = (int)Param[1];
	Ncycles_show = (int)Param[2];
	Nmu1 = (int)Param[3];
	Nmu2 = (int)Param[4];
	Nnu = (int)Param[5];
	Nphi = (int)Param[6];	
	Nr = (int)Param[7];
	// Time steps
	dtMax = Param[8];
	dtRecord = Param[9];
	// density
	rho = Param[10];
	// Time parameters (s) 
	Tc = Param[11];
	Ta = Param[12];
	Tca = Param[13];
	Tca_shift = Param[14];

	// Reference configuration parameters (LV)
	a0 = Param[15];
	a40 = Param[16];
	a60 = Param[17];
	muin0 = Param[18];
	mumid0 = Param[19];
	muout0 = Param[20];
	phi_d = Param[21];
	phi_s = Param[22];
	phi_t = Param[23];
	alpha = Param[24];
	nu_up = Param[25];

	//Alternative geometric parameters
	RefVolLV = Param[26];
	WallVolLV = Param[27];
	WallVolSep = Param[28];
	LAtoSA = Param[29];

	// Muscle fiber angles
	psi_in_b0 = Param[30];
	psi_out_b0 = Param[31];

	// Muscle fiber sarcomere parameters
	Ls0 = Param[32];
	Lsmax = Param[33];
	Lsw = Param[34];

	// LV active force generation, activation-dependent viscous resistance
	km = Param[35];
	kav = Param[36];
	// LV activation with sarcomere length
	kp = Param[37];	//not used
	// LV, RV steepen activation onset (if kdfac < 1)
	kdfac = Param[38];

	// LV passive stress parameters: viscous, elastic
	kv = Param[39];
	c1 = Param[40];
	bff = Param[41];
	bxx = Param[42];
	bfx = Param[43];

	// Right ventricle parameters
	rin0_rv = Param[44];
	rout0_rv = Param[45];
	f_rv = Param[46];
	RefVolRV = Param[47];
	WallVolRV = Param[48];
	km_rv = Param[49];
	kav_rv = Param[50];
	kv_rv = Param[51];
	c1_rv = Param[52];
	bff_rv = Param[53];
	bxx_rv = Param[54];
	bfx_rv = Param[55];
	br_rv = Param[56];
	bperp_rv = Param[57];

	// Left atrium parameters
	rin0_la = Param[58];
	rout0_la = Param[59];
	RefVolLA = Param[60];
	WallVolLA = Param[61];
	km_la = Param[62];
	kav_la = Param[63];
	kv_la = Param[64];
	c1_la = Param[65];
	br_la = Param[66];
	bperp_la = Param[67];

	// Right atrium parameters
	rin0_ra = Param[68];
	rout0_ra = Param[69];
	RefVolRA = Param[70];
	WallVolRA = Param[71];
	km_ra = Param[72];
	kav_ra = Param[73];
	kv_ra = Param[74];
	c1_ra = Param[75];
	br_ra = Param[76];
	bperp_ra = Param[77];

	// Valve parameters (Mynard dissertation) area VA in cm2, length VL in cm
	VA_aovo = Param[78];
	VL_aov = Param[79];
	VA_mivo = Param[80];
	VL_miv = Param[81];
	VA_puvo = Param[82];
	VL_puv = Param[83];
	VA_tcvo = Param[84];
	VL_tcv = Param[85];
	VA_aovc = Param[86];
	VA_mivc = Param[87];
	VA_puvc = Param[88];
	VA_tcvc = Param[89];

	//Valve rate constants,1/(kPa s) from Mynard in cm2/dyn/s
	Kaovo = Param[90];
	Kmivo = Param[91];
	Kpuvo = Param[92];
	Ktcvo = Param[93];
	Kaovc = Param[94];
	Kmivc = Param[95];
	Kpuvc = Param[96];
	Ktcvc = Param[97];

	// Vein inertia parameters
	VA_pv = Param[98];
	VL_pv = Param[99];
	VA_sv = Param[100];
	VL_sv = Param[101];

	// Systemic resistance
	Rsp = Param[109];
	Rsv = Param[110];

	// Pulmonary resistance
	Rpa = Param[111];
	Rpp = Param[112];
	Rpv = Param[113];

	// Compliances
	Csa = Param[114];
	Csv = Param[115];
	Cpa = Param[116];
	Cpv = Param[117];

	// External tissue pressure
	Pext = Param[118];

	// Initial pressures
	Psp0 = Param[119];
	Psv0 = Param[120];
	Ppp0 = Param[121];
	Ppv0 = Param[122];

	//pericardium parameters
	V0_peri = Param[123];
	s_peri = Param[124];

	//LVAD parameters: max flow in cm3/s, stall pressure in mmHg, start time in cycles
	q0_LVAD = Param[125];
	Ps_LVAD = Param[126];
	t0_LVAD = Param[127];

	L_pv = 1.e-4 * rho * VL_pv / VA_pv;
	L_sv = 1.e-4 * rho * VL_sv / VA_sv;

	//Alternative geometric parameters
	if (RefVolLV > 0.) {
		if (odeSystem == "LVRVmodel2024") {
			// Fraction of volume between mumid0 and muout0 belonging to RV free wall
			// May need adjustment if alpha and phi values are altered
			RVfract = 0.417;
			Vmuout0 = RefVolLV + WallVolLV / (1. - RVfract);
		}
		else Vmuout0 = RefVolLV + WallVolLV;
		muin0 = atanh((cos(nu_up) + 1.) / 2. / LAtoSA);
		Param[18] = muin0;
		CSS = (cos(nu_up) + 1.) * cosh(muin0) * DSQR(sinh(muin0))
			+ (cosh(muin0) - 1.) * cos(nu_up) * DSQR(sin(nu_up));
		a0 = pow(RefVolLV * 1.5 / PI / CSS, third);
		Param[15] = a0;
		c0 = cos(nu_up);
		acube = DCUB(a0);
		pp = -(1. - c0 + c0 * c0);
		qq = -c0 * (1. - c0) - 3. * Vmuout0 / (2. * PI * acube * (1. + c0));
		sdel = sqrt(DSQR(qq / 2.) + DCUB(pp / 3.));
		root = pow(-qq / 2. + sdel, third) + pow(-qq / 2. - sdel, third);
		chout0 = real(root); // The imaginary part is 0, but still need to take only real part
		muout0 = acosh(chout0);
		Param[20] = muout0;
		if (odeSystem == "LVRVmodel2024") {
			Vmumid0 = RefVolLV + WallVolSep / RVfract;
			qq = -c0 * (1. - c0) - 3. * Vmumid0 / (2. * PI * acube * (1. + c0));
			sdel = sqrt(DSQR(qq / 2.) + DCUB(pp / 3.));
			root = pow(-qq / 2. + sdel, third) + pow(-qq / 2. - sdel, third);
			chout0 = real(root);
			mumid0 = acosh(chout0);
			Param[19] = mumid0;
		}
	}
	if (RefVolRV > 0. && odeSystem != "LVRVmodel2024") {
		rin0_rv = pow(RefVolRV * 0.75 / PI / f_rv, third);
		rout0_rv = pow((RefVolRV + WallVolRV) * 0.75 / PI / f_rv, third);
		Param[44] = rin0_rv;
		Param[45] = rout0_rv;
	}
	if (RefVolLA > 0.) {
		rin0_la = pow(RefVolLA * 0.75 / PI, third);
		rout0_la = pow((RefVolLA + WallVolLA) * 0.75 / PI, third);
		Param[58] = rin0_la;
		Param[59] = rout0_la;
	}
	if (RefVolRA > 0.) {
		rin0_ra = pow(RefVolRA * 0.75 / PI, third);
		rout0_ra = pow((RefVolRA + WallVolRA) * 0.75 / PI, third);
		Param[68] = rin0_ra;
		Param[69] = rout0_ra;
	}
}

