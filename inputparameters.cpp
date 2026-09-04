// inputparameters.cpp
// edited by TWS, November 2019
#define _CRT_SECURE_NO_WARNINGS	//TWS
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include "nrutil.h"

typedef std::complex<double> CX;
using namespace std;

// This finds a vector of values from a line that contains the key "string"
void findArray(ifstream &fileID, string key, int N, double *values) {
	double success = 0;
	string str;
	size_t foundKey;
	stringstream numStr;
	fileID.clear();
	fileID.seekg(0, ios::beg);
	int k = 0;
	while (std::getline(fileID, str)) {
		// look through file for the key
		foundKey = str.find(key);
		char * checkChar = new char[str.length() + 1];
		std::strcpy(checkChar, str.c_str());
		size_t foundEscape = str.find("//");

		if (foundKey != string::npos
			&& (foundKey == 0 || isspace(checkChar[foundKey - 1]))
			&& (isspace(checkChar[foundKey + key.length()]) || (checkChar[foundKey + key.length()] == '='))
			&& foundEscape == string::npos) {

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// add something on the end so it doesnt miss the last number
			str.append(" end line ");

			// once the key is found, search for the number
			char * strChar = new char[str.length() + 1];
			std::strcpy(strChar, str.c_str());
			success = 0.5;

			// once the key is found, search for the values
			char c;
			char cprev = 'e';
			int i = 0;

			while (strChar[i]) {
				c = strChar[i];
				if (isdigit(c) || c == '.' || (c == 'e' && isdigit(cprev)) || c == '+' || c == '-') numStr << c;
				else if (numStr.str().length() > 0) {
					numStr >> values[k];
					k = k + 1;
					if (k > N) {
						success = -1;
						break;
					}
					numStr.str(std::string());
					numStr.clear();
				}
				else {
					// if its not a number or a . and no number is stored then just skip
				}
				cprev = c;
				i++;
			}
			delete[] strChar;
		}
		delete[] checkChar;
	}
	// success if k got to n
	if (k == N && success != -1) {
		//cout << "Vector " << key << " successfully copied with " << N << " elements." << endl;
		success = 1;
	}
	if (success == 0) cout << "!!! Vector " << key << " not found. Vector arbitrarily set to 0. !!!" << endl;
	else if (success == 0.5)
		cout << "!!! Vector " << key << " did not have enough elements! Remaining elements arbitrarily set to 0. !!!" << endl;
	else if (success == -1)
		cout << "!!! Vector " << key << " had too many elements! Returning vector with only " << N << " elements !!!" << endl;
}

// This finds a single value from a line that contains the key "string"
// Doesnt read lines that contain "//"
double findValue(ifstream &fileID, string key) {
	double value;
	bool success = 0;
	string str;
	size_t foundKey;
	stringstream numStr;
	fileID.clear();
	fileID.seekg(0, ios::beg);

	while (std::getline(fileID, str)) {
		// look through file for the key
		foundKey = str.find(key);
		char * checkChar = new char[str.length() + 1];
		std::strcpy(checkChar, str.c_str());
		size_t foundEscape = str.find("//");

		if (foundKey != string::npos
			&& (foundKey == 0 || isspace(checkChar[foundKey - 1]))
			&& (isspace(checkChar[foundKey + key.length()]) || (checkChar[foundKey + key.length()] == '='))
			&& foundEscape == string::npos) {

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// once the key is found, search for the number
			char * strChar = new char[str.length() + 1];
			std::strcpy(strChar, str.c_str());
			int i = 0;
			char c;
			char cprev = 'e';
			while (strChar[i]) {
				c = strChar[i];
				if (isdigit(c) || c == '.' || (c == 'e' && isdigit(cprev)) || c == '+' || c == '-') numStr << c;
				else if (numStr.str().length() > 0) break;
				cprev = c;
				i++;
			}
			numStr >> value;
			success = 1;
			delete[] strChar;
		}
		delete[] checkChar;
	}
	if (success == 1) printf(""); // printf("Set %s = %6f\n", key.c_str(), value);
	else {
		cout << "*** Parameter " << key << " not found. Parameter arbitrarily set to 0. ***" << endl;
		value = 0;
	}
	return value;
}

// This finds a string that follows the key
string findString(ifstream &fileID, string key) {
	string outString;
	bool success = 0;
	string str;
	size_t foundKey;
	stringstream valueStr;
	fileID.clear();
	fileID.seekg(0, ios::beg);

	while (std::getline(fileID, str)) {
		// look through file for the key
		foundKey = str.find(key);
		char * checkChar = new char[str.length() + 1];
		std::strcpy(checkChar, str.c_str());
		size_t foundEscape = str.find("//");
		if (foundKey != string::npos
			&& (foundKey == 0 || isspace(checkChar[foundKey - 1]))
			&& (isspace(checkChar[foundKey + key.length()]) || (checkChar[foundKey + key.length()] == '='))
			&& foundEscape == string::npos) {

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// once the key is found, search for the number
			char * strChar = new char[str.length() + 1];
			std::strcpy(strChar, str.c_str());
			int i = 0;
			char c;
			char cprev = 'e';
			while (strChar[i]) {
				c = strChar[i];
				// skip '=' or spaces
				if (c != ' ' && c != '=' && c != '\r') valueStr << c;
				else if (valueStr.str().length() > 0) 	break;
				cprev = c;
				i++;
			}
			outString = valueStr.str();
			success = 1;
			delete[] strChar;
		}
		delete[] checkChar;
	}
	if (success == 1)	printf(""); //cout << "Set string " << key << " to \"" << outString << "\"" << endl;
	else {
		cout << "*** String " << key << " not found. String is left empty. ***" << endl;
		outString = "";
	}
	return outString;
}

void inputparameters(int run)
{
	extern string endPause, odeSystem, septalDisp, LVAD, RVhypertrophy, pericardium;

	// LV parameters
	extern int Nvars, Nmu1, Nmu2, Nmu, Nnu, Nphi, Ncycles, Ncycles_show, Npar;
	extern double nu_up, Tc, Ta;
	extern double dtMax, dtRecord;
	extern double a0, a40, a60, muin0, muin0, mumid0, muout0, phi_t, phi_s, phi_d, alpha;
	extern double psi_in_b0, psi_out_b0;
	extern double Ls0, Lsmax, Lsw;
	extern double km, kav, kdfac;
	extern double kv, c1, bff, bfx, bxx;

	// Circulatory parameters
	extern double Csa, Csv, Cpa, Cpv;
	extern double Rsp, Rsv, Rpa, Rpp, Rpv, Rra;

	// Other chamber parameters
	extern int Nr;
	extern double Tca, Tca_shift;
	extern double km_la, kav_la, kv_la, br_la, c1_la, bperp_la, rin0_la, rout0_la;
	extern double km_rv, kav_rv, kv_rv, br_rv, c1_rv, bperp_rv, rin0_rv, rout0_rv;
	extern double bff_rv, bxx_rv, bfx_rv;
	extern double km_ra, kav_ra, kv_ra, br_ra, c1_ra, bperp_ra, rin0_ra, rout0_ra;
	extern double *y_init;
	extern double f_rv;

	//valve parameters
	extern double B_lv, B_la, B_rv, B_ra, L_lv, L_la, L_rv, L_ra;
	extern double VA_aovo, VA_aovc, VA_mivo, VA_mivc, VA_puvo, VA_puvc, VA_tcvo, VA_tcvc, VL_aov, VL_miv, VL_puv, VL_tcv, rho;
	extern double Kaovo, Kaovc, Kmivo, Kmivc, Kpuvo, Kpuvc, Ktcvo, Ktcvc;

	//vein parameters
	extern double L_pv, L_sv;
	double VA_pv, VA_sv, VL_pv, VL_sv;

	//aorta parameters
	extern double Laorta, A00, A01, Gaorta0, Gaorta1, tpulse, impedancefac;

	//External tissue pressure, inital pressures
	extern double Pext, Psp0, Psv0, Ppp0, Ppv0;

	//pericardium parameters
	extern double V0_peri, s_peri;

	//LVAD parameters
	extern double q0_LVAD, Ps_LVAD, t0_LVAD;

///////////////////////////////////////////////////////////////////////////////////////////
	double RefVolLV, RefVolLA, RefVolRV, RefVolRA, WallVolLV, WallVolLA, WallVolRV, WallVolRA, WallVolSep;
	double Vmuout0 = 0., Vmumid0, CSS, LAtoSA, RVfract, third = 1. / 3.;

	// used for cubic solving
	double c0 = 0., pp, qq, acube = 0., chout0;
	CX sdel, root;
	
	string InputFile;
	char fname[80], fname1[80];
	sprintf(fname, "parameters%03i.txt", run);
	sprintf(fname1, "initialValues%03i.txt", run);
	InputFile = fname;
	
	ifstream paramFile(InputFile.c_str(), ios::in);

	// Import parameter values from settings file
	// Use LVRV model?
	odeSystem = findString(paramFile, "odeSystem");
	// Pause with console open after run?
	endPause = findString(paramFile, "endPause");
	// include septal displacement (a5)?
	septalDisp = findString(paramFile, "septalDisp");
	// include RVhypertrophy (a60)?
	RVhypertrophy = findString(paramFile, "RVhypertrophy");
	// include LVAD? Choose on or off
	LVAD = findString(paramFile, "LVAD");
	// include pericardium? Choose on or off
	pericardium = findString(paramFile, "pericardium");

	//density
	rho = findValue(paramFile, "rho");

	// Physical domain parameters
	Nmu1 = (int)findValue(paramFile, "Nmu1");
	Nmu2 = (int)findValue(paramFile, "Nmu2");
	Nmu = Nmu1 + Nmu2;
	Nnu = (int)findValue(paramFile, "Nnu");
	Nphi = (int)findValue(paramFile, "Nphi");
	nu_up = findValue(paramFile, "nu_up");

	// Time(/ cycle) parameters
	Ncycles = (int)findValue(paramFile, "Ncycles");
	Ncycles_show = (int)findValue(paramFile, "Ncycles_show");
	Tc = findValue(paramFile, "Tc");
	Ta = findValue(paramFile, "Ta");
	Tca = findValue(paramFile, "Tca");
	Tca_shift = findValue(paramFile, "Tca_shift");
	
	// Time steps
	dtMax = findValue(paramFile, "dtMax");
	dtRecord = findValue(paramFile, "dtRecord");

	// Initial structure parameters
	a0 = findValue(paramFile, "a0");
	a40 = findValue(paramFile, "a40");
	a60 = findValue(paramFile, "a60");
	muin0 = findValue(paramFile, "muin0");
	mumid0 = findValue(paramFile, "mumid0");
	muout0 = findValue(paramFile, "muout0");
	phi_d = findValue(paramFile, "phi_d");
	phi_s = findValue(paramFile, "phi_s");
	phi_t = findValue(paramFile, "phi_t");
	alpha = findValue(paramFile, "alpha");

	// Muscle fiber parameterization angles
	psi_in_b0 = findValue(paramFile, "psi_in_b0");
	psi_out_b0 = findValue(paramFile, "psi_out_b0");

	// Muscle fiber sarcomere parameters
	Ls0 = findValue(paramFile, "Ls0");
	Lsmax = findValue(paramFile, "Lsmax");
	Lsw = findValue(paramFile, "Lsw");

	// Active stress parameters
	km = findValue(paramFile, "km");
	kav = findValue(paramFile, "kav");
	kdfac = findValue(paramFile, "kdfac");

	// Passive stress parameters
	kv = findValue(paramFile, "kv");
	c1 = findValue(paramFile, "c1");
	bff = findValue(paramFile, "bff");
	bxx = findValue(paramFile, "bxx");
	bfx = findValue(paramFile, "bfx");

	// Compliances
	Csa = findValue(paramFile, "Csa");
	Csv = findValue(paramFile, "Csv");
	Cpa = findValue(paramFile, "Cpa");
	Cpv = findValue(paramFile, "Cpv");

	// Hydraulic resistance
	Rsp = findValue(paramFile, "Rsp");
	Rsv = findValue(paramFile, "Rsv");
	Rpa = findValue(paramFile, "Rpa");
	Rpp = findValue(paramFile, "Rpp");
	Rpv = findValue(paramFile, "Rpv");

	// Left atrium parameters	
	Nr = (int)findValue(paramFile, "Nr");
	km_la = findValue(paramFile, "km_la");
	kav_la = findValue(paramFile, "kav_la");
	kv_la = findValue(paramFile, "kv_la");
	c1_la = findValue(paramFile, "c1_la");
	br_la = findValue(paramFile, "br_la");
	bperp_la = findValue(paramFile, "bperp_la");
	rin0_la = findValue(paramFile, "rin0_la");
	rout0_la = findValue(paramFile, "rout0_la");

	// Right ventricle parameters
	km_rv = findValue(paramFile, "km_rv");
	f_rv = findValue(paramFile, "f_rv");
	kav_rv = findValue(paramFile, "kav_rv");
	kv_rv = findValue(paramFile, "kv_rv");
	c1_rv = findValue(paramFile, "c1_rv");
	bff_rv = findValue(paramFile, "bff_rv");
	bxx_rv = findValue(paramFile, "bxx_rv");
	bfx_rv = findValue(paramFile, "bfx_rv");
	br_rv = findValue(paramFile, "br_rv");
	bperp_rv = findValue(paramFile, "bperp_rv");
	rin0_rv = findValue(paramFile, "rin0_rv");
	rout0_rv = findValue(paramFile, "rout0_rv");

	// Right atrium parameters
	km_ra = findValue(paramFile, "km_ra");
	kav_ra = findValue(paramFile, "kav_ra");
	kv_ra = findValue(paramFile, "kv_ra");
	c1_ra = findValue(paramFile, "c1_ra");
	br_ra = findValue(paramFile, "br_ra");
	bperp_ra = findValue(paramFile, "bperp_ra");
	rin0_ra = findValue(paramFile, "rin0_ra");
	rout0_ra = findValue(paramFile, "rout0_ra");

	// Alternative geometric parameters
	RefVolLV = findValue(paramFile, "RefVolLV");
	WallVolLV = findValue(paramFile, "WallVolLV");
	WallVolSep = findValue(paramFile, "WallVolSep");
	LAtoSA = findValue(paramFile, "LAtoSA");
	RefVolRV = findValue(paramFile, "RefVolRV");
	WallVolRV = findValue(paramFile, "WallVolRV");

	// Initial conditions
	Psp0 = findValue(paramFile, "Psp0");
	Psv0 = findValue(paramFile, "Psv0");
	Ppp0 = findValue(paramFile, "Ppp0");
	Ppv0 = findValue(paramFile, "Ppv0");

	if (RefVolLV > 0.) {
		if (odeSystem == "LVRVmodel2024") {
			// Fraction of volume between mumid0 and muout0 belonging to RV free wall
			// May need adjustment if alpha and phi values are altered
			RVfract = 0.417;
			Vmuout0 = RefVolLV + WallVolLV / (1. - RVfract);
		}
		else Vmuout0 = RefVolLV + WallVolLV;
		muin0 = atanh((cos(nu_up) + 1.) / 2. / LAtoSA);
		CSS = (cos(nu_up) + 1.) * cosh(muin0) * DSQR(sinh(muin0))
			+ (cosh(muin0) - 1.) * cos(nu_up) * DSQR(sin(nu_up));
		a0 = pow(RefVolLV * 1.5 / PI / CSS, third);
		c0 = cos(nu_up);
		acube = DCUB(a0);
		pp = -(1. - c0 + c0 * c0);
		qq = -c0 * (1. - c0) - 3. * Vmuout0 / (2. * PI * acube * (1. + c0));
		sdel = sqrt(DSQR(qq / 2.) + DCUB(pp / 3.));
		root = pow(-qq / 2. + sdel, third) + pow(-qq / 2. - sdel, third);
		chout0 = real(root); // The imaginary part is 0, but still need to take only real part
		muout0 = acosh(chout0);
		if (odeSystem == "LVRVmodel2024") {
			Vmumid0 = RefVolLV + WallVolSep / RVfract;
			qq = -c0 * (1. - c0) - 3. * Vmumid0 / (2. * PI * acube * (1. + c0));
			sdel = sqrt(DSQR(qq / 2.) + DCUB(pp / 3.));
			root = pow(-qq / 2. + sdel, third) + pow(-qq / 2. - sdel, third);
			chout0 = real(root);
			mumid0 = acosh(chout0);
		}
	}
	if (RefVolRV > 0. && odeSystem != "LVRVmodel2024") {
		rin0_rv = pow(RefVolRV * 0.75 / PI / f_rv, third);
		rout0_rv = pow((RefVolRV + WallVolRV) * 0.75 / PI / f_rv, third);
	}
	RefVolLA = findValue(paramFile, "RefVolLA");
	WallVolLA = findValue(paramFile, "WallVolLA");
	if (RefVolLA > 0.) {
		rin0_la = pow(RefVolLA * 0.75 / PI, third);
		rout0_la = pow((RefVolLA + WallVolLA) * 0.75 / PI, third);
	}
	RefVolRA = findValue(paramFile, "RefVolRA");
	WallVolRA = findValue(paramFile, "WallVolRA");
	if (RefVolRA > 0.) {
		rin0_ra = pow(RefVolRA * 0.75 / PI, third);
		rout0_ra = pow((RefVolRA + WallVolRA) * 0.75 / PI, third);
	}

	// Valve parameters
	VA_aovo = findValue(paramFile, "VA_aovo");
	VL_aov = findValue(paramFile, "VL_aov");
	VA_mivo = findValue(paramFile, "VA_mivo");
	VL_miv = findValue(paramFile, "VL_miv");
	VA_puvo = findValue(paramFile, "VA_puvo");
	VL_puv = findValue(paramFile, "VL_puv");
	VA_tcvo = findValue(paramFile, "VA_tcvo");
	VL_tcv = findValue(paramFile, "VL_tcv");
	VA_pv = findValue(paramFile, "VA_pv");
	VL_pv = findValue(paramFile, "VL_pv");
	VA_sv = findValue(paramFile, "VA_sv");
	VL_sv = findValue(paramFile, "VL_sv");
	L_pv = 1.e-4 * rho * VL_pv / VA_pv;
	L_sv = 1.e-4 * rho * VL_sv / VA_sv;
	VA_aovc = findValue(paramFile, "VA_aovc");
	VA_mivc = findValue(paramFile, "VA_mivc");
	VA_puvc = findValue(paramFile, "VA_puvc");
	VA_tcvc = findValue(paramFile, "VA_tcvc");
	//Kaovo, Kaovc, Kmivo, Kmivc, Kpuvo, Kpuvc, Ktcvo, Ktcvc;
	Kaovo = findValue(paramFile, "Kaovo");
	Kmivo = findValue(paramFile, "Kmivo");
	Kpuvo = findValue(paramFile, "Kpuvo");
	Ktcvo = findValue(paramFile, "Ktcvo");
	Kaovc = findValue(paramFile, "Kaovc");
	Kmivc = findValue(paramFile, "Kmivc");
	Kpuvc = findValue(paramFile, "Kpuvc");
	Ktcvc = findValue(paramFile, "Ktcvc");

	// External tissue pressure
	Pext = findValue(paramFile, "Pext");

	//pericardium parameters
	V0_peri = findValue(paramFile, "V0_peri");
	s_peri = findValue(paramFile, "s_peri");

	//LVAD parameters: max flow in cm3/s, stall pressure in mmHg, start time in cycles
	q0_LVAD = findValue(paramFile, "q0_LVAD");
	Ps_LVAD = findValue(paramFile, "Ps_LVAD");
	t0_LVAD = findValue(paramFile, "t0_LVAD");

	//aorta parameters
	Laorta = findValue(paramFile, "Laorta");
	A00 = findValue(paramFile, "A00");
	A01 = findValue(paramFile, "A01");
	Gaorta0 = findValue(paramFile, "Gaorta0"); // entrance compliance in (dyn/cm2)^-1
	Gaorta1 = findValue(paramFile, "Gaorta1"); // decrease in compliance with distance
	tpulse = findValue(paramFile, "tpulse"); //length of pulse
	impedancefac = findValue(paramFile, "impedancefac");	//distal impedance factor, matched = 1
	paramFile.close();
}
