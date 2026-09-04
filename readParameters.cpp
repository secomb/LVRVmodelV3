// readParameters.cpp
// Reads parameters.txt
// TWS, July 2025

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "nrutil.h"

using namespace std;

// This finds a single value from a line that contains the key "string"
// Doesn't read lines that contain "//"
double findValue(ifstream& fileID, string key) {
	double value;
	bool success = 0;
	string str;
	size_t foundKey;
	stringstream numStr;
	fileID.clear();
	fileID.seekg(0, ios::beg);

	while (std::getline(fileID, str)) {				// look through file for the key
		foundKey = str.find(key);
		char* checkChar = new char[str.length() + 1];
		std::strcpy(checkChar, str.c_str());
		size_t foundEscape = str.find("//");

		if (foundKey != string::npos
			&& (foundKey == 0 || isspace(checkChar[foundKey - 1]))
			&& (isspace(checkChar[foundKey + key.length()]) || (checkChar[foundKey + key.length()] == '='))
			&& foundEscape == string::npos) {

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// once the key is found, search for the number
			char* strChar = new char[str.length() + 1];
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
	if (success == 1) printf(""); //printf("Set %s = %6f\n", key.c_str(), value);
	else {
		//cout << "*** Parameter " << key << " not found. Set to 0. ***" << endl;
		value = 0;
	}
	return value;
}

// This finds a vector of values from a line that contains the key "string"
void findArray(ifstream& fileID, string key, int N, double* values) {
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
		char* checkChar = new char[str.length() + 1];
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
			char* strChar = new char[str.length() + 1];
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
		cout << "Vector " << key << " successfully copied with " << N << " elements." << endl;
		success = 1;
	}
	if (success == 0) cout << "Vector " << key << " not found. Vector arbitrarily set to 0." << endl;
	else if (success == 0.5)
		cout << "Vector " << key << " did not have enough elements. Remaining elements arbitrarily set to 0." << endl;
	else if (success == -1)
		cout << "Vector " << key << " had too many elements. Returning vector with only " << N << " elements." << endl;
}

// This finds a string that follows the key
string findString(ifstream& fileID, string key) {
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
		char* checkChar = new char[str.length() + 1];
		std::strcpy(checkChar, str.c_str());
		size_t foundEscape = str.find("//");
		if (foundKey != string::npos
			&& (foundKey == 0 || isspace(checkChar[foundKey - 1]))
			&& (isspace(checkChar[foundKey + key.length()]) || (checkChar[foundKey + key.length()] == '='))
			&& foundEscape == string::npos) {

			// erase the key, as it might contain a digit
			str.erase(foundKey, key.length());

			// once the key is found, search for the number
			char* strChar = new char[str.length() + 1];
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
	if (success == 1) printf(""); //cout << "Set string " << key << " to \"" << outString << "\"" << endl;
	else {
		cout << "*** String " << key << " not found. String is left empty. ***" << endl;
		outString = "";
	}
	return outString;
}

void readParameters(int run)
{
	extern int Nparam, Nopt, aortaVar;		//test if aorta parameters are being varied
	extern int* ParamVar, * ia;
	extern double* Param0, *xx, *xxlb, *xxub;
	extern std::string ParamName[128];	//note: zero-based array, include extra element
	extern string odeSystem, septalDisp, LVAD, RVhypertrophy, pericardium, Optimize, Skip;

	int i, j;
	double* xx00;

	Nparam = 127;
	Param0 = dvector(1, Nparam);
	ParamVar = ivector(1, Nparam);

	// Number of cycles
	ParamName[1] = "Ncycles";
	ParamName[2] = "Ncycles_show";

	// Domain discretization parameters. LV, RV: Nmu1, Nmu2, Nphi, Nnu, other chambers: Nr
	ParamName[3] = "Nmu1";
	ParamName[4] = "Nmu2";
	ParamName[5] = "Nnu";
	ParamName[6] = "Nphi";
	ParamName[7] = "Nr";

	// Time steps
	ParamName[8] = "dtMax";
	ParamName[9] = "dtRecord";

	// density
	ParamName[10] = "rho";

	// Time parameters (s) 
	ParamName[11] = "Tc";

	// *** The following parameters can be used as optimization variables ***
	// Time parameters (s) 
	ParamName[12] = "Ta";
	ParamName[13] = "Tca";
	ParamName[14] = "Tca_shift";

	// Reference configuration parameters (LV)
	ParamName[15] = "a0";
	ParamName[16] = "a40";
	ParamName[17] = "a60";
	ParamName[18] = "muin0";
	ParamName[19] = "mumid0";
	ParamName[20] = "muout0";
	ParamName[21] = "phi_d";
	ParamName[22] = "phi_s";
	ParamName[23] = "phi_t";
	ParamName[24] = "alpha";
	ParamName[25] = "nu_up";

	// Alternative geometric parameters
	ParamName[26] = "RefVolLV";
	ParamName[27] = "WallVolLV";
	ParamName[28] = "WallVolSep";
	ParamName[29] = "LAtoSA";
	
	// Muscle fiber parameterization angles
	ParamName[30] = "psi_in_b0";
	ParamName[31] = "psi_out_b0";

	// Muscle fiber sarcomere parameters
	ParamName[32] = "Ls0";
	ParamName[33] = "Lsmax";
	ParamName[34] = "Lsw";

	// LV active force generation, activation-dependent viscous resistance
	ParamName[35] = "km";
	ParamName[36] = "kav";
	// LV activation with sarcomere length
	ParamName[37] = "kp";	//not used
	// LV, RV steepen activation onset (if kdfac < 1)
	ParamName[38] = "kdfac";

	// LV passive stress parameters: viscous, elastic
	ParamName[39] = "kv";
	ParamName[40] = "c1";
	ParamName[41] = "bff";
	ParamName[42] = "bxx";
	ParamName[43] = "bfx";

	// Right ventricle parameters
	ParamName[44] = "rin0_rv";
	ParamName[45] = "rout0_rv";
	ParamName[46] = "f_rv";
	ParamName[47] = "RefVolRV";
	ParamName[48] = "WallVolRV";
	ParamName[49] = "km_rv";
	ParamName[50] = "kav_rv";
	ParamName[51] = "kv_rv";
	ParamName[52] = "c1_rv";
	ParamName[53] = "bff_rv";
	ParamName[54] = "bxx_rv";
	ParamName[55] = "bfx_rv";
	ParamName[56] = "br_rv";
	ParamName[57] = "bperp_rv";

	// Left atrium parameters
	ParamName[58] = "rin0_la";
	ParamName[59] = "rout0_la";
	ParamName[60] = "RefVolLA";
	ParamName[61] = "WallVolLA";
	ParamName[62] = "km_la";
	ParamName[63] = "kav_la";
	ParamName[64] = "kv_la";
	ParamName[65] = "c1_la";
	ParamName[66] = "br_la";
	ParamName[67] = "bperp_la";

	// Right atrium parameters
	ParamName[68] = "rin0_ra";
	ParamName[69] = "rout0_ra";
	ParamName[70] = "RefVolRA";
	ParamName[71] = "WallVolRA";
	ParamName[72] = "km_ra";
	ParamName[73] = "kav_ra";
	ParamName[74] = "kv_ra";
	ParamName[75] = "c1_ra";
	ParamName[76] = "br_ra";
	ParamName[77] = "bperp_ra";

	// Valve parameters (Mynard dissertation) area VA in cm2, length VL in cm
	ParamName[78] = "VA_aovo";
	ParamName[79] = "VL_aov";
	ParamName[80] = "VA_mivo";
	ParamName[81] = "VL_miv";
	ParamName[82] = "VA_puvo";
	ParamName[83] = "VL_puv";
	ParamName[84] = "VA_tcvo";
	ParamName[85] = "VL_tcv";
	ParamName[86] = "VA_aovc";
	ParamName[87] = "VA_mivc";
	ParamName[88] = "VA_puvc";
	ParamName[89] = "VA_tcvc";

	//Valve rate constants,1/(kPa s) from Mynard in cm2/dyn/s
	ParamName[90] = "Kaovo";
	ParamName[91] = "Kmivo";
	ParamName[92] = "Kpuvo";
	ParamName[93] = "Ktcvo";
	ParamName[94] = "Kaovc";
	ParamName[95] = "Kmivc";
	ParamName[96] = "Kpuvc";
	ParamName[97] = "Ktcvc";

	// Vein inertia parameters
	ParamName[98] = "VA_pv";
	ParamName[99] = "VL_pv";
	ParamName[100] = "VA_sv";
	ParamName[101] = "VL_sv";

	//aorta parameters: length cm, area at root cm2, exponential taper per cm, compliance (dyn/cm2)^-1, length of flow pulse (s)
	//distal impedance factor, matched = 1, Gaorta0, compliance at entrance Gaorta1, increase of compliance with distance
	ParamName[102] = "Laorta";
	ParamName[103] = "A00";
	ParamName[104] = "A01";
	ParamName[105] = "Gaorta0"; // entrance compliance in (dyn/cm2)^-1
	ParamName[106] = "Gaorta1"; // decrease in compliance with distance
	ParamName[107] = "tpulse"; //length of pulse
	ParamName[108] = "impedancefac";	//distal impedance factor, matched = 1

	// Systemic resistance
	ParamName[109] = "Rsp";
	ParamName[110] = "Rsv";

	// Pulmonary resistance
	ParamName[111] = "Rpa";
	ParamName[112] = "Rpp";
	ParamName[113] = "Rpv";

	// Compliances
	ParamName[114] = "Csa";
	ParamName[115] = "Csv";
	ParamName[116] = "Cpa";
	ParamName[117] = "Cpv";

	// External tissue pressure
	ParamName[118] = "Pext";

	// Initial pressures (kPa) - determine fluid loading
	ParamName[119] = "Psp0";
	ParamName[120] = "Psv0";
	ParamName[121] = "Ppp0";
	ParamName[122] = "Ppv0";

	//pericardium parameters
	ParamName[123] = "V0_peri";
	ParamName[124] = "s_peri";

	//LVAD parameters: max flow in cm3/s, stall pressure in mmHg, start time in cycles
	ParamName[125] = "q0_LVAD";
	ParamName[126] = "Ps_LVAD";
	ParamName[127] = "t0_LVAD";

	char fname[80];
	sprintf(fname, "parameters%03i.txt", run);
	string InputFile = fname;
	ifstream paramFile(InputFile.c_str(), ios::in);

	Skip = findString(paramFile, "Skip");		// Skip?
	 
	if (Skip != "on") {
		for (i = 1; i <= Nparam; i++) {
			Param0[i] = findValue(paramFile, ParamName[i]);
			ParamVar[i] = (int)findValue(paramFile, "ipar_" + ParamName[i]);
		}

		aortaVar = 0;
		for (i = 99; i <= 105; i++) if (ParamVar[i]) aortaVar = 1;	//set if any aorta parameter is being varied

		// Use LVRV model?
		odeSystem = findString(paramFile, "odeSystem");
		// include septal displacement (a5)?
		septalDisp = findString(paramFile, "septalDisp");
		// include RVhypertrophy (a60)?
		RVhypertrophy = findString(paramFile, "RVhypertrophy");
		// include LVAD? Choose on or off
		LVAD = findString(paramFile, "LVAD");
		// include pericardium? Choose on or off
		pericardium = findString(paramFile, "pericardium");
		// Optimize? Choose on or off
		Optimize = findString(paramFile, "Optimize");

		if (Optimize == "on") {			// Determine Nopt
			Nopt = 0;
			for (i = 1; i <= Nparam; i++) if (ParamVar[i] > Nopt) Nopt = ParamVar[i];
			printf("Nopt = %i\n", Nopt);
			if (Nopt > 0) {		//check which optimization parameters are used
				printf("If parameter is used, ia = 1 ... ");
				ia = ivector(1, Nopt);
				for (j = 1; j <= Nopt; j++) {
					ia[j] = 0;
					for (i = 1; i <= Nparam; i++) if (ParamVar[i] == j) ia[j] = 1;
					printf("ia[%i] = %i ", j, ia[j]);
				}
				printf("\n");
				// Read initial values and bounds for multiplicative factors
				xx00 = dvector(0, Nopt - 1);
				xx = dvector(1, Nopt);				//optimization parameters
				xxlb = dvector(1, Nopt);			//lower bounds
				xxub = dvector(1, Nopt);			//upper bounds
				findArray(paramFile, "xx", Nopt, xx00);
				for (i = 1; i <= Nopt; i++) xx[i] = xx00[i - 1];
				findArray(paramFile, "xxlb", Nopt, xx00);
				for (i = 1; i <= Nopt; i++) xxlb[i] = xx00[i - 1];
				findArray(paramFile, "xxub", Nopt, xx00);
				for (i = 1; i <= Nopt; i++) xxub[i] = xx00[i - 1];
				free_dvector(xx00, 0, Nopt - 1);
			}
		}
	}
	paramFile.close();
}
