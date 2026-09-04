// updateActivation.cpp
// edited by TWS, November 2019
#include <math.h>
#include "nrutil.h"
#include <string>
using namespace std;

void updateActivation(double t) {
	extern int Nmu, Nnu;
	extern double At, At_la, Tc, Ta, kdfac, Tca_shift, Tca;
	
	double truncT = t - Tc + Ta - (int)((t - Tc + Ta )/ Tc) * Tc ;	//truncT = 0 at the start of systole

	// Left and right ventricle activation function
	At = 0.;
	if (truncT < Ta && truncT >= 0)  // eliminated kd, no need to do multiple calculations. TWS 6/2025
		At = sin(PI * pow(truncT / Ta, kdfac));	// Left ventricle activation function, modified 6/25

	// Left and right atrium activation function
	At_la = 0.;
	if (fabs(truncT - Tc + Tca_shift) <= Tca / 2.) At_la = cos(PI / Tca * (truncT - Tc + Tca_shift));
	else if (fabs(truncT + Tca_shift) <= Tca / 2.) At_la = cos(PI / Tca * (truncT + Tca_shift));
}
