/***********************************************************************
Numerical Recipes routine: gaussj.c for Gauss-Jordan Elimination
Version 3.0, May 17, 2011. 1-based version, Feb. 2024
************************************************************************/
#include <math.h>
#include "nrutil.h"
#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}

void gaussj(double **a, int n, double *b)
{
	int *indxc, *indxr, *ipiv;
	int i, icol, irow, j, k, l, ll;
	int ii, jj;
	double big, dum, pivinv;
	indxc = ivector(1,n);
	indxr = ivector(1,n);
	ipiv = ivector(1,n);
	for (j = 1; j <= n; j++) ipiv[j] = 0;
	for (i = 1; i <= n; i++) {
		big = 0.0;
		for (j = 1; j <= n; j++) if (ipiv[j] != 1) for (k = 1; k <= n; k++) if (ipiv[k] == 0) {
			if (fabs(a[j][k]) >= big) {
				big = fabs(a[j][k]);
				irow = j;
				icol = k;
			}
		}
		else if (ipiv[k] > 1)
			nrerror("GAUSSJ: Singular Matrix-1");
		++(ipiv[icol]);
		if (irow != icol) {
			for (l = 1; l <= n; l++) SWAP(a[irow][l], a[icol][l]);
			SWAP(b[irow], b[icol]);
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0)
			nrerror("GAUSSJ: Singular Matrix-2");
		pivinv = 1.0 / a[icol][icol];
		a[icol][icol] = 1.0;
		for (l = 1; l <= n; l++) a[icol][l] *= pivinv;
		b[icol] *= pivinv;
		for (ll = 1; ll <= n; ll++) if (ll != icol) {
			dum = a[ll][icol];
			a[ll][icol] = 0.0;
			for (l = 1; l <= n; l++) a[ll][l] -= a[icol][l] * dum;
			b[ll] -= b[icol] * dum;
		}
	}
	for (l = n; l > 0; l--) if (indxr[l] != indxc[l]) for (k = 1; k <= n; k++) SWAP(a[k][indxr[l]], a[k][indxc[l]]);
	free_ivector(ipiv, 1, n);
	free_ivector(indxr, 1, n);
	free_ivector(indxc, 1, n);
}
#undef SWAP