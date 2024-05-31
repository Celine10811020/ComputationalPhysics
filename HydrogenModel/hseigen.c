#include <stdio.h>
#include <math.h>
#include "nrutil.h"

void hseigen(int n, double **A, double **B, double *eigenvalue)
{
	int i, j, m, n2, i2;
	double d, **L, **LT, **L1, **LT1, *p, *q;

	L = dmatrix(1, n, 1, n);
	LT = dmatrix(1, n, 1, n);
	L1 = dmatrix(1, n, 1, n);
	LT1 = dmatrix(1, n, 1, n);
	p = dvector(1, n);
	q = dvector(1, n);

	cholesky(B, L, n);
	inverse(L, L1, n);
	transpose(L, LT, n);
	inverse(LT, LT1, n);
	multiply(L1, A, n);
	multiply(A, LT1, n);

	tred2(LT1, n, p, q);
	tqli(p, q, n, LT1);

	for(i=1; i<=n; i++)
	{
		for(j=1; j<=n; j++)
		{
			A[i][j] = LT1[i][j];
		}
	}

	for(i=1; i<=n; i++)
	{
		eigenvalue[i] = p[i];
	}

	free_dvector(q, 1, n);
	free_dvector(p, 1, n);
	free_dmatrix(L, 1, n, 1, n);
	free_dmatrix(LT, 1, n, 1, n);
	free_dmatrix(L1, 1, n, 1, n);
	free_dmatrix(LT1, 1, n, 1, n);
}
