//gcc -o main main.c tqli.c tred2.c nrutil.c pythag.c hseigen.c matrix.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

#define PI 3.1415926536

double xi[16] = {0.048307665687738, 0.144471961582796, 0.239287362252137, 0.331868602282127, 0.421351276130635, 0.506899908932229, 0.58771575724076, 0.663044266930215, 0.732182118740289, 0.794483795967942, 0.849367613732570, 0.896321155766052, 0.9349060759377390, 0.964762255587506, 0.985611511545268, 0.997263861849481};
double wi[16] = {0.096540088514727, 0.095638720079274, 0.093844399080804, 0.091173878695763, 0.087652093004403, 0.083311924226946, 0.07819389578707, 0.072345794108848, 0.065822222776361, 0.058684093478535, 0.050998059262374, 0.042835898022218, 0.0342738629130210, 0.025392065309262, 0.016274394730905, 0.007018610009469};

double IntK(double a, double b, int n, int m);
double IntV1(double a, double b, int n, int m);
double IntV2(double a, double b, int n, int m);
double IntS(double a, double b, int n, int m);

int main()
{
	int i, j, n, m, l;
	double L;
	double ansOne, ansTwo;
	double ansK, ansV1, ansV2, ansS;
	double **H, **S, *eigenvalue, *q;

	n = 3;
	L = 10;
	l = 1;

	H = dmatrix(1, n, 1, n);
	S = dmatrix(1, n, 1, n);
	eigenvalue = dvector(1, n);
	q = dvector(1, n);

	for(i=1; i<=n; i++)
	{
			for(j=1; j<=n; j++)
			{
					ansK = IntK(0, L, i, j);
					ansV1 = IntV1(0, L, i, j);
					ansV2 = IntV2(0, L, i, j);
					H[i][j] = (-2*PI*PI*i*j/L/L)*ansK + ansV1 - (l*(l+1)/2)*ansV2;

					ansS = IntS(0, L, i, j);
					S[i][j] = ansS;
			}
	}

	for(i=1; i<=n; i++)
	{
			for(j=1; j<=n; j++)
			{
					printf("%d, %d\t%25.16f\n", i, j, H[i][j]);
			}
			printf("\n");
	}
	printf("\n");

	for(i=1; i<=n; i++)
	{
			for(j=1; j<=n; j++)
			{
					printf("%d, %d\t%25.16f\n", i, j, S[i][j]);
			}
			printf("\n");
	}
	printf("\n");

	//tred2(H, n, eigenvalue, q);
	//tqli(eigenvalue, q, n, H);

	hseigen(n, H, S, eigenvalue);

	for(i=1; i<=n; i++)
	{
		printf("value %d: %25.16lf\n", i, eigenvalue[i]);
	}

	printf("\n");
  for(i=1; i<=n; i++)
	{
    printf("Eigenvector %d: ", i);

    for(j=1; j<=n; j++)
		{
      printf("%25.16lf ", H[j][i]);
    }
    printf("\n");
	}

	free_dvector(q, 1, n);
	free_dvector(eigenvalue, 1, n);
	free_dmatrix(H, 1, n, 1, n);
	free_dmatrix(S, 1, n, 1, n);

	return 0;
}

//integral of kinetic from a to b
double IntK(double a, double b, int n, int m)
{
    int i;
    double L = b-a;
    double one, two;
    double sum = 0.0;

    for(i=0; i<16; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

        sum = sum + wi[i] * cos(2*PI*n*one / L) * one * one * cos(2*PI*m*one / L);
        sum = sum + wi[i] * cos(2*PI*n*two / L) * two * two * cos(2*PI*m*two / L);
    }

    return (b-a)/2.0 * sum;
}

//integral of potential of Coulomb from a to b
double IntV1(double a, double b, int n, int m)
{
    int i;
    double L = b-a;
    double one, two;
    double sum = 0.0;

    for(i=0; i<16; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

        sum = sum + wi[i] * sin(2*PI*n*one / L) * (one) * sin(2*PI*m*one / L);
        sum = sum + wi[i] * sin(2*PI*n*two / L) * (two) * sin(2*PI*m*two / L);
    }

    return (b-a)/2.0 * sum;
}

//integral of potential of angular momentum from a to b
double IntV2(double a, double b, int n, int m)
{
    int i;
    double L = b-a;
    double one, two;
    double sum = 0.0;

    for(i=0; i<16; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

        sum = sum + wi[i] * sin(2*PI*n*one / L) * sin(2*PI*m*one / L);
        sum = sum + wi[i] * sin(2*PI*n*two / L) * sin(2*PI*m*two / L);
    }

    return (b-a)/2.0 * sum;
}

//integral of Ju Zhen S from a to b
double IntS(double a, double b, int n, int m)
{
    int i;
    double L = b-a;
    double one, two;
    double sum = 0.0;

    for(i=0; i<16; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

        sum = sum + wi[i] * sin(2*PI*n*one / L) * one * one * sin(2*PI*m*one / L);
        sum = sum + wi[i] * sin(2*PI*n*two / L) * two * two * sin(2*PI*m*two / L);
    }

    return (b-a)/2.0 * sum;
}
