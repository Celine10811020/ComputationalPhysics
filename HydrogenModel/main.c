//gcc -o main main.c tqli.c tred2.c nrutil.c pythag.c hseigen.c matrix.c sort.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

#define PI 3.1415926535897932384626433832795

double xi[48] = {0.016276744849602969579,0.048812985136049731112,
                 0.081297495464425558994,0.113695850110665920911,
                 0.145973714654896941989,0.178096882367618602759,
                 0.210031310460567203603,0.241743156163840012328,
                 0.273198812591049141487,0.304364944354496353024,
                 0.335208522892625422616,0.365696861472313635031,
                 0.395797649828908603285,0.425478988407300545365,
                 0.454709422167743008636,0.483457973920596359768,
                 0.511694177154667673586,0.539388108324357436227,
                 0.566510418561397168404,0.593032364777572080684,
                 0.618925840125468570386,0.644163403784967106798,
                 0.668718310043916153953,0.692564536642171561344,
                 0.715676812348967626225,0.738030643744400132851,
                 0.759602341176647498703,0.780369043867433217604,
                 0.800308744139140817229,0.819400310737931675539,
                 0.837623511228187121494,0.854959033434601455463,
                 0.871388505909296502874,0.886894517402420416057,
                 0.901460635315852341319,0.915071423120898074206,
                 0.927712456722308690965,0.939370339752755216932,
                 0.950032717784437635756,0.959688291448742539300,
                 0.968326828463264212174,0.975939174585136466453,
                 0.982517263563014677447,0.988054126329623799481,
                 0.992543900323762624572,0.995981842987209290650,
                 0.998364375863181677724,0.999689503883230766828};

double wi[48] = {0.032550614492363166242,0.032516118713868835987,
                 0.032447163714064269364,0.032343822568575928429,
                 0.032206204794030250669,0.032034456231992663218,
                 0.031828758894411006535,0.031589330770727168558,
                 0.031316425596861355813,0.031010332586313837423,
                 0.030671376123669149014,0.030299915420827593794,
                 0.029896344136328385984,0.029461089958167905970,
                 0.028994614150555236543,0.028497411065085385646,
                 0.027970007616848334440,0.027412962726029242823,
                 0.026826866725591762198,0.026212340735672413913,
                 0.025570036005349361499,0.024900633222483610288,
                 0.024204841792364691282,0.023483399085926219842,
                 0.022737069658329374001,0.021966644438744349195,
                 0.021172939892191298988,0.020356797154333324595,
                 0.019519081140145022410,0.018660679627411467385,
                 0.017782502316045260838,0.016885479864245172450,
                 0.015970562902562291381,0.015038721026994938006,
                 0.014090941772314860916,0.013128229566961572637,
                 0.012151604671088319635,0.011162102099838498591,
                 0.010160770535008415758,0.009148671230783386633,
                 0.008126876925698759217,0.007096470791153865269,
                 0.006058545504235961683,0.005014202742927517693,
                 0.003964554338444686674,0.002910731817934946408,
                 0.001853960788946921732,0.000796792065552012429};

double IntK(double a, double b, int n, int m);
double IntV1(double a, double b, int n, int m);
double IntV2(double a, double b, int n, int m);
double IntS(double a, double b, int n, int m);

int main()
{
	int i, j, n, m;
	double l, L;
	double ansOne, ansTwo;
	double ansK, ansV1, ansV2, ansS;
	double **H, **S, *eigenvalue, *q;
	double **K, **V1, **V2;

	n = 50;
	L = 50.0;
	l = 1.0;

	H = dmatrix(1, n, 1, n);
	K = dmatrix(1, n, 1, n);
	V1 = dmatrix(1, n, 1, n);
	V2 = dmatrix(1, n, 1, n);
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
					H[i][j] = ansK + ansV1 + ansV2;

					K[i][j] = ansK;
					V1[i][j] = ansV1;
					V2[i][j] = ansV2;

					ansS = IntS(0, L, i, j);
					S[i][j] = ansS;
			}
	}
/*
	printf("H\n");
	for(i=1; i<=n; i++)
	{
			for(j=1; j<=n; j++)
			{
					printf("%d, %d\t%25.16f\n", i, j, H[i][j]);
			}
			printf("\n");
	}
	printf("\n");

	printf("K\n");
	for(i=1; i<=n; i++)
	{
			for(j=1; j<=n; j++)
			{
					printf("%d, %d\t%25.16f\n", i, j, K[i][j]);
			}
			printf("\n");
	}
	printf("\n");

	printf("V1\n");
	for(i=1; i<=n; i++)
	{
			for(j=1; j<=n; j++)
			{
					printf("%d, %d\t%25.16f\n", i, j, V1[i][j]);
			}
			printf("\n");
	}
	printf("\n");

	printf("V2\n");
	for(i=1; i<=n; i++)
	{
			for(j=1; j<=n; j++)
			{
					printf("%d, %d\t%25.16f\n", i, j, V2[i][j]);
			}
			printf("\n");
	}
	printf("\n");
*/
	/*printf("S\n");
	for(i=1; i<=n; i++)
	{
			for(j=1; j<=n; j++)
			{
					printf("%d, %d\t%25.16f\n", i, j, S[i][j]);
			}
			printf("\n");
	}
	printf("\n");*/

	hseigen(n, H, S, eigenvalue);

  mergeSort(eigenvalue, H, 1, n);

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

    for(i=0; i<48; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

				sum = sum + wi[i] * cos(PI*n*one / L) * one*one * (1.0/2.0)*(PI*n/L)*(PI*m/L) * cos(PI*m*one / L);
        sum = sum + wi[i] * cos(PI*n*two / L) * two*two * (1.0/2.0)*(PI*n/L)*(PI*m/L) * cos(PI*m*two / L);
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

    for(i=0; i<48; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

				sum = sum + wi[i] * sin(PI*n*one / L) * (-1.0*one) * sin(PI*m*one / L);
        sum = sum + wi[i] * sin(PI*n*two / L) * (-1.0*two) * sin(PI*m*two / L);
    }

    return (b-a)/2.0 * sum;
}

//integral of potential of angular momentum from a to b
double IntV2(double a, double b, int n, int m)
{
    int i;
    double L = b-a;
		int l = 1;
    double one, two;
    double sum = 0.0;

    for(i=0; i<48; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

				sum = sum + wi[i] * sin(PI*n*one / L) * l*(l+1)/2 * sin(PI*m*one / L);
        sum = sum + wi[i] * sin(PI*n*two / L) * l*(l+1)/2 * sin(PI*m*two / L);
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

    for(i=0; i<48; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

				sum = sum + wi[i] * sin(PI*n*one / L) * one*one * sin(PI*m*one / L);
        sum = sum + wi[i] * sin(PI*n*two / L) * two*two * sin(PI*m*two / L);
    }

    return (b-a)/2.0 * sum;
}
