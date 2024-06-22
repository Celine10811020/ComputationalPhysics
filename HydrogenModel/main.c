//gcc -o main main.c tqli.c tred2.c nrutil.c pythag.c hseigen.c matrix.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

#define PI 3.1415926535897932384626433832795

double xi[32] = {0.024350292663424432509,0.072993121787799039450,
                 0.121462819296120554470,0.169644420423992818037,
                 0.217423643740007084150,0.264687162208767416374,
                 0.311322871990210956158,0.357220158337668115950,
                 0.402270157963991603696,0.446366017253464087985,
                 0.489403145707052957479,0.531279464019894545658,
                 0.571895646202634034284,0.611155355172393250249,
                 0.648965471254657339858,0.685236313054233242564,
                 0.719881850171610826849,0.752819907260531896612,
                 0.783972358943341407610,0.813265315122797559742,
                 0.840629296252580362752,0.865999398154092819761,
                 0.889315445995114105853,0.910522137078502805756,
                 0.929569172131939575821,0.946411374858402816062,
                 0.961008799652053718919,0.973326827789910963742,
                 0.983336253884625956931,0.991013371476744320739,
		 					 	 0.996340116771955279347,0.999305041735772139457};

double wi[32] = {0.048690957009139720383,0.048575467441503426935,
                 0.048344762234802957170,0.047999388596458307728,
                 0.047540165714830308662,0.046968182816210017325,
                 0.046284796581314417296,0.045491627927418144480,
                 0.044590558163756563060,0.043583724529323453377,
                 0.042473515123653589007,0.041262563242623528610,
                 0.039953741132720341387,0.038550153178615629129,
                 0.037055128540240046040,0.035472213256882383811,
                 0.033805161837141609392,0.032057928354851553585,
                 0.030234657072402478868,0.028339672614259483228,
                 0.026377469715054658672,0.024352702568710873338,
                 0.022270173808383254159,0.020134823153530209372,
                 0.017951715775697343085,0.015726030476024719322,
                 0.013463047896718642598,0.011168139460131128819,
                 0.008846759826363947723,0.006504457968978362856,
		 					 	 0.004147033260562467635,0.001783280721696432947};

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

	n = 20;
	L = 20.0;
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
					//H[i][j] = (PI*PI*i*j/L/L/-2.0)*ansK + ansV1 - (l*(l+1)/2)*ansV2;
					H[i][j] = ansK + ansV1 + ansV2;

					K[i][j] = ansK;
					V1[i][j] = ansV1;
					V2[i][j] = ansV2;

					ansS = IntS(0, L, i, j);
					S[i][j] = ansS;

					//printf("%d, %d\t%25.16f\n", i, j, (-2*PI*PI*i*j/L/L)*ansK);
			}
	}

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

    for(i=0; i<32; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

        //sum = sum + wi[i] * cos(PI*n*one / L) * one * one * cos(PI*m*one / L); // * one * one;
        //sum = sum + wi[i] * cos(PI*n*two / L) * two * two * cos(PI*m*two / L); // * two * two;

				sum = sum + wi[i] * cos(PI*n*one / L) * one*one * (1.0/2.0)*(PI*n/L)*(PI*m/L) * cos(PI*m*one / L); // * one * one;
        sum = sum + wi[i] * cos(PI*n*two / L) * two*two * (1.0/2.0)*(PI*n/L)*(PI*m/L) * cos(PI*m*two / L); // * two * two;
    }

    return (b-a)/2.0 * sum;
		//return sum;
}

//integral of potential of Coulomb from a to b
double IntV1(double a, double b, int n, int m)
{
    int i;
    double L = b-a;
    double one, two;
    double sum = 0.0;

    for(i=0; i<32; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

        //sum = sum + wi[i] * sin(PI*n*one / L) * (one) * sin(PI*m*one / L); // * one * one;
        //sum = sum + wi[i] * sin(PI*n*two / L) * (two) * sin(PI*m*two / L); // * two * two;

				sum = sum + wi[i] * sin(PI*n*one / L) * (-1.0*one) * sin(PI*m*one / L); // * one * one;
        sum = sum + wi[i] * sin(PI*n*two / L) * (-1.0*two) * sin(PI*m*two / L); // * two * two;
    }

    return (b-a)/2.0 * sum;
		//return sum;
}

//integral of potential of angular momentum from a to b
double IntV2(double a, double b, int n, int m)
{
    int i;
    double L = b-a;
		int l = 1;
    double one, two;
    double sum = 0.0;

    for(i=0; i<32; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

        //sum = sum + wi[i] * sin(PI*n*one / L) * sin(PI*m*one / L); // * one * one;
        //sum = sum + wi[i] * sin(PI*n*two / L) * sin(PI*m*two / L); // * two * two;

				sum = sum + wi[i] * sin(PI*n*one / L) * l*(l+1)/2 * sin(PI*m*one / L); // * one * one;
        sum = sum + wi[i] * sin(PI*n*two / L) * l*(l+1)/2 * sin(PI*m*two / L); // * two * two
    }

    return (b-a)/2.0 * sum;
		//return sum;
}

//integral of Ju Zhen S from a to b
double IntS(double a, double b, int n, int m)
{
    int i;
    double L = b-a;
    double one, two;
    double sum = 0.0;

    for(i=0; i<32; i++)
    {
        one = (b-a)/2.0 * (xi[i]+1) + a;
        two = (b-a)/2.0 * (xi[i]*-1+1) + a;

        //sum = sum + wi[i] * sin(PI*n*one / L) * one * one * sin(PI*m*one / L); // * one * one;
        //sum = sum + wi[i] * sin(PI*n*two / L) * two * two * sin(PI*m*two / L); // * two * two;

				sum = sum + wi[i] * sin(PI*n*one / L) * one*one * sin(PI*m*one / L); // * one * one;
        sum = sum + wi[i] * sin(PI*n*two / L) * two*two * sin(PI*m*two / L);
    }

    return (b-a)/2.0 * sum;
		//return sum;
}
