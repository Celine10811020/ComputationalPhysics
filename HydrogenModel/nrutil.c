#include <stdlib.h>
#include <malloc.h>
#ifndef STDIO
#include <stdio.h>
#define STDIO 1
#endif

void nrerror(char error_text[])
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *vector(int nl, int nh)
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

int *ivector(int nl, int nh)
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

double *dvector(int nl, int nh)
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

long double **ldmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	long double **m;

	m=(long double **) malloc((unsigned) (nrh-nrl+1)*sizeof(long double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(long double *) malloc((unsigned) (nch-ncl+1)*sizeof(long double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_vector(double *v, int nl, int nh)
{ free((char*) (v+nl)); }

void free_ivector(int *v, int nl, int nh)
{ free((char*) (v+nl)); }

void free_dvector(double *v, int nl, int nh)
{ free((char*) (v+nl)); }

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh; i>=nrl; i--)
	{
		free((char*) (m[i]+ncl));
	}
	free((char*) (m+nrl));
}

void free_ldmatrix(long double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh; i>=nrl; i--)
	{
		free((char*) (m[i]+ncl));
	}
	free((char*) (m+nrl));
}
