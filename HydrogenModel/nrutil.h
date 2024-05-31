#ifndef NRUTIL_H
#define NRUTIL_H

double *vector(int nl, int nh);
double *dvector(int nl, int nh);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
long double **ldmatrix(int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);

void free_vector(double *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_ldmatrix(long double **m, int nrl, int nrh, int ncl, int nch);
void nrerror(char error_text[]);
void free_ivector(int *v, int nl, int nh);

void tred2(double **a, int n, double d[], double e[]);
void tqli(double d[], double e[], int n, double **z);
double pythag(double a,double b);
void hseigen(int n, double **A, double **B, double *eigenvector);
void cholesky(double **B, double **L, int n);
void inverse(double **A, double **B, int n);
void multiply(double **A, double **B, int n);
void transpose(double **A, double **B, int n);

#endif
