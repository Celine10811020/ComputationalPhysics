#include <math.h>
#include "nrutil.h"

void cholesky(double **B, double **L, int n)
{
  int i, j, k;
  double s;

  for(i=1; i<=n; i++)
  {
    for(j=1; j<(i+1); j++)
    {
      s = 0;

      for(k=1; k<j; k++)
      {
        s += L[i][k] * L[j][k];
      }

      L[i][j] = (i==j) ? sqrt(B[i][i] - s) : (1.0 / L[j][j] * (B[i][j] - s));
    }
  }
}

void inverse(double **A, double **B, int n)
{
  int i, j, k;
  double temp, **C;

  C = dmatrix(1, n, 1, n);

  for(i=1; i<=n; i++)
  {
    for(j=1; j<=n; j++)
    {
      C[i][j] = A[i][j];

      if(i == j)
      {
        B[i][j] = 1.0;
      }else
        B[i][j] = 0.0;
      }
  }

  for(i=1; i<=n; i++)
  {
    temp = C[i][i];

    for(j=1; j<=n; j++)
    {
      C[i][j] /= temp;
      B[i][j] /= temp;
    }
    for(k=1; k<=n; k++)
    {
      if(k != i)
      {
        temp = C[k][i];
        for(j=1; j<=n; j++)
        {
          C[k][j] -= C[i][j] * temp;
          B[k][j] -= B[i][j] * temp;
        }
      }
    }
  }

  free_dmatrix(C, 1, n, 1, n);
}

void transpose(double **A, double **B, int n)
{
	int i, j;

	for(i=1; i<=n; i++)
	{
		for(j=1; j<=n; j++)
		{
			B[j][i] = A[i][j];
		}
	}
}

void multiply(double **A, double **B, int n)
{
  int i, j, k;
  double sum, **C;

  C = dmatrix(1, n, 1, n);

  for(i=1; i<=n; i++)
  {
    for(j=1; j<=n; j++)
    {
      sum = 0;
      for(k=1; k<=n; k++)
      {
        sum += A[i][k] * B[k][j];
      }

      C[i][j] = sum;
    }
  }

  for(i=1; i<=n; i++)
  {
    for(j=1; j<=n; j++)
    {
      B[i][j] = C[i][j];
    }
  }

  free_dmatrix(C, 1, n, 1, n);
}
