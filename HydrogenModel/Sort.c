#include <stdlib.h>
#include "nrutil.h"

void mergeSort(double *eigenvalues, double **H, int head, int tail)
{
  if(head < tail)
  {
    int mid = head+(tail-head) / 2;

    mergeSort(eigenvalues, H, head, mid);
    mergeSort(eigenvalues, H, mid+1, tail);

    merge(eigenvalues, H, head, mid, tail);
  }
}

void merge(double *eigenvalue, double **H, int head, int mid, int tail)
{
  int lenA = mid - head + 1;
  int lenB = tail - mid;
  double *headSub = (double*)malloc(sizeof(double)*lenA);
  double *rightSub = (double*)malloc(sizeof(double)*lenB);
  double **Lvec = dmatrix(1, lenA, 1, tail);
  double **Rvec = dmatrix(1, lenB, 1, tail);

  int headIndex = 0;
  int rightIndex = 0;
  int j;

  for(headIndex=0; headIndex<lenA; headIndex++)
  {
    headSub[headIndex] = eigenvalue[head + headIndex];
    for(j=1; j<=tail; j++)
    {
      Lvec[headIndex+1][j] = H[j][head+headIndex];
    }
  }

  for(rightIndex=0; rightIndex<lenB; rightIndex++)
  {
    rightSub[rightIndex] = eigenvalue[mid + 1 + rightIndex];
    for(j=1; j<=tail; j++)
    {
      Rvec[rightIndex+1][j] = H[j][mid+1+rightIndex];
    }
  }

  headIndex = 0;
  rightIndex = 0;

  int writePointer = head;
  while(headIndex<lenA && rightIndex<lenB)
  {
    if (headSub[headIndex] <= rightSub[rightIndex])
    {
      eigenvalue[writePointer] = headSub[headIndex];
      for(j=1; j<=tail; j++)
      {
        H[j][writePointer] = Lvec[headIndex+1][j];
      }
      headIndex++;
    }else
    {
      eigenvalue[writePointer] = rightSub[rightIndex];
      for(j=1; j<=tail; j++)
      {
        H[j][writePointer] = Rvec[rightIndex+1][j];
      }
        rightIndex++;
    }
    writePointer++;
  }

  while(headIndex < lenA)
  {
    eigenvalue[writePointer] = headSub[headIndex];
    for(j=1; j<=tail; j++)
    {
      H[j][writePointer] = Lvec[headIndex+1][j];
    }
    headIndex++;
    writePointer++;
  }

  while (rightIndex < lenB)
  {
    eigenvalue[writePointer] = rightSub[rightIndex];
    for(j=1; j<=tail; j++)
    {
      H[j][writePointer] = Rvec[rightIndex+1][j];
    }
    rightIndex++;
    writePointer++;
  }

  free(headSub);
  free(rightSub);
  free_dmatrix(Lvec, 1, lenA, 1, tail);
  free_dmatrix(Rvec, 1, lenB, 1, tail);
}
