//以2^n邊形計算pi值

#include <stdio.h>
#include <math.h>

double recursive(double an)
{
	return sqrt( pow(an, 2) + 1) + an;
}

int main()
{
	double an, an_1, pi;

	an = 1;
	
	for(int i=1; i<=30; i++)
	{
		an_1 = recursive(an);
		an = an_1;
		
		pi = 4 * pow(2, i) / sqrt( pow(an_1, 2) + 1 );

	    printf("%25.16lf\n", pi);
	}

	return 0;
}
