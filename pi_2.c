#include <stdio.h>

int main()
{
	double ans;

	ans = 1;
	
	for(int i=2; i<=30000; i++)
	{
		if(i%2 == 0) //even
		{
		    ans = ans - 1.0 / (2*i - 1);
		}else //odd
		{
		    ans = ans + 1.0 / (2*i - 1);
		}
		
		printf("%25.16lf\n", ans*4);
	}
	
	return 0;
}