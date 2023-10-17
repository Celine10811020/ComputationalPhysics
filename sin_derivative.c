#include <stdio.h>
#include <math.h>

int main()
{
	double h = 0.01;
	double x = 0.1;
	double ans_1, ans_2, ans_3;

	ans_1 = ( sin(x+h) - sin(x) ) / h;
	ans_2 = ( sin(x+h) - sin(x-h) ) / (2*h);
	ans_3 = ( -sin(x+2*h) + 8*sin(x+h) - 8*sin(x-h) + sin(x-2*h) ) / (12*h);

	printf("ans_1 = %20.16lf\n", ans_1);
	printf("ans_2 = %20.16lf\n", ans_2);
	printf("ans_3 = %20.16lf\n", ans_3);
	printf("ans_4 = %20.16lf\n", cos(x));
	
	return 0;
}