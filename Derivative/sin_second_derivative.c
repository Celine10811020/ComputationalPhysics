#include <stdio.h>
#include <math.h>

int main()
{
	double h = 0.01;
	double x = 0.1;
	double ans_1, ans_2, ans_3;

	ans_1 = ( sin(x) - 2*sin(x-h) + sin(x-2*h) ) / (h*h); //向後差
	ans_2 = ( sin(x+h) - 2*sin(x) + sin(x-h) ) / (h*h); //兩點中央差
	ans_3 = ( -sin(x+2*h) + 16*sin(x+h) -30*sin(x) + 16*sin(x-h) - sin(x-2*h) ) / (12*h*h); //四點中央差

	printf("ans_1 = %20.16lf\n", ans_1);
	printf("ans_2 = %20.16lf\n", ans_2);
	printf("ans_3 = %20.16lf\n", ans_3);
	printf("ans_4 = %20.16lf\n", -sin(x));

	return 0;
}

