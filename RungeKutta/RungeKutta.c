/*
    y' + y = 0
    ans: y = e^(-x)

    let y' = -y = f(y)
*/

#include <stdio.h>
#include <math.h>

int main()
{
    int i, num;
    double h, k_1, k_2, k_3, k_4;
    double  x[1001], y_1[1001]/*1st order*/, y_2[1001]/*2nd order*/, y_4[1001]/*4th order*/;

    num = 1000;
    h = 1.0/num;

    for(i=0; i<=num; i++)
    {
        x[i] = i * h;
    }

    y_1[0] = 1.0;
    y_2[0] = 1.0;
    y_4[0] = 1.0;

    for(i=0; i<=num; i++)
    {
        k_1 = -(y_1[i]);
        y_1[i+1] = y_1[i] + h*k_1;

        k_1 = -(y_2[i]);
        k_2 = -(y_2[i] + h*k_1);
        y_2[i+1] = y_2[i] + (h/2)*(k_1+k_2);


        k_1 = -(y_4[i]);
        k_2 = -(y_4[i] + (h/2)*k_1);
        k_3 = -(y_4[i] + (h/2)*k_2);
        k_4 = -(y_4[i] + h*k_3);
        y_4[i+1] = y_4[i] + (h/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    }

    for(i=0; i<=num; i+=50)
    {
        printf("%.3lf\n  RK1:\t%.16lf\n  RK2:\t%.16lf\n  RK4:\t%.16lf\ne^(-x):\t%.16lf\n\n", x[i], y_1[i], y_2[i], y_4[i], exp(-x[i]));
    }

    return 0;
}
