//integral sin(x)

#include <stdio.h>
#include <math.h>

#define PI 3.1415926536

double xi[8] = {0.095012509837637, 0.281603550779259, 0.458016777657227, 0.617876244402644, 0.755404408355003, 0.865631202387832, 0.944575023073233, 0.989400934991650};
double wi[8] = {0.189450610455069, 0.182603415044924, 0.169156519395003, 0.149595988816577, 0.124628971255534, 0.095158511682493, 0.062253523938648, 0.027152459411754};

double GaussLegendreIntegral(double a, double b);

int main()
{
    double L;
    double ans;

    L = 20;

    ans = GaussLegendreIntegral(0, L);

    printf("Result of integration: %lf\n", ans);

    return 0;
}

//integral from a to b
double GaussLegendreIntegral(double a, double b)
{
    int i;
    double one, two;
    double sum;

    sum = 0.0;

    for(i=0; i<8; i++)
    {
        one = (b-a)/2 * (xi[i]+1) + a;
        two = (b-a)/2 * (xi[i]*-1+1) + a;

        sum = sum + wi[i] * sin(one);
        sum = sum + wi[i] * sin(two);
    }

    return (b-a)/2 * sum;
}
