#include <stdio.h>
#include <math.h>

#define PI 3.1415926536

// int sin(2*pi*i*x/L) * sin(2*pi*j*x/L) dx
// int sin(2*pi*i*x/L) * x * x * sin(2*pi*j*x/L) dx

double xi[16] = {0.048307665687738, 0.144471961582796, 0.239287362252137, 0.331868602282127, 0.421351276130635, 0.506899908932229, 0.58771575724076, 0.663044266930215, 0.732182118740289, 0.794483795967942, 0.849367613732570, 0.896321155766052, 0.9349060759377390, 0.964762255587506, 0.985611511545268, 0.997263861849481};
double wi[16] = {0.096540088514727, 0.095638720079274, 0.093844399080804, 0.091173878695763, 0.087652093004403, 0.083311924226946, 0.07819389578707, 0.072345794108848, 0.065822222776361, 0.058684093478535, 0.050998059262374, 0.042835898022218, 0.0342738629130210, 0.025392065309262, 0.016274394730905, 0.007018610009469};

double IntSin(double a, double b, int n, int m);

int main()
{
    int i, j ,n;
    double L;
    double ansOne, ansTwo;
    double ans;

    n = 10; //Ju Zhen nxn
    L = 20;

    for(i=1; i<n; i++)
    {
        //ansOne = IntSin(0, L, i);

        for(j=1; j<n; j++)
        {
            //ansTwo = IntSin(-1, 1, j);
            ans = IntSin(0, L, i, j);

            /*
            if(ansOne!=0 && ansTwo!=0)
            {
                printf("\t%d %d\n", i, j);
            }
            printf("%d, %d\t%f, %f\n", i, j, ansOne, ansTwo, ansOne);
            */

            printf("%d, %d\t%25.16f\n", i, j, ans);
        }

        printf("\n");
    }

    //printf("Result of integration: %lf\n", ans);

    return 0;
}

//integral from a to b
double IntSin(double a, double b, int n, int m)
{
    int i;
    double L = b-a;
    double oneN, twoN, oneM, twoM;
    //double one, two;
    double sum = 0.0;

    for(i=0; i<16; i++)
    {
        oneN = (b-a)/2 * (xi[i]+1) + a;
        twoN = (b-a)/2 * (xi[i]*-1+1) + a;
        oneM = (b-a)/2 * (xi[i]+1) + a;
        twoM = (b-a)/2 * (xi[i]*-1+1) + a;

        //one = (b-a)/2 * (xi[i]+1) + a;
        //two = (b-a)/2 * (xi[i]*-1+1) + a;

        sum = sum + wi[i] * sin(n*oneN) * sin(m*oneM);
        sum = sum + wi[i] * sin(n*twoN) * sin(m*twoM);

        //sum = sum + wi[i] * sin(2*PI*n*one / L) * sin(2*PI*m*one / L);
        //sum = sum + wi[i] * sin(2*PI*n*two / L) * sin(2*PI*m*two / L);

        //sum = sum + wi[i] * sin(2*PI*n*one / L) * one * one * sin(2*PI*m*one / L);
        //sum = sum + wi[i] * sin(2*PI*n*two / L) * two * two * sin(2*PI*m*two / L);
    }

    return L/(2*PI) * L/(2*PI) * (b-a)/2 * sum;
    //return (b-a)/2 * sum;
}

