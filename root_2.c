#include <stdio.h>
#include <math.h>

int main()
{
    int i;
    double num;
    double x;

    printf("Please enter a number: ");
    scanf("%lf", &num);

    for(i=0; i<num; i++)
    {
        x = num - i*i;

        if(x < 0)
        {
            x = i - 1;
            break;
        }
    }

    //printf("\t%25.16lf\n", x);

    x = num;

    for(i=0; i<10; i++)
    {
        x = (x + num/x) / 2;

        //printf("\t%d. %25.16lf\n", i, x);
    }

    printf("Calculate Ans: %25.16lf\n", x);
    printf("Standard  Ans: %25.16lf\n", sqrt(num));

    return 0;
}
