//taylor

#include <stdio.h>
#include <math.h>

int main()
{
    int i, sign;
    double temp;
    double num;
    double ans;

    printf("Please enter a number: ");
    scanf("%lf", &num);

    ans = 1;
    temp = 1;
    sign = -1;

    //printf("\t%d. %25.16lf\n", 0, ans);

    for(i=2; i<30; i=i+2)
    {
        temp = temp * (num*num) / (i*(i-1));
        ans = ans + sign*temp;

        sign = sign * -1;

        printf("\t%d. %25.16lf\n", i/2, ans);
    }

    printf("Calculate Ans: %25.16lf\n", ans);
    printf("Standard  Ans: %25.16lf\n", cos(num));

    return 0;
}
