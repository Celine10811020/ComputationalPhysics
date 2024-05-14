//勘根定理

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    int correct;
    double number;
    double x_1, x_2, y_1, y_2;
    double gap;
    double ans;

    printf("Please enter a number: ");
    scanf("%lf", &number);
    printf("Accurate to _ decimal places. ");
    scanf("%d", &correct);

    x_1 = 0;
    gap = 1;
    ans = 0;

    while(1)
    {
        //printf("\t\tgap: %lf, correct:%d\n", gap, correct);
        if(correct < 0)
        {
            ans = (x_1 + x_2) / 2;
            break;
        }
        if(ans != 0)
        {
            break;
        }

        y_1 = x_1*x_1 - number;

        while(1)
        {
            x_2 = x_1 + gap;

            y_2 = x_2*x_2 - number;

            //printf("\tx: %16.25lf, %16.25lf\n", x_1, x_2);
            //printf("\ty: %16.25lf, %16.25lf\n", y_1, y_2);
            //system("PAUSE");

            if(y_2 > 0)
            {
                break;
            }else if(y_2 == 0)
            {
                ans = x_2;
                break;
            }else
            {
                y_1 = y_2;
                x_1 = x_2;
            }
        }

        gap = gap * 0.1;
        correct--;
    }

    printf("Calculate Ans: %25.16lf\n", ans);
    printf("Standard  Ans: %25.16lf\n", sqrt(number));
    //printf("%16.25lf\n%16.25lf", ans_1, ans_2);
    return 0;
}
