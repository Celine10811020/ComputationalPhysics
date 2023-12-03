//y = x^2;

#include <stdio.h>

int main()
{
    int i;
    int n;
    double det_x;
    double sum_1, sum_2, ans_1, ans_2, ans_3;
    double y[1001];

    n = 1000;
    det_x = 1.0/n;
    sum_1 = 0;
    sum_2 = 0;

    for(i=0; i<=n; i++)
    {
        y[i] = i*det_x*i*det_x;
    }

    for(i=1; i<n; i+=2)
    {
        sum_1 = sum_1 + y[i]*det_x;
    }
    for(i=2; i<n; i+=2)
    {
        sum_2 = sum_2 + y[i]*det_x;
    }

    ans_1 = sum_1 + sum_2 + y[n]*det_x;
    ans_2 = (1.0/2) * (double)(2*sum_1 + 2*sum_2 + y[n]*det_x);
    ans_3 = (1.0/3) * (double)(4*sum_1 + 2*sum_2 + y[n]*det_x);

    printf("黎　曼　和：%25.16lf\n", ans_1);
    printf("梯 形 公式：%25.16lf\n", ans_2);
    printf("辛普森方法：%25.16lf\n", ans_3);

    return 0;
}
