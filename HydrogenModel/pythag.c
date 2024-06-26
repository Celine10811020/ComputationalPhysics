#include <math.h>
#include "nrutil.h"
#define SQR(x) ((x)*(x))

double pythag(double a, double b)
{
	double absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if(absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
