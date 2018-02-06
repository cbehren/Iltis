//declares routines for the rejection method.

double rtsafe(void (*funcd)(double, double *, double *), double x1, double x2,
	double xacc);
void erf_min(double u_naught,double *fn,double *df);
double rejection_method(double x,double a);//following Laursen 2009
