#include "stdio.h"
#include <cmath>
#include "Utilities.H"
#include "RejectionMethod.H"
#include "RandomNumbers.H"
#ifdef _OPENMP
#include "omp.h"
#endif
double erf_const;
#ifdef _OPENMP
#pragma omp threadprivate (erf_const)
#endif
double rejection_method(double x,double a)//following Laursen 2009
{
 int tid = 0;
#ifdef _OPENMP
 tid = omp_get_thread_num();
#endif
  
  if(x!=x)
    {
    printf("Error: x is nan.\n");
    std::abort();
    }
  int sign=1;
  double u0=4.3;
  double x_cw=log(a);
  x_cw=1.59-0.60*x_cw-0.03*x_cw*x_cw;
  double right,left;
  double p;
  double u;
  double dummy;
  double fraction;
  int trycount=1;
  
  double xi=x;
  
  if(x<0.0)
  {
    sign=-1;
    xi*=-1;
    
  }
  if(xi>100*x_cw)
  {
    erf_const = 2.0*(RNG::uniform(tid)-0.5);
    u=rtsafe(&erf_min,-1.0e20,1.0e20,1.0e-4);
    return u*sign;
  }
  if(xi<0.2)
    u0=0.0;
   else if(xi<x_cw)
     u0=xi-0.01*pow(a,1.0/6.0)*exp(1.2*xi);

  double R;
  double theta0;
  while(1)//Following Zheng-Zheng 2002, p17
  {
    R=RNG::uniform(tid);
    theta0=atan((u0-xi)/a);
    p = ((theta0+M_PI/2.0)/((1.0-exp(-(u0*u0)))*
            theta0+(1.0+exp(-(u0*u0)))*M_PI/2.0));
    if(R<=p)
    {
        right=theta0;
        left=-M_PI/2.0;
    }
    else
    {
        right=M_PI/2.0;
        left=theta0;   
    }
    u = (right-left)*(RNG::uniform(tid))+left;
    u = a*tan(u)+xi;
    dummy = RNG::uniform(tid);
    fraction = exp(-(u*u));
    if(u>u0)
        fraction = fraction/(exp(-(u0*u0)));
    
    if(dummy < fraction)
    {
        u*=sign;
        return u;
    }
    trycount++;
    
    if(trycount==1000000)
    {
            printf("WARNING: Rejection says x %le u0 %le x_cw %le a %le u %le xi %le\n",x,u0,x_cw,a,u,xi);
    }
 } 
  
}

void erf_min(double u_naught,double *fn,double *df)
{

*df=1.128379167095513*exp(-u_naught*u_naught);
*fn=erf(u_naught)-erf_const;
}
