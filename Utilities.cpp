#include "Utilities.H"
#include <cmath>
void crossproduct(const double x[3],const double y[3],double result[])
{
  
  result[0]=x[1]*y[2]-x[2]*y[1];
  result[1]=x[2]*y[0]-x[0]*y[2];
  result[2]=x[0]*y[1]-x[1]*y[0];
  
  
}
void normalize(double* vec)
{
  int i;
  
  double abs = norm(vec);
  
  for(i=0;i<3;i++)
    vec[i]=vec[i]/abs;
  return;
  
}

double norm(const double a[3])
{
    double r = 0.0;
    for(int i=0;i<3;i++)
        r += a[i]*a[i];
    return sqrt(r);
}
double scalar(const double a[3],const double b[3])
{
    double r = 0.0;
    for(int i=0;i<3;i++)
        r += a[i]*b[i];
    return r;
}
