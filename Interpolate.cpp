#include "Interpolate.H"
#include "stdlib.h"


void Interpolate::init()
{
  X=NULL;
  Y=NULL;
  N=0;
  bc_left=ABORT;
  bc_right=ABORT;
}
double Interpolate::left_boundary(double x)
{
  switch(bc_left)
  {
    case ABORT:
      std::cout << "Interpolate::right_boundary(): x = " << x << std::endl;
      abort();
      break;
    case PERIODIC:
      return y(x+(X[N-1]-X[0]));
      break;
    case CONSTANT:
      return Y[0];
      break;
    case LINEAR:
      //get the slope.
      double m=(Y[1]-Y[0])/(X[1]-X[0]);
      return m*(x-X[0])+Y[0];
      break;
    
  }
  return 0;
  
}
double Interpolate::right_boundary(double x)
{
    switch(bc_right)
  {
    case ABORT: 
      std::cout << "Interpolate::right_boundary(): x = " << x << std::endl;
      abort();
      break;
    case PERIODIC:
      return y(x-(X[N-1]-X[0]));
      break;
    case CONSTANT:
      return Y[N-1];
      break;
    case LINEAR:
      //get the slope.
      double m=(Y[N-1]-Y[N-2])/(X[N-1]-X[N-2]);
      return m*(x-X[N-1])+Y[N-1];
      break;
    
  }
  return 0;
  
  
}


Interpolate::Interpolate(double *x,double *y,long n, BOUNDARY_COND leftedge_bc, BOUNDARY_COND rightedge_bc)
{
  init();
  Init(x,y,n,leftedge_bc,rightedge_bc);
  
}
Interpolate::Interpolate()
{
  init();
}
void Interpolate::Init(double *x,double *y,long n, BOUNDARY_COND leftedge_bc, BOUNDARY_COND rightedge_bc)
{
  X = new double[n];
  Y = new double[n];
  N = n;


  for(int i=0;i<N;i++)
  {

    X[i] = x[i];
    Y[i] = y[i];
  }

  //make sure order is ascending.
  double temp=X[0];
  for(int i=1;i<N;i++)
    if(X[i]<temp)
    {
     std::cout << "Interpolate.cpp: X order must be ascending. (" <<X[i] << "<"<< temp << ")" << std::endl;
     abort();
    }
    else
      temp=X[i];
    
  //set up the boundary conditions.
  bc_left = leftedge_bc;
  bc_right = rightedge_bc;
  
}

double Interpolate::y(double x)
{
  if(x<X[0])
    return left_boundary(x);
  if(x>=X[N-1])
    return right_boundary(x);
  
  long right_index=0;
  while(x>X[right_index])
  { 
    right_index++;
  }
  long left_index = right_index-1;
//   std::cout << "index: " <<  left_index << " " << right_index << std::endl;
  if(left_index<0)
    return Y[0];
  
  //get the fractions.
  double lfrac=1.0-(x-X[left_index])/(X[right_index]-X[left_index]);
  if(lfrac>1.0)
   abort();
  double rfrac=1.0-lfrac;

  return lfrac*Y[left_index]+rfrac*Y[right_index];
  
 
}


Interpolate::~Interpolate()
{
 if(X)
   delete[] X;
 if(Y)
   delete[] Y;
  
}
